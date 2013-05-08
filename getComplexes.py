# STAT 245E, Spring 2013
# Final Project

# gets the list of TF's which interact with the TF's that are part of complexes
# (validation data for clustering results)

import re  # regular expression
import urllib2
import webbrowser          
from subprocess import call 
import sys

# appends module location to library path (!this depends on your computer)
sys.path.append("/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages")   

# for parsing html
from bs4 import BeautifulSoup


       
## function def:

'''Input: string of the url of the page of the given TF 
on http://www.yeastgenome.org/

Returns a list of three lists, which are the intersection of the set of 
interacting TF/genes with the three subsets of our list, TF's using standard
name, TF's which cannot be accessed directly url with their names, TF's
using alias, respectively
Note that the three subset (form a partition) are represented by different data
structures: stfset ("s": standard name), ptflist ("p": problem), and atfdict
("a": alias)
If there is no link to the "interactions" tab on the web page, a "None" object
will be returned and a message will be printed to console 
'''                                                      

def getInteractingTF(contenturl):
	# finds interaction link
	soup = BeautifulSoup(urllib2.urlopen(contenturl).read())
	spans = soup.findAll('span', attrs={'class': 'right small'})
	# len(spans)  # 5   
	ilink = ''
	for s in spans:
		match = re.findall('href=\"(http://.*?interactions.pl\?.*?)\">', str(s))  
		if len(match) > 0:
			ilink = match[0]
	# finds out lists of interacting TF's for this given TF
	if len(ilink) > 0: 
		ilisturl = ''.join([ilink, "&rm=no_js"])
		soup = BeautifulSoup(urllib2.urlopen(ilisturl).read())      
		table = soup.find('table', attrs={'class': 'result_table_border'})
		rows = soup.findAll('tr')
		# uses set to store the interacting TF's to eliminate duplicates automatically
		iset = set([])                                                                
		# omits the first two rows
		for index in range(2, len(rows)):
			tf = rows[index].findAll('td')[1]  # python is 0-based
			iset.add(str(tf.text).strip()) 
		return([iset.intersection(stfset),iset.intersection(set(ptflist)), 
			iset.intersection(set(atfdict.keys()))])
	else:
		return(None)





### MAIN

# reads in negation words
# "complexes.txt" is generated in getComplexes.R
fc = open("./complexes.txt", "r")
clines = fc.readlines()
fc.close()

cdict = {}
for line in clines:
	out = line.strip() # removes leading and trailing whitespace
	out = out.split(" ")
	if cdict.has_key(out[0].upper()):
		cdict.get(out[0].upper()).append(out[1].upper())
	else:
		cdict[out[0].upper()] = [out[1].upper()]

		
# len(cdict)  # 19, check
for key, tflist in cdict.items():
	print key, ":", len(tflist)



# check if they are the standard names  
ctfset = set()
for key, tflist in cdict.items():
	ctfset |= set(tflist)


# len(ctfset)  # 120, check
# there are 120 TF's


# deals with alias
# if the given TF name is alias, then we still get to the page with the 
# standard name, e.g,:  
# http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=SKN7
# http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=BRY1
# http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=POS9

# pattern to find standard name
pat = '''Standard Name.+?<span class=['\"]i['\"]>(.+?)</span>''' 
# list of TF's using alias, keys are the standard names and values are alias
catfdict = {}              
# list of TF's with problem opening url directly
cptflist = []

	
for tf in ctfset:
	contenturl = ''.join(["http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=", tf])
	soup = BeautifulSoup(urllib2.urlopen(contenturl).read())
	table = soup.find('table', attrs={'class': 'summary'})  # the first one  
	if (str(table) == 'None'):
		print "str(table) == 'None': ", tf      
		cptflist.append(tf)  
	else:                                    
		# finds standard name
		sname = re.findall(pat, str(table))  
		sname = sname[0].strip()
		if sname != tf:
			print tf, sname
			catfdict[tf] = sname



#TAF90 TAF5
#TAF17 TAF9
#TAF19 TAF13
#ARGR3 ARG82
#ARGR2 ARG81
#ARGR1 ARG80
#TAF145 TAF1
#RPB6 RPO26
#RPB1 RPO21
#SIG1 MOT2
#TBP1 SPT15
#TAF47 TAF3
#A1 MATA1
#str(table) == 'None':  HAP3
#MED3 PGD1
#TAF30 TAF14
#ADA4 GCN5
#ANC1 TAF14
#MED9 CSE2
#SRB11 SSN8
#TAF61 TAF12
#SRB10 SSN3
#MED10 NUT2
#SRB9 SSN2
#RAD25 SSL2
#ADA1 HFI1
#ADA3 NGG1
#TFG3 TAF14
#TAF150 TAF2
#TAF25 TAF10
#TAF60 TAF6
#TAF67 TAF7
#SSN6 CYC8
#CDC46 MCM5
#str(table) == 'None':  HAP2


cstfset = ctfset.difference(set(cptflist).union(set(catfdict.keys())))    
len(cstfset)  # 85
len(cptflist)  # 2

# Note that catfdict is mapping alias to standard name
len(catfdict)  # 33
len(set(catfdict.values())) # 31
# there are alias for the SAME standard nmae
#TAF30 TAF14
#ANC1 TAF14 
#TFG3 TAF14


cptflist
# ['HAP3', 'HAP2']
# have been checked in getInteractions.py, they are OK



# constructs a dictionary from TF to its url
urldict = {
	'HAP3': "http://www.yeastgenome.org/cgi-bin/locus.fpl?rm=display_result&dbid=S000000117",
	'HAP2': "http://www.yeastgenome.org/cgi-bin/locus.fpl?rm=display_result&dbid=S000003206"	
}  

for tf in cstfset:
	urldict[tf] = ''.join(["http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=", tf]) 


for tf in catfdict.keys():
	urldict[tf] = ''.join(["http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=", catfdict.get(tf)]) 


len(urldict)  #  120
                        

###########################################

tfset = set([])      

# gets the list of TF's of interest and stores the names in a set
contenturl = "http://oreganno.org/tfview/cgi-bin/tflist.pl?speciesname=Saccharomyces%20cerevisiae"
soup = BeautifulSoup(urllib2.urlopen(contenturl).read())
rows = soup.findAll('tr')                          

# omits the first row, which contains column names 
for index in range(1, len(rows)):
	tf = rows[index].findAll('td')[0]  # python is 0-based
	tfname = tf.contents[0].string.strip()
	tfset.add(tfname)

    
# excluds "UNKNOWN" 
tfset = tfset.difference(set([u'UNKNOWN']))  
# len(tfset)  # 118

# deals with alias
# if the given TF name is alias, then we still get to the page with the 
# standard name, e.g,:  
# http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=SKN7
# http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=BRY1
# http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=POS9

# pattern to find standard name
pat = '''Standard Name.+?<span class=['\"]i['\"]>(.+?)</span>''' 
# list of TF's using alias, keys are the standard names and values are alias
atfdict = {}              
# list of TF's with problem opening url directly
ptflist = []

# check if they are the standard names  
for tf in tfset:
	contenturl = ''.join(["http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=", tf])
	soup = BeautifulSoup(urllib2.urlopen(contenturl).read())
	table = soup.find('table', attrs={'class': 'summary'})  # the first one  
	if (str(table) == 'None'):
		print "str(table) == 'None': ", tf      
		ptflist.append(tf)  
	else:                                    
		# finds standard name
		sname = re.findall(pat, str(table))  
		sname = sname[0].strip()
		if sname != tf:
			print tf, sname
			atfdict[sname] = tf



# partitions the set of TF's, tfset, into three subsets:
# TF's using standard name, TF's which cannot be accessed directly 
# url with their names, TF's using alias  
# Note that the three subset (form a partition) are represented by different data
# structures: stfset ("s": standard name), ptflist ("p": problem), and atfdict
# ("a": alias)  
# Recall that atfdict uses standard names as keys
stfset = tfset.difference(set(ptflist).union(set(atfdict.values())))    
len(stfset)  # 107
len(ptflist)  # 7
len(atfdict)  # 4


		
##################################      
   

# finds intersections
results = {}

for key, tflist in cdict.items():
	for ctf in tflist:
		out = getInteractingTF(urldict.get(ctf)) 
		if out is None:
			print "check", ctf, ": there is no link to \"interaction\" tab on the webpage"
		else:
			for j in range(2): # stfset & ptflist  
				if len(out[j]) > 0:  
					for tf in out[j]:
						if results.has_key(key):
							results.get(key).add(tf)
						else:
							results[key] = set([tf])
			if len(out[2]) > 0:  # atfdict
				for tf in out[2]:
					if results.has_key(key):
						results.get(key).add(atfdict.get(tf))
					else:
						results[key] = set([atfdict.get(tf)])


                                                                                           
# check A1 : there is no link to "interaction" tab on the webpage
# check TAF145 : there is no link to "interaction" tab on the webpage     		

for key, tflist in results.items():
	print tflist
	print set(tflist).difference(tfset), "\n"
	




# opens output file
try:
	fout = open("complexesList.txt", "w")
except IOError:
	print "output file could not be opened"


# writes pairset into a plain text file
for key, tflist in results.items():
	for tf in tflist:
		print >> fout, ' '.join([key, tf])


fout.close() 

# checks the number of lines in the output file
# note that "call" is from module "subprocess"
call(["wc", "-l", "complexesList.txt"])
#  109







 