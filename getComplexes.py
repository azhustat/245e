# STAT 245E, Spring 2013
# Final Project



import re  # regular expression
import urllib2
import webbrowser          
from subprocess import call 
import sys

# appends module location to library path (!this depends on your computer)
sys.path.append("/opt/local/Library/Frameworks/Python.framework/Versions/2.7/lib/python2.7/site-packages")   

# for parsing html
from bs4 import BeautifulSoup


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

cstfset = ctfset.difference(set(cptflist).union(set(catfdict.values())))    
len(cstfset)  # 87
len(cptflist)  # 2

# Note that catfdict is mapping alias to standard name
len(catfdict)  # 31

cptflist
# ['HAP3', 'HAP2']
# have been checked in getInteractions.py, they are OK


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
	for tf in tflist:
		out = ['', '']
		if tf in catfdict.keys():
			sname = catfdict.get(tf)
		else:
			sname = tf		
		if sname in atfdict.keys():
			out[0] = key
			out[1] = atfdict.get(sname)
		elif sname in stfset.union(set(ptflist)):
			out[0] = key
			out[1] = sname			
		if out[0] != '':
			if results.has_key(out[0]):
				results.get(out[0]).append(out[1])
			else:
				results[out[0]] = [out[1]]
			

for key, tflist in results.items():
	print tflist
	print tfset.intersection(tflist)




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
#  13












		