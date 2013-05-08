# STAT 245E, Spring 2013
# Final Project    

# gets verified pairwise interactions of TF's from http://www.yeastgenome.org/


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
			#print str(tf.text).strip()
		return([iset.intersection(stfset), iset.intersection(set(ptflist)), 
			iset.intersection(set(atfdict.keys()))])
	else:
         return(None)




### MAIN   

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


# console outputs:
#str(table) == 'None':  YAP3
#str(table) == 'None':  GAT1
#YDR520C URC2
#str(table) == 'None':  SUT1
#str(table) == 'None':  RDS1
#RLR1 THO2
#RCS1 AFT1
#YML081W TDA9
#str(table) == 'None':  HAP1
#str(table) == 'None':  HAP3
#str(table) == 'None':  HAP2
     


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


# examines the list of TF's with problem
for tf in ptflist:
	contenturl = ''.join(["http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=", tf])
	webbrowser.open_new_tab(contenturl)

for tf in ptflist:
	print ''.join(["'", tf, "': \"\","])
	

# constructs a dictionary from TF to its url
urldict = {
	'YAP3': 'http://www.yeastgenome.org/cgi-bin/locus.fpl?rm=display_result&dbid=S000001001',
	'GAT1': "http://www.yeastgenome.org/cgi-bin/locus.fpl?rm=display_result&dbid=S000001873",
	'SUT1': "http://www.yeastgenome.org/cgi-bin/locus.fpl?rm=display_result&dbid=S000003130",
	'RDS1': "http://www.yeastgenome.org/cgi-bin/locus.fpl?rm=display_result&dbid=S000000703",
	'HAP1': "http://www.yeastgenome.org/cgi-bin/locus.fpl?rm=display_result&dbid=S000004246",
	'HAP3': "http://www.yeastgenome.org/cgi-bin/locus.fpl?rm=display_result&dbid=S000000117",
	'HAP2': "http://www.yeastgenome.org/cgi-bin/locus.fpl?rm=display_result&dbid=S000003206"	
}  

for tf in stfset:
	urldict[tf] = ''.join(["http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=", tf]) 

for tf in atfdict.keys():
	urldict[atfdict.get(tf)] = ''.join(["http://www.yeastgenome.org/cgi-bin/locus.fpl?locus=", tf]) 


# len(urldict)  # 118
                        

# gets interacting pairs and uses set to eliminate duplicates   
pairset = set([])
for tf in urldict.keys():  
	out = getInteractingTF(urldict.get(tf)) 
	if out is None:
		print "check", tf, ": there is no link to \"interaction\" tab on the webpage"
	else:
		for j in range(2): # stfset & ptflist  
			if len(out[j]) > 0:  
				for k in out[j]:
					pair = [tf.encode('ascii', 'ignore'), k]
					pair.sort()
					pairset.add(tuple(pair))  
		if len(out[2]) > 0:  # atfdict
			for k in out[2]:
				pair = [tf.encode('ascii', 'ignore'), 
					atfdict.get(k).encode('ascii', 'ignore')]
				pair.sort()
				pairset.add(tuple(pair)) 
   

# console outputs:
# check MATA1 : there is no link to "interaction" tab on the webpage

len(pairset)  # 265


# opens output file
try:
	fout = open("interactingPairs.txt", "w")
except IOError:
	print "output file could not be opened"
	

# writes pairset into a plain text file
for pair in pairset:	
	print >> fout, ' '.join(pair)

fout.close() 
                                     
# checks the number of lines in the output file
# note that "call" is from module "subprocess"
call(["wc", "-l", "interactingPairs.txt"])
#  265 ../output/interactingPairs.txt





