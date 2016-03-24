''' Author: Victor C Mason '''


import os
import re

fileextension = '.v4.tab' # change this to specify which files are accessed.
extension = re.compile('.+%s$' % (fileextension))

def MakeDictFromTwoEntryTabFile(f):
	d = {}
	FILE = open(f, 'r')
	line = FILE.readline()
	while line:
		d[line.strip().split('\t')[0]] = line.strip().split('\t')[1]
		line = FILE.readline()
	FILE.close()
	return(d)
	
def One2OneOrthologByBitValue():

	filelist = []
	files = filter(os.path.isfile, os.listdir('.'))
	for filename in files:
		if extension.match(filename) != None:
			filelist.append(filename)
	dlist = []
	for f in sorted(filelist):
		d = MakeDictFromTwoEntryTabFile(f)
		dlist.append(d)
	
	dlistintersect = list(set.intersection(*(set(d.keys()) for d in dlist))) # finds shared keys in list of dictionaries and returns list of keys. Note: replace .keys() to make key and value the same and change list() to dict() to return a dictionary.
	#dsets = (set(d.keys) for d in dlist) # find shared keys in list of dictionaries python
	#dintersect = list(set.itersection(*dsets))

	output = 'Query_Sequence_Name\t' + '\t'.join(sorted(filelist)) + '\n'
	for k in dlistintersect:
		templist = [k]
		for d in dlist:
			templist.append(d[k])
		output += '\t'.join(templist) + '\n'
    
	OUT = open('SharedQuerySeqsGVACVOTCH_BLASTBitFilter.ThreeWayHMD_%d.tab' % (len(dlistintersect)), 'w')
	OUT.write(output)
	OUT.close()
    
One2OneOrthologByBitValue()
