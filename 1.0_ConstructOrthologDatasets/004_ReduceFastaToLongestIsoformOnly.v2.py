import os
import re

f2 = 'TCH_SyntesizedCDSSeqsUsingTCHGffFile_32729_AA.fas'


def RetrieveEnsemblAAByEnsemblGeneID():
	allprot = {} # Reads in whole file
	F2 = open(f2, 'r')
	line = F2.readline()
	while line:
		if line[0] == '>':
			name = line.strip()
			allprot[name] = allprot.get(name, '')
		else:
			allprot[name] += line.strip() # Fucking interleaved...
		line = F2.readline()
	F2.close()
    
	allgene = {}
	allgenename = {}
	for key in allprot.keys():
		geneObj = re.search( r'\|(product=.+) isoform X\d+$', key) # raw strings as r'expression'
		if geneObj: # if there is a match object... continue
			geneid = geneObj.group(1) # collapsing all 'isoform X#' to one gene and keeping the longest isoform of that gene.
		else:
			geneid = key.split('|')[5] # this should be accessing the product= .+ field
		try:
			allgene[geneid]
		except:
			allgene[geneid] = allprot[key] # if no allgene[geneid] then initiate one with value of allprot[key]
			allgenename[geneid] = key
		else:
			if len(allprot[key]) > len(allgene[geneid]): # keep the longest isoform of each geneid
				allgene[geneid] = allprot[key]
				allgenename[geneid] = key
        
	output = ''
	for gid in allgene.keys():
		output += '%s\n%s\n' % (allgenename[gid], allgene[gid])
	
	print 'Number of Total number of protein AA isoforms: %d' % (len(allprot.keys()))
	print 'Number of Unique genes with longest AA Isoforms: %d' % (len(allgene.keys()))
	
	OUT = open('.'.join(f2.split('.')[:-1]) + '.LongestIsoforms_%d.fas' % (len(allgene.keys())), 'w')
	OUT.write(output)
	OUT.close()
	
	
RetrieveEnsemblAAByEnsemblGeneID()
