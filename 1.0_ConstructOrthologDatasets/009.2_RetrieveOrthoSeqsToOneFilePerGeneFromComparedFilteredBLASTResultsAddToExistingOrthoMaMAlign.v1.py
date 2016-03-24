''' Author: Victor C Mason '''


import os
import re

indir = 'OrthoMaMv9_HMGMC_NTraw'
orthofile = 'SharedQuerySeqsGVACVOTCH_BLASTBitFilter.ThreeWayHMD_5466.tab' # This is the file produced after running CompareAllFilteredBLASTResultsRetrieveCommonQuerySequences.v1.py
filelist = ['CVO_SyntesizedCDSSeqsUsingGVAGffFile_27054.fas','GVA_SyntesizedCDSSeqsUsingGVAGffFile_32137.fas','TCH_SyntesizedCDSSeqsUsingTCHGffFile_32729.fas'] # This list must be in order of the columns in file produced from CompareAllFilteredBLASTResultsRetrieveCommonQuerySequences.v1.py. They should be alphabetical

#fileextension = '.fas' # change this to specify which files are accessed.
#extension = re.compile('.+%s$' % (fileextension))

def AppendStringToFile(string, f):
	with open(f, 'a') as FILE:
		FILE.write(string)
	return()

def RetrieveSeqsFromFileWithDictionaryList(dlist, orthofile):
	newd = {}
	FILE = open(orthofile, 'r')
	line = FILE.readline() # header line
	line = FILE.readline()
	while line:
		output = '' # need one file per gene
		s = line.strip().split('\t')
		n = s[0]
		for d,o in map(None, dlist, s[1:]):
			output += d[o.split('|')[0]] + '\n'
		newd[n] = output
		#WriteOUT('%s.OrthoUnalign.NT.fa' % (n), output)
		line = FILE.readline()
	FILE.close()
	return(newd)

def ParseFasta(fle):
	seqs = {}
	FILE = open(fle, 'r')
	line = FILE.readline()
	while line:
		if line[0] == '>':
			n = line[1:].strip().split('|')[0] # limit to only protein ID
			seqs[n] = line # for other programs seqs[n] = ''
		else:
			seqs[n] += line.strip()
		line = FILE.readline()
	FILE.close()
	return(seqs)

def Main(filelist, orthofile):

	#filelist = []
	#files = filter(os.path.isfile, os.listdir('.'))
	#for filename in files:
	#	if extension.match(filename) != None:
	#		filelist.append(filename)

	count = 1
	dlist = []
	for f in filelist:
		seqs = ParseFasta(f)
		dlist.append(seqs)
		count += 1
	
	newd = RetrieveSeqsFromFileWithDictionaryList(dlist, orthofile)
	print type(newd)
	for n in newd.keys():
		AppendStringToFile(newd[n], './%s/%s' % (indir, '_'.join(n.split('_')[:-1]) + '_NT.fasta'))

Main(filelist, orthofile)
