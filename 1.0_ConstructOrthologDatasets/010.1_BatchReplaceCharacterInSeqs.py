''' Author: Victor C Mason '''
''' Date: April 12, 2015 '''
''' Convert Phylip To Fasta'''

import os
import re
import multiprocessing as mp
import operator

fileextension = '.fa' # change this to specify which files are accessed.
extension = re.compile('.+%s$' % (fileextension))

#stops = ['TGA', 'TAG', 'TAA']
old = '_'
new = 'X'

def ReplaceOldCharWithNewCharInSeqs(d, old, new):
	out = {}
	for k in d.keys():
		seq = d[k].replace(old, new) # replaces old char with new char in all values of dictionary d
		out[k] = seq
	return(out)

def ReplaceStopCodonsWithMissDataInSeqs(d):  ### NOT USED ###
	new = {}
	for k in d.keys():
		seq = d[k]
		for stop in stops: # replace all stop codons in sequence string with missingdatachar
			seq = seq.replace(stop, 'NNN')
		new[k] = seq
	return(new)

def WriteOUT(outfile, output):
	OUT = open(outfile, 'w')
	OUT.write(output)
	OUT.close()

def FormatDictionaryOfNucleotideSeqsToPhylip(seqs): ### NOT USED ###
	longestseq = max(seqs.values(), key=len) # returns longest string in list
	longestname = max(seqs.keys(), key=len)
	o = '%d %d\n' % (len(seqs.keys()), len(longestseq))
	for n in seqs.keys():
		nameout = n + (len(longestname)-len(n))*' '
		seqout = seqs[n] + (len(longestseq)-len(seqs[n]))*missdatachar
		o += '%s %s\n' % (nameout, seqout)
	return(o)

def FormatDictionaryOfNucleotideSeqsToFasta(d):
	o = ''
	for k in d.keys():
		o += '>%s\n%s\n' % (k, d[k])
	return(o)

def ParseFasta(fle):
	seqs = {}
	FILE = open(fle, 'r')
	line = FILE.readline()
	while line:
		if line[0] == '>':
			n = line[1:].strip()
			seqs[n] = ''
		else:
			seqs[n] += line.strip()
		line = FILE.readline()
	FILE.close()
	return(seqs)

def ParsePhylip(fle): ### NOT USED ###
	seqs = {}
	FILE = open(fle, 'r')
	first = FILE.readline()
	line = FILE.readline()
	while line:
		n, s = line.strip().split()[0], line.strip().split()[1]
		seqs[n] = s
		line = FILE.readline()
	FILE.close()
	return(seqs)

def main(old, new):
	filelist = []
	files = filter(os.path.isfile, os.listdir('.'))
	for filename in files:
		if extension.match(filename) != None:
			filelist.append(filename)	

	for fle in filelist:
		outfile = fle #   [:-len(fileextension)] + '.fasta'
		#seqs = ParsePhylip(fle)
		print fle
		seqs = ParseFasta(fle)
		#seqs = ReplaceStopCodonsWithMissDataInSeqs(seqs)
		seqs = ReplaceOldCharWithNewCharInSeqs(seqs, old, new)
		output = FormatDictionaryOfNucleotideSeqsToFasta(seqs)
		WriteOUT(outfile, output)

main(old, new)