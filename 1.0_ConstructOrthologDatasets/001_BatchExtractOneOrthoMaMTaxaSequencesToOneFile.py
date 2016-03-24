''' Author: Victor C Mason '''
''' Date: April 12, 2015 '''
''' Convert Phylip To Fasta'''

import os
import re
import multiprocessing as mp
import operator

name = 'Canis' # Only records fasta seqs if Homo is name
gapcharacter = '-'
outfile = 'ExtractOrthoMaM_CanisOne2OneWithHMMDP_AA.fa'
fileextension = '.fasta' # change this to specify which files are accessed.
extension = re.compile('.+%s$' % (fileextension))

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

def FormatDictionaryOfNucleotideSeqsToFastaUseFilename(d, f):
	o = ''
	for k in d.keys():
		o += '>%s_%s\n%s\n' % (f[:-len(fileextension)], k, d[k])
	return(o)

def FormatDictionaryOfNucleotideSeqsToFasta(d):
	o = ''
	for k in d.keys():
		o += '>%s\n%s\n' % (k, d[k])
	return(o)

def RemoveCharFromSeq(d, gap):
	newseqs = {}
	for k in d.keys():
		out = ''
		for nucl in d[k]:
			if nucl != gap:
				out += nucl
		newseqs[k] = out
	return(newseqs)

def ParseFastaRecordIfNameIsKey(fle, name): ### NOT USED ###
	seqs = {}
	FILE = open(fle, 'r')
	line = FILE.readline()
	while line:
		if line[0] == '>':
			n = line[1:].strip()
			if n == name:
				seqs[n] = ''
		else:
			if n == name:
				seqs[n] += line.strip()
		line = FILE.readline()
	FILE.close()
	return(seqs)

def ParseFasta(fle): ### NOT USED ###
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

def main():
	filelist = []
	files = filter(os.path.isfile, os.listdir('.'))
	for filename in files:
		if extension.match(filename) != None:
			filelist.append(filename)	
	output = ''
	for fle in filelist:
		seqs = ParseFastaRecordIfNameIsKey(fle, name)
		seqs = RemoveCharFromSeq(seqs, gapcharacter)
		temp = FormatDictionaryOfNucleotideSeqsToFastaUseFilename(seqs, fle)
		output += temp
	WriteOUT(outfile, output)

main()
