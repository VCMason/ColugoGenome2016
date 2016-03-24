''' Author: Victor C Mason '''
''' Date: April 12, 2015 '''
''' Will concatenate fasta files that have the same individuals, and all sequences are the same length. '''

import os
import re
import multiprocessing as mp
import operator

reqtaxa = ['Homo', 'Pan', 'Gorilla', 'Pongo', 'Nomascus', 'Macaca', 'Callithrix', 'Otolemur', 'GVA', 'CVO', 'TCH', 'Ptilocercus', 'Mus', 'Rattus', 'Cavia', 'Ictidomys', 'Oryctolagus', 'Bos', 'Canis', 'Felis', 'Loxodonta'] # required number of taxa in each alignment

fileextension = '.fa' # change this to specify which files are accessed.
extension = re.compile('.+%s$' % (fileextension))

def WriteOUT(outfile, output):
	OUT = open(outfile, 'w')
	OUT.write(output)
	OUT.close()

def FormatDictionaryOfNucleotideSeqsToPhylip(seqs):
	o = ''
	longestseq = max(seqs.values(), key=len) # returns longest string in list
	longestname = max(seqs.keys(), key=len)
	o = '%d %d\n' % (len(seqs.keys()), len(longestseq))
	for n in seqs.keys():
		nameout = n + (len(longestname)-len(n))*' '
		seqout = seqs[n] + (len(longestseq)-len(seqs[n]))*'N'
		o += '%s %s\n' % (nameout, seqout)
	return(o)

def MakeSeqsInDictSameLength(seqs):
	d = {}
	longestseq = max(seqs.values(), key=len) # returns longest string in list
	for n in seqs.keys():
		seqout = seqs[n] + (len(longestseq)-len(seqs[n]))*'-'
		d[n] = seqout
	return(d)

def CheckIfSeqsAreSameLength(seqs):
	length = len(seqs.values()[0])
	gate = 1
	for v in seqs.values():
		if length != v:
			gate == 0
	return(gate)

def CheckIfIndividualsMatch(seqs, reqtaxa):
	gate = 1
	for n in seqs.keys():
		if n not in reqtaxa:
			gate = 0
	return(gate)

def ParseFasta(fle):
	seqs = {}
	FILE = open(fle, 'r')
	line = FILE.readline()
	while line:
		if line[0] == '>':
			n = line[1:].strip().split()[0] # .split() b/c of trimal
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
	
	allout = ''
	physeqs = {}
	names = [] # name of each partition window.
	genenum = []
	start = [] # starting position of each partition window.
	end = [] # ending position of each partition window.
	pos = 0	
	count = 0
	countgenes = 0
	s = 1
	concatseqs = {}
	l = []
	for fle in filelist:
		seqs = ParseFasta(fle)
		gate = CheckIfIndividualsMatch(seqs, reqtaxa)
		if gate == 1:
			temp = CheckIfSeqsAreSameLength(seqs)
			if temp == 0:
				l.append(fle)
			seqs = MakeSeqsInDictSameLength(seqs) # make sure all sequences are padded to be equal length before concatenating
			for n in seqs.keys():
				concatseqs[n] = concatseqs.get(n, '') + seqs[n]
				length = len(seqs[n])
			names.append(fle)
			start.append(count+1)
			count += length
			end.append(count)
			countgenes += 1
			genenum.append(countgenes)
		else:
			pass
		
	print 'Number of files detected: %d' % (countgenes)
	print l
	outfile = 'Concatenated%dOrthoGenes.21TaxawGVAwTCHwCVOwPLO.gb.Frame.PC50.Filter0.4.AllSites.NU.phy' % (countgenes) # name of concatenated alignment
	for n in concatseqs.keys():
		print len(concatseqs[n])
	
	output = FormatDictionaryOfNucleotideSeqsToPhylip(concatseqs)
	WriteOUT(outfile, output)
	
	output = ''
	for n, s, e, c in map(None, names, start, end, genenum):
		output += '%d\t%s\t%d\t%d\n' % (c, n, s, e)
	WriteOUT('GenePartitionCoordinates.txt', output)
	
	
main()
