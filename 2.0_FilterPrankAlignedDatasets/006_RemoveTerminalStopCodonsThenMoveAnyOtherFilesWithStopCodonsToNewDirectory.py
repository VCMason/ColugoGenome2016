''' Author: Victor C Mason '''
''' Date: April 12, 2015 '''
''' Convert Phylip To Fasta'''

import os
import re
import multiprocessing as mp

stopcodons = ['TAG', 'TAA', 'TGA'] # all uppercase letters required
stopcodondir = 'AlignmentsWithStopCodons'
misschar = 'N' # missing data character
fileextension = '.fa' # change this to specify which files are accessed.
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

def FormatDictionaryOfNucleotideSeqsToFasta(d):
	o = ''
	for k in d.keys():
		o += '>%s\n%s\n' % (k, d[k])
	return(o)

def ScanForStopCodonInDictValues(seqs, stopcodons):
	gate = 1
	for k in seqs.keys():
		for stop in stopcodons:
			for i in range(0, len(seqs[k]), 3):
				if seqs[k][i:i+3].upper() == stop:
					gate = 0
					break
	return(gate)

def ReplaceStopCodonsFromLastThreeNucleotidesWithMissChar(seqs, stopcodons, misschar):
	trimseqs = {}
	for k in seqs.keys():
		revseq = seqs[k][::-1]
		revlast3 = ''
		count = 0
		gate = 0
		for base in revseq:
			if base != '-' and base.upper() != misschar:
				revlast3 += base
				gate = 1
			elif gate == 1:
				revlast3 += base
			count += 1
			if len(revlast3) == 3:
				break
		last3 = revlast3[::-1].upper()
		if last3 in stopcodons:
			if count == 3:
				trimseqs[k] = seqs[k][:-count] + misschar*3
			elif count > 3:
				trimseqs[k] = seqs[k][:-count] + misschar*3 + seqs[k][(-count+3):]
		else:
			#print 'NotStopCodon %s' % (last3)
			trimseqs[k] = seqs[k]
	return(trimseqs)

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

def main():
	filelist = []
	files = filter(os.path.isfile, os.listdir('.'))
	for filename in files:
		if extension.match(filename) != None:
			filelist.append(filename)	
	countnostop = 0
	countstop = 0
	for fle in filelist:
		seqs = ParseFasta(fle)
		newseqs = ReplaceStopCodonsFromLastThreeNucleotidesWithMissChar(seqs, stopcodons, misschar)
		gate = ScanForStopCodonInDictValues(newseqs, stopcodons)
		output = FormatDictionaryOfNucleotideSeqsToFasta(newseqs)
		if gate == 1:
			outfile = fle[:-len(fileextension)] + '.nostop.fasta'
			WriteOUT(outfile, output)
			countnostop += 1
		elif gate == 0:
			outfile = fle[:-len(fileextension)] + '.stop.fasta'
			if not os.path.exists(stopcodondir):
				os.makedirs(stopcodondir)
			WriteOUT('./%s/%s' % (stopcodondir, outfile), output)
			countstop += 1
	print 'Terminal stop codons replaced with missing data character: %s' % (misschar*3)
	print 'Number of Files with Zero Stop Codons: %d' % (countnostop)
	print 'Number of Files with One or More Stop Codons: %d' % (countstop)

main()