''' Author: Victor C Mason '''
''' Date: September 07, 2015 '''
''' Input: aligned fasta '''
''' Output: '''
''' Summary: Removes every third nucleotide in each sequence'''

import os
import re
import timeit

start = timeit.default_timer()

outdir = 'FirstAndSecondSites_NU'
outdir1 = 'FirstSitesOnly_NU'
outdir2 = 'SecondSitesOnly_NU'
outdir3 = 'ThirdSitesOnly_NU'

missdatachar = 'N'
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


def RemoveEveryThirdCharacterInDictSeqs(seqs):
	removethirdposition = {}
	for n in seqs.keys():
		count = 0
		newseq = ''
		for char in seqs[n]:
			if count%3 == 2:
				pass
			else:
				newseq += char
			count += 1			
				
		removethirdposition[n] = newseq
	return(removethirdposition)

def Isolate3rdSites(seqs):
	removethirdposition = {}
	for n in seqs.keys():
		count = 0
		newseq = ''
		for char in seqs[n]:
			if count%3 == 0 or count%3 == 1:
				pass
			else:
				newseq += char
			count += 1			
				
		removethirdposition[n] = newseq
	return(removethirdposition)

def Isolate2ndSites(seqs):
	removethirdposition = {}
	for n in seqs.keys():
		count = 0
		newseq = ''
		for char in seqs[n]:
			if count%3 == 0 or count%3 == 2:
				pass
			else:
				newseq += char
			count += 1			
				
		removethirdposition[n] = newseq
	return(removethirdposition)

def Isolate1stSites(seqs):
	removethirdposition = {}
	for n in seqs.keys():
		count = 0
		newseq = ''
		for char in seqs[n]:
			if count%3 == 1 or count%3 == 2:
				pass
			else:
				newseq += char
			count += 1			
				
		removethirdposition[n] = newseq
	return(removethirdposition)

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
	os.system('mkdir ./%s' % (outdir))
	os.system('mkdir ./%s' % (outdir1))
	os.system('mkdir ./%s' % (outdir2))
	os.system('mkdir ./%s' % (outdir3))
	
	filelist = []
	files = filter(os.path.isfile, os.listdir('.'))
	for filename in files:
		if extension.match(filename) != None:
			filelist.append(filename)
	
	for fle in sorted(filelist): # sort files by filename because of the number prefix to keep them in order [00001_*, 00002_*, ... 07118_*]		
		seqs = ParseFasta(fle)
		twothree = RemoveEveryThirdCharacterInDictSeqs(seqs)
		output = FormatDictionaryOfNucleotideSeqsToFasta(twothree)
		outfile = fle[:-len(fileextension)] + '.red.nostop.twothree.fasta'
		outpath = '%s/%s/%s' % (os.getcwd(), outdir, outfile)
		WriteOUT(outpath, output)
		
		one = Isolate1stSites(seqs)
		output = FormatDictionaryOfNucleotideSeqsToFasta(one)
		outfile = fle[:-len(fileextension)] + '.red.nostop.one.fasta'
		outpath = '%s/%s/%s' % (os.getcwd(), outdir1, outfile)
		WriteOUT(outpath, output)

		two = Isolate2ndSites(seqs)
		output = FormatDictionaryOfNucleotideSeqsToFasta(two)
		outfile = fle[:-len(fileextension)] + '.red.nostop.two.fasta'
		outpath = '%s/%s/%s' % (os.getcwd(), outdir2, outfile)
		WriteOUT(outpath, output)
		
		three = Isolate3rdSites(seqs)
		output = FormatDictionaryOfNucleotideSeqsToFasta(three)
		outfile = fle[:-len(fileextension)] + '.red.nostop.three.fasta'
		outpath = '%s/%s/%s' % (os.getcwd(), outdir3, outfile)
		WriteOUT(outpath, output)
		
main()
		
stop = timeit.default_timer()
		
print stop-start
