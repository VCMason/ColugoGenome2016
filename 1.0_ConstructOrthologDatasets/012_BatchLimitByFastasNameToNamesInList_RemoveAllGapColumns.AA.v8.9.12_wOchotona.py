''' Author: Victor C Mason '''
''' Date: May 14, 2015 '''
''' Output: '''
''' Summary: Counts all INDELs shared between all individuals in alignment. '''

import os
import re
import timeit

start = timeit.default_timer()

namelist = ['Homo', 'Pan', 'Macaca', 'Callithrix', 'Otolemur', 'Oryctolagus', 'Ochotona', 'Canis', 'Felis', 'GVA', 'CVO', 'TCH']
outdir = 'ReducedNames_12TaxawGVACVOTCH_wOchotona_AAAlignments_6'
missdatachar = 'X' # use capital letter
fileextension = '_AA.seedaligned.fa.fixname' # change this to specify which files are accessed.
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
		o += '>%s\n' % (k)
		step = 60
		for i in range(0, len(d[k]), step):
			o += '%s\n' % (d[k][i:i+step])
	return(o)

def ReplaceGapsInFrontAndBackOfEachSequenceWithMissDataCharFromDict(seqs): ### NOT USED ###
	Xseqs = {}
	for key in seqs.keys():
		
		Obj = re.search(r'(^-+[^-]).+', seqs[key])
		try:
			Obj.group(1)
		except:
			frontfix = seqs[key]
		else:
			front = len(Obj.group(1)) - 1
			frontfix = re.sub(r'^-+', missdatachar*front, seqs[key]) # minus 1 because matching NON dash character too in regex.
			
		Obj = re.search(r'.+([^-]-+$)', seqs[key])
		try:
			Obj.group(1)
		except:
			Xseqs[key] = frontfix
		else:
			back = len(Obj.group(1)) - 1
			Xseqs[key] = re.sub(r'-+$', missdatachar*back, frontfix)
	return(Xseqs)

def RemoveAllGapColumnsFromTransposedDict(dtrans, names, seqs, fle):
	redseqs = {}
	goodcolumns = [] #record columns that pass complete deletion criterion
	count = 0
	for col in dtrans.keys():
		gate = 1
		if len(set(dtrans[col])) == 1 and dtrans[col][0] == '-': # for each nucl in list of nucl for that column
			gate = 0
		elif len(set(dtrans[col])) == 1 and dtrans[col][0] == ' ': # for each nucl in list of nucl for that column
			gate = 0
		if gate == 1: # keep column
			goodcolumns.append(count)
		count += 1
	for name in names: # construct reduced seqs excluding columns of only gaps
		s = ''
		for col in goodcolumns:
			#print '%s, %s, %s\n' % (fle, name, col) # use this print for error testing
			s += seqs[name][col]
		redseqs[name] = s
	return(redseqs)

def TransposeValuesOfDictionary(d):
	names = []
	dtrans = {}
	names = sorted(d.keys())
	for k in names:
		count = 0
		for char in d[k]:
			dtrans[count] = dtrans.get(count, []) + [char] # makes dict with column as key (starting @ 0), and list of each char in column for each individual (in the same order as names) as the value.
			count += 1
	return(dtrans, names)

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

def main():
	os.system('mkdir ./%s' % (outdir))
	
	filelist = []
	files = filter(os.path.isfile, os.listdir('.'))
	for filename in files:
		if extension.match(filename) != None:
			filelist.append(filename)
	
	count = 0
	for fle in sorted(filelist): # sort files by filename because of the number prefix to keep them in order [00001_*, 00002_*, ... 07118_*]		
		seqs = ParseFasta(fle)
		newseqs = {}
		for n in namelist: # limit taxa before limit by columns with only gaps
			try:
				seqs[n]
			except:
				pass
			else:
				newseqs[n] = seqs[n]
		if len(newseqs.keys()) != len(namelist):
			pass
		else:
			dtrans, names = TransposeValuesOfDictionary(newseqs)
			redseqs = RemoveAllGapColumnsFromTransposedDict(dtrans, names, newseqs, fle) # send transposed seqs, ordered list of names that represent order of individuals' bases in value of dtrans, and original sequences. 		
			#outseqs = ReplaceGapsInFrontAndBackOfEachSequenceWithMissDataCharFromDict(redseqs)
		
			output = FormatDictionaryOfNucleotideSeqsToFasta(redseqs) # in v7.2 and below this was outseqs not redseqs
			WriteOUT('%s/%s/%s' % (os.getcwd(), outdir, fle[:-len(fileextension)] + '.AA.Align.red.fa'), output)
			count += 1
			
	print 'Number alignments queried: %d' % (len(filelist))
	print 'Number of alignments with all desired taxa: %d' % (count)

main()

stop = timeit.default_timer()

print stop-start
