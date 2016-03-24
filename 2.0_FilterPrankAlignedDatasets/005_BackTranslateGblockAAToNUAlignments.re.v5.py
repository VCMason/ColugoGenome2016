import re
import os
import glob
import multiprocessing as mp

outdir = 'GblockAAGuidedNUAlignments'
inNTdir = 'OrthoMaM_wMus_wGVACVOTCH_wPLO_NTUnaligned'
inAADir = '21TaxaAA_PrankAligned'

def WriteOUT(outfile, output):
	OUT = open(outfile, 'w')
	OUT.write(output)
	OUT.close()

def FormatDictionaryOfNucleotideSeqsToPhylip(seqs): ### Not Used ###
	longestseq = max(seqs.values(), key=len) # returns longest string in list
	longestname = max(seqs.keys(), key=len)
	o = '%d %d\n' % (len(seqs.keys()), len(longestseq))
	for n in seqs.keys():
		nameout = n + (len(longestname)-len(n))*' '
		seqout = seqs[n] + (len(longestseq)-len(seqs[n]))*'N'
		o += '%s %s\n' % (nameout, seqout)
	return(o)

def FormatDictionaryOfNucleotideSeqsToFasta(d, names):
	o = ''
	for k in names:
		o += '>%s\n%s\n' % (k, d[k])
	return(o)

def TrimAlignedNucleotidesUsingAAGblockCoordinate(start, end, alignNU):
	aligntrimNU = {}
	for s,e in map(None,start,end):
		for name in alignNU.keys():
			aligntrimNU[name] = aligntrimNU.get(name, '') + alignNU[name][(s-1)*3:e*3] # have to subtract 1 from s to include the first base. # Multiply by three because s and e are AA coordinates.
	return(aligntrimNU)

def AlignOriginalNucleotideSeqsWithOriginalAAAlignment(nucseqs, AAalign):
	alignNU = {}
	for name in AAalign.keys():
		pos = 1 # position in nucleotide alignment
		alignednuc = ''
		for AA in AAalign[name]:
			if AA == '-':
				alignednuc += '---' # add three dashes to the aligned nucleotide sequence if dash is in AA original (Non-Gblock trimmed) alignment
			else:
				alignednuc += nucseqs[name][pos-1:pos+2] # if AA present then add the representing nucleotide bases for that AA from original nucleotide sequence
				pos += 3
		alignNU[name] = alignednuc
	return(alignNU)

def IsolateStartAndEndCoordinatesOfGblocksAATxtFile(f):
	start = []
	end = []
	FILE = open(f, 'r')
	line = FILE.readline()
	while line:
		if line[:8] == 'Flanks: ':
			s = re.split(r'[\[\]]', line.strip()) # returns anything between two [] as element in list
			#print f
			#print s
			for elem in s:
				if elem.strip() != 'Flanks:' and elem.strip() != '':
					#print elem
					start.append(int(elem.split()[0]))
					end.append(int(elem.split()[1]))
			#print start
			#print end
		line = FILE.readline()
	FILE.close()
	return(start, end)

def RemoveCharFromSeq(d, char):
	newseqs = {}
	for k in d.keys():
		out = ''
		for nucl in d[k]:
			if nucl != char:
				out += nucl
		newseqs[k] = out
	return(newseqs)

def ParseFasta(fle):
	seqs = {}
	names = []
	FILE = open(fle, 'r')
	line = FILE.readline()
	while line:
		if line[0] == '>':
			n = line[1:].strip()
			seqs[n] = ''
			names.append(n)
		else:
			seqs[n] += line.strip()
		line = FILE.readline()
	FILE.close()
	
	return(seqs, names)

def main():
	os.system('mkdir %s' % (outdir))
	
	cwd = os.getcwd()
	path = cwd + '/GblocksFiltered_AA_Alignments/*.txt'
	filelist = glob.glob(path)

	cmdlist = []
	for f in filelist: # f is complete path to file with filename
		stump = f.split('/')[-1][:-len('.AA.Align.red.prank.fas-gb.txt')]
		a = '%s/%s/%s' % (cwd, inAADir, stump + '.AA.Align.red.prank.fas')
		nucl = '%s/%s/%s' % (cwd, inNTdir, stump + '_NT.fasta.fixname')
		o = stump + '_AAguidedNTAlign.gb.fa' # extension .fas is added by prank

		AAalign, aanames = ParseFasta(a)
		start, end = IsolateStartAndEndCoordinatesOfGblocksAATxtFile(f) # start and end contain coordinates of original alignment KEPT by Gblocks
		
		if len(start) > 0:
			nucseqs, nunames = ParseFasta(nucl)
			nucseqs = RemoveCharFromSeq(nucseqs, '-')
			alignNU = AlignOriginalNucleotideSeqsWithOriginalAAAlignment(nucseqs, AAalign)
			aligntrimNU = TrimAlignedNucleotidesUsingAAGblockCoordinate(start, end, alignNU)
			output = FormatDictionaryOfNucleotideSeqsToFasta(aligntrimNU, aanames)
			outfile = './%s/%s' % (outdir, o)
			WriteOUT(outfile, output)			
			
	
main()
