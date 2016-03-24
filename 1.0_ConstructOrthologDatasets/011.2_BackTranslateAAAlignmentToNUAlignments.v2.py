import re
import os
import glob
import multiprocessing as mp

gap = '-'
nucldir = 'OrthoMaMv9_12Taxa_NTraw_wMus_wGVACVOTCHUnaligned'
outdir = 'AAGuidedNTAlignments'
fileextension = '.seedaligned.fa.fixname' # change this to specify which files are accessed.
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
		seqout = seqs[n] + (len(longestseq)-len(seqs[n]))*'N'
		o += '%s %s\n' % (nameout, seqout)
	return(o)

def FormatDictionaryOfNucleotideSeqsToFasta(d):
	longestseq = max(d.values(), key=len)
	o = ''
	for k in d.keys():
		o += '>%s\n%s\n' % (k, d[k] + (len(longestseq)-len(d[k]))*gap) # should be d[k] second but this should fix a small problem where 1 dash is missing from seemingly random sequences.
	return(o)

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
	os.system('mkdir %s' % (outdir))
	
	filelist = []
	files = filter(os.path.isfile, os.listdir('.'))
	for filename in files:
		if extension.match(filename) != None:
			filelist.append(filename)	

	cmdlist = []
	for f in filelist:
		AAalign = ParseFasta(f)
		nucl = './%s/%s' % (nucldir, '_'.join(f.split('_')[:-1]) + '_NT.fasta.fixname') # making path to daughter directory and nucleotide filename that corresponds to this AA alignment
		nucseqs = ParseFasta(nucl)
		nucseqs = RemoveCharFromSeq(nucseqs, gap) # removes gap character from nucleotide alignment
		alignNU = AlignOriginalNucleotideSeqsWithOriginalAAAlignment(nucseqs, AAalign)
		output = FormatDictionaryOfNucleotideSeqsToFasta(alignNU)
		#output = FormatDictionaryOfNucleotideSeqsToPhylip(alignNU)
		outfile = './%s/%s' % (outdir, f[:-len(fileextension)] + '.AAGuided.NT.Alignment.fa')
		WriteOUT(outfile, output)			
		
main()
