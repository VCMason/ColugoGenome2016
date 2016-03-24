import re
import os
import glob
import multiprocessing as mp

fileextension = '_AAalinged.gb.frame.filter.fa'
missingcharAA = 'X'
missingcharNU = 'N'
filteredgblocksAAdir = 'GblocksAAFastaRemFramAndFilteredAlignments_MedianSig2Min0.40Win10' # input AA alignments
gblocksnucldir = 'GblockAAGuidedNUAlignments' # input directory containing AA guided gblocks reduced NUCLEOTIDE phylip alignments
outputdir = 'GblocksNTFastaRemFramAndFilteredAlignments_MedianSig2Min0.40Win10' # output directory

def WriteOUT(outfile, output):
	OUT = open(outfile, 'w')
	OUT.write(output)
	OUT.close()

def FormatDictionaryOfNucleotideSeqsToPhylip(seqs, Ns):
	longestseq = max(seqs.values(), key=len) # returns longest string in list
	longestname = max(seqs.keys(), key=len)
	o = '%d %d\n' % (len(seqs.keys()), len(longestseq))
	for n in Ns:
		nameout = n + (len(longestname)-len(n))*' '
		seqout = seqs[n] + (len(longestseq)-len(seqs[n]))*'N'
		o += '%s %s\n' % (nameout, seqout)
	return(o)

def FormatDictionaryOfNucleotideSeqsToFasta(d, Ns):  
	o = ''
	for k in Ns:
		o += '>%s\n%s\n' % (k, d[k])
	return(o)

def MaskMissingCharFromAAOntoNU(nucseqs, AAalign, missingcharAA, missingcharNU):
	for name in AAalign.keys():
		pos = 1 # position in nucleotide alignment
		for AA in AAalign[name]:
			if AA == missingcharAA:
				nucseqs[name] = nucseqs[name][:pos-1] + missingcharNU*3 + nucseqs[name][pos+2:] # if AA present then add the representing nucleotide bases for that AA from original nucleotide sequence
				pos += 3
			else:
				pos += 3
	return(nucseqs)

def ParsePhylip(fle):
	seqs = {} # names keys, sequences value
	Ns = []	# names
	Ss = [] # Sequences
	FILE = open(fle, 'r')
	first = FILE.readline()
	line = FILE.readline()
	while line:
		n, s = line.strip().split()[0], line.strip().split()[1]
		seqs[n] = s
		Ns.append(n)
		Ss.append(s)
		line = FILE.readline()
	FILE.close()
	return(seqs, Ns) # , Ns, Ss

def ParseFasta(fle):
	seqs = {}
	Ns = []	# names
	FILE = open(fle, 'r')
	line = FILE.readline()
	while line:
		if line[0] == '>':
			n = line[1:].strip()
			seqs[n] = ''
			Ns.append(n)
		else:
			seqs[n] += line.strip()
		line = FILE.readline()
	FILE.close()
	return(seqs, Ns)

def main():
	os.system('mkdir %s' % (outputdir))
	
	cwd = os.getcwd()
	path = '%s/%s/*%s' % (cwd, filteredgblocksAAdir, fileextension)
	filelist = glob.glob(path)

	count = 0
	cmdlist = []
	for f in filelist: # f is complete path to file with filename # f = Filtered for high divergence windows Gblocks trimmed prank alinged AA seqs in phylip
		stump = f.split('/')[-1][:-len(fileextension)] 
		nucl = '%s/%s/%s' % (cwd, gblocksnucldir, stump + '_AAguidedNTAlign.gb.nostop.fasta') # nucleotide gblocks reduced AA guided alignment (has not yet been filtered for high divergence windows)
		outfilename = '%s_NT_AAguidedNTAlignment.gb.frame.filter.fa' % (stump)
		print nucl
		if os.path.isfile(nucl):
			AAalign, Ns = ParseFasta(f) # get seqs of High divergence filtered AA gblock alignments
			nucseqs, notusing = ParseFasta(nucl)
			maskNU = MaskMissingCharFromAAOntoNU(nucseqs, AAalign, missingcharAA, missingcharNU)
			output = FormatDictionaryOfNucleotideSeqsToFasta(maskNU, Ns)
			outfile = './%s/%s' % (outputdir, outfilename)
			WriteOUT(outfile, output)
		else:
			count += 1
			
	print 'Number of Nucleotide Files not present: %d' % (count)
main()
