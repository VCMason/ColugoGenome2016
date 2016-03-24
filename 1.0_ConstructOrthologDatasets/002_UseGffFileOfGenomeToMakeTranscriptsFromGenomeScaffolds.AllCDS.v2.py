import os
import re

organismabv = 'TCH' # abreviation for organism
f2 = '/WMLab/victor/Orthologs/DetermineGVAOrthologousCDSSeqs/AllCDSSeqsFromAnnotationGffFile/TCH/GCF_000334495.1_TupChi_1.0_genomic.gff' # reference genome .gff file
f3 = '/WMLab/victor/Orthologs/DetermineGVAOrthologousCDSSeqs/AllCDSSeqsFromAnnotationGffFile/TCH/GCF_000334495.1_TupChi_1.0_genomic.fna' # Consensus sequences of aligned reads to reference
cutoff = 50.0 # percent base coverage required in assembled mRNA molecules


def PercentNotN(seq):
	count = 0.0
	for nucl in seq:
		if nucl != 'N':
			count += 1.0
	pcov = (count/len(seq))*100
	return(pcov)

def RevComp(seq):
	revcomp = ''
	rev = seq[::-1]
	#print rev
	for nucl in rev:
		if nucl.upper() == 'A':
			revcomp += 'T'
		elif nucl.upper() == 'T':
			revcomp += 'A'
		elif nucl.upper() == 'G':
			revcomp += 'C'
		elif nucl.upper() == 'C':
			revcomp += 'G'
		elif nucl.isalpha():
			revcomp += 'N'
		else:
			print 'Illegal Character: ' + nucl
	return(revcomp)

def RetrieveEnsemblAAByEnsemblGeneID():
	#F1 = open(f1, 'r')
	#g = {}
	#line = F1.readline()
	#while line:
	#	g[line.strip().split('\t')[2].split('|')[3]] = line.strip() # orthologous protein ID for reference seqeunce is key
	#	line = F1.readline()
	#F1.close()
	#print '1'
	allscaf = {} # Reads in whole consensus sequence file
	F3 = open(f3, 'r')
	line = F3.readline()
	while line:
		if line[0] == '>':
			name = line[1:].strip().split()[0] # isolate scaffold name
			allscaf[name] = allscaf.get(name, '')
		else:
			allscaf[name] += line.strip() # Fucking interleaved...
		line = F3.readline()
	F3.close()
	print '2'
	F2 = open(f2, 'r') # .gff annotation file
	g = {}
	pscaf = {} # scaffolds in the GVA genome, however by calling consensus seqs for CVO, CVO CNS seqs now have the same names as the GVA scaffolds.
	ps = {}
	pe = {}
	pstrand = {}
	porf = {}
	pgene = {}
	pdesc = {}
	line = F2.readline()
	while line:
		if line[0] == '#':
			pass
		elif line[0] != '#':
			s = line.strip().split('\t')
			scaffold, cds, start, end, strand, orf, notes = s[0], s[2], s[3], s[4], s[6], s[7], s[8]
			if cds == 'CDS':
				#print line.strip()
				notesObj = re.search( r'.+(gene=.+);(product=.+);protein_id=(.._.+\.1)', notes) # search raw string notes for name of protein
				try:
					pid = notesObj.group(3) # some lines are labeled as CDS sequences in column 3 but have no protein ID... which makes notesObj.group(3) to fail
				except:
					pass
				else:
					pid = notesObj.group(3)
					gene = notesObj.group(1)
					desc = notesObj.group(2)
					g[pid] = line.strip()
					pscaf[pid] = scaffold # always the same
					ps[pid] = ps.get(pid, []) + [int(start)] # need ordered list b/c first elem is for first exon
					pe[pid] = pe.get(pid, []) + [int(end)] # need ordered list b/c first elem is for first exon
					pstrand[pid] = strand # always the same
					porf[pid] = porf.get(pid, []) + [int(orf)] # need ordered list b/c first elem is for first exon
					pgene[pid] = gene # always the same
					pdesc[pid] = desc # always the same
		line = F2.readline()
	F2.close()
	print '3'
	allcds = {}
	removedBCcutoff = 0
	ERRORZ = ''
	for pid in g.keys(): # iterate through every protein id
		gate = 1 # gate open: stays open if we have seqs for all exons. If 1 exon missing, then gate == 0
		try:
			allscaf[pscaf[pid]]
		except:
			pass
		else:
			totcdsseq = ''
			c = 0
			for s, e, o in map(None, ps[pid], pe[pid], porf[pid]): # Iterate through ordered list of values in order of 1stexon, 2ndexon...	
				c += 1
				if pstrand[pid] == '+':
					seq = allscaf[pscaf[pid]][(s-1):(e)] # o (ORF) is not used, Only indicates which base 1st, 2nd, 3rd base (0,1,2) is the first base of a codon, NOT which is the first base used in an mRNA molecule for that exon.
					if len(seq) != 0 and len(allscaf[pscaf[pid]]) >= e:
						totcdsseq += seq
					else:
						gate = 0
						print 'Missing or incomplete Exon, Protein Excluded! PID:%s Exon:%d/%d has %d nucleotides! ExpectedLength:%d Start:%d End:%d ORF:%d Scaffold:%s' % (pid, c, len(ps[pid]), len(seq), e-s+1, s, e, o, pscaf[pid])
						ERRORZ += '>Missing or incomplete Exon, Protein Excluded! PID:%s Exon:%d/%d has %d nucleotides! ExpectedLength:%d Start:%d End:%d ORF:%d Scaffold:%s\n%s\n' % (pid, c, len(ps[pid]), len(seq), e-s+1, s, e, o, pscaf[pid], allscaf[pscaf[pid]])		
				elif pstrand[pid] == '-':
					seq = allscaf[pscaf[pid]][(s-1):(e)]
					if len(seq) != 0 and len(allscaf[pscaf[pid]]) >= e:
						revcomp = RevComp(seq)
						totcdsseq += revcomp
					else:
						gate = 0
						print 'Missing or incomplete Exon, Protein Excluded! RevComp PID:%s Exon:%d/%d has %d nucleotides! ExpectedLength:%d Start:%d End:%d ORF:%d Scaffold:%s' % (pid, c, len(ps[pid]), len(seq), e-s+1, s, e, o, pscaf[pid])
						ERRORZ += '>Missing or incomplete Exon, Protein Excluded! RevComp PID:%s Exon:%d/%d has %d nucleotides! ExpectedLength:%d Start:%d End:%d ORF:%d Scaffold:%s\n%s\n' % (pid, c, len(ps[pid]), len(seq), e-s+1, s, e, o, pscaf[pid], allscaf[pscaf[pid]])
			if totcdsseq != '':
				pcov = PercentNotN(totcdsseq)
				if gate == 1 and pcov >= cutoff: # gate == 1: then all exons were included in protein
					allcds['>%s|%s|%s|%s|Exons:%d|%s' % (pid, organismabv, pgene[pid], pstrand[pid], len(ps[pid]), pdesc[pid])] = totcdsseq
				elif pcov < cutoff:
					removedBCcutoff += 1
	print '4'			
	exclude = 'Number of Proteins excluded because < %.0f percent coverage: %d' % (cutoff, removedBCcutoff)
	keep = 'Number of CDS sequences synthesized from CVO map to GVA, GVA.gff, and CVO consensus sequences: %d' % (len(allcds.keys()))

	print exclude
	print keep

	outs = '%s\n%s\n##\n%s' % (keep, exclude, ERRORZ)
	ERR = open('%sSummaryExcludedProtiens.stat.fa' % (organismabv), 'w')
	ERR.write(outs)
	ERR.close()
	
	output = ''
	pnames = ''
	for key in allcds.keys():
		output += '%s\n%s\n' % (key, allcds[key])
		pnames += '%s\n' % (key[1:])
	
	OUT = open('%s_SyntesizedCDSSeqsUsingTCHGffFile_%d.fa' % (organismabv, len(allcds.keys())), 'w')
	OUT.write(output)
	OUT.close()

	POUT = open('%s_SynthmRNAProteinIDs.tab' % (organismabv), 'w')
	POUT.write(pnames)
	POUT.close()
	
	
RetrieveEnsemblAAByEnsemblGeneID()
