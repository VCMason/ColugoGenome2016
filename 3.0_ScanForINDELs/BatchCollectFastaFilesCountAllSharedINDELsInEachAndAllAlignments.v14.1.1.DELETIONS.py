''' Author: Victor C Mason '''
''' Date: May 14, 2015 '''
''' Output: '''
''' Summary: v14.1 only requires the required individuals from the hypothesis file to be present to report the INDEL. '''
''' Summary: v13.7 requires that the number of individuals with the INDEL must be >= number of individuals from hypothesis minus number of individuals with missing data from hypothesis that overlaps the INDEL. '''


import os
import re
import shutil
import multiprocessing as mp
import timeit

start = timeit.default_timer()

hypotheses = 'IndelPhylogeneticHypotheses_AllTaxa.v13.1.txt'
missdatachar = 'X' # use capital letter
min = 3 # lowest number of taxa allowed to grant support to a hypothesis
notallowedalone = ['TCH', 'Tupaia'] # not allowed to represent a shared indel with only these two taxa
outdir = 'GenesWithINDELsSupportingAHypothesis_v14_1_1_Primatomorpha'
ordertaxaoutput = ['Homo', 'Pan', 'Gorilla', 'Pongo', 'Nomascus', 'Papio', 'Macaca', 'Callithrix', 'Tarsius', 'Otolemur', 'Microcebus', 'GVA', 'TCH', 'Tupaia', 'Mus', 'Rattus', 'Dipodomys', 'Cavia', 'Ictidomys', 'Oryctolagus', 'Ochotona', 'Ovis', 'Bos', 'Tursiops', 'Sus', 'Vicugna', 'Mustela', 'Ailuropoda', 'Canis', 'Felis', 'Equus', 'Pteropus', 'Myotis', 'Sorex', 'Erinaceous', 'Dasypus', 'Choloepus', 'Procavia', 'Loxodonta', 'Echinops', 'Macropus', 'Sarcophilus', 'Monodelphis', 'Ornithorhynchus'] # NOTE: !!! ONLY TAXA IN THIS LIST WILL BE PRINTED OUT TO HIGHLIGHTED ALIGNED FASTA FILES FOR HIGHLIGHTED INDELs !!!
fileextension = '.red.fa' # change this to specify which files are accessed.
extension = re.compile('.+%s$' % (fileextension))

def WriteOUT(outfile, output):
	OUT = open(outfile, 'w')
	OUT.write(output)
	OUT.close()
	
def FormatDictionaryOfNucleotideSeqsToFasta(d):
	o = ''
	for k in ordertaxaoutput: # d.keys()
		try:
			d[k]
		except:
			pass
		else:
			o += '>%s\n%s\n' % (k, d[k])
	return(o)

def FormatOutput(d, f): #list of dictionaries, output, fasta(ish) name (in this case gene name)
	o = ''
	o += '>%s\n' % (f)
	for key in sorted(d.keys(), key=len)[::-1]: # sort keys by length and reverse so decending in length
		o += '%d\t%s\n' % (d[key], key)
	return(o)

def FormatOutputLists(listdict, listnames): #dictionary, output, fasta(ish) name (in this case gene name)
	o = ''
	count = 0
	for d, g in map(None, listdict, listnames):
		o += '>%s\n' % (g)
		for key in sorted(d.keys(), key=len)[::-1]:
			o += '%s\t%s\n' % (key, ','.join(d[key]))
			count += 1
	return(o)

def FormatGeneFastasToHighlightINDELs(seqs, indelspergene):
	listoutseqs = []
	
	starts = []
	ends = []
	#print indelspergene.values()
	for listBS in indelspergene.values():
		for span in listBS:
			#print span
			#print type(span)
			starts.append(int(span.split('_')[0]))
			ends.append(int(span.split('_')[1]))
		            
	ssort = sorted(starts)
	esort = [x for (y,x) in sorted(zip(starts,ends))] # sort ends based on starts	
	
	count = 0
	for s,e in map(None, ssort, esort): # access sorted list of start and end coordinates for this gene
		outseqs = {} # only make this shit around one INDEL at a time, so if there are three indels in a gene then you will have three files of the same gene with one of the three indels highlighted in each gene file.
		for n in seqs.keys():
			outseqs[n] = seqs[n][:s-10] + missdatachar*10 + seqs[n][s-10:e+10] + missdatachar*10 + seqs[n][e+10:]
		listoutseqs.append(outseqs)
		count += 1
	return(listoutseqs)		
	#ssort = sorted(s)[::-1]
	#esort = [x for (y,x) in sorted(zip(s,e))][::-1] # sort e based on s, then reverse list
	
def RecordINDELsThatSupportEachHypothesis(h, hreq, sharedspans, overlapmiss, oneindeldict, gene, seqs):
	indelspergene = {} # key is taxa combinations, value is spans (start_end) coordinates of INDELs			
	for s in sharedspans.keys():
		for hyp in h.keys():
			nummiss = 0
			try:
				overlapmiss[s]
			except:
				pass
			else:
				nummiss = len(list(set(overlapmiss[s]) & set(h[hyp])))
			
			if len(sharedspans[s]) >= min and len(sharedspans[s]) < len(seqs) and len(sharedspans[s]) >= (len(hreq[hyp])) and set(sharedspans[s]) <= set(h[hyp]) and set(hreq[hyp]) <= set(sharedspans[s]) and len(list(set(sharedspans[s]) ^ set(notallowedalone))) >= 1: # this is basically the gateway to hypothesis heaven. # requires the shared DEL be present more than min taxa. requires that the shared DEL be present in less than numtaxa - 1. and len(list(set(taxlist) & set(sharedspans[s]))) >= 1 ##### len(sharedspans[s]) >= (len(hreq[hyp])) should be replaced with len(sharedspans[s]) >= (len(h[hyp]-nummiss)) for more stringency #####
				oneindeldict[hyp] = oneindeldict.get(hyp, 0) + 1 # INDEL supports hypothesis so record it.
				k = '_'.join(sorted(list(sharedspans[s]))) # change our set back into list to sort alphabetically before
				indelspergene[k] = indelspergene.get(k, []) + [s]
				
	return(oneindeldict, indelspergene)

def CollectAllContinuousStringsOfCharFromDict_Span(char, d):
	spans = {}
	for n in d.keys(): # record spans of all indels from all taxa
		matchObjs = [m for m in re.finditer(r'%s+' % (char), d[n])] # record all indel events
		for m in matchObjs:
			s = '%d_%d' % (m.span()[0], m.span()[1])
			spans[s] = spans.get(s, []) + [n]
	return(spans)
	
def RecordIndelsFromDict(d):
	
	spans = CollectAllContinuousStringsOfCharFromDict_Span('-', d) # records all spans(INDELs) ('start_end') as key and all taxa that have that span(INDEL) as list in value.
	spanmiss = CollectAllContinuousStringsOfCharFromDict_Span(missdatachar, d) # records all spans(of missdatachar) ('start_end') as key and all taxa that have that span(of missdatachar) as list in value.
	
	sharedspans = {}
	for s in spans.keys():
		spans[s] = set(spans[s]) # make list of taxa a set of unique taxa (as value) that has INDEL (span) as key, also variable naming is completely ridic.
		if len(spans[s]) >= 2: # then it means the span(INDEL) is shared with at least two taxa.
			sharedspans[s] = spans[s] # save set of taxa to sharedspans if there is more than one taxa in value spans[s]
			
	overlapmiss = {} # indel span as key, and value is set()tuple of all taxa that have a missingdatachar inside of indel span.
	for mspan in spanmiss.keys(): # identify which taxa have missing data within the range(span) of the INDELs identified above.
		ms, me = int(mspan.split('_')[0]), int(mspan.split('_')[1])
		rmiss = range(ms,me)
		for spanindel in sharedspans.keys():
			ins, ine = int(spanindel.split('_')[0]), int(spanindel.split('_')[1])
			rindel = range(ins,ine)
			intersectlist = list(set(rmiss) & set(rindel))
			if len(intersectlist) > 0:
				overlapmiss[spanindel] = spanmiss[mspan] # record list of all taxa with missing data for this span
	for k in overlapmiss.keys():
		overlapmiss[k] = set(overlapmiss[k]) # just making sure that the list of taxa is unique.
		
	return(sharedspans, overlapmiss)
### NOT USED, but correct ###
def ReplaceGapsInFrontAndBackOfEachSequenceWithMissDataCharFromDict(seqs): ### NOT USED, but correct ###
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
	return(Xseqs) ### NOT USED, but correct ###

def ParseFasta(fle):
	seqs = {}
	FILE = open(fle, 'r')
	line = FILE.readline()
	while line:
		if line[0] == '>':
			n = line[1:].strip()
			seqs[n] = ''
		elif line.strip() == '#####': # for hypotheses file
			break
		else:
			seqs[n] += line.strip()
		line = FILE.readline()
	FILE.close()
	return(seqs)

def main():
	
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	
	filelist = []
	files = filter(os.path.isfile, os.listdir('.'))
	for filename in files:
		if extension.match(filename) != None:
			filelist.append(filename)
	
	hyps = ParseFasta(hypotheses)
	h = {} # will taxa that represent various hypotheses
	hreq = {}
	for key in hyps.keys():
		h[key] = ''.join(hyps[key].split('|')[0].split()).split(',') # Reads in hypotheses from .fasta(ish) file and has all taxa as list. i.e. >Sundatheria\nGVA, TCH, Tupaia\n
		hreq[key] = ''.join(hyps[key].split('|')[1].split()).split(',') # Reads in Taxa defined to be required for that hypothesis to be valid.
		
	oneindeldict = {} # counts all INDELs that support each hypothesis
	allgene = []
	allindelspergene = []
	for fle in sorted(filelist): # sort files by filename because of the number prefix to keep them in order [00001_*, 00002_*, ... 07118_*]
		gene = '_'.join(fle.split('_')[:2]) # makes gene name from filename		
		seqs = ParseFasta(fle)
		sharedspans, overlapmiss= RecordIndelsFromDict(seqs)
		oneindeldict, indelspergene = RecordINDELsThatSupportEachHypothesis(h, hreq, sharedspans, overlapmiss, oneindeldict, gene, seqs)
		allgene.append(gene)
		allindelspergene.append(indelspergene)
		if indelspergene != {}: # if gene has an indel that supports a hypothesis then copy it to directory outdir
			listofDicts = FormatGeneFastasToHighlightINDELs(seqs, indelspergene)
			count = 0
			for d in listofDicts:
				output = FormatDictionaryOfNucleotideSeqsToFasta(d)
				dest = os.path.join(os.getcwd(), outdir, '%s_%d.fa' % (fle[:-len(fileextension)], count))
				WriteOUT(dest, output)
				count += 1
			#src = os.path.join(os.getcwd(), fle)
			#shutil.copy(src, dest)
		
	print 'Number Alignments Queried: %d' % (len(filelist))
	
	output = FormatOutput(oneindeldict, 'NumberOfINDELsThatSupportEachHypothesis')
	WriteOUT('INDELsSupportingHypothesesCounts.v14.1.1.txt', output)
	output = FormatOutputLists(allindelspergene, allgene)
	WriteOUT('INDELsSupportingHypotheses_PerGene.v14.1.1.txt', output)
	
main()

stop = timeit.default_timer()

print stop-start