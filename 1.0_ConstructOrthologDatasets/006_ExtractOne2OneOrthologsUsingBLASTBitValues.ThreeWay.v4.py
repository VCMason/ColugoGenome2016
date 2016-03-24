''' All files must be BLAST against the same reference database where query taxa are used to determine orthology form X number of taxa. '''
''' If 5 files present, then all five taxa must have best BLAST hit to same DB sequence, and only once per taxa '''


import os
import re

fileextension = '.fa.out' # change this to specify which files are accessed.
qcovcutoff = 50.0 # Percent of query coverage cutoff required
bitcutoff = 0.99 # this value 1.0-0.01 times the bit value of the best BLAST hit must be greater than the bit value of the second best BLAST hit. bitcutoff*(BBH) > 2ndBBH

extension = re.compile('.+%s$' % (fileextension))

def GetINFOFromBLAST(f, count, one2one, countone2one):
	bits, qcoverage, querygenes, numberhits = {}, {}, {}, {}

	BLAST = open(f, 'r')
	line = BLAST.readline()
	while line:
		querygeneid, dbgeneid, bitval, qcov = line.strip().split('\t')[0], line.strip().split('\t')[1], float(line.strip().split('\t')[11]), float(line.strip().split('\t')[13])
		querygeneid = '_'.join(querygeneid.split('_')[:-1]) # remove this if all query sequences from different taxa have the same name. Modify if they have part of the sequence different.
		numberhits[querygeneid] = numberhits.get(querygeneid, 0) + 1 # counts total number of hits to each Homogene
		if numberhits[querygeneid] <= 2:
			bits[querygeneid] = bits.get(querygeneid, []) + [bitval] # records only the best 2 BLAST hit's bitvalues
			querygenes[querygeneid] = querygenes.get(querygeneid, []) + [dbgeneid] # key is query name and value is list of top two BLAST database names.
			qcoverage[querygeneid] = qcoverage.get(querygeneid, []) + [qcov] # key is query name and value is list of top two BLAST percent query coverages.
		line = BLAST.readline()
	BLAST.close()
	
	red = {} # key is query gene name and val is database name. Reduced set of BBH that have 75% of best blast hit bitvalue is > than second best blast hit and query coverage is greater than qcovcutoff.
	uniquequery = {} # 1 if query genes only had best BLAST hit to one database gene
	for query in bits.keys(): # key of bits is query name
		if len(bits[query]) == 2:
			if bits[query][0]*bitcutoff > bits[query][1] and qcoverage[query][0] > qcovcutoff: # if 75% of best blast hit bitvalue is > than second best blast hit => one2one ortholog 
				red[query] = querygenes[query][0]
		elif len(bits[query]) == 1 and qcoverage[query][0] > qcovcutoff:
			red[query] = querygenes[query][0]
		uniquequery[querygenes[query][0]] = uniquequery.get(querygenes[query][0], 0) + 1 #count number of occurances of each database sequence as Best BLAST target for query sequences. If = 1 then only one BBH to this database sequence in file.

	for query in red.keys(): # checks if one2one vs one2many. Looks if a database sequence has one or more than count Best BLAST hits.
		if uniquequery[querygenes[query][0]] == 1:
			one2one[query] = querygenes[query][0]
			countone2one[query] = countone2one.get(query, 0) + 1
		else:
			pass
	
	return(one2one, countone2one)

def One2OneOrthologByBitValue():

	filelist = []
	files = filter(os.path.isfile, os.listdir('.'))
	for filename in files:
		if extension.match(filename) != None:
			filelist.append(filename)
	one2one = {}
	countone2one = {}
	count = 1
	for f in filelist:
		one2one, countone2one = GetINFOFromBLAST(f, count, one2one, countone2one)
		count += 1
	
	realone2one = {}
	for query in countone2one.keys():
		#print countone2one[query]
		if countone2one[query] == len(filelist):
			realone2one[query] = one2one[query]
	print 'Number of One2One Orthologs identified in at least one of input BLAST files: %d' % (len(countone2one))
	print 'Number of One2One Orthologs in all input BLAST files: %d' % (len(realone2one.keys()))
    
	output = ''
	for k in realone2one.keys():
		output += '%s\t%s\n' % (k, realone2one[k])
    
	OUT = open('TCH_HMGMC_BLASTBitFilter.ThreeWayHMD_%d_qcov%.0f_bit%.2f.v4.tab' % (len(realone2one.keys()), qcovcutoff, bitcutoff), 'w')
	OUT.write(output)
	OUT.close()
    
    
One2OneOrthologByBitValue()
