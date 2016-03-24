''' Author: Victor C Mason '''
''' Date: September 07, 2015 '''
''' Input: aligned fasta '''
''' Output: '''
''' Summary: Identifies outlier branch lengths in newick trees'''

import os
import re
import glob
import shutil
import timeit
import math
from ete2 import Tree
import itertools
import operator

start = timeit.default_timer()

zscore = 3.75 # specify number of standard deviations away from the median. # set as 0 if do not want to exclude individuals based on the number of standard deviations away from the mean.
cutofftaxa = ['Homo', 'Mus'] #  # specify two reliable taxa that span basal nodes in the phylogeny (i.e. human and the outgroup taxa) that will represent a minimum branch length cutoff (multiplied by two) # if no pairs of taxa are reliable just leave empty i.e. [] and the program will only exclude outliers only when branch lenth > median+(Zscore*stddev)
multiplier = 2.50 # value to multiply the branch length between cutoff taxa. # if multiplier = 3.0 then a branch length is an outlier if it is >= 3.0*BLofCutofftaxa
AndOr = 'or' # specify as 'and' or 'or'. # determines if you want to determine outliers by zscore AND cutofftaxaBL or by zscore OR cutofftaxaBL

indir = '21TaxaNT_Trees'
NoOutlierDir = '21Species_NoOutlierTrees_z%.2f_c%.2f' % (zscore, multiplier)
OutlierDir = '21Species_OutlierTrees_z%.2f_c%.2f' % (zscore, multiplier)
fileextension = '.tre' # change this to specify which input files are accessed.
extension = re.compile('.+%s$' % (fileextension))

def WriteOUT(outfile, output):
	OUT = open(outfile, 'w')
	OUT.write(output)
	OUT.close()

def RecordValueIfGreaterThanNumberOfStandardDevsAwayFromMedian(dist, zscore, cutoff):
	outliers = {}
	median = Median(dist.values())
	stddev = StandardDeviationUsingMedian(dist.values()) # average of squared differences of values from median of values
	countmincut = 0.0
	countsdev = 0.0
	countboth = 0.0
	for k in dist.keys():
		if AndOr == 'or':
			if dist[k] > (median+(zscore*stddev)) or dist[k] >= multiplier*cutoff: # value must be significantly different according to Zscore and also greater than pairwise distance cutoff.
				outliers[k] = dist[k] # records the pair of taxa as key and distance that is identified as an outlier as value
		else:
			if dist[k] > (median+(zscore*stddev)) and dist[k] >= multiplier*cutoff: # value must be significantly different according to Zscore and also greater than pairwise distance cutoff.
				outliers[k] = dist[k] # records the pair of taxa as key and distance that is identified as an outlier as value
		
		if dist[k] > (median+(zscore*stddev)) and dist[k] >= multiplier*cutoff:
			countboth += 1.0
		elif dist[k] >= multiplier*cutoff:
			countmincut += 1.0
		elif dist[k] > (median+(zscore*stddev)):
			countsdev += 1.0
	return(outliers, median, stddev, countboth, countmincut, countsdev)

def StandardDeviationUsingMedian(mylist):
	mn = Median(mylist)
	variance = sum([(e-mn)**2 for e in mylist]) / (len(mylist)-1) # average of the squared differences of values from mylist from their median value. # minus 1 in denominator because window is sample of larger population (Bessel's correction)
	return(math.sqrt(variance)) # stddev = sigma = sqrt of variance

def StandardDeviation(mylist):
	mn = Mean(mylist)
	variance = sum([(e-mn)**2 for e in mylist]) / (len(mylist)-1) # average of the squared differences of values from mylist from their mean value. # minus 1 in denominator because window is sample of larger population (Bessel's correction)
	return(math.sqrt(variance)) # stddev = sigma = sqrt of variance

def Median(mylist):
	sorts = sorted(mylist)
	length = len(sorts)
	if not length % 2: # if not odd (if even)
		return ((sorts[length / 2] + sorts[length / 2 - 1]) / 2.0)
	return sorts[length / 2] # if not even then it is odd

def Mean(mylist):
	return((sum(mylist)*1.0)/len(mylist))

def GetDistanceBetweenAllLeafNodes(allpairs, t):
	dist = {} # pairs of taxa are key, distance between them is value
	alldistperleaf = {}
	for pair in allpairs:
		d = t.get_distance(pair[0],pair[1])
		alldistperleaf[pair[0]] = alldistperleaf.get(pair[0], []) + [d]
		alldistperleaf[pair[1]] = alldistperleaf.get(pair[1], []) + [d]
		dist['_'.join(pair)] = d
	return(dist, alldistperleaf)
	
def ReadNewick(f):
	FILE = open(f, 'r')
	nw = FILE.readline() # assume one tree per file
	FILE.close()
	return(nw)

def main():
	
	cwd = os.getcwd()
	path = '%s/%s/*%s' % (cwd, indir, fileextension)
	filelist = glob.glob(path)	
			
	summary = ''
	log = ''
	countoffenders = {}
	totalboth = 0.0
	totalmincut = 0.0
	totalsdev = 0.0
	gate = 1
	nooutliers = []
	alloutliers = []
	for fle in sorted(filelist): # sort files by filename because of the number prefix to keep them in order [00001_*, 00002_*, ... 07118_*]		
		nw = ReadNewick(fle)
		t = Tree(nw)
		if gate == 1: # only calculate this block once, assumes same taxa in all trees.
			s = t.get_leaf_names() # list of species in tree
			allpairs = [list(map(str,comb)) for comb in itertools.combinations(s, 2)] # all unique combinations of elements of specified length # makes list of lists
		dist, alldistperleaf = GetDistanceBetweenAllLeafNodes(allpairs, t)
		if len(cutofftaxa) == 2:
			cutoff = t.get_distance(cutofftaxa[0], cutofftaxa[1])
		else:
			cutoff = 0 # if no cutofftaxa were identified
			print 'No valid cutoff taxa were identified, so cutoff = 0'
		outliers, median, stddev, countboth, countmincut, countsdev = RecordValueIfGreaterThanNumberOfStandardDevsAwayFromMedian(dist, zscore, cutoff)
		totalboth += countboth
		totalmincut += countmincut
		totalsdev += countsdev
		if len(outliers) >= 1: # if there are outliers
			alloutliers.append(os.path.split(fle)[1]) # required
			summary += '>%s_MEDIAN:%f_SD:%f_MINCUT:%fMEAN:%f\n' % (os.path.split(fle)[1], median, stddev, multiplier*cutoff, Mean(dist.values())) # below is all output formating in this block
			for k in outliers:
				summary += '%s\t%s\n' % (k, outliers[k])
			summary += '@\n'
			mediandistperleaf = {}
			for n in t.get_leaf_names():
				mediandistperleaf[n] = Median(alldistperleaf[n])
				summary += '%s\tMEDIAN:%f\n' % (n, Median(alldistperleaf[n]))
			log += '>%s_MEDIAN:%f_SD:%f_MINCUT:%fMEAN:%f\n' % (os.path.split(fle)[1], median, stddev, multiplier*cutoff, Mean(dist.values()))
			indivoutliers, md, mn, countboth, countmincut, countsdev = RecordValueIfGreaterThanNumberOfStandardDevsAwayFromMedian(mediandistperleaf, zscore, cutoff)
			l1,l2 = [], []
			for k in indivoutliers.keys():
				l1.append(k)
				l2.append(indivoutliers[k])
				countoffenders[k] = countoffenders.get(k, 0) + 1
			log += '%s\n%s\n' % ('\t'.join(l1), '\t'.join([str(f) for f in l2]))
		else:
			nooutliers.append(os.path.split(fle)[1]) # required
			#print '>%s_MEDIAN:%f_SD:%f_MINCUT:%f\n' % (os.path.split(fle)[1], median, stddev, cutoff)
		gate = 0

	outfile = 'DetailsOfOutliers.%dTaxa_z%.2f_c%.2f.log.fq' % (len(s), zscore, multiplier) # s is list of taxa
	outpath = os.path.join(os.getcwd(), outfile)
	WriteOUT(outpath, summary) # yes summary is the detailed output
	outfile = 'SummaryOfOutliers.%dTaxa_z%.2f_c%.2f.log.fa' % (len(s), zscore, multiplier) # s is list of taxa
	outpath = os.path.join(os.getcwd(), outfile)
	WriteOUT(outpath, log)
	
	if not os.path.exists(NoOutlierDir):
		os.mkdir(NoOutlierDir)
	for f in nooutliers:
		src = os.path.join(os.getcwd(), indir, f)
		dst = os.path.join(os.getcwd(), NoOutlierDir, f)
		shutil.copy(src, dst)
	if not os.path.exists(OutlierDir):
		os.mkdir(OutlierDir)
	for f in alloutliers:
		src = os.path.join(os.getcwd(), indir, f)
		dst = os.path.join(os.getcwd(), OutlierDir, f)
		shutil.copy(src, dst)
	out = ''
	out += 'Total number of branches identified as outliers because >= %f * branch length of %s AND >= %f standard deviations away from median: %f\n' % (multiplier, '_'.join(cutofftaxa), zscore, totalboth)
	out += 'Total number of branches identified as outliers because >= %f * branch length of %s: %f\n' % (multiplier, '_'.join(cutofftaxa), totalmincut)
	out += 'Total number of branches identified as outliers because >= %f standard deviations away from median: %f\n' % (zscore, totalsdev)
	out += 'Total number of branches identified as outliers because >= %f * branch length of %s AND >= %f standard deviations away from median: %f\n' % (multiplier, '_'.join(cutofftaxa), zscore, (totalboth/(totalboth+totalsdev+totalmincut)*100.0))	
	out += 'Percent of Branches excluded because >= %f * branch length of %s: %f\n' % (multiplier, '_'.join(cutofftaxa), (totalmincut/(totalboth+totalsdev+totalmincut)*100.0))
	out += 'Percent of Branches excluded because >= %f standard deviations away from median: %f\n' % (zscore, (totalsdev/(totalboth+totalsdev+totalmincut)*100.0))
	out += 'Total outlier trees: %d\n' % (len(alloutliers))
	out += 'Total Non-outlier trees: %d\n' % (len(nooutliers))
	for tup in sorted(countoffenders.items(), key=operator.itemgetter(1), reverse=True):
		print 'Taxa:%s was identified as an outlier:%d times\n' % (tup[0], tup[1])
		out += 'Taxa:%s was identified as an outlier:%d times\n' % (tup[0], tup[1])
	outfile = 'NumberOfOutliers.%dTaxa_z%.2f_c%.2f.log.txt' % (len(s), zscore, multiplier) # s is list of taxa
	outpath = os.path.join(os.getcwd(), outfile)
	WriteOUT(outpath, out)	

main()
		
stop = timeit.default_timer()
		
print 'Program execution Length: %f' % (stop-start)
