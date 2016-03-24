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

start = timeit.default_timer()

indir = 'GblocksNTFastaRemFramAndFilteredAlignments_MedianSig2Min0.40Win10' # Sequence Alignments
NoOutlierDir = '21Species_NoOutlierTrees_z3.75_c2.50' # input no outlier trees directory
outdir = '21TaxaNT_BranchLengthFilteredAlignments'
fileextension = '.fa' # change this to specify which sequence alignments are accessed.
extension = re.compile('.+%s$' % (fileextension))
fileextension2 = '.tre' # change this to specify which tree files are accessed.
extension2 = re.compile('.+%s$' % (fileextension))

def WriteOUT(outfile, output):
	OUT = open(outfile, 'w')
	OUT.write(output)
	OUT.close()

def main():
	
	cwd = os.getcwd()
	path = '%s/%s/*%s' % (cwd, NoOutlierDir, fileextension2)
	filelist = glob.glob(path)
	
	if not os.path.exists(outdir):
		os.mkdir(outdir)			
	
	for fle in filelist: # list of tree file paths
		f = '.'.join(os.path.split(fle)[1].split('.')[1:-2]) + fileextension # converts tree filename to align filename
		src = os.path.join(os.getcwd(), indir, f)
		dst = os.path.join(os.getcwd(), outdir, f)
		shutil.copy(src, dst)	

main()
		
stop = timeit.default_timer()
		
print 'Program execution Length: %f' % (stop-start)
