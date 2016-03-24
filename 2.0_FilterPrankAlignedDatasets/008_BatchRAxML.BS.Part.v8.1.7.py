''' Author: Victor C Mason '''
''' Note: to exclude columns from an alignment execute a command with the following structure: '''
''' raxml -s infile.trim.phy -x 1234 -p 123 -N 100 -f a -k -m GTRCAT -T 32 -n outfile.phy -E in.excludefile '''


import os
import re

fileextension = '.phy' # change this to specify which files are accessed.

extension = re.compile('.+%s$' % (fileextension))

N = 1000 # num bootstraps
t = 32 # num threads
part = 'Concatenated3926OrthoGenes.NU.part' # partition filename

def main():
	filelist = []
	files = filter(os.path.isfile, os.listdir('.'))
	for filename in files:
		if extension.match(filename) != None:
			filelist.append(filename)
 
	#os.system('mkdir GblocksRAxMLNUTrees')
	for fle in filelist:
		out = '%s.N%d.tre' % (fle[:-len(fileextension)], N)
		
		#command = 'raxmlHPC-PTHREADS-SSE3 -f a -m GTRGAMMA -x 1234 -p 12345 -N %d -T %d -s %s -n %s' % (N, t, fle, out) # Finds max-likelihood tree and does boostrap
		command = 'raxml -f d -m GTRGAMMA -p 12345 -T %d -s %s -n %s' % (t, fle, out) # -f d is default rapid-hill climb alg. Finds max-likelihood tree. # raxmlHPC-PTHREADS-SSE3
		#out = '%sPart.tre' % (fle[:-len(fileextension)])
		#command = 'raxmlHPC-PTHREADS-SSE3 -f a -m GTRGAMMA -x 1234 -p 12345 -N %d -T %d -s %s -q %s -n %s' % (N, t, fle, part, out) # Finds max-likelihood tree and does boostrap
		#command = 'raxml -f d -m GTRGAMMA -p 12345 -T %d -s %s -q %s -n %s' % (t, fle, part, out) # -f d is default rapid-hill climb alg. Finds max-likelihood tree. With Partition File # raxmlHPC-PTHREADS-SSE3
		os.system(command)


main()
