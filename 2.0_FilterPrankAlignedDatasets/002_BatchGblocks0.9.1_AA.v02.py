import os
import re
import glob
import multiprocessing as mp

#b = 1000 # num bootstraps
t = 5 # num threads

def ExecuteCMD(cmd):
	os.system(cmd)

def main():
	
	cwd = os.getcwd()
	os.system('mkdir %s/12TaxaAA_PrankAligned/GblocksFiltered_AA_Alignments' % (cwd))
	path = cwd + '/12TaxaAA_PrankAligned/*.AA.Align.red.prank.fas'
	filelist = glob.glob(path)
 
	cmdlist = []
	for fle in filelist:
		#stump = f.split('/')[-1][:-len('_AA_aligned.fasta.1.fas')]
		#o = stump + '_AA_aligned.GB' # extension .fas is added by prank
		command = '/home/vmason/bin/Gblocks %s -t=p -p=t' % (fle) # -t=p (protein) -t=d (DNA)
		cmdlist.append(command)

	p = mp.Pool(t)
	p.map(ExecuteCMD, cmdlist)
	p.close()
	p.join()


main()
