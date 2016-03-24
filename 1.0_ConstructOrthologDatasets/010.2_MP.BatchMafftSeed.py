'''Note: MAFFT v7.127b (2013/10/29) '''

import os
import re
import multiprocessing as mp

fileextension = '.OrthoUnalign.fa' # change this to specify which files are accessed.
extension = re.compile('.+%s$' % (fileextension))

mafft = 'mafft-linsi' # options: 'mafft', 'mafft-linsi'
msa = '--seed' # options: '', '--seed', '--add', '--addfragments'
direction = '' # options: '', '--adjustdirection'
order = '--inputorder' # options: --inputorder, or --reorder
thrd = 2
t = 8

def ExecuteCMD(cmd):
	os.system(cmd)

def RetrieveEnsemblAAByEnsemblGeneID():
	filelist = []
	files = filter(os.path.isfile, os.listdir('.'))
	for filename in files:
		if extension.match(filename) != None:
			filelist.append(filename)	
	cmdlist = []
	for f in filelist:
		s = f[:-len(fileextension)] + '.fasta'
		o = '%s.seedaligned.fa' % (f[:-len(fileextension)])

		cmd = '/usr/bin/mafft-linsi --thread %d --seed %s %s > %s' % (thrd, s, f, o)
		#cmd = '%s %s %s %s %s %s > %s' % (mafft, msa, direction, t, order, f, o)
		cmdlist.append(cmd)
	
	p = mp.Pool(t)
	p.map(ExecuteCMD, cmdlist)
	p.close()
	p.join()

RetrieveEnsemblAAByEnsemblGeneID()
