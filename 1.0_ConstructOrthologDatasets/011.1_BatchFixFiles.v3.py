import os
import re

fileextension = '.seedaligned.fa' # change this to specify which files are accessed.
extension = re.compile('.+%s$' % (fileextension))

def main():
	filelist = []
	files = filter(os.path.isfile, os.listdir('.'))
	for filename in files:
		if extension.match(filename) != None:
			filelist.append(filename)
	for f in filelist:
		FILE = open(f, 'r')
		line = FILE.readline()
		output = ''
		while line:
			if line[0] == '>' and '|' in line:
				s = line.strip().split('|')
				output += '>%s\n' % (s[1])
			elif line[0] == '>' and line[:13] == '>_seed__seed_':
				s = line.strip().split('_')
				output += '>%s\n' % (s[4])
			elif line[0] == '>' and line[:7] == '>_seed_':
				s = line.strip().split('_')
				output += '>%s\n' % (s[2])
			elif line[0] == '>':
				n = line[1:].strip()
				output += '>%s\n' % (n)
			else:
				output += line	
			line = FILE.readline()
		FILE.close()
		OUT = open('%s.fixname' % (f), 'w')
		OUT.write(output)
		OUT.close()

main()
