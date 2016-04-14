import os
import re
import shutil

fileextension = '.fq' # change this to specify which files are accessed.
outextension = '.DAMtrim.fq'

extension = re.compile('.+%s$' % (fileextension))


def main():
    filelist = []
    files = filter(os.path.isfile, os.listdir('.'))
    for filename in files:
        if extension.match(filename) != None:
            filelist.append(filename)
 
    uniquefle = {}
    for fle in filelist:
	OUT = open('%s%s' % (fle[:-len(fileextension)], outextension), 'a')
	FILE = open(fle, 'r')
	line = FILE.readline()
	count = 0
	while line:
	    modulus = count%4
	    if modulus == 0:
		OUT.write(line)
	    elif modulus == 1:
		line = line.strip()
		output = '%s\n' % (line[3:-3])
		OUT.write(output)
	    elif modulus == 2:
		OUT.write(line)
	    elif modulus == 3:
		line = line.strip()
		output = '%s\n' % (line[3:-3])
		OUT.write(output)
	    else:
		pass
	    line = FILE.readline()
	    count += 1
	FILE.close()
	OUT.close()

	
main()
