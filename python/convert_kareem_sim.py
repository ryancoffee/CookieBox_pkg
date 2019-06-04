#!/reg/g/psdm/sw/conda/inst/miniconda2-prod-rhel7/envs/ana-1.3.59-py3/bin/python3

from numpy import savetxt,loadtxt,ndarray
import sys
import re

def main():
	if len(sys.argv)<2:
		print('Massive fail... give a filename as the arguement')
		print('syntax = execname filename [offset scale [outlen] ]')
		return '' ;
	filename = sys.argv[1]
	print(filename)
	outlen = 256
	offset = 330
	scale = 0.35
	if len(sys.argv)>3:
		offset = sys.argv[2]
		scale = sys.argv[3]
		if len(sys.argv)>4:
			outlen = sys.argv[4]
	data = loadtxt(filename).T
	outdata = scale * data[:,offset:offset+outlen]
	outname = re.sub(r'\.\w+$','.out',filename)
	savetxt(outname,outdata,fmt='%.6e')
	return outname;

if __name__ == '__main__':
	main();
