#!/usr/bin/env python

from sys import argv, stdout
from math import fabs, log10, sqrt, sin, cos, pi
from json import loads
from subprocess import Popen, PIPE

def main():
	if len(argv) < 2:
		raise Exception('>>> ERROR! Please supply a data file name <<<')
	dataFile = open(argv[1], 'r')
	interval = int(argv[2])
	line = dataFile.readline()
	n = 0
	while line:
		p = loads(line)
		if (n % interval == 0):
			print >> stdout, line
			command = str(p['tau']) + ' 2 ' + str(p['r']) + ' ' + str(cos(p['th'])) + ' ' + str(p['t']) + ' ' + str(p['ph']) + ' ' + str(sqrt(p['R'])) + ' ' + str(- sin(p['th']) * sqrt(p['THETA'])) + ' ' + str(p['tDot']) + ' ' + str(p['phDot']) + ' ' + str(-1) + ' ' + str(0) + ' ' + str(0) + ' ' + str(0) + ' ' + str(0) + ' ' + str(0) + ' ' + str(0) + ' ' + str(1) + ' ' + str(-1) + ' ' + str(1) + ' ' + str(0) + ' ' + str(0)
			print >> stdout, command
			print >> stdout, ''
			f = open('out.' + str(n) + '.ppm','w')
			Popen(["/home/ian/projects/c/kerr-image/kerr-image"], stdin=PIPE, stdout=f).communicate(input=command)[0]
			f.close()
		line = dataFile.readline()
		n += 1

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
