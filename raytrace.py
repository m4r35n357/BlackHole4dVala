#!/usr/bin/env python

from sys import argv, stdout
from math import fabs, log10, sqrt, sin, cos, pi
from json import loads

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
			print >> stdout, str(p['tau']) + ' 2 ' + str(p['r']) + ' ' + str(cos(p['th'])) + ' ' + str(p['t']) + ' ' + str(p['ph']) + ' ' + str(sqrt(p['R'])) + ' ' + str(- sin(p['th']) * sqrt(p['THETA'])) + ' ' + str(p['tDot']) + ' ' + str(p['phDot']) + ' ' + str(-1) + ' ' + str(0) + ' ' + str(0) + ' ' + str(0) + ' ' + str(0) + ' ' + str(0) + ' ' + str(0) + ' ' + str(1) + ' ' + str(-1) + ' ' + str(1) + ' ' + str(0) + ' ' + str(0)
		line = dataFile.readline()
		n += 1

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
