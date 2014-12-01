#!/usr/bin/env python

from sys import argv, stdin, stdout
from math import fabs, log10, sqrt, sin, cos, pi
from json import loads

def main():
	if len(argv) < 1:
		raise Exception('>>> ERROR! Please supply a plotting poing interval <<<')
	interval = int(argv[1])
	line = stdin.readline()
	n = 0
	while line:
		p = loads(line)
		if (n % interval == 0):
			print >> stdout, str(p['tau']) + ' 2 ' + str(p['r']) + ' ' + str(cos(p['th'])) + ' ' + str(p['t']) + ' ' + str(p['ph']) + ' ' + str(sqrt(p['R'])) + ' ' + str(- sin(p['th']) * sqrt(p['THETA'])) + ' ' + str(p['tDot']) + ' ' + str(p['phDot']) + ' ' + str(-1) + ' ' + str(0) + ' ' + str(0) + ' ' + str(0) + ' ' + str(0) + ' ' + str(0) + ' ' + str(0) + ' ' + str(1) + ' ' + str(0) + ' ' + str(-1) + ' ' + str(0) + ' ' + str(0)
		line = stdin.readline()
		n += 1

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
