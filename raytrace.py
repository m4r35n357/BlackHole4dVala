#!/usr/bin/env python

from sys import argv, stdin, stdout, stderr
from math import sqrt, sin, cos
from json import loads

def main():
	if len(argv) < 1:
		raise Exception('>>> ERROR! Please supply a plotting interval <<<')
	interval = int(argv[1])
	line = stdin.readline()
	n = 0
	while line:
		p = loads(line)
		if (n % interval == 0):
			print >> stdout, str(p['tau']) + ' 2 ' \
				+ str(p['r']) + ' ' + str(cos(p['th'])) + ' ' + str(p['t']) + ' ' + str(p['ph']) + ' ' \
				+ str(p['rDot']) + ' ' + str(- sin(p['th']) * p['thDot']) + ' ' + str(p['tDot']) + ' ' + str(p['phDot']) + ' ' \
				+ str(-1) + ' ' + str(0) + ' ' + str(0) + ' ' + str(0) + ' ' \
				+ str(0) + ' ' + str(0) + ' ' + str(0) + ' ' + str(1) + ' ' \
				+ str(0) + ' ' + str(1) + ' ' + str(0) + ' ' + str(0)
		line = stdin.readline()
		n += 1

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
