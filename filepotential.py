#!/usr/bin/env python

from sys import argv
from matplotlib import pyplot
from json import loads

def main():
	if len(argv) < 2:
		raise Exception('>>> ERROR! Please supply a potential file name <<<')
	dataFile = open(argv[1], 'r')
	line = dataFile.readline()
	ax1 = pyplot.figure().add_subplot(111)
        pyplot.grid(b=True, which='major', color='k', linestyle='-')
	ax1.set_xlabel('r, theta')
	ax1.set_ylabel('r, R(r)', color='b')
        ax1.set_ylim(-50, 30)
	ax2 = ax1.twinx()
	ax2.set_ylabel('theta, THETA(theta)', color='r')
        ax2.set_ylim(-50, 30)
	n = 0
	while line:
		p = loads(line)
		ax1.plot(p['r'], p['R'], 'b.', markersize=2)
		ax2.plot(p['theta'], p['THETA'], 'r.', markersize=2)
		line = dataFile.readline()
		n += 1
	pyplot.show()

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
