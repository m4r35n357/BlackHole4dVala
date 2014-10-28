#!/usr/bin/env python

from sys import argv
from matplotlib import pyplot
from json import loads

def main():
	if len(argv) < 3:
		raise Exception('>>> ERROR! Please supply a data file name and a plotting interval <<<')
	dataFile = open(argv[1], 'r')
	interval = int(argv[2])
	line = dataFile.readline()
	ax1 = pyplot.figure().add_subplot(111)
	ax1.set_xlabel('tau')
	ax1.set_ylabel('r, R(r)', color='b')
	ax2 = ax1.twinx()
	ax2.set_ylabel('theta, THETA(theta)', color='r')
	n = 0
	while line:
		p = loads(line)
		if (n % interval == 0):
			ax1.plot(p['tau'], p['r'], 'b.', markersize=2)
			ax2.plot(p['tau'], p['th'], 'r.', markersize=2)
			ax1.plot(p['tau'], p['R'], 'g.', markersize=1)
			ax2.plot(p['tau'], p['THETA'], 'm.', markersize=1)
		line = dataFile.readline()
		n += 1
	pyplot.show()

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
