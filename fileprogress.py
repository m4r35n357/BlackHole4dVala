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
	ax1.set_xlabel('Proper Time, tau')
	ax1.set_ylabel('Error, dB', color='b')
#	ax2 = ax1.twinx()
#	ax2.set_ylabel('Cumulative Error', color='b')
	n = 0
	while line:
		p = loads(line)
		if (n % interval == 0):
			ax1.plot(p['tau'], p['ETh'], 'r.', markersize=1)
			ax1.plot(p['tau'], p['ER'], 'g.', markersize=1)
			ax1.plot(p['tau'], p['E'], 'b.', markersize=2)
			ax1.plot(p['tau'], p['EC'], 'k.', markersize=2)
		line = dataFile.readline()
		n += 1
	pyplot.show()

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
