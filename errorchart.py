#!/usr/bin/env python

from sys import argv
from math import fabs, log10
from matplotlib import pyplot
from json import loads

def main():
	if len(argv) < 3:
		raise Exception('>>> ERROR! Please supply a data file name and a plotting interval <<<')
	dataFile = open(argv[1], 'r')
	interval = int(argv[2])
	line = dataFile.readline()
	ax1 = pyplot.figure().add_subplot(111)
	ax1.set_xlabel('Mino Time, lambda')
	ax1.set_ylabel('Radial and Latitudinal Error, dB', color='b')
        ax1.set_ylim(-180.0, 0.0)
	ax2 = ax1.twinx()
	ax2.set_ylabel('4-Velocity Norm Error', color='g')
        delta = 1.0e-6
#        ax2.set_ylim(1.0 - delta, 1.0 + delta)
        ax2.set_ylim(-180.0, 0.0)
	n = 0
	while line:
		p = loads(line)
		if (n % interval == 0):
			ax1.plot(p['mino'], p['ETh'], 'r.', markersize=1)
			ax1.plot(p['mino'], p['ER'], 'b.', markersize=1)
			ax1.plot(p['mino'], p['E'], 'k.', markersize=2)
			ax2.plot(p['mino'], p['v4e'], 'g.', markersize=2)
		line = dataFile.readline()
		n += 1
#        pyplot.legend(['E', 'Er', 'Eth', 'EC'], loc='best')
        try:
            pyplot.show()
        except AttributeError as e:
            print('ATTRIBUTE ERROR: ' + str(argv[0]) + ': ' + str(e))
            exit(-1)		

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
