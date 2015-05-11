#!/usr/bin/env python

from sys import argv, stdin, stdout
from matplotlib import pyplot
from json import loads
from array import array
from scipy.interpolate import interp1d
import numpy as np

def main():
	if len(argv) < 3:
		raise Exception('>>> ERROR! Please supply a data file name, the number of points to plot, and a coordinate to plot <<<')
	dataFile = open(argv[1], 'r')
	line = dataFile.readline()
	nPoints = int(argv[2])
	coordinate = argv[3]
	tau = array('d')
	c = array('d')
	cDot = array('d')
	tauMax = 0.0
	while line:
		p = loads(line)
		tauValue = p['tau']
		tauMax = tauValue if tauValue > tauMax else tauMax
		tau.append(tauValue)
		c.append(p[coordinate])
		cDot.append(p[coordinate + 'Dot'])
		line = dataFile.readline()
	# interpolate here
	cI = interp1d(tau, c)
	cDotI = interp1d(tau, cDot)
	ax1 = pyplot.figure().add_subplot(111)
	ax1.set_xlabel('tau', color='k')
	ax1.set_ylabel(coordinate, color='b')
	ax2 = ax1.twinx()
	ax2.set_ylabel(coordinate + 'Dot', color='r')
	tauI = np.linspace(0, tauMax, num = nPoints)
	for i in range(len(tauI)):
		ax1.plot(tauI[i], cI(tauI[i]), 'b.', markersize=2)
		ax2.plot(tauI[i], cDotI(tauI[i]), 'r.', markersize=2)
	pyplot.show()

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"

