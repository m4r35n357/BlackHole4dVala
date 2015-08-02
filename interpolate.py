#!/usr/bin/env python

from sys import argv, stdin, stdout
from matplotlib import pyplot
from json import loads
from array import array
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy as np

def main():
	if len(argv) < 3:
		raise Exception('>>> ERROR! Please supply a data file name, a time variable, the number of points to plot, and a coordinate to plot <<<')
	dataFile = open(argv[1], 'r')
	timeCoordinate = str(argv[2])
	nPoints = int(argv[3])
	coordinate = argv[4]
	tau = array('d')
	c = array('d')
	cDot = array('d')
	tauMax = 0.0
	line = dataFile.readline()
	while line:  # build raw data arrays
		p = loads(line)
		tauValue = p[timeCoordinate]
		tauMax = tauValue if tauValue > tauMax else tauMax
		tau.append(tauValue)
		c.append(p[coordinate])
		cDot.append(p[coordinate + 'Dot'])
		line = dataFile.readline()
        try:  # interpolate here
        	cI = InterpolatedUnivariateSpline(tau, c, k = 1)
		cDotI = InterpolatedUnivariateSpline(tau, cDot, k = 1)
        except ValueError as e:
            print('DATA ERROR: ' + str(argv[0]) + ': ' + str(e))
            exit(-2)		
	ax1 = pyplot.figure().add_subplot(111)
        pyplot.grid(b=True, which='major', color='0.25', linestyle='-')
	ax1.set_xlabel('Time: ' + timeCoordinate, color='0.20')
	ax1.set_ylabel(coordinate, color='b')
	ax2 = ax1.twinx()
	ax2.set_ylabel(coordinate + 'Dot', color='r')
	tauI = np.linspace(0, tauMax, num = nPoints)
	for i in range(len(tauI)):
		ax1.plot(tauI[i], cI(tauI[i]), 'b.', markersize=2)
		ax2.plot(tauI[i], cDotI(tauI[i]), 'r.', markersize=2)
		#ax1.plot(cI(tauI[i]), cDotI(tauI[i]), 'm.', markersize=2)
        try:
            pyplot.show()
        except AttributeError as e:
            print('ATTRIBUTE ERROR: ' + str(argv[0]) + ':' + str(coordinate) + ': ' + str(e))
            exit(-3)		

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"

