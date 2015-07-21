#!/usr/bin/env python

from sys import argv
from math import fabs, log10
from matplotlib import pyplot
from json import loads
from array import array
from scipy.interpolate import interp1d
import numpy as np

def main():
	if len(argv) < 3:
		raise Exception('>>> ERROR! Please supply a data file name, a time variable, and the number of points to plot <<<')
	dataFile = open(argv[1], 'r')
	timeCoordinate = str(argv[2])
	nPoints = int(argv[3])
	line = dataFile.readline()
	tau = array('d')
	eR = array('d')
	eTh = array('d')
	ev4 = array('d')
	tauMax = 0.0
	while line:  # build raw data arrays
		p = loads(line)
		tauValue = p[timeCoordinate]
		tauMax = tauValue if tauValue > tauMax else tauMax
		tau.append(tauValue)
		eR.append(p['ER'])
		eTh.append(p['ETh'])
		ev4.append(p['v4e'])
		line = dataFile.readline()
        try:  # interpolate here
        	eRI = interp1d(tau, eR, kind='linear', copy=False)
		eThI = interp1d(tau, eTh, kind='linear', copy=False)
		ev4I = interp1d(tau, ev4, kind='linear', copy=False)
        except ValueError as e:
            print('DATA ERROR: ' + str(argv[0]) + ': ' + str(e))
            exit(-2)		
	ax1 = pyplot.figure().add_subplot(111)
        pyplot.grid(b=True, which='major', color='0.25', linestyle='-')
	ax1.set_xlabel('Time: ' + timeCoordinate)
	ax1.set_ylabel('Radial and Latitudinal Errors, dB', color='k')
        ax1.set_ylim(-180.0, 0.0)
	ax2 = ax1.twinx()
	ax2.set_ylabel('4-Velocity Norm Error', color='g')
        ax2.set_ylim(-180.0, 0.0)
	tauI = np.linspace(0, tauMax, num = nPoints)
	for i in range(len(tauI)):
		ax1.plot(tauI[i], eRI(tauI[i]), 'b.', markersize=2)
		ax1.plot(tauI[i], eThI(tauI[i]), 'r.', markersize=2)
		ax2.plot(tauI[i], ev4I(tauI[i]), 'g.', markersize=2)
        try:
            pyplot.show()
        except AttributeError as e:
            print('ATTRIBUTE ERROR: ' + str(argv[0]) + ':' + str(coordinate) + ': ' + str(e))
            exit(-3)		

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
