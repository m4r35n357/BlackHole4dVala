#!/usr/bin/env python

from sys import argv, stdin, stdout
from matplotlib import pyplot
from json import loads
from array import array
from scipy.interpolate import interp1d
import numpy as np

def main():
	if len(argv) < 2:
		raise Exception('>>> ERROR! Please supply a data file name, a plotting interval, and a coordinate to plot <<<')
	dataFile = stdin
	line = dataFile.readline()
	nData = int(argv[1])
	tau = array('d')
	x = array('d')
	y = array('d')
	z = array('d')
	e = array('d')
	tauMax = 0.0
	while line:
		p = loads(line)
		tauValue = p['tau']
		tauMax = tauValue if tauValue > tauMax else tauMax
		tau.append(tauValue)
		x.append(p['x'])
		y.append(p['y'])
		z.append(p['z'])
		e.append(p['E'])
		line = dataFile.readline()
	# interpolate here
	xI = interp1d(tau, x)
	yI = interp1d(tau, y)
	zI = interp1d(tau, z)
	eI = interp1d(tau, e)
	tauI = np.linspace(0, tauMax, num = nData)
	for i in range(len(tauI)):
		print >> stdout, '{"E":%.1f, "x":%.9e, "y":%.9e, "z":%.9e}' % (eI(tauI[i]), xI(tauI[i]), yI(tauI[i]), zI(tauI[i]))  # Log data


if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"

