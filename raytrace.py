#!/usr/bin/env python

from sys import argv, stdin, stdout, stderr
from math import sqrt, sin, cos, pi
from json import loads
from array import array
from scipy.interpolate import interp1d
import numpy as np

def main():
	if len(argv) < 2:
		raise Exception('>>> ERROR! Please supply number of points to plot and a viewing angle <<<')
	nData = int(argv[1])
	line = stdin.readline()
	tau = array('d')
	t = array('d')
	r = array('d')
	th = array('d')
	ph = array('d')
	tDot = array('d')
	rDot = array('d')
	thDot = array('d')
	phDot = array('d')
	tauMax = 0.0
	while line:  # build raw data arrays line by line
		p = loads(line)
		tauValue = p['tau']
		tauMax = tauValue if tauValue > tauMax else tauMax
		tau.append(tauValue)
		t.append(p['t'])
		r.append(p['r'])
		th.append(p['th'])
		ph.append(p['ph'])
		tDot.append(p['tDot'])
		rDot.append(p['rDot'])
		thDot.append(p['thDot'])
		phDot.append(p['phDot'])
		line = stdin.readline()	
        try: # interpolate raw data arrays
		tI = interp1d(tau, t, kind='linear')
		rI = interp1d(tau, r, kind='linear')
		thI = interp1d(tau, th, kind='linear')
		phI = interp1d(tau, ph, kind='linear')
		tDotI = interp1d(tau, tDot, kind='linear')
		rDotI = interp1d(tau, rDot, kind='linear')
		thDotI = interp1d(tau, thDot, kind='linear')
		phDotI = interp1d(tau, phDot, kind='linear')
        except ValueError as e:
            print('DATA ERROR: ' + str(argv[0]) + ': ' + str(e))
            exit(-2)		
	x = np.linspace(0, tauMax, num = nData)  # create evenly spaced time values
	angle = long(argv[2]) * pi / 180.0
	c = cos(angle)
	s = sin(angle)
	for i in range(len(x)):
		print >> stdout, str(x[i]) + ' 2 ' \
			+ str(rI(x[i])) + ' ' + str(cos(thI(x[i]))) + ' ' + str(tI(x[i])) + ' ' + str(phI(x[i])) + ' ' \
			+ str(rDotI(x[i])) + ' ' + str(- sin(thI(x[i])) * thDotI(x[i])) + ' ' + str(tDotI(x[i])) + ' ' + str(phDotI(x[i])) + ' ' \
			+ str(-c) + ' ' + str(0) + ' ' + str(0) + ' ' + str(s) + ' ' \
			+ str(s) + ' ' + str(0) + ' ' + str(0) + ' ' + str(c) + ' ' \
			+ str(0) + ' ' + str(1) + ' ' + str(0) + ' ' + str(0)

if __name__ == "__main__":
	main()
else:
	print >> sys.stderr, __name__ + " module loaded"
