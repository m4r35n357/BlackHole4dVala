#!/usr/bin/env python

from sys import argv, stderr, exit
from math import sqrt, sin, fabs
from visual import scene, sphere, curve, points, rate, ellipsoid, ring, color
from json import loads
from array import array
from scipy.interpolate import interp1d
import numpy as np

def isco (a):
	z1 = 1.0 + pow(1.0 - a * a, 1.0 / 3.0) * (pow(1.0 + a, 1.0 / 3.0) + pow(1.0 - a, 1.0 / 3.0))
	z2 = sqrt(3.0 * a * a + z1 * z1)
	if a >= 0.0:
		return 3.0 + z2 - sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))
	else:
		return 3.0 + z2 + sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))

def main():
	if len(argv) < 5:
		raise Exception('>>> ERROR! Please supply a data file name, a parameter file name, a time variable, and the number of points to use <<<')
	dataFile = open(argv[1], 'r')
	parameterFile = open(argv[2], 'r')
	timeCoordinate = str(argv[3])
	nData = int(argv[4])
        try:  # get parameters
            parameters = loads(parameterFile.read())
        except ValueError as e:
            print('PARAMETER ERROR: ' + str(e))
            exit(-1)		
	a = parameters['a']
	m = parameters['M']
	horizon = m * (1.0 + sqrt(1.0 - a * a));
	cauchy = m * (1.0 - sqrt(1.0 - a * a));
	# get raw data
	timeMax = 0.0
	time = array('d')
	x = array('d')
	y = array('d')
	z = array('d')
	e = array('d')
        dataLine = dataFile.readline()
	while dataLine:  # build raw data arrays
		data = loads(dataLine)
		timeValue = data[timeCoordinate]
		timeMax = timeValue if timeValue > timeMax else timeMax
		time.append(timeValue)
		x.append(data['x'])
		y.append(data['y'])
		z.append(data['z'])
		e.append(data['v4e'])
		dataLine = dataFile.readline()
        try:  # interpolate here
		xI = interp1d(time, x, kind='linear', copy=False)
		yI = interp1d(time, y, kind='linear', copy=False)
		zI = interp1d(time, z, kind='linear', copy=False)
		eI = interp1d(time, e, kind='linear', copy=False)
        except ValueError as e:
            print('DATA ERROR: ' + str(argv[0]) + ': ' + str(e))
            exit(-2)		
	#  set up the scene
	scene.center = (0.0, 0.0, 0.0)
	scene.width = scene.height = 1000.0
	scene.range = (20.0, 20.0, 20.0)
        inner = 2.0 * sqrt(cauchy**2 + a**2)
	ellipsoid(pos = scene.center, length = inner, height = inner, width = 2.0 * cauchy, color = color.blue, opacity = 0.4)  # Inner Horizon
        outer = 2.0 * sqrt(horizon**2 + a**2)
	ellipsoid(pos = scene.center, length = outer, height = outer, width = 2.0 * horizon, color = color.blue, opacity = 0.3)  # Outer Horizon
        ergo = 2.0 * sqrt(4.0 + a**2)
	ellipsoid(pos = scene.center, length = ergo, height = ergo, width = 2.0 * horizon, color = color.gray(0.7), opacity = 0.2)  # Ergosphere
        if fabs(a) > 0.0:
	    ring(pos=scene.center, axis=(0, 0, 1), radius = a, color = color.white, thickness=0.01)  # Singularity
        else:
            sphere(pos=scene.center, radius = 0.05, color = color.white)  # Singularity
	ring(pos=scene.center, axis=(0, 0, 1), radius = sqrt(isco(a)**2 + a**2), color = color.magenta, thickness=0.01)  # ISCO
	# animate!
	ball = sphere()  # Particle
	counter = 0
	timeI = np.linspace(0, timeMax, num = nData)
	for i in range(len(timeI)):
		if counter % 1000 == 0:
			ball.visible = False
			ball = sphere(radius = 0.2)  # Particle
			ball.trail = curve(size = 1)  #  trail
		rate(60)
		error = eI(timeI[i])
		if error < -120.0:
			ball.color = color.green
		elif error < -90.0:
			ball.color = color.cyan
		elif error < -60.0:
			ball.color = color.yellow
		elif error < -30.0:
			ball.color = color.orange
		else:
			ball.color = color.red
		ball.pos = (xI(timeI[i]), yI(timeI[i]), zI(timeI[i]))
		ball.trail.append(pos = ball.pos, color = ball.color)
		counter += 1

if __name__ == "__main__":
	main()
else:
	print >> stderr, __name__ + " module loaded"

