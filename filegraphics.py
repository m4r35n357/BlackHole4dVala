#!/usr/bin/env python

from sys import argv, stderr, exit, stdout
from math import sqrt, sin, cos, fabs
from visual import scene, sphere, curve, points, rate, ellipsoid, ring, color
from json import loads
from array import array
from scipy.interpolate import interp1d, InterpolatedUnivariateSpline
import numpy as np

def isco (a):
    z1 = 1.0 + pow(1.0 - a * a, 1.0 / 3.0) * (pow(1.0 + a, 1.0 / 3.0) + pow(1.0 - a, 1.0 / 3.0))
    z2 = sqrt(3.0 * a * a + z1 * z1)
    if a >= 0.0:
        return 3.0 + z2 - sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))
    else:
        return 3.0 + z2 + sqrt((3.0 - z1) * (3.0 + z1 + 2.0 * z2))

def main():
    if len(argv) < 4:
        raise Exception('>>> ERROR! Please supply a data file name, a parameter file name, and a time variable <<<')
    dataFile = open(argv[1], 'r')
    parameterFile = open(argv[2], 'r')
    timeCoordinate = str(argv[3])
    try:  # get parameters
        parameters = loads(parameterFile.read())
    except ValueError as e:
        print('PARAMETER ERROR: ' + str(e))
        exit(-1)
    a = parameters['a']
    m = parameters['M']
    horizon = m * (1.0 + sqrt(1.0 - a * a));
    cauchy = m * (1.0 - sqrt(1.0 - a * a));
    #  set up the scene
    scene.center = (0.0, 0.0, 0.0)
    scene.width = 1024
    scene.height = 1024
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
#        for j in range(2, 21, 2):
#            ring(pos=scene.center, axis=(0, 0, 1), radius = j, color = color.gray(0.3), thickness=0.01)
#            ring(pos=scene.center, axis=(0, 1, 0), radius = j, color = color.gray(0.3), thickness=0.01)
#            ring(pos=scene.center, axis=(1, 0, 0), radius = j, color = color.gray(0.3), thickness=0.01)
    # animate!
    ball = sphere()  # Particle
    counter = 0
    dataLine = dataFile.readline()
    while dataLine:  # build raw data arrays
        data = loads(dataLine)
        r = float(data['r'])
        th = float(data['th'])
        ph = float(data['ph'])
        ra = sqrt(r**2 + a**2)
        sth = sin(th)
        x = ra * sth * cos(ph)
        y = ra * sth * sin(ph)
        z = r * cos(th)
        e = float(data['v4e'])
        if counter % 1000 == 0:
            ball.visible = False
            ball = sphere(radius = 0.2)  # Particle
            ball.trail = curve(size = 1)  #  trail
        rate(60)
        if e < -120.0:
            ball.color = color.green
        elif e < -90.0:
            ball.color = color.cyan
        elif e < -60.0:
            ball.color = color.yellow
        elif e < -30.0:
            ball.color = color.orange
        else:
            ball.color = color.red
        ball.pos = (x, y, z)
        ball.trail.append(pos = ball.pos, color = ball.color)
        counter += 1
        dataLine = dataFile.readline()

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"

