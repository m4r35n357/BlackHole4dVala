#!/usr/bin/env python

from sys import argv, stdin
from math import log10, sqrt
from visual import scene, sphere, points, rate
from json import loads

def main():
    # scene basics
    scene.center = (0,0,0)
    scene.width = 1800
    scene.height = 1000
#   scene.range = (10.0, 10.0, 10.0)
    scene.range = (1.0e9, 1.0e9, 1.0e9)
    # get data dimensions
    line = stdin.readline()
    bodies = loads(line)
    pRange = range(len(bodies))
    #  set up the balls
    colours = [ (1.0, 1.0, 0.0), (1.0, 1.0, 1.0), (0.0, 1.0, 0.0), (0.0, 0.5, 0.5), (1.0, 0.0, 0.0), (0.5, 1.0, 0.0), (1.0, 0.0, 1.0), (1.0, 0.5, 0.0), (0.0, 1.0, 1.0), (1.0, 1.0, 1.0) ]
    spheres = []
    for j in pRange:
        p = bodies[j]
        r = 1000000.0 * log10(p['mass']**(1.0 / 3.0))
#       ball = sphere(pos = (p['qX'], p['qY'], p['qZ']), radius = 0.1 * p['mass']**(1.0 / 3.0), color = colours[j])
        ball = sphere(pos = (p['qX'], p['qY'], p['qZ']), radius = r, color = colours[j])
#       ball = sphere(pos = (p['qX'], p['qY'], p['qZ']), radius = 10000000.0, color = colours[j])
        ball.trail = points(color = ball.color, size = 1)
        spheres.append(ball)
    while line:
        rate(60)
        bodies = loads(line)
        X = Y = Z = mT = 0.0;
        for j in pRange:  # COG correction
            p = bodies[j]
            X += p['qX'] * p['mass']
            Y += p['qY'] * p['mass']
            Z += p['qZ'] * p['mass']
            mT += p['mass']
        for j in pRange:
            p = bodies[j]
            ball = spheres[j]
            position = (p['qX'] - X / mT, p['qY'] - Y / mT, p['qZ'] - Z / mT)
            ball.pos = position
            ball.trail.append(pos = position, retain = 1000)
        line = stdin.readline()

if __name__ == "__main__":
    main()
else:
    print >> sys.stderr, __name__ + " module loaded"

