#!/usr/bin/env python
'''
Copyright (c) 2014, 2015, 2016, 2017 Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''
from sys import argv, stdin, stderr
from math import log10, sqrt
from visual import scene, sphere, curve, points, rate
from json import loads
from os import popen

def main():
    print "N-Body Plotter: {}".format(argv)
    # scene basics
    scene.center = (0,0,0)
    scene.width = scene.height = 1024
    #scene.range = (10.0, 10.0, 10.0)
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
        #ball = sphere(pos = (p['qX'], p['qY'], p['qZ']), radius = 0.1 * p['mass']**(1.0 / 3.0), color = colours[j])
        ball = sphere(pos = (p['qX'], p['qY'], p['qZ']), radius = r, color = colours[j])
        #ball = sphere(pos = (p['qX'], p['qY'], p['qZ']), radius = 10000000.0, color = colours[j])
        ball.trail = points(color = ball.color, size = 1)
        spheres.append(ball)
    counter = 0
    while line:
        rate(60)
        #if counter % 1000 == 0:
        #    ball.visible = False
        #    ball = sphere(radius = 0.2)  # Particle
        #    ball.trail = curve(size = 1)  #  trail
        bodies = loads(line)
        X = Y = Z = mT = 0.0
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
        #popen('import -window 0x3200003 -compress None VPythonOutput/' + str(counter).zfill(4) + '.png')
        counter += 1
        line = stdin.readline()

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"

