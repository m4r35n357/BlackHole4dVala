#!/usr/bin/env pypy
'''
Copyright (c) 2014, 2015, 2016, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''
from sys import stdin, stdout, stderr
from math import fabs, log10, sqrt
from json import loads
from array import array
from symplectic import Integrator, logError

class Particle(object):

    def __init__(self, qX, qY, qZ, pX, pY, pZ, mass):
        self.qX = qX
        self.qY = qY
        self.qZ = qZ
        self.pX = pX
        self.pY = pY
        self.pZ = pZ
        self.mass = mass

    def __str__(self):
        return "{\"qX\":%.6e,\"qY\":%.6e,\"qZ\":%.6e,\"pX\":%.6e,\"pY\":%.6e,\"pZ\":%.6e,\"mass\":%.3e}" % (self.qX, self.qY, self.qZ, self.pX, self.pY, self.pZ, self.mass)

class NBody(object):
    def __init__(self, g, runTime, timeStep, errorLimit, bodies, order):
        self.bodies = bodies
        self.pRange = range(len(bodies))
        self.g = g
        self.ts = timeStep
        self.eMax = errorLimit
        self.T = runTime
        self.integrator = Integrator(self, order)

    @staticmethod
    def dist (xA, yA, zA, xB, yB, zB):  # Euclidean distance between point A and point B
        return sqrt((xB - xA)**2 + (yB - yA)**2 + (zB - zA)**2)

    def h (self):  # Conserved energy
        energy = 0.0
        for i in self.pRange:
            a = self.bodies[i]
            energy += 0.5 * (a.pX**2 + a.pY**2 + a.pZ**2) / a.mass
            for j in self.pRange:
                if i > j:
                    b = self.bodies[j]
                    energy -= self.g * a.mass * b.mass / self.dist(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ)
        return energy

    def qUp (self, d):  # Update Positions
        for i in self.pRange:
            a = self.bodies[i]
            tmp = d * self.ts / a.mass
            a.qX += a.pX * tmp
            a.qY += a.pY * tmp
            a.qZ += a.pZ * tmp

    def pUp (self, c):  # Update Momenta
        for i in self.pRange:
            a = self.bodies[i]
            for j in self.pRange:
                if i > j:
                    b = self.bodies[j]
                    tmp = - c * self.ts * self.g * a.mass * b.mass / self.dist(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ)**3
                    dPx = tmp * (b.qX - a.qX)
                    dPy = tmp * (b.qY - a.qY)
                    dPz = tmp * (b.qZ - a.qZ)
                    a.pX -= dPx
                    a.pY -= dPy
                    a.pZ -= dPz
                    b.pX += dPx
                    b.pY += dPy
                    b.pZ += dPz

    def print_out (self, time, hNow, h0, dbValue):
        data = []
        for i in self.pRange:
            data.append(str(self.bodies[i]))
        print >> stdout, "[" + ','.join(data) + "]"  # Log data

def icJson ():
    ic = loads(stdin.read())
    bodies = []
    for a in ic['bodies']:
        if 'pX' in a and 'pY' in a and 'pZ' in a:  # momenta specified
            bodies.append(Particle(a['qX'], a['qY'], a['qZ'], a['pX'], a['pY'], a['pZ'], a['mass']))
        elif 'vX' in a and 'vY' in a and 'vZ' in a:  # velocities specified, convert to momenta
            mass = a['mass']
            bodies.append(Particle(a['qX'], a['qY'], a['qZ'], mass * a['vX'], mass * a['vY'], mass * a['vZ'], mass))
        else:
            raise Exception('>>> ERROR! Specify either momenta or velocites consistently <<<')
    return NBody(ic['g'], ic['simulationTime'], ic['timeStep'], ic['errorLimit'], bodies, ic['integratorOrder'])

def main ():  # Need to be inside a function to return . . .
    s = icJson()  # Create a symplectic integrator object from JSON input
    h0 = s.h()  # Set up error reporting
    t = 0.0
    while True:
        hNow = s.h()
        dbValue = logError(fabs(hNow / h0 - 1.0))
        s.print_out(t, hNow, h0, dbValue)
        if fabs(t) > s.T or dbValue > s.eMax:
            return
        s.integrator.compose()
        t += s.ts

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"

