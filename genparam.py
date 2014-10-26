#!/usr/bin/env python

from sys import stdin, stdout, stderr
from math import fabs, log10, sqrt, sin, cos, pi
from json import loads
from array import array
import pprint
import numpy as np
from scipy import linalg

class InitialConditions(object):

    def __init__(self, particle, rMin, rMax, thetaMin, a, factorL, integrator):
        self.M = 1.0        
	self.mu = 1.0 if (particle == True) else 0.0
        tolerance = 1.0e-6
	singular = fabs(rMax - rMin) > 2.0 * tolerance;
	self.r0 = rMin - tolerance if singular else rMin
	self.r1 = rMax - tolerance if singular else rMax
	self.th0 = thetaMin if thetaMin > 0.01 else 0.01
        self.a = a
	self.factorL = factorL
	self.integrator = integrator
        self.duration = 20.0
        self.timestep = 0.001
        self.E = 1.0
        self.L = 2.0
        self.Q = 0.0

    def rDot (self, r):
	return ((r * r + self.a * self.a) * self.E - self.a * self.L) * ((r * r + self.a * self.a) * self.E - self.a * self.L) - (r * r - 2.0 * self.M * r + self.a * self.a) * (self.mu * self.mu * r * r + (self.L - self.a * self.E) * (self.L - self.a * self.E) + self.Q)

    def thDot (self, theta):
	return self.Q - cos(theta) * cos(theta) * (self.a * self.a * (self.mu * self.mu - self.E * self.E) + self.L * self.L / (sin(theta) * sin(theta)))

    def qDot (self):
	return np.array([self.rDot(self.r0), self.rDot(self.r1), self.thDot(self.th0)])

    def solve (self):
	while np.dot(self.qDot(), self.qDot()) > 1.0e-18:
            p0 = self.E * (self.r0 * self.r0 + self.a * self.a) - self.a * self.L
	    p1 = self.E * (self.r1 * self.r1 + self.a * self.a) - self.a * self.L
	    delta0 = self.r0 * self.r0 - 2.0 * self.M * self.r0 + self.a * self.a
	    delta1 = self.r1 * self.r1 - 2.0 * self.M * self.r1 + self.a * self.a
	    l_ae = (self.L - self.a * self.E)
	    j = np.array([[2.0 * (self.r0 * self.r0 + self.a * self.a) * p0 + 2.0 * self.a * l_ae * delta0, - 2.0 * self.a * p0 - 2.0 * l_ae * delta0, - delta0],
                          [2.0 * (self.r1 * self.r1 + self.a * self.a) * p1 + 2.0 * self.a * l_ae * delta1, - 2.0 * self.a * p1 - 2.0 * l_ae * delta1, - delta1],
                          [2.0 * cos(self.th0) * cos(self.th0) * self.a * self.a * self.E, - 2.0 * cos(self.th0) * cos(self.th0) * self.L / (sin(self.th0) * sin(self.th0)), 1.0]])
            print "J:"
            pprint.pprint(j)
            correction = np.linalg.solve(j, self.qDot())
	    self.E -= correction[0]
	    self.L -= correction[1]
	    self.Q -= correction[2]

def main ():
	ic = InitialConditions(True, 9.0, 12.0, pi / 2.0, -1.0, 1.0, 8)
	ic.solve()
	print >> stdout, ""
	print >> stdout, "{ \"M\" : " + str(ic.M) + ","
	print >> stdout, "  \"a\" : " + str(ic.a) + ","
	print >> stdout, "  \"mu\" : " + str(ic.mu) + ","
	print >> stdout, "  \"E\" : " + str(ic.E) + ","
	print >> stdout, "  \"Lz\" : " + str(ic.L * ic.factorL) + ","
	print >> stdout, "  \"C\" : " + str(ic.Q) + ","
	print >> stdout, "  \"r\" : " + str(ic.r1) + ","
	print >> stdout, "  \"theta\" : " + str(ic.th0) + ","
	print >> stdout, "  \"time\" : " + str(ic.duration) + ","
	print >> stdout, "  \"step\" : " + str(ic.timestep) + ","
	print >> stdout, "  \"integratorOrder\" : " + str(ic.integrator)
	print >> stdout, "}"

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"


