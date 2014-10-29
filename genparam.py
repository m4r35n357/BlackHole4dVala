#!/usr/bin/env python

from sys import argv, stdout, stderr
from math import fabs, sin, cos, pi, sqrt
import numpy as np

class InitialConditions(object):
    def __init__(self, particle, rMin, rMax, thetaMin, a, factorL, integrator):
        self.M = 1.0        
	self.mu = 1.0 if (particle == True) else 0.0
        tolerance = 1.0e-6
	nonsingular = fabs(rMax - rMin) > 2.0 * tolerance;
	self.r0 = rMin if nonsingular else rMin - tolerance
	self.r1 = rMax if nonsingular else rMax + tolerance
	self.th0 = thetaMin if thetaMin > 0.01 else 0.01
        self.a = a
	self.factorL = factorL
	self.integrator = integrator
        self.duration = 20.0
        self.timestep = 0.01
        self.E = 1.0
        self.L = 2.0
        self.Q = 0.0

    def rDot (self, r):
	return ((r**2 + self.a**2) * self.E - self.a * self.L)**2 - (r**2 - 2.0 * self.M * r + self.a**2) * (self.mu**2 * r**2 + (self.L - self.a * self.E)**2 + self.Q)

    def thDot (self, theta):
	return self.Q - cos(theta)**2 * (self.a**2 * (self.mu**2 - self.E**2) + self.L**2 / sin(theta)**2)

    def qDot (self):
	return np.array([self.rDot(self.r0), self.rDot(self.r1), self.thDot(self.th0)])

    def solve (self):
        a2 = self.a**2
	while np.dot(self.qDot(), self.qDot()) > 1.0e-21:
            p0 = self.E * (self.r0**2 + a2) - self.a * self.L
	    p1 = self.E * (self.r1**2 + a2) - self.a * self.L
	    delta0 = self.r0**2 - 2.0 * self.M * self.r0 + a2
	    delta1 = self.r1**2 - 2.0 * self.M * self.r1 + a2
	    l_ae = self.L - self.a * self.E
	    j = np.array([[2.0 * (self.r0**2 + a2) * p0 + 2.0 * self.a * l_ae * delta0, - 2.0 * self.a * p0 - 2.0 * l_ae * delta0, - delta0],
                          [2.0 * (self.r1**2 + a2) * p1 + 2.0 * self.a * l_ae * delta1, - 2.0 * self.a * p1 - 2.0 * l_ae * delta1, - delta1],
                          [2.0 * a2 * self.E * cos(self.th0)**2, - 2.0 * self.L * (cos(self.th0) / sin(self.th0))**2, 1.0]])
            correction = np.linalg.solve(j, self.qDot())
	    self.E -= correction[0]
	    self.L -= correction[1]
	    self.Q -= correction[2]

    def circular (self):
        sqrtR = sqrt(self.r1)
        tmp = sqrt(self.r1**2 - 3.0 * self.r1 + 2.0 * self.a * sqrtR)
        self.E = (self.r1**2 - 2.0 * self.r1 + self.a * sqrtR) / (self.r1 * tmp)
        self.L = (self.r1**2 - 2.0 * self.a * sqrtR + self.a**2) / (sqrtR * tmp)
        self.Q = 0.0

    def plummet (self):
        self.L = 0.0
        self.Q = 0.0
        a2 = self.a**2
	while self.rDot(self.r1)**2 > 1.0e-12:
            p1 = self.E * (self.r1**2 + a2)
	    delta1 = self.r1**2 - 2.0 * self.M * self.r1 + a2
	    self.E -= self.rDot(self.r1) / (2.0 * (self.r1**2 + a2) * p1 - 2.0 * self.a * self.a * self.E * delta1)

def main ():
    if len(argv) == 6:
        ic = InitialConditions(True, float(argv[1]), float(argv[2]), float(argv[3]) * pi, float(argv[4]), float(argv[5]), 8)
	ic.solve()
    elif len(argv) == 4:
        ic = InitialConditions(True, 0.0, float(argv[1]), 0.5 * pi, float(argv[2]), float(argv[3]), 8)
        ic.circular()
    elif len(argv) == 3:
        ic = InitialConditions(True, 0.0, float(argv[1]), 0.5 * pi, float(argv[2]), 1.0, 8)
        ic.plummet()
    else:
        print >> stderr, "Bad input data!"
        return
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


