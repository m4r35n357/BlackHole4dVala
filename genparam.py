#!/usr/bin/env python

from sys import argv, stdout, stderr
from math import fabs, sin, cos, pi, sqrt
import numpy as np
from scipy.optimize import minimize

class InitialConditions(object):
    def __init__(self, particle, rMin, rMax, thetaMin, a, factorL, integrator):
        self.M = 1.0        
	self.mu = 1.0 if (particle == True) else 0.0
        tolerance = 1.0e-6
	nonsingular = fabs(rMax - rMin) > 2.0 * tolerance;
	self.r0 = rMin if nonsingular else rMin - tolerance
	self.r1 = rMax if nonsingular else rMax + tolerance
	self.th0 = thetaMin if thetaMin > 0.01 else 0.01
#	self.th0 = pi - thetaMin 
        self.a = a
	self.factorL = factorL
	self.integrator = integrator
        self.duration = 20.0
        self.timestep = 0.01
        self.E = 1.0
        self.L = 2.0
        self.Q = 0.0

    def rDot2 (self, r):
	return ((r**2 + self.a**2) * self.E - self.a * self.L)**2 - (r**2 - 2.0 * self.M * r + self.a**2) * (self.mu**2 * r**2 + (self.L - self.a * self.E)**2 + self.Q)

    def thDot2 (self, theta):
	return self.Q - cos(theta)**2 * (self.a**2 * (self.mu**2 - self.E**2) + self.L**2 / sin(theta)**2)

    def delta (self, r):
        return r**2 - 2.0 * self.M * r + self.a**2

    def PA (self, r, E, L):
        return (r**2 + self.a**2) * E - self.a * L

    def PB (self, r, E, L, Q):
        return Q + (L - self.a * E)**2 + self.mu**2 * r**2

    def constantR (self, x):
        E = x[0]
        L = x[1]
        Q = x[2]
        return (self.PA(self.r0, E, L)**2 - self.delta(self.r0) * self.PB(self.r0, E, L, Q))**2 + \
               (4.0 * self.r0 * E * self.PA(self.r0, E, L) - 2.0 * (self.r0 - self.M) * self.PB(self.r0, E, L, Q) - 2.0 * self.mu**2 * self.r0 * self.delta(self.r0))**2 + \
               (Q - cos(self.th0)**2 * (self.a**2 * (self.mu**2 - E**2) + L**2 / sin(self.th0)**2))**2

    def variableR (self, x):
        E = x[0]
        L = x[1]
        Q = x[2]
        return (self.PA(self.r0, E, L)**2 - self.delta(self.r0) * self.PB(self.r0, E, L, Q))**2 + \
               (self.PA(self.r1, E, L)**2 - self.delta(self.r1) * self.PB(self.r1, E, L, Q))**2 + \
               (Q - cos(self.th0)**2 * (self.a**2 * (self.mu**2 - E**2) + L**2 / sin(self.th0)**2))**2

    def solve (self, function):
        res = minimize(function, np.array([0.0, 0.0, 0.0]), method='Nelder-Mead', options={'xtol': 1e-12, 'ftol': 1e-12, 'maxiter': 1.0e6, 'maxfev': 1.0e6, 'disp': False})
        #print(res.x)
        self.E = res.x[0]
        self.L = res.x[1]
        self.Q = res.x[2]
        
    def plummet (self):
        self.L = 0.0
        self.Q = 0.0
        a2 = self.a**2
	while self.rDot2(self.r1)**2 > 1.0e-12:
	    self.E -= self.rDot2(self.r1) / (2.0 * self.E * (self.r1**2 + a2)**2 - 2.0 * a2 * self.E * (self.r1**2 - 2.0 * self.M * self.r1 + a2))

def main ():
    if len(argv) == 6:
        ic = InitialConditions(True, float(argv[1]), float(argv[2]), float(argv[3]) * pi, float(argv[4]), float(argv[5]), 8)
	ic.solve(ic.variableR)
        rValue = 0.5 * (ic.r0 + ic.r1)
        thValue = 0.5 * pi
    elif len(argv) == 5:
        ic = InitialConditions(True, float(argv[1]), 0.0, float(argv[2]) * pi, float(argv[3]), float(argv[4]), 8)
	ic.solve(ic.constantR)
        rValue = ic.r0
        thValue = 0.5 * pi
    elif len(argv) == 4:
        ic = InitialConditions(True, 0.0, float(argv[1]), float(argv[2]) * pi, float(argv[3]), 1.0, 8)
        ic.plummet()
        rValue = ic.r1
        thValue = ic.th0
    else:
        print >> stderr, "Bad input data!"
        return
    print >> stdout, "{ \"M\" : " + str(ic.M) + ","
    print >> stdout, "  \"a\" : " + str(ic.a) + ","
    print >> stdout, "  \"mu\" : " + str(ic.mu) + ","
    print >> stdout, "  \"E\" : " + str(ic.E) + ","
    print >> stdout, "  \"Lz\" : " + str(ic.L * ic.factorL) + ","
    print >> stdout, "  \"C\" : " + str(ic.Q) + ","
    print >> stdout, "  \"r\" : " + str(rValue) + ","
    print >> stdout, "  \"theta\" : " + str(thValue) + ","
    print >> stdout, "  \"time\" : " + str(ic.duration) + ","
    print >> stdout, "  \"step\" : " + str(ic.timestep) + ","
    print >> stdout, "  \"integratorOrder\" : " + str(ic.integrator)
    print >> stdout, "}"
    rscale = rValue + 10
    thscale = 0.5 * pi
    for x in range(0, 1000):
        print >> stderr, "{ \"r\":" + str(0.001 * x * rscale) + ", \"R\":" + str(ic.rDot2(0.001 * x * rscale)) + ", \"theta\":" + str(0.001 * x * rscale) + ", \"THETA\":" + str(ic.thDot2(0.5 * pi + 0.001 * x * thscale)) + " }"

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"


