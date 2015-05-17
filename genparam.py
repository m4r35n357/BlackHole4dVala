#!/usr/bin/env python
'''
Copyright (c) 2014, 2015, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''
from sys import argv, stdout, stderr
from math import fabs, sin, cos, pi, sqrt
import numpy as np
from scipy.optimize import minimize

class InitialConditions(object):
    def __init__(self, particle, rMin, rMax, thetaMin, a, factorL, integrator):
        self.M = 1.0        
	self.mu = 1.0 if (particle == True) else 0.0
        self.mu2 = self.mu**2
	self.r0 = rMin
	self.r1 = rMax
	self.th0 = pi * (1.0 - thetaMin) 
        self.a = - a  # convention
        self.a2 = self.a**2
	self.factorL = factorL
	self.integrator = integrator
        self.duration = 20.0
        self.timestep = 0.01

    def delta (self, r):
        return r**2 - 2.0 * self.M * r + self.a2

    def PA (self, r, E, L):
        return (r**2 + self.a2) * E - self.a * L

    def PB (self, r, E, L, Q):
        return Q + (L - self.a * E)**2 + self.mu2 * r**2

    def THETA (self, th, E, L, Q):
        return Q - cos(th)**2 * (self.a2 * (self.mu2 - E**2) + (L / sin(th))**2)

    def constantR (self, x):
        E = x[0]
        L = x[1]
        Q = x[2]
        return (self.PA(self.r0, E, L)**2 - self.delta(self.r0) * self.PB(self.r0, E, L, Q))**2 + \
               (4.0 * self.r0 * E * self.PA(self.r0, E, L) - \
                2.0 * (self.r0 - self.M) * self.PB(self.r0, E, L, Q) - \
                2.0 * self.mu2 * self.r0 * self.delta(self.r0))**2 + \
               (self.THETA(self.th0, E, L, Q))**2

    def variableR (self, x):
        E = x[0]
        L = x[1]
        Q = x[2]
        return (self.PA(self.r0, E, L)**2 - self.delta(self.r0) * self.PB(self.r0, E, L, Q))**2 + \
               (self.PA(self.r1, E, L)**2 - self.delta(self.r1) * self.PB(self.r1, E, L, Q))**2 + \
               (self.THETA(self.th0, E, L, Q))**2

    def solve (self, function):
        res = minimize(function, np.array([0.0, 0.0, 0.0]), method='Nelder-Mead', \
                       options={'xtol': 1e-12, 'ftol': 1e-12, 'maxiter': 1.0e6, 'maxfev': 1.0e6, 'disp': False})
        #print(res.x)
        self.E = res.x[0]
        self.L = res.x[1]
        self.Q = res.x[2]
        
def main ():
    if len(argv) == 6:
        ic = InitialConditions(True, float(argv[1]), float(argv[2]), float(argv[3]), float(argv[4]), float(argv[5]), 8)
	ic.solve(ic.variableR)
        rValue = 0.5 * (ic.r0 + ic.r1)
        thValue = 0.5 * pi
    elif len(argv) == 5:
        ic = InitialConditions(True, float(argv[1]), 0.0, float(argv[2]), float(argv[3]), float(argv[4]), 8)
	ic.solve(ic.constantR)
        rValue = ic.r0
        thValue = 0.5 * pi
    elif len(argv) == 4:
        ic = InitialConditions(True, float(argv[1]), 0.0, float(argv[2]), float(argv[3]), 0.0, 8)
	ic.solve(ic.constantR)
        rValue = ic.r0
        thValue = ic.th0
    else:
        print >> stderr, "Bad input data!"
        return
    print >> stdout, "{ \"M\" : " + str(ic.M) + ","
    print >> stdout, "  \"a\" : " + str(- ic.a) + ","  # convention
    print >> stdout, "  \"mu\" : " + str(ic.mu) + ","
    print >> stdout, "  \"E\" : " + str(ic.E) + ","
    print >> stdout, "  \"Lz\" : " + str(- ic.L * ic.factorL) + ","  # convention
    print >> stdout, "  \"C\" : " + str(ic.Q) + ","
    print >> stdout, "  \"r\" : " + str(rValue) + ","
    print >> stdout, "  \"theta\" : " + str(thValue) + ","
    print >> stdout, "  \"time\" : " + str(ic.duration) + ","
    print >> stdout, "  \"step\" : " + str(ic.timestep) + ","
    print >> stdout, "  \"integratorOrder\" : " + str(ic.integrator)
    print >> stdout, "}"
    rscale = rValue + 5.0
    thscale = 0.5 * pi
    nPoints = 1000 + 1
    for x in range(0, nPoints):
        scaledX = 1.0 * x / nPoints
        print >> stderr, "{ \"x\":" + str(scaledX * rscale) \
                       + ", \"R\":" + str(ic.PA(scaledX * rscale, ic.E, ic.L)**2 - ic.delta(scaledX * rscale) * ic.PB(scaledX * rscale, ic.E, ic.L, ic.Q)) \
                       + ", \"THETA\":" + str(ic.THETA(0.5 * pi + scaledX * thscale, ic.E, ic.L, ic.Q)) + " }"

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"


