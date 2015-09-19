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
from sys import argv, stdout, stderr, exit
from math import fabs, sin, cos, pi, sqrt, copysign
from array import array
import numpy as np
from scipy.optimize import minimize

class InitialConditions(object):
    def __init__(self, particle, rMin, rMax, thetaMin, a, factorL):
        self.M = 1.0        
        self.mu2 = 1.0 if (particle == True) else 0.0
	self.r0 = rMin
	self.r1 = rMax
	self.th0 = pi * (1.0 - thetaMin) 
        self.a = a
        self.a2 = self.a**2
	self.factorL = factorL
	self.integrator = 'sb2'
        self.starttime = 0.0
        self.duration = 50.0
        self.timestep = 0.001
        self.ic = np.array([1.0, copysign(5.0, a), 0.0]) if a >= 0.0 else np.array([1.0, - copysign(5.0, a), 0.0])
        self.ic = np.array([1.0, 0.0, 5.0]) if thetaMin < 0.01 else self.ic

    def coefficients (self, E, L, Q):
        E2_1 = E**2 - self.mu2
        return array('d', [E2_1, 2.0 * self.mu2, self.a2 * E2_1 - L**2 - Q, 2.0 * ((self.a * E - L)**2 + Q), - self.a2 * Q])

    def R (self, r, c):
        return (((c[0] * r + c[1]) * r + c[2]) * r + c[3]) * r + c[4]

    def dR (self, r, c):
        return ((4.0 * c[0] * r + 3.0 * c[1]) * r + 2.0 * c[2]) * r + c[3]

    def THETA (self, th, E, L, Q):
        return Q - cos(th)**2 * (self.a2 * (self.mu2 - E**2) + (L / sin(th))**2)

    def constantR (self, x):
        c = self.coefficients(x[0], x[1], x[2])
        return self.R(self.r0, c)**2 + self.dR(self.r0, c)**2 + self.THETA(self.th0, x[0], x[1], x[2])**2

    def variableR (self, x):
        c = self.coefficients(x[0], x[1], x[2])
        return self.R(self.r0, c)**2 + self.R(self.r1, c)**2 + self.THETA(self.th0, x[0], x[1], x[2])**2

    def solve (self, function):
        res = minimize(function, self.ic, method='Nelder-Mead', options={'xtol': 1e-12, 'ftol': 1e-12, 'maxiter': 1.0e6, 'maxfev': 1.0e6, 'disp': False})
        self.fun = res.fun
        self.message = res.message
        self.success = res.success
        self.E = res.x[0]
        self.L = res.x[1]
        self.Q = res.x[2]
        
def main ():
    if len(argv) == 6:
        ic = InitialConditions(True, float(argv[1]), float(argv[2]), float(argv[3]), float(argv[4]), float(argv[5]))
        ic.solve(ic.variableR)
    elif len(argv) == 5:
        ic = InitialConditions(True, float(argv[1]), 0.0, float(argv[2]), float(argv[3]), float(argv[4]))
        ic.solve(ic.constantR)
    elif len(argv) == 4:
        ic = InitialConditions(True, float(argv[1]), 0.0, float(argv[2]), float(argv[3]), 0.0)
        ic.solve(ic.constantR)
    else:
        print >> stderr, "Bad input data!"
        return
    # Initial conditions file
    print >> stdout, "{ \"M\" : " + str(ic.M) + ","
    print >> stdout, "  \"a\" : " + str(ic.a) + ","
    print >> stdout, "  \"mu\" : " + str(ic.mu2) + ","
    print >> stdout, "  \"E\" : " + repr(ic.E) + ","
    print >> stdout, "  \"Lz\" : " + repr(ic.factorL * ic.L) + ","
    print >> stdout, "  \"C\" : " + repr(ic.Q) + ","
    print >> stdout, "  \"r\" : " + repr(ic.r0) + ","
    print >> stdout, "  \"theta\" : " + repr(0.5 * pi) + ","
    print >> stdout, "  \"start\" : " + str(ic.starttime) + ","
    print >> stdout, "  \"duration\" : " + str(ic.duration) + ","
    print >> stdout, "  \"step\" : " + str(ic.timestep) + ","
    print >> stdout, "  \"integrator\" : \"" + ic.integrator + "\","
    print >> stdout, "  \"error\" : " + str(ic.fun) + ","
    print >> stdout, "  \"success\" : \"" + str(ic.success) + "\","
    print >> stdout, "  \"message\" : \"" + str(ic.message) + "\""
    print >> stdout, "}"
    # Potential plot data
    c = ic.coefficients(ic.E, ic.L, ic.Q)
    rMax = ic.r0 if (ic.r0 > ic.r1) else ic.r1
    nPoints = 1001
    for x in range(1, nPoints - 1):
        xValue = 1.0 * x / nPoints
        print >> stderr, "{ \"x\":" + str(xValue * rMax) \
                       + ", \"R\":" + str(ic.R(xValue * rMax, c)) \
                       + ", \"THETA\":" + str(ic.THETA(xValue * pi, ic.E, ic.L, ic.Q)) + " }"
    # Return value
    if not ic.success or ic.fun > 1.0e-6:
        exit(-1)

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"

