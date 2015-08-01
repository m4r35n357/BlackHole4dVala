#!/usr/bin/env pypy
'''
Copyright (c) 2014, 2015, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''
from sys import argv, stdin, stdout, stderr
from math import fabs, log10, sqrt, sin, cos, pi
from json import loads
from array import array

class BL(object):   # Boyer-Lindquist coordinates on the Kerr le2
    def __init__(self, bhMass, spin, pMass2, energy, momentum, carter, r0, thetaMin, simtime, timestep, order):
    	self.a = spin
        self.mu2 = pMass2
    	self.E = energy
    	self.L = momentum
    	self.Q = carter
    	self.r = r0
    	self.th = thetaMin
    	self.h = timestep
        self.cbrt2 = 2.0**(1.0 / 3.0)
        self.cbrt2one = 1.0 - self.cbrt2
        self.cbrt2two = 2.0 - self.cbrt2
	if order == '2b':  # Second order base
            self.base = self.stormerVerlet2;
            self.coeff = array('d', [1.0])
	elif order == '4c':  # Fourth order, composed from Second order
            self.base = self.stormerVerlet2;
            self.coeff = array('d', [1.0 / self.cbrt2two, - self.cbrt2 / self.cbrt2two, 1.0 / self.cbrt2two])
	elif order == '4b':  # Fourth order base
            self.base = self.stormerVerlet4;
            self.coeff = array('d', [1.0])
	elif order == '6c':  # Sixth order, composed from Fourth order
            self.base = self.stormerVerlet4;
            fthrt2 = 2.0**(1.0 / 5.0)
            self.coeff = array('d', [1.0 / (2.0 - fthrt2), - fthrt2 / (2.0 - fthrt2), 1.0 / (2.0 - fthrt2)])
	else:  # Wrong value for integrator order
            raise Exception('>>> ERROR! Integrator order must be 2b, 4c, 4b or 6c <<<')
        self.simtime = abs(simtime)
        self.t = self.ph = 0.0
    	self.a2 = self.a**2
        self.aE = self.a * self.E
        self.a2E = self.a2 * self.E
        self.L2 = self.L**2
        self.aL = self.a * self.L
        E2_mu2 = self.E**2 - self.mu2
        self.c = array('d', [E2_mu2, 2.0 * self.mu2, self.a2 * E2_mu2 - self.L2 - self.Q, 2.0 * ((self.aE - self.L)**2 + self.Q), - self.a2 * self.Q])
        self.a2xE2_mu2 = - self.a2 * E2_mu2
        self.coefficients = range(len(self.coeff))

    def errors (self, R, THETA, tDot, rDot, thDot, phDot):  # Error analysis
        def logError (e):
            return 10.0 * log10(e if e > 1.0e-15 else 1.0e-15) 
        def modH (xDot, X):
            return 0.5 * fabs(xDot**2 - X)
        def v4Error (tDot, rDot, thDot, phDot):  # norm squared, xDot means dx/dTau !!!
            return fabs(self.mu2 + self.sth2 / self.S * (self.a * tDot - self.ra2 * phDot)**2 + self.S / self.D * rDot**2 + self.S * thDot**2 - self.D / self.S * (tDot - self.a * self.sth2 * phDot)**2)
        self.eR = logError(modH(rDot, R))
        self.eTh = logError(modH(thDot, THETA))
        self.v4e = logError(v4Error(tDot / self.S, rDot / self.S, thDot / self.S, phDot / self.S))  # d/dTau = 1/sigma * d/dLambda !!!

    def refresh (self, r, th):  # Update quantities that depend on r or theta
        self.sth = sin(th)
        self.cth = cos(th)
        self.sth2 = self.sth**2
        cth2 = 1.0 - self.sth2
        self.ra2 = r**2 + self.a2
	self.D = (r - 2.0) * r + self.a2
	self.S = r**2 + self.a2 * cth2
        self.R = (((self.c[0] * r + self.c[1]) * r + self.c[2]) * r + self.c[3]) * r + self.c[4]
	self.TH = self.a2xE2_mu2 + self.L2 / self.sth2
	self.THETA = self.Q - cth2 * self.TH
	P = self.ra2 * self.E - self.aL
        self.tDot = self.ra2 * P / self.D + self.aL - self.a2E * self.sth2
        self.phDot = self.a * P / self.D - self.aE + self.L / self.sth2
	
    def qUp (self, c):  # Coordinate updates, x += dH/dxDot (i.e. dT/dxDot).  N.B. x = r or theta; t and phi are just along for the ride . . .
        self.t += c * self.tDot
        self.r += c * self.rDot
        self.th += c * self.thDot
        self.ph += c * self.phDot
        self.refresh(self.r, self.th)

    def qDotUp (self, c):  # Velocity (momentum) updates, xDot += -dH/dx (i.e. dV/dx, minus sign cancels with the one in the pseudo-Hamiltonian)
        self.rDot += c * (((4.0 * self.c[0] * self.r + 3.0 * self.c[1]) * self.r + 2.0 * self.c[2]) * self.r + self.c[3]) * 0.5
        self.thDot += c * (self.cth * self.sth * self.TH + self.L2 * (self.cth / self.sth)**3)

    def stormerVerlet2 (self, y):  # Compose higher orders from this second-order symplectic base
        self.qUp(0.5 * y)
        self.qDotUp(y)
        self.qUp(0.5 * y)

    def stormerVerlet4 (self, y):  # Compose higher orders from this fourth-order symplectic base
        self.qUp(0.5 / self.cbrt2two * y)
        self.qDotUp(1.0 / self.cbrt2two * y)
        self.qUp(0.5 * self.cbrt2one / self.cbrt2two * y)
        self.qDotUp(- self.cbrt2 / self.cbrt2two * y)
        self.qUp(0.5 * self.cbrt2one / self.cbrt2two * y)
        self.qDotUp(1.0 / self.cbrt2two * y)
        self.qUp(0.5 / self.cbrt2two * y)

def main ():  # Need to be inside a function to return . . .
    ic = loads(stdin.read())
    bl = BL(ic['M'], ic['a'], ic['mu'], ic['E'], ic['Lz'], ic['C'], ic['r'], ic['theta'], ic['time'], ic['step'], ic['integratorOrder'])
    bl.refresh(bl.r, bl.th)
    bl.rDot = - sqrt(bl.R if bl.R > 0.0 else 0.0)
    bl.thDot = - sqrt(bl.THETA if bl.THETA > 0.0 else 0.0)
    mino = tau = 0.0
    while not abs(mino) > bl.simtime:
        bl.errors(bl.R, bl.THETA, bl.tDot, bl.rDot, bl.thDot, bl.phDot)
	print >> stdout, '{"mino":%.9e, "tau":%.9e, "v4e":%.1f, "ER":%.1f, "ETh":%.1f, "t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "tDot":%.9e, "rDot":%.9e, "thDot":%.9e, "phDot":%.9e}' % (mino, tau, bl.v4e, bl.eR, bl.eTh, bl.t, bl.r, bl.th, bl.ph, bl.tDot / bl.S, bl.rDot / bl.S, bl.thDot / bl.S, bl.phDot / bl.S)  # Log data,  d/dTau = 1/sigma * d/dLambda !!!
        for i in bl.coefficients:  # Composition happens in this loop
            bl.base(bl.coeff[i] * bl.h)
        mino += bl.h
        tau += bl.h * bl.S  # dTau = sigma * dLambda !!!

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"

