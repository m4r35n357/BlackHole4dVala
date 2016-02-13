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
    def __init__(self, bhMass, spin, pMass2, energy, momentum, carter, r0, thetaMin, starttime, duration, timestep, order):
        self.a = spin
        self.mu2 = pMass2
        self.E = energy
        self.L = momentum
        self.Q = carter
        self.r = r0
        self.th = thetaMin
        self.starttime = abs(starttime)
        self.duration = abs(duration)
        self.endtime = self.starttime + self.duration
        self.h = timestep
        self.cbrt2 = 2.0**(1.0 / 3.0)
        self.f2 = 1.0 / (2.0 - self.cbrt2)
        self.coefficients = array('d', [0.5 * self.f2, self.f2, 0.5 * (1.0 - self.cbrt2) * self.f2, - self.cbrt2 * self.f2])
        if order == 'sb2':  # Second order base
            self.base = self.base2;
            self.w = array('d', [1.0])
        elif order == 'sc4':  # Fourth order, composed from Second order
            self.base = self.base2;
            self.w = array('d', [self.coefficients[1], self.coefficients[3], self.coefficients[1]])
        elif order == 'sb4':  # Fourth order base
            self.base = self.base4;
            self.w = array('d', [1.0])
        elif order == 'sc6':  # Sixth order, composed from Fourth order
            self.base = self.base4;
            fthrt2 = 2.0**(1.0 / 5.0)
            self.w = array('d', [1.0 / (2.0 - fthrt2), - fthrt2 / (2.0 - fthrt2), 1.0 / (2.0 - fthrt2)])
        else:  # Wrong value for integrator order
            raise Exception('>>> ERROR! Integrator order must be sb2, sc4, sb4, or sc6 <<<')
        self.wRange = range(len(self.w))
        self.t = self.ph = self.v4cum = 0.0
        self.count = 0
        self.a2 = self.a**2
        self.aE = self.a * self.E
        self.a2E = self.a2 * self.E
        self.L2 = self.L**2
        self.aL = self.a * self.L
        E2_mu2 = self.E**2 - self.mu2
        self.c = array('d', [E2_mu2, 2.0 * self.mu2, self.a2 * E2_mu2 - self.L2 - self.Q, 2.0 * ((self.aE - self.L)**2 + self.Q), - self.a2 * self.Q])
        self.a2xE2_mu2 = - self.a2 * E2_mu2

    def errors (self, R, THETA, tP, rP, thP, phP):  # Error analysis
        def logError (e):
            return 10.0 * log10(e if e > 1.0e-18 else 1.0e-18)
        def modH (xDot, X):
            return 0.5 * fabs(xDot**2 - X)
        def v4Error (tP, rP, thP, phP):  # norm squared, xDot means dx/dTau !!!
            return fabs(self.mu2 + self.sth2 / self.S * (self.a * tP - self.ra2 * phP)**2 + self.S / self.D * rP**2 + self.S * thP**2 - self.D / self.S * (tP - self.a * self.sth2 * phP)**2)
        self.eR = logError(modH(rP, R))
        self.eTh = logError(modH(thP, THETA))
        error = v4Error(tP / self.S, rP / self.S, thP / self.S, phP / self.S)
        self.v4cum += error
        self.v4c = logError(self.v4cum / self.count)
        self.v4e = logError(error)  # d/dTau = 1/sigma * d/dLambda !!!

    def refresh (self, r, th):  # Update quantities that depend on current values of r or theta
        r2 = r * r
        self.sth = sin(th)
        self.cth = cos(th)
        self.sth2 = self.sth**2
        self.cth2 = 1.0 - self.sth2
        self.ra2 = r2 + self.a2
        self.D = self.ra2 - 2.0 * r
        self.S = r2 + self.a2 * self.cth2
        self.R = (((self.c[0] * r + self.c[1]) * r + self.c[2]) * r + self.c[3]) * r + self.c[4]
        self.TH = self.a2xE2_mu2 + self.L2 / self.sth2
        self.THETA = self.Q - self.cth2 * self.TH
        P_D = (self.ra2 * self.E - self.aL) / self.D
        self.tP = self.ra2 * P_D + self.aL - self.a2E * self.sth2
        self.phP = self.a * P_D - self.aE + self.L / self.sth2

    def qUp (self, d):  # q += d * dq/dTau, where dq/dTau = dH/dp (i.e. dT/dp).  N.B. here q = r or theta; t and phi are just along for the ride . . .
        self.t += d * self.tP
        self.r += d * self.rP
        self.th += d * self.thP
        self.ph += d * self.phP
        self.refresh(self.r, self.th)

    def pUp (self, c):  # p += c * dp/dTau, where dp/dTau = -dH/dq (i.e. dV/dq, minus sign cancels with the one in the pseudo-Hamiltonian)
        self.rP += c * (((4.0 * self.c[0] * self.r + 3.0 * self.c[1]) * self.r + 2.0 * self.c[2]) * self.r + self.c[3]) * 0.5
        self.thP += c * (self.cth * self.sth * self.TH + self.L2 * (self.cth / self.sth)**3)

    def base2 (self, w):  # Compose higher orders from this second-order symplectic base (d2 = 0.0)
        self.pUp(w * 0.5)  # c1 = 0.5
        self.qUp(w)        # d1 = 1.0
        self.pUp(w * 0.5)  # c2 = 0.5

    def base4 (self, w):  # Compose higher orders from this fourth-order symplectic base (d4 = 0.0)
        self.pUp(w * self.coefficients[0])  # w * c1
        self.qUp(w * self.coefficients[1])  # w * d1
        self.pUp(w * self.coefficients[2])  # w * c2
        self.qUp(w * self.coefficients[3])  # w * d2
        self.pUp(w * self.coefficients[2])  # w * c3
        self.qUp(w * self.coefficients[1])  # w * d3
        self.pUp(w * self.coefficients[0])  # w * c4

def main ():  # Need to be inside a function to return . . .
    ic = loads(stdin.read())
    bl = BL(ic['M'], ic['a'], ic['mu'], ic['E'], ic['Lz'], ic['C'], ic['r'], ic['theta'], ic['start'], ic['duration'], ic['step'], ic['integrator'])
    bl.refresh(bl.r, bl.th)
    bl.rP = sqrt(fabs(bl.R))
    bl.thP = sqrt(fabs(bl.THETA))
    mino = tau = 0.0
    while not abs(mino) > bl.endtime:
        bl.count += 1
        bl.errors(bl.R, bl.THETA, bl.tP, bl.rP, bl.thP, bl.phP)
        if abs(mino) > bl.starttime:
            print >> stdout, '{"mino":%.9e, "tau":%.9e, "v4e":%.1f, "v4c":%.1f, "ER":%.1f, "ETh":%.1f, "t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "tP":%.9e, "rP":%.9e, "thP":%.9e, "phP":%.9e}' % (mino, tau, bl.v4e, bl.v4c, bl.eR, bl.eTh, bl.t, bl.r, bl.th, bl.ph, bl.tP / bl.S, bl.rP / bl.S, bl.thP / bl.S, bl.phP / bl.S)  # Log data,  d/dTau = 1/sigma * d/dLambda !!!
        for i in bl.wRange:  # Composition happens in this loop
            bl.base(bl.w[i] * bl.h)
        mino += bl.h
        tau += bl.h * bl.S  # dTau = sigma * dLambda !!!

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"

