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
from math import log10, sqrt, sin
from json import loads

class BL(object):
    def __init__(self, bhMass, spin, pMass2, energy, momentum, carter, r0, thetaMin, starttime, duration, timestep, tratio):
        self.a = spin
        self.mu2 = pMass2
        self.E = energy
        self.L = momentum
        self.Q = carter
        self.r = r0
        self.th = thetaMin
        self.starttime = starttime
        self.endtime = starttime + duration
        self.h = timestep
        self.tr = tratio
        self.kt = [0.0, 0.0, 0.0, 0.0]
        self.kr = [0.0, 0.0, 0.0, 0.0]
        self.kth = [0.0, 0.0, 0.0, 0.0]
        self.kph = [0.0, 0.0, 0.0, 0.0]
        self.sgnR = self.sgnTH = 1.0
        self.t = self.ph = self.v4cum = 0.0
        self.a2 = self.a**2
        self.aE = self.a * self.E
        self.a2E = self.a2 * self.E
        self.L2 = self.L**2
        self.aL = self.a * self.L
        self.L_aE2 = (self.L - self.aE)**2
        self.a2xE2_mu2 = - self.a2 * (self.E**2 - self.mu2)
        self.derivatives(self.r, self.th, 0)

    def derivatives (self, radius, theta, i):
        r2 = radius**2
        self.sth2 = sin(theta)**2
        self.cth2 = 1.0 - self.sth2
        self.ra2 = r2 + self.a2
        self.D = self.ra2 - 2.0 * radius
        self.S = r2 + self.a2 * self.cth2
        self.R = (self.ra2 * self.E - self.aL)**2 - self.D * (self.Q + self.L_aE2 + self.mu2 * r2)
        self.THETA = self.Q - self.cth2 * (self.a2xE2_mu2 + self.L2 / self.sth2)
        P_D = (self.ra2 * self.E - self.aL) / self.D
        self.kt[i] = (self.ra2 * P_D + self.aL - self.a2E * self.sth2) / self.S
        self.kr[i] = sqrt(self.R if self.R >= 0.0 else -self.R) / self.S
        self.kth[i] = sqrt(self.THETA if self.THETA >= 0.0 else -self.THETA) / self.S
        self.kph[i] = (self.a * P_D - self.aE + self.L / self.sth2) / self.S

    def rk4Step (self):
        def update (kx):
            return (kx[0] + 3.0 * (kx[1] + kx[2]) + kx[3]) * self.h / 8.0
        self.sgnR = self.sgnR if self.R > 0.0 else - self.sgnR
        self.sgnTH = self.sgnTH if self.THETA > 0.0 else - self.sgnTH
        self.derivatives(self.r + 1.0 / 3.0 * self.sgnR * self.kr[0] * self.h, self.th + 1.0 / 3.0 * self.sgnTH * self.kth[0] * self.h, 1)
        self.derivatives(self.r + 2.0 / 3.0 * self.sgnR * self.kr[1] * self.h, self.th + 2.0 / 3.0 * self.sgnTH * self.kth[1] * self.h, 2)
        self.derivatives(self.r + self.sgnR * self.kr[2] * self.h, self.th + self.sgnTH * self.kth[2] * self.h, 3)
        self.t += update(self.kt)
        self.r += update(self.kr) * self.sgnR
        self.th += update(self.kth) * self.sgnTH
        self.ph += update(self.kph)
        self.derivatives(self.r, self.th, 0)

    def output (self, tau):
        e = self.mu2 + self.sth2 / self.S * (self.a * self.kt[0] - self.ra2 * self.kph[0])**2 + self.S / self.D * self.kr[0]**2 + self.S * self.kth[0]**2 - self.D / self.S * (self.kt[0] - self.a * self.sth2 * self.kph[0])**2
        e = e if e >= 0.0 else -e
        print >> stdout, '{"tau":%.9e, "v4e":%.1f, "t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "tP":%.9e, "rP":%.9e, "thP":%.9e, "phP":%.9e}' % (tau, 10.0 * log10(e if e > 1.0e-18 else 1.0e-18), self.t, self.r, self.th, self.ph, self.kt[0], self.kr[0], self.kth[0], self.kph[0])  # Log data

def main ():
    ic = loads(stdin.read())
    bl = BL(ic['M'], ic['a'], ic['mu'], ic['E'], ic['Lz'], ic['C'], ic['r'], ic['theta'], ic['start'], ic['duration'], ic['step'], ic['plotratio'])
    count = 0
    tau = 0.0
    while tau <= bl.endtime:
        if tau >= bl.starttime and count % bl.tr == 0:
            bl.output(tau)
        bl.rk4Step()
        count += 1
        tau += bl.h
    bl.output(tau)

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"

