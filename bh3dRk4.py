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
from math import log10, sqrt, sin, pi
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
        self.sgnR = self.sgnTH = -1.0
        self.t = self.ph = self.v4cum = 0.0
        self.a2 = self.a**2
        self.horizon = 1.0 + sqrt(1.0 - self.a2)
        self.aE = self.a * self.E
        self.a2E = self.a2 * self.E
        self.L2 = self.L**2
        self.aL = self.a * self.L
        E2_mu2 = self.E**2 - self.mu2
        self.a2xE2_mu2 = self.a2 * E2_mu2
        self.c = [E2_mu2, 2.0 * self.mu2, self.a2xE2_mu2 - self.L2 - self.Q, 2.0 * ((self.aE - self.L)**2 + self.Q), - self.a2 * self.Q]
        self.f(self.r, self.th, 0)

    def f (self, radius, theta, stage):
        r2 = radius**2
        self.sth2 = sin(theta)**2
        self.cth2 = 1.0 - self.sth2
        self.ra2 = r2 + self.a2
        self.D = self.ra2 - 2.0 * radius
        self.S = r2 + self.a2 * self.cth2
        self.R = (((self.c[0] * radius + self.c[1]) * radius + self.c[2]) * radius + self.c[3]) * radius + self.c[4]
        self.THETA = self.Q - self.cth2 * (self.L2 / self.sth2 - self.a2xE2_mu2)
        P_D = (self.ra2 * self.E - self.aL) / self.D
        self.tDot = (self.ra2 * P_D + self.aL - self.a2E * self.sth2) / self.S
        self.rDot = sqrt(self.R if self.R >= 0.0 else -self.R) / self.S
        self.thDot = sqrt(self.THETA if self.THETA >= 0.0 else -self.THETA) / self.S
        self.phDot = (self.a * P_D - self.aE + self.L / self.sth2) / self.S
        self.kt[stage] = self.h * self.tDot
        self.kr[stage] = self.h * self.rDot
        self.kth[stage] = self.h * self.thDot
        self.kph[stage] = self.h * self.phDot

    def rk4Step (self):
        def sumK (kx):
            return (kx[0] + 2.0 * (kx[1] + kx[2]) + kx[3]) / 6.0
        self.sgnR = self.sgnR if self.R > 0.0 else - self.sgnR
        self.sgnTH = self.sgnTH if self.THETA > 0.0 else - self.sgnTH
        self.f(self.r + 0.5 * self.kr[0], self.th + 0.5 * self.kth[0], 1)
        self.f(self.r + 0.5 * self.kr[1], self.th + 0.5 * self.kth[1], 2)
        self.f(self.r + self.kr[2], self.th + self.kth[2], 3)
        self.t += sumK(self.kt)
        self.r += sumK(self.kr) * self.sgnR
        self.th += sumK(self.kth) * self.sgnTH
        self.ph += sumK(self.kph)
        self.f(self.r, self.th, 0)

    def output (self, tau):
        e = self.mu2 + self.sth2 / self.S * (self.a * self.tDot - self.ra2 * self.phDot)**2 + self.S / self.D * self.rDot**2 + self.S * self.thDot**2 - self.D / self.S * (self.tDot - self.a * self.sth2 * self.phDot)**2
        e = e if e >= 0.0 else -e
        print >> stdout, '{"tau":%.9e, "v4e":%.1f, "v4c":%.1f, "ER":%.1f, "ETh":%.1f, "t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "tP":%.9e, "rP":%.9e, "thP":%.9e, "phP":%.9e}' % (tau, 10.0 * log10(e if e > 1.0e-18 else 1.0e-18), -180.0, -180.0, -180.0, self.t, self.r, self.th, self.ph, self.tDot, self.rDot, self.thDot, self.phDot)  # Log data

def main ():
    ic = loads(stdin.read())['IC']
    bl = BL(ic['M'], ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], (90.0 - ic['th0']) * pi / 180.0, ic['start'], ic['duration'], ic['step'], ic['plotratio'])
    count = 0
    tau = 0.0
    while tau <= bl.endtime and bl.r >= bl.horizon:
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

