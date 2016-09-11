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
    def __init__(self, Lambda, spin, pMass2, energy, momentum, carter, r0, thetaMin, starttime, duration, timestep, tratio):
        self.l_3 = Lambda / 3.0
        self.a = spin
        self.mu2 = pMass2
        self.E = energy
        self.L = momentum
        self.a2 = self.a**2
        self.a2l_3 = self.a2 * self.l_3
        self.a2mu2 = self.a2 * self.mu2
        self.horizon = 1.0 + sqrt(1.0 - self.a2)
        self.aE = self.a * self.E
        self.aL = self.a * self.L
        self.X2 = (1.0 + self.a2l_3)**2
        self.K = carter + self.X2 * (self.L - self.aE)**2
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
        self.r = r0
        self.th = (90.0 - thetaMin) * pi / 180.0
        self.f(self.r, self.th, 0)

    def f (self, radius, theta, stage):
        r2 = radius**2
        self.sth2 = sin(theta)**2
        cth2 = 1.0 - self.sth2
        self.ra2 = r2 + self.a2
        P = self.ra2 * self.E - self.aL
        self.D_r = (1.0 - self.l_3 * r2) * self.ra2 - 2.0 * radius
        self.R = self.X2 * P**2 - self.D_r * (self.mu2 * r2 + self.K)
        T = self.aE * self.sth2 - self.L
        self.D_th = 1.0 + self.a2l_3 * cth2
        self.TH = self.D_th * (self.K - self.a2mu2 * cth2) - self.X2 * T**2 / self.sth2
        P_Dr = P / self.D_r
        T_Dth = T / self.D_th
        self.S = r2 + self.a2 * cth2
        X2_S = self.X2 / self.S
        self.tDot = (P_Dr * self.ra2 - T_Dth * self.a) * X2_S
        self.rDot = sqrt(self.R if self.R >= 0.0 else -self.R) / self.S
        self.thDot = sqrt(self.TH if self.TH >= 0.0 else -self.TH) / self.S
        self.phDot = (P_Dr * self.a - T_Dth / self.sth2) * X2_S
        self.kt[stage] = self.h * self.tDot
        self.kr[stage] = self.h * self.rDot
        self.kth[stage] = self.h * self.thDot
        self.kph[stage] = self.h * self.phDot

    def rk4Step (self):
        def sumK (kx):
            return (kx[0] + 2.0 * (kx[1] + kx[2]) + kx[3]) / 6.0
        self.sgnR = self.sgnR if self.R > 0.0 else - self.sgnR
        self.sgnTH = self.sgnTH if self.TH > 0.0 else - self.sgnTH
        self.f(self.r + 0.5 * self.kr[0], self.th + 0.5 * self.kth[0], 1)
        self.f(self.r + 0.5 * self.kr[1], self.th + 0.5 * self.kth[1], 2)
        self.f(self.r + self.kr[2], self.th + self.kth[2], 3)
        self.t += sumK(self.kt)
        self.r += sumK(self.kr) * self.sgnR
        self.th += sumK(self.kth) * self.sgnTH
        self.ph += sumK(self.kph)
        self.f(self.r, self.th, 0)

    def output (self, tau):
        e = self.mu2 + self.sth2 * self.D_th / (self.S * self.X2) * (self.a * self.tDot - self.ra2 * self.phDot)**2 \
                     + self.S / self.D_r * self.rDot**2 + self.S / self.D_th * self.thDot**2 \
                     - self.D_r / (self.S * self.X2) * (self.tDot - self.a * self.sth2 * self.phDot)**2
        e = e if e >= 0.0 else -e
        print >> stdout, '{"tau":%.9e, "v4e":%.1f, "D_r":%.1f, "D_th":%.1f, "S":%.1f,' % \
                          (tau, 10.0 * log10(e if e > 1.0e-18 else 1.0e-18), self.D_r, self.D_th, self.S),
        print >> stdout, '"t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e,' % (self.t, self.r, self.th, self.ph),
        print >> stdout, '"tP":%.9e, "rP":%.9e, "thP":%.9e, "phP":%.9e}' % (self.tDot, self.rDot, self.thDot, self.phDot)  # Log data

def main ():
    ic = loads(stdin.read())['IC']
    bl = BL(ic['lambda'], ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'], ic['start'], ic['duration'], ic['step'], ic['plotratio'])
    count = 0
    tau = 0.0
    while tau <= bl.endtime and bl.r >= bl.horizon and bl.D_r >= 0.0:
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

