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
    def __init__(self, Lambda, spin, pMass2, energy, momentum, carter, r0, thetaMin, starttime, duration, timestep, tratio, integrator):
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
        self.tau = self.t = self.ph = self.v4cum = 0.0
        self.r = r0
        self.th = (90.0 - thetaMin) * pi / 180.0
        if integrator == 'rk4':
            self.evaluator = self.evaluator4
            self.updater = self.updater4
        elif integrator == 'rk438':
            self.evaluator = self.evaluator438
            self.updater = self.updater438
        else:
            print >> stderr("Bad integrator type, valid choices are: [ rk4 | rk438 ]")
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
        self.Ut = (P_Dr * self.ra2 - T_Dth * self.a) * X2_S
        self.Ur = sqrt(self.R if self.R >= 0.0 else -self.R) / self.S
        self.Uth = sqrt(self.TH if self.TH >= 0.0 else -self.TH) / self.S
        self.Uph = (P_Dr * self.a - T_Dth / self.sth2) * X2_S
        self.kt[stage] = self.h * self.Ut
        self.kr[stage] = self.h * self.Ur
        self.kth[stage] = self.h * self.Uth
        self.kph[stage] = self.h * self.Uph

    @staticmethod
    def updater4 (kx):
        return (kx[0] + 2.0 * (kx[1] + kx[2]) + kx[3]) / 6.0

    @staticmethod
    def updater438 (kx):
        return (kx[0] + 3.0 * (kx[1] + kx[2]) + kx[3]) / 8.0

    def evaluator4 (self):
        self.f(self.r + 0.5 * self.kr[0], self.th + 0.5 * self.kth[0], 1)
        self.f(self.r + 0.5 * self.kr[1], self.th + 0.5 * self.kth[1], 2)
        self.f(self.r + self.kr[2], self.th + self.kth[2], 3)

    def evaluator438 (self):
        self.f(self.r + 1.0 / 3.0 * self.kr[0], self.th + 1.0 / 3.0 * self.kth[0], 1)
        self.f(self.r - 1.0 / 3.0 * self.kr[0] + self.kr[1], self.th - 1.0 / 3.0 * self.kth[0] + self.kth[1], 2)
        self.f(self.r + self.kr[0] - self.kr[1] + self.kr[2], self.th + self.kth[0] - self.kth[1] + self.kth[2], 3)

    def iterate (self):
        self.sgnR = self.sgnR if self.R > 0.0 else - self.sgnR
        self.sgnTH = self.sgnTH if self.TH > 0.0 else - self.sgnTH
        self.evaluator()
        self.t += self.updater(self.kt)
        self.r += self.updater(self.kr) * self.sgnR
        self.th += self.updater(self.kth) * self.sgnTH
        self.ph += self.updater(self.kph)
        self.f(self.r, self.th, 0)

    def output (self):
        SX2 = self.S * self.X2
        e = self.mu2 + self.sth2 * self.D_th / SX2 * (self.a * self.Ut - self.ra2 * self.Uph)**2 \
                        + self.S / self.D_r * self.Ur**2 + self.S / self.D_th * self.Uth**2 \
                        - self.D_r / SX2 * (self.Ut - self.a * self.sth2 * self.Uph)**2
        e = e if e >= 0.0 else -e
        print >> stdout, '{"tau":%.9e, "v4e":%.1f, "D_r":%.1f, "D_th":%.1f, "S":%.1f,' \
                         % (self.tau, 10.0 * log10(e if e > 1.0e-18 else 1.0e-18), self.D_r, self.D_th, self.S),
        print >> stdout, '"t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "tP":%.9e, "rP":%.9e, "thP":%.9e, "phP":%.9e}' \
                         % (self.t, self.r, self.th, self.ph, self.Ut, self.Ur, self.Uth, self.Uph)  # Log data

def main ():
    ic = loads(stdin.read())['IC']
    bl = BL(ic['lambda'], ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'], \
            ic['start'], ic['duration'], ic['step'], ic['plotratio'], ic['integrator'])
    count = 0
    while bl.tau <= bl.endtime and bl.r >= bl.horizon and bl.D_r >= 0.0:
        if bl.tau >= bl.starttime and count % bl.tr == 0:
            bl.output()
        bl.iterate()
        count += 1
        bl.tau += bl.h
    bl.output()

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"

