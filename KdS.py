#!/usr/bin/env pypy
"""
Copyright (c) 2014, 2015, 2016, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
from json import loads
from math import sqrt, sin, pi, fabs, log10, cos
from sys import stdin, stderr, stdout, argv


class KdSBase(object):
    def __init__(self, Lambda, a, mu2, E, L, C, r0, thetaMin, start, end, timestep, tratio):
        self.l_3 = Lambda / 3.0
        self.a = a
        self.mu2 = mu2
        self.E = E
        self.L = L
        self.a2 = a**2
        self.a2l_3 = self.a2 * self.l_3
        self.a2mu2 = self.a2 * mu2
        self.aE = a * E
        self.aL = a * L
        self.X2 = (1.0 + self.a2l_3)**2
        self.K = C + self.X2 * (L - self.aE)**2
        self.start_time = start
        self.end_time = end
        self.h = timestep
        self.tr = tratio
        self.sgnR = self.sgnTH = -1
        self.t = self.ph = 0.0
        self.r = r0
        self.th = (90.0 - thetaMin) * pi / 180.0

    def refresh(self, radius, theta):
        self.r2 = radius**2
        self.sth = sin(theta)
        self.cth = cos(theta)
        self.sth2 = self.sth**2
        self.cth2 = 1.0 - self.sth2
        self.ra2 = self.r2 + self.a2
        self.P = self.ra2 * self.E - self.aL
        self.D_r = (1.0 - self.l_3 * self.r2) * self.ra2 - 2.0 * radius
        self.R = self.X2 * self.P**2 - self.D_r * (self.mu2 * self.r2 + self.K)
        self.T = self.aE * self.sth2 - self.L
        self.D_th = 1.0 + self.a2l_3 * self.cth2
        self.TH = self.D_th * (self.K - self.a2mu2 * self.cth2) - self.X2 * self.T**2 / self.sth2
        P_Dr = self.P / self.D_r
        T_Dth = self.T / self.D_th
        self.S = self.r2 + self.a2 * self.cth2
        self.Ut = (P_Dr * self.ra2 - T_Dth * self.a) * self.X2
        self.Uph = (P_Dr * self.a - T_Dth / self.sth2) * self.X2

    def v4_error(self, Ut, Ur, Uth, Uph):  # norm squared, xDot means dx/dTau !!!
        SX2 = self.S * self.X2
        return fabs(self.mu2 + self.sth2 * self.D_th / SX2 * (self.a * Ut - self.ra2 * Uph)**2
                    + self.S / self.D_r * Ur**2 + self.S / self.D_th * Uth**2
                    - self.D_r / SX2 * (Ut - self.a * self.sth2 * Uph)**2)

    @staticmethod
    def log_error(e):
        return 10.0 * log10(e if e > 1.0e-18 else 1.0e-18)


class BhRk4(KdSBase):
    def __init__(self, Lambda, spin, pMass2, energy, momentum, carter, r0, thetaMin, starttime, duration, timestep, tratio, integrator):
        super(BhRk4, self).__init__(Lambda, spin, pMass2, energy, momentum, carter, r0, thetaMin, starttime, duration, timestep, tratio)
        self.kt = [0.0, 0.0, 0.0, 0.0]
        self.kr = [0.0, 0.0, 0.0, 0.0]
        self.kth = [0.0, 0.0, 0.0, 0.0]
        self.kph = [0.0, 0.0, 0.0, 0.0]
        self.sgnR = self.sgnTH = -1.0
        if integrator == 'rk4':
            self.evaluator = self.evaluator4
            self.updater = self.updater4
        elif integrator == 'rk438':
            self.evaluator = self.evaluator438
            self.updater = self.updater438
        else:
            raise Exception('>>> ERROR! Integrator order must be rk4 or rk438, was {} <<<'.format(integrator))
        self.max_iterations = round(duration / self.h)
        self.f(self.r, self.th, 0)

    def f(self, radius, theta, stage):
        self.refresh(radius, theta)
        self.Ut /= self.S
        self.Ur = sqrt(self.R if self.R >= 0.0 else -self.R) / self.S
        self.Uth = sqrt(self.TH if self.TH >= 0.0 else -self.TH) / self.S
        self.Uph /= self.S
        self.kt[stage] = self.h * self.Ut
        self.kr[stage] = self.h * self.Ur
        self.kth[stage] = self.h * self.Uth
        self.kph[stage] = self.h * self.Uph

    @staticmethod
    def updater4(kx):
        return (kx[0] + 2.0 * (kx[1] + kx[2]) + kx[3]) / 6.0

    @staticmethod
    def updater438(kx):
        return (kx[0] + 3.0 * (kx[1] + kx[2]) + kx[3]) / 8.0

    def evaluator4(self):
        self.f(self.r + 0.5 * self.kr[0], self.th + 0.5 * self.kth[0], 1)
        self.f(self.r + 0.5 * self.kr[1], self.th + 0.5 * self.kth[1], 2)
        self.f(self.r + self.kr[2], self.th + self.kth[2], 3)

    def evaluator438(self):
        self.f(self.r + 1.0 / 3.0 * self.kr[0], self.th + 1.0 / 3.0 * self.kth[0], 1)
        self.f(self.r - 1.0 / 3.0 * self.kr[0] + self.kr[1], self.th - 1.0 / 3.0 * self.kth[0] + self.kth[1], 2)
        self.f(self.r + self.kr[0] - self.kr[1] + self.kr[2], self.th + self.kth[0] - self.kth[1] + self.kth[2], 3)

    def iterate(self):
        self.sgnR = self.sgnR if self.R > 0.0 else - self.sgnR
        self.sgnTH = self.sgnTH if self.TH > 0.0 else - self.sgnTH
        self.evaluator()
        self.t += self.updater(self.kt)
        self.r += self.updater(self.kr) * self.sgnR
        self.th += self.updater(self.kth) * self.sgnTH
        self.ph += self.updater(self.kph)
        self.f(self.r, self.th, 0)

    def solve(self):
        tau = 0.0
        iterationCount = plotCount = 0
        while tau < self.end_time:
            if tau >= self.start_time and iterationCount % self.tr == 0:
                self.plot(tau)
                plotCount += 1
            self.iterate()
            iterationCount += 1
            tau = iterationCount * self.h
        self.plot(tau)
        return iterationCount, plotCount

    def plot(self, tau):
        print >> stdout, '{"tau":%.9e, "v4e":%.1f, "D_r":%.9e, "D_th":%.9e, "S":%.9e,' \
                         % (tau, self.log_error(self.v4_error(self.Ut, self.Ur, self.Uth, self.Uph)), self.D_r, self.D_th, self.S),
        print >> stdout, '"t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "tP":%.9e, "rP":%.9e, "thP":%.9e, "phP":%.9e}' \
                         % (self.t, self.r, self.th, self.ph, self.Ut, self.Ur, self.Uth, self.Uph)  # Log data


class BhSymp(KdSBase):
    def __init__(self, Lambda, spin, pMass2, energy, momentum, carter, r0, thetaMin, starttime, duration, timestep, tratio, integrator):
        super(BhSymp, self).__init__(Lambda, spin, pMass2, energy, momentum, carter, r0, thetaMin, starttime, duration, timestep, tratio)
        self.integrator = Symplectic(self, integrator)
        self.refresh(self.r, self.th)
        self.Ur = - sqrt(fabs(self.R))
        self.Uth = - sqrt(fabs(self.TH))

    @staticmethod
    def modH(xdot, x):
        return 0.5 * fabs(xdot**2 - x)

    def qUp(self, d):
        self.t += d * self.h * self.Ut
        self.r += d * self.h * self.Ur
        self.th += d * self.h * self.Uth
        self.ph += d * self.h * self.Uph
        self.refresh(self.r, self.th)

    def pUp(self, c):
        self.Ur += c * self.h * (2.0 * self.r * self.E * self.P * self.X2
                                 - (self.r * (1.0 - self.l_3 * self.r2) - 1.0 - self.l_3 * self.r * self.ra2) * (self.K + self.mu2 * self.r2)
                                 - self.mu2 * self.r * self.D_r)
        self.Uth += c * self.h * (self.cth * self.sth * self.a2 * (self.mu2 * self.D_th - self.l_3 * (self.K - self.a2mu2 * self.cth2))
                                  + self.cth * self.X2 * self.T / self.sth * (self.T / self.sth2 - 2.0 * self.aE))

    def plot(self, mino, tau):
        eR = self.log_error(self.modH(self.Ur, self.R))
        eTh = self.log_error(self.modH(self.Uth, self.TH))
        v4e = self.log_error(self.v4_error(self.Ut / self.S, self.Ur / self.S, self.Uth / self.S, self.Uph / self.S))  # d/dTau = 1/sigma * d/dLambda !!!
        print >> stdout, '{"mino":%.9e, "tau":%.9e, "v4e":%.1f, "v4c":%.1f, "ER":%.1f, "ETh":%.1f, "t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "tP":%.9e, "rP":%.9e, "thP":%.9e, "phP":%.9e}' \
                % (mino, tau, v4e, -180.0, eR, eTh, self.t, self.r, self.th, self.ph, self.Ut / self.S, self.Ur / self.S, self.Uth / self.S, self.Uph / self.S)  # Log data,  d/dTau = 1/sigma * d/dLambda !!!

    def solve(self):
        mino = tau = 0.0
        iterationCount = plotCount = 0
        while tau <= self.end_time:
            if tau >= self.start_time and iterationCount % self.tr == 0:
                self.plot(mino, tau)
                plotCount += 1
            self.integrator.compose()
            iterationCount += 1
            mino = iterationCount * self.h
            tau += self.h * self.S  # dTau = sigma * dlambda  - NB lambda is affine parameter here, not the cc !!!
        self.plot(mino, tau)
        return iterationCount, plotCount


class Symplectic(object):
    def __init__(self, model, order):
        self.cbrt2 = 2.0**(1.0 / 3.0)
        self.f2 = 1.0 / (2.0 - self.cbrt2)
        self.coefficients = [0.5 * self.f2, self.f2, 0.5 * (1.0 - self.cbrt2) * self.f2, - self.cbrt2 * self.f2]
        if order == 'sb2':  # Second order base
            self.base = self.base2
            self.w = [1.0]
        elif order == 'sc4':  # Fourth order, composed from Second order
            self.base = self.base2
            self.w = [self.coefficients[1], self.coefficients[3], self.coefficients[1]]
        elif order == 'sb4':  # Fourth order base
            self.base = self.base4
            self.w = [1.0]
        elif order == 'sc6':  # Sixth order, composed from Fourth order
            self.base = self.base4
            fthrt2 = 2.0**(1.0 / 5.0)
            self.w = [1.0 / (2.0 - fthrt2), - fthrt2 / (2.0 - fthrt2), 1.0 / (2.0 - fthrt2)]
        else:
            raise Exception('>>> ERROR! Integrator order must be sb2, sc4, sb4, or sc6, was "{}" <<<'.format(order))
        self.wRange = range(len(self.w))
        self.model = model

    def base2(self, w):  # Compose higher orders from this second-order symplectic base (d2 = 0.0)
        self.model.pUp(w * 0.5)  # c1 = 0.5
        self.model.qUp(w)        # d1 = 1.0
        self.model.pUp(w * 0.5)  # c2 = 0.5

    def base4(self, w):  # Compose higher orders from this fourth-order symplectic base (d4 = 0.0)
        self.model.pUp(w * self.coefficients[0])  # w * c1
        self.model.qUp(w * self.coefficients[1])  # w * d1
        self.model.pUp(w * self.coefficients[2])  # w * c2
        self.model.qUp(w * self.coefficients[3])  # w * d2
        self.model.pUp(w * self.coefficients[2])  # w * c3
        self.model.qUp(w * self.coefficients[1])  # w * d3
        self.model.pUp(w * self.coefficients[0])  # w * c4

    def compose(self):
        for i in self.wRange:  # Composition happens in this loop
            self.base(self.w[i])


def main():
    print >> stderr, "Executable: {}".format(argv[0])
    input_data = stdin.read()
    ic = loads(input_data)['IC']
    print >> stderr, input_data
    if not ic.get("integrator") or "rk4" in ic['integrator']:
        BhRk4(ic['lambda'], ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'],
                ic['start'], ic['end'], ic['step'], ic['plotratio'], ic['integrator']).solve()
    else:
        BhSymp(ic['lambda'], ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'],
                ic['start'], ic['end'], ic['step'], ic['plotratio'], ic['integrator']).solve()

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"
