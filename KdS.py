#!/usr/bin/env pypy
"""
Copyright (c) 2014, 2015, 2016, 2017, Ian Smith (m4r35n357)
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


class Symplectic(object):
    def __init__(self, model, order):
        if order == 'sb1':
            self.step = self.symplectic_euler
        elif order == 'sb2':
            self.step = self.stormer_verlet
        elif order == 'sb4':
            self.cbrt2 = 2.0 ** (1.0 / 3.0)
            self.f2 = 1.0 / (2.0 - self.cbrt2)
            self.coefficients = [0.5 * self.f2, self.f2, 0.5 * (1.0 - self.cbrt2) * self.f2, - self.cbrt2 * self.f2]
            self.step = self.forest_ruth
        else:
            raise Exception('>>> Integrator must be sb2 or sb4, was "{}" <<<'.format(order))
        self.model = model

    def symplectic_euler(self, w):
        self.model.qUp(w)
        self.model.pUp(w)

    def stormer_verlet(self, w):
        self.model.qUp(w * 0.5)
        self.model.pUp(w)
        self.model.qUp(w * 0.5)

    def forest_ruth(self, w):
        self.model.qUp(w * self.coefficients[0])
        self.model.pUp(w * self.coefficients[1])
        self.model.qUp(w * self.coefficients[2])
        self.model.pUp(w * self.coefficients[3])
        self.model.qUp(w * self.coefficients[2])
        self.model.pUp(w * self.coefficients[1])
        self.model.qUp(w * self.coefficients[0])


class BhSymp(object):
    def __init__(self, Lambda, a, mu2, E, L, C, r0, th0, order):
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
        self.sgnR = self.sgnTH = -1
        self.t = self.ph = 0.0
        self.r = r0
        self.th = (90.0 - th0) * pi / 180.0
        self.integrator = Symplectic(self, order)

    def refresh(self):
        self.r2 = self.r**2
        self.sth = sin(self.th)
        self.sth2 = self.sth**2
        self.cth2 = 1.0 - self.sth2
        self.ra2 = self.r2 + self.a2
        self.P = self.ra2 * self.E - self.aL
        self.D_r = (1.0 - self.l_3 * self.r2) * self.ra2 - 2.0 * self.r
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
        return fabs(self.mu2 + self.sth2 * self.D_th / SX2 * (self.a * Ut - self.ra2 * Uph)**2 + self.S / self.D_r * Ur**2
                    + self.S / self.D_th * Uth**2 - self.D_r / SX2 * (Ut - self.a * self.sth2 * Uph)**2)

    @staticmethod
    def log_error(e):
        return 10.0 * log10(e if e > 1.0e-18 else 1.0e-18)

    @staticmethod
    def modH(xdot, x):
        return 0.5 * fabs(xdot**2 - x)

    def qUp(self, h):
        self.t += h * self.Ut
        self.r += h * self.Ur
        self.th += h * self.Uth
        self.ph += h * self.Uph
        self.refresh()

    def pUp(self, h):
        self.Ur += h * (self.r * (2.0 * self.E * self.P * self.X2 - self.mu2 * self.D_r)
                        - (self.r * (1.0 - self.l_3 * self.r2) - self.l_3 * self.r * self.ra2 - 1.0) * (self.K + self.mu2 * self.r2))
        self.Uth += h * (cos(self.th) * (self.sth * self.a2 * (self.mu2 * self.D_th - self.l_3 * (self.K - self.a2mu2 * self.cth2))
                        + self.X2 * self.T / self.sth * (self.T / self.sth2 - 2.0 * self.aE)))

    def plot(self, mino, tau):
        eR = self.log_error(self.modH(self.Ur, self.R))
        eTh = self.log_error(self.modH(self.Uth, self.TH))
        v4e = self.log_error(self.v4_error(self.Ut / self.S, self.Ur / self.S, self.Uth / self.S, self.Uph / self.S))  # d/dTau = 1/sigma * d/dLambda !!!
        print >> stdout, '{"mino":%.9e, "tau":%.9e, "v4e":%.1f, "v4c":%.1f, "ER":%.1f, "ETh":%.1f, "t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "tP":%.9e, "rP":%.9e, "thP":%.9e, "phP":%.9e}' \
                % (mino, tau, v4e, -180.0, eR, eTh, self.t, self.r, self.th, self.ph, self.Ut / self.S, self.Ur / self.S, self.Uth / self.S, self.Uph / self.S)  # Log data,  d/dTau = 1/sigma * d/dLambda !!!

    def solve(self, start, end, h, tr):
        mino = tau = 0.0
        i = plotCount = 0
        self.refresh()
        self.Ur = - sqrt(fabs(self.R))
        self.Uth = - sqrt(fabs(self.TH))
        while tau <= end:
            if tau >= start and i % tr == 0:
                self.plot(mino, tau)
                plotCount += 1
            self.integrator.step(h)
            i += 1
            mino = h * i
            tau += h * self.S  # dTau = sigma * dlambda  - NB lambda is affine parameter here, not the cc !!!
        self.plot(mino, tau)
        return i, plotCount


print >> stderr, "Executable: {}".format(argv[0])
input_data = stdin.read()
ic = loads(input_data)['IC']
print >> stderr, input_data
BhSymp(ic['lambda'], ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'], ic['integrator']).solve(
    ic['start'], ic['end'], ic['step'], ic['plotratio'])
