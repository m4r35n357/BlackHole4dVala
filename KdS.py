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


class KdSBase(object):
    def __init__(self, Lambda, a, mu2, E, L, C, r0, thetaMin):
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
    def __init__(self, Lambda, spin, pMass2, energy, momentum, carter, r0, thetaMin, integrator):
        super(BhRk4, self).__init__(Lambda, spin, pMass2, energy, momentum, carter, r0, thetaMin)
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

    def f(self, radius, theta, h, stage):
        self.refresh(radius, theta)
        self.Ut /= self.S
        self.Ur = sqrt(self.R if self.R >= 0.0 else -self.R) / self.S
        self.Uth = sqrt(self.TH if self.TH >= 0.0 else -self.TH) / self.S
        self.Uph /= self.S
        self.kt[stage] = h * self.Ut
        self.kr[stage] = h * self.Ur
        self.kth[stage] = h * self.Uth
        self.kph[stage] = h * self.Uph

    @staticmethod
    def updater4(kx):
        return (kx[0] + 2.0 * (kx[1] + kx[2]) + kx[3]) / 6.0

    @staticmethod
    def updater438(kx):
        return (kx[0] + 3.0 * (kx[1] + kx[2]) + kx[3]) / 8.0

    def evaluator4(self, h):
        self.f(self.r + 0.5 * self.kr[0], self.th + 0.5 * self.kth[0], h, 1)
        self.f(self.r + 0.5 * self.kr[1], self.th + 0.5 * self.kth[1], h, 2)
        self.f(self.r + self.kr[2], self.th + self.kth[2], h, 3)

    def evaluator438(self, h):
        self.f(self.r + self.kr[0] / 3.0, self.th + self.kth[0] / 3.0, h, 1)
        self.f(self.r - self.kr[0] / 3.0 + self.kr[1], self.th - self.kth[0] / 3.0 + self.kth[1], h, 2)
        self.f(self.r + self.kr[0] - self.kr[1] + self.kr[2], self.th + self.kth[0] - self.kth[1] + self.kth[2], h, 3)

    def iterate(self, h):
        self.sgnR = self.sgnR if self.R > 0.0 else - self.sgnR
        self.sgnTH = self.sgnTH if self.TH > 0.0 else - self.sgnTH
        self.evaluator(h)
        self.t += self.updater(self.kt)
        self.r += self.updater(self.kr) * self.sgnR
        self.th += self.updater(self.kth) * self.sgnTH
        self.ph += self.updater(self.kph)
        self.f(self.r, self.th, h, 0)

    def solve(self, start, end, h, tr):
        tau = 0.0
        i = plotCount = 0
        self.f(self.r, self.th, h, 0)
        while tau < end:
            if tau >= start and i % tr == 0:
                self.plot(tau)
                plotCount += 1
            self.iterate(h)
            i += 1
            tau = i * h
        self.plot(tau)
        return i, plotCount

    def plot(self, tau):
        print >> stdout, '{"tau":%.9e, "v4e":%.1f, "D_r":%.9e, "D_th":%.9e, "S":%.9e,' \
                         % (tau, self.log_error(self.v4_error(self.Ut, self.Ur, self.Uth, self.Uph)), self.D_r, self.D_th, self.S),
        print >> stdout, '"t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "tP":%.9e, "rP":%.9e, "thP":%.9e, "phP":%.9e}' \
                         % (self.t, self.r, self.th, self.ph, self.Ut, self.Ur, self.Uth, self.Uph)  # Log data


class BhSymp(KdSBase):
    def __init__(self, Lambda, spin, pMass2, energy, momentum, carter, r0, thetaMin, integrator):
        super(BhSymp, self).__init__(Lambda, spin, pMass2, energy, momentum, carter, r0, thetaMin)
        self.integrator = Symplectic(self, integrator)

    @staticmethod
    def modH(xdot, x):
        return 0.5 * fabs(xdot**2 - x)

    def qUp(self, h):
        self.t += h * self.Ut
        self.r += h * self.Ur
        self.th += h * self.Uth
        self.ph += h * self.Uph
        self.refresh(self.r, self.th)

    def pUp(self, h):
        self.Ur += h * (self.r * (2.0 * self.E * self.P * self.X2 - self.mu2 * self.D_r)
                        - (self.r * (1.0 - self.l_3 * self.r2) - self.l_3 * self.r * self.ra2 - 1.0) * (self.K + self.mu2 * self.r2))
        self.Uth += h * (self.cth * (self.sth * self.a2 * (self.mu2 * self.D_th - self.l_3 * (self.K - self.a2mu2 * self.cth2))
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
        self.refresh(self.r, self.th)
        self.Ur = - sqrt(fabs(self.R))
        self.Uth = - sqrt(fabs(self.TH))
        while tau <= end:
            if tau >= start and i % tr == 0:
                self.plot(mino, tau)
                plotCount += 1
            self.integrator.compose(h)
            i += 1
            mino = i * h
            tau += h * self.S  # dTau = sigma * dlambda  - NB lambda is affine parameter here, not the cc !!!
        self.plot(mino, tau)
        return i, plotCount


class Symplectic(object):
    def __init__(self, model, order):
        self.cbrt2 = 2.0**(1.0 / 3.0)
        self.f2 = 1.0 / (2.0 - self.cbrt2)
        self.coefficients = [0.5 * self.f2, self.f2, 0.5 * (1.0 - self.cbrt2) * self.f2, - self.cbrt2 * self.f2]
        if order == 'sb1':  # Second order base
            self.base = self.base1
            self.w = [1.0]
        elif order == 'sb2':  # Second order base
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
        elif order == 'sh6':  # Sixth order, composed from Second order
            self.base = self.base2
            self.w = [
                0.78451361047755726381949763,
                0.23557321335935813368479318,
                -1.17767998417887100694641568,
                1.31518632068391121888424973,
                -1.17767998417887100694641568,
                0.23557321335935813368479318,
                0.78451361047755726381949763
            ]
        elif order == 'sh8':  # Eighth order, composed from Second order
            self.base = self.base2
            self.w = [
                0.74167036435061295344822780,
                -0.40910082580003159399730010,
                0.19075471029623837995387626,
                -0.57386247111608226665638773,
                0.29906418130365592384446354,
                0.33462491824529818378495798,
                0.31529309239676659663205666,
                -0.79688793935291635401978884,
                0.31529309239676659663205666,
                0.33462491824529818378495798,
                0.29906418130365592384446354,
                -0.57386247111608226665638773,
                0.19075471029623837995387626,
                -0.40910082580003159399730010,
                0.74167036435061295344822780
            ]
        elif order == 'sh10':  # Tenth order, composed from Second order
            self.base = self.base2
            self.w = [
                0.09040619368607278492161150,
                0.53591815953030120213784983,
                0.35123257547493978187517736,
                -0.31116802097815835426086544,
                -0.52556314194263510431065549,
                0.14447909410225247647345695,
                0.02983588609748235818064083,
                0.17786179923739805133592238,
                0.09826906939341637652532377,
                0.46179986210411860873242126,
                -0.33377845599881851314531820,
                0.07095684836524793621031152,
                0.23666960070126868771909819,
                -0.49725977950660985445028388,
                -0.30399616617237257346546356,
                0.05246957188100069574521612,
                0.44373380805019087955111365,
                0.05246957188100069574521612,
                -0.30399616617237257346546356,
                -0.49725977950660985445028388,
                0.23666960070126868771909819,
                0.07095684836524793621031152,
                -0.33377845599881851314531820,
                0.46179986210411860873242126,
                0.09826906939341637652532377,
                0.17786179923739805133592238,
                0.02983588609748235818064083,
                0.14447909410225247647345695,
                -0.52556314194263510431065549,
                -0.31116802097815835426086544,
                0.35123257547493978187517736,
                0.53591815953030120213784983,
                0.09040619368607278492161150
            ]
        else:
            raise Exception('>>> Integrator must be sb2, sc4, sb4, sc6, sh6, sh8, or sh10, was "{}" <<<'.format(order))
        self.wRange = range(len(self.w))
        self.model = model

    def base1(self, w):  # Symplectic Euler
        self.model.qUp(w)
        self.model.pUp(w)

    def base2(self, w):  # Compose higher orders from this second-order symplectic base (d2 = 0.0)
        self.model.qUp(w * 0.5)  # c1 = 0.5
        self.model.pUp(w)        # d1 = 1.0
        self.model.qUp(w * 0.5)  # c2 = 0.5

    def base4(self, w):  # Compose higher orders from this fourth-order symplectic base (d4 = 0.0)
        self.model.qUp(w * self.coefficients[0])  # w * c1
        self.model.pUp(w * self.coefficients[1])  # w * d1
        self.model.qUp(w * self.coefficients[2])  # w * c2
        self.model.pUp(w * self.coefficients[3])  # w * d2
        self.model.qUp(w * self.coefficients[2])  # w * c3
        self.model.pUp(w * self.coefficients[1])  # w * d3
        self.model.qUp(w * self.coefficients[0])  # w * c4

    def compose(self, h):
        for i in self.wRange:  # Composition happens in this loop
            self.base(self.w[i] * h)


def main():
    print >> stderr, "Executable: {}".format(argv[0])
    input_data = stdin.read()
    ic = loads(input_data)['IC']
    print >> stderr, input_data
    if not ic.get("integrator") or "rk4" in ic['integrator']:
        BhRk4(ic['lambda'], ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'],
                ic['integrator']).solve(ic['start'], ic['end'], ic['step'], ic['plotratio'])
    else:
        BhSymp(ic['lambda'], ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'],
                ic['integrator']).solve(ic['start'], ic['end'], ic['step'], ic['plotratio'], )

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"
