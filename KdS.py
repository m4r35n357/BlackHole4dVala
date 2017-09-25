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
from math import sin, pi, cos, sqrt
from sys import stdin, stderr, argv


class Symplectic(object):
    def __init__(self, model, h, order):
        self.model = model
        self.h = h
        self.x1 = 1.0 / (4.0 - 4.0**(1.0 / 7.0))
        self.y1 = 1.0 / (4.0 - 4.0**(1.0 / 5.0))
        self.z1 = 1.0 / (4.0 - 4.0**(1.0 / 3.0))
        self.x3 = 1.0 - 4.0 * self.x1
        self.y3 = 1.0 - 4.0 * self.y1
        self.z3 = 1.0 - 4.0 * self.z1
        if order == 'sb1':
            print >> stderr, "Python First order"
            self.method = self.first_order
        elif order == 'sb2':
            print >> stderr, "Python Second order"
            self.method = self.second_order
        elif order == 'sb4':
            print >> stderr, "Python Fourth order"
            self.method = self.fourth_order
        elif order == 'sb6':
            print >> stderr, "Python Sixth order"
            self.method = self.sixth_order
        elif order == 'sb8':
            print >> stderr, "Python Eightth order"
            self.method = self.eightth_order
        else:
            raise Exception('>>> Integrator must be sb1, sb2 or sb4, was "{}" <<<'.format(order))

    def first_order(self):
        self.model.qUpdate(self.h)
        self.model.pUpdate(self.h)

    def base2(self, s):
        self.model.qUpdate(self.h * s * 0.5)
        self.model.pUpdate(self.h * s)
        self.model.qUpdate(self.h * s * 0.5)

    def second_order(self):
        self.base2(1.0)

    def base4(self, s):
        self.base2(s * self.z1)
        self.base2(s * self.z1)
        self.base2(s * self.z3)
        self.base2(s * self.z1)
        self.base2(s * self.z1)

    def fourth_order(self):
        self.base4(1.0)

    def base6(self, s):
        self.base4(s * self.y1)
        self.base4(s * self.y1)
        self.base4(s * self.y3)
        self.base4(s * self.y1)
        self.base4(s * self.y1)

    def sixth_order(self):
        self.base6(1.0)

    def base8(self, s):
        self.base6(s * self.x1)
        self.base6(s * self.x1)
        self.base6(s * self.x3)
        self.base6(s * self.x1)
        self.base6(s * self.x1)

    def eightth_order(self):
        self.base8(1.0)


class BhSymp(object):
    def __init__(self, Lambda, a, mu2, E, L, C, r0, th0, xh):
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
        self.two_EX2 = 2.0 * E * self.X2
        self.two_aE = 2.0 * a * E
        self.K = C + self.X2 * (L - self.aE)**2
        self.t = self.ph = 0.0
        self.r = r0
        self.th = (90.0 - th0) * pi / 180.0
        self.cross = xh

    def refresh(self):
        self.r2 = self.r**2
        self.sth = sin(self.th)
        self.cth = cos(self.th)
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

    def v4_error(self, Ut, Ur, Uth, Uph):
        SX2 = self.S * self.X2
        return self.mu2 + self.sth2 * self.D_th / SX2 * (self.a * Ut - self.ra2 * Uph)**2 + self.S / self.D_r * Ur**2 \
               + self.S / self.D_th * Uth**2 - self.D_r / SX2 * (Ut - self.a * self.sth2 * Uph)**2

    def qUpdate(self, c):
        self.t += c * self.Ut
        self.r += c * self.Ur
        self.th += c * self.Uth
        self.ph += c * self.Uph
        self.refresh()

    def pUpdate(self, d):
        self.Ur += d * (self.r * (self.two_EX2 * self.P - self.mu2 * self.D_r) - (self.r * (1.0 - self.l_3 * (self.r2 + self.ra2)) - 1.0) * (self.K + self.mu2 * self.r2))
        self.Uth += d * (self.cth * (self.sth * self.a2 * (self.mu2 * self.D_th - self.l_3 * (self.K - self.a2mu2 * self.cth2)) + self.X2 * self.T / self.sth * (self.T / self.sth2 - self.two_aE)))

    def solve(self, method, h, start, end, tr):
        mino = tau = 0.0
        i = plotCount = 0
        self.refresh()
        self.Ur = - sqrt(self.R if self.R >= 0.0 else -self.R)
        self.Uth = - sqrt(self.TH if self.TH >= 0.0 else -self.TH)
        while (tau < end) and (self.cross or self.D_r > 0.0):
            if tau >= start and i % tr == 0:
                self.plot(mino, tau, self.Ut / self.S, self.Ur / self.S, self.Uth / self.S, self.Uph / self.S)
                plotCount += 1
            method()
            i += 1
            mino = h * i
            tau += h * self.S
        self.plot(mino, tau, self.Ut / self.S, self.Ur / self.S, self.Uth / self.S, self.Uph / self.S)
        return i, plotCount

    def plot(self, mino, tau, Ut ,Ur, Uth, Uph):
        eR = Ur**2 - self.R / self.S**2
        eTh = Uth**2 - self.TH / self.S**2
        v4e = self.v4_error(Ut, Ur, Uth, Uph)
        print '{{"mino":{:.9e},"tau":{:.9e},"v4e":{:.9e},"ER":{:.9e},"ETh":{:.9e},"t":{:.9e},"r":{:.9e},"th":{:.9e},"ph":{:.9e},"tP":{:.9e},"rP":{:.9e},"thP":{:.9e},"phP":{:.9e}}}'.format(
            mino, tau, v4e, eR, eTh, self.t, self.r, self.th, self.ph, Ut, Ur, Uth, Uph)


if __name__ == "__main__":
    print >> stderr, "Simulator: {}".format(argv[0])
    input_data = stdin.read()
    ic = loads(input_data)['IC']
    print >> stderr, input_data
    bh = BhSymp(ic['lambda'], ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'], ic['cross'])
    integrator = Symplectic(bh, ic['step'], ic['integrator']).method
    bh.solve(integrator, ic['step'], ic['start'], ic['end'], ic['plotratio'])
else:
    print >> stderr, __name__ + " module loaded"
