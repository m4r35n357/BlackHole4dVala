#!/usr/bin/env python3
"""
Copyright (c) 2014-2018, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from json import loads
from sys import stdin, stderr, argv
from gmpy2 import get_context, mpfr, sqrt, sin, cos, acos
get_context().precision = 236  # Set this BEFORE importing any Taylor Series stuff!
from Symplectic import Symplectic, D1, D2, D0, D3
from taylor import to_mpfr


class BhSymp(object):
    def __init__(self, cc, a, mu2, e, l, q, r0, th0, xh):
        self.l_3 = cc / D3
        self.a = a
        self.mu2 = mu2
        self.E = e
        self.L = l
        self.a2 = a**2
        self.a2l_3 = self.a2 * self.l_3
        self.a2mu2 = self.a2 * mu2
        self.aE = a * e
        self.aL = a * l
        self.X2 = (D1 + self.a2l_3)**2
        self.two_EX2 = D2 * e * self.X2
        self.two_aE = D2 * a * e
        self.K = q + self.X2 * (l - self.aE)**2
        self.t = self.ph = D0
        self.r = r0
        self.th = (to_mpfr(90.0) - th0) * acos(to_mpfr(-1)) / to_mpfr(180.0)
        self.cross = xh
        self.refresh()
        self.Ur = - sqrt(self.r_potential if self.r_potential >= D0 else -self.r_potential)
        self.Uth = - sqrt(self.th_potential if self.th_potential >= D0 else -self.th_potential)

    def refresh(self):
        self.r2 = self.r**2
        self.sth = sin(self.th)
        self.cth = cos(self.th)
        self.sth2 = self.sth**2
        self.cth2 = D1 - self.sth2
        self.ra2 = self.r2 + self.a2
        self.P = self.ra2 * self.E - self.aL
        self.d_r = (D1 - self.l_3 * self.r2) * self.ra2 - D2 * self.r
        self.r_potential = self.X2 * self.P**2 - self.d_r * (self.mu2 * self.r2 + self.K)
        self.T = self.aE * self.sth2 - self.L
        self.d_th = D1 + self.a2l_3 * self.cth2
        self.th_potential = self.d_th * (self.K - self.a2mu2 * self.cth2) - self.X2 * self.T**2 / self.sth2
        p_dr = self.P / self.d_r
        t_dth = self.T / self.d_th
        self.S = self.r2 + self.a2 * self.cth2
        self.Ut = (p_dr * self.ra2 - t_dth * self.a) * self.X2
        self.Uph = (p_dr * self.a - t_dth / self.sth2) * self.X2

    def v4_error(self, ut, ur, uth, uph):
        sx2 = self.S * self.X2
        return self.mu2 + self.sth2 * self.d_th / sx2 * (self.a * ut - self.ra2 * uph)**2 + self.S / self.d_r * ur**2 \
               + self.S / self.d_th * uth**2 - self.d_r / sx2 * (ut - self.a * self.sth2 * uph)**2

    def q_update(self, c):
        self.t += c * self.Ut
        self.r += c * self.Ur
        self.th += c * self.Uth
        self.ph += c * self.Uph
        self.refresh()

    def p_update(self, d):
        self.Ur += d * (self.r * (self.two_EX2 * self.P - self.mu2 * self.d_r) - (self.r * (D1 - self.l_3 * (self.r2 + self.ra2)) - D1) * (self.K + self.mu2 * self.r2))
        self.Uth += d * (self.cth * (self.sth * self.a2 * (self.mu2 * self.d_th - self.l_3 * (self.K - self.a2mu2 * self.cth2)) + self.X2 * self.T / self.sth * (self.T / self.sth2 - self.two_aE)))

    def solve(self, method, h, start, end, tr):
        mino = tau = 0.0
        i = plot_count = 0
        while (tau < end) and (self.cross or self.d_r > D0):
            if tau >= start and i % tr == 0:
                self.plot(mino, tau, self.Ut / self.S, self.Ur / self.S, self.Uth / self.S, self.Uph / self.S)
                plot_count += 1
            method()
            i += 1
            mino = h * i
            tau += h * self.S
        self.plot(mino, tau, self.Ut / self.S, self.Ur / self.S, self.Uth / self.S, self.Uph / self.S)
        return i, plot_count

    def plot(self, mino, tau, ut, ur, uth, uph):
        er = ur**2 - self.r_potential / self.S**2
        eth = uth**2 - self.th_potential / self.S**2
        v4e = self.v4_error(ut, ur, uth, uph)
        print('{{"mino":{:.9e},"tau":{:.9e},"v4e":{:.9e},"ER":{:.9e},"ETh":{:.9e},"t":{:.9e},"r":{:.9e},"th":{:.9e},"ph":{:.9e},"tP":{:.9e},"rP":{:.9e},"thP":{:.9e},"phP":{:.9e}}}'.format(
            mino, tau, v4e, er, eth, self.t, self.r, self.th, self.ph, ut, ur, uth, uph))


if __name__ == "__main__":
    print("Simulator: {}".format(argv[0]), file=stderr)
    input_data = stdin.read()
    ic = loads(input_data, parse_float=mpfr)['IC']
    print(input_data, file=stderr)
    bh = BhSymp(ic['lambda'], ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'], ic['cross'])
    step = ic['step']
    bh.solve(Symplectic(bh, step, ic['integrator'], ic['scheme']).method, step, ic['start'], ic['end'], ic['plotratio'])
else:
    print(__name__ + " module loaded", file=stderr)
