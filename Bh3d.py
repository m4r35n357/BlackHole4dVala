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
from gmpy2 import get_context, mpfr, acos, sqrt
get_context().precision = 113  # Set this BEFORE importing any Taylor Series stuff!
from Symplectic import Symplectic, D1, D2, D0
from dual import Dual, make_mpfr


class BhSymp(object):
    def __init__(self, a, μ2, e, lz, cc, r_0, θ_0, xh):
        self.a = a
        self.μ2 = μ2
        self.E = e
        self.L = lz
        self.a2 = a**2
        self.a2μ2 = self.a2 * μ2
        self.aE = a * e
        self.aL = a * lz
        self.K = cc + (lz - self.aE)**2
        self.t = make_mpfr(0)
        self.r = Dual.from_number(r_0, variable=True)
        self.θ = Dual.from_number((make_mpfr(90.0) - θ_0) * acos(make_mpfr(-1)) / make_mpfr(180.0), variable=True)
        self.φ = make_mpfr(0)
        self.cross = xh
        self.refresh()
        self.u_r = - sqrt(self.R.val if self.R.val >= D0 else - self.R.val)
        self.u_θ = - sqrt(self.Θ.val if self.Θ.val >= D0 else - self.Θ.val)

    def refresh(self):
        r2 = self.r.sqr
        self.ra2 = r2 + self.a2
        P = self.ra2 * self.E - self.aL
        self.Δ = self.ra2 - D2 * self.r
        self.R = P.sqr - self.Δ * (self.μ2 * r2 + self.K)
        self.sin2θ = self.θ.sin.sqr
        cos2θ = D1 - self.sin2θ
        T = self.aE * self.sin2θ - self.L
        self.Θ = self.K - self.a2μ2 * cos2θ - T.sqr / self.sin2θ
        P_Δ = P.val / self.Δ.val
        self.Σ = r2.val + self.a2 * cos2θ.val
        self.u_t = P_Δ * self.ra2.val - T.val * self.a
        self.u_φ = P_Δ * self.a - T.val / self.sin2θ.val

    def p4_error(self, u_t, u_r, u_θ, u_φ):
        return (self.μ2 + self.sin2θ.val / self.Σ * (self.a * u_t - self.ra2.val * u_φ)**2
                + self.Σ / self.Δ.val * u_r**2 + self.Σ * u_θ**2
                - self.Δ.val / self.Σ * (u_t - self.a * self.sin2θ.val * u_φ)**2)

    def q_update(self, c):
        self.t += c * self.u_t
        self.r.val += c * self.u_r
        self.θ.val += c * self.u_θ
        self.φ += c * self.u_φ
        self.refresh()

    def p_update(self, d):
        self.u_r += 0.5 * d * self.R.der
        self.u_θ += 0.5 * d * self.Θ.der

    def solve(self, method, h, start, end, tr):
        mino = τ = 0.0
        i = 0
        while (τ < end) and (self.cross or self.Δ.val > D0):
            if τ >= start and i % tr == 0:
                self.plot(mino, τ, self.u_t / self.Σ, self.u_r / self.Σ, self.u_θ / self.Σ, self.u_φ / self.Σ)
            method()
            i += 1
            mino = h * i
            τ += h * self.Σ
        self.plot(mino, τ, self.u_t / self.Σ, self.u_r / self.Σ, self.u_θ / self.Σ, self.u_φ / self.Σ)

    def plot(self, mino, τ, u_t, u_r, u_θ, u_φ):
        er = u_r**2 - self.R.val / self.Σ**2
        eθ = u_θ**2 - self.Θ.val / self.Σ**2
        v4e = self.p4_error(u_t, u_r, u_θ, u_φ)
        print('{{"mino":{:.9e},"tau":{:.9e},"v4e":{:.9e},"ER":{:.9e},"ETh":{:.9e},"t":{:.9e},"r":{:.9e},"th":{:.9e},"ph":{:.9e}}}'.format(
            mino, τ, v4e, er, eθ, self.t, self.r.val, self.θ.val, self.φ))


if __name__ == "__main__":
    print("Simulator: {}".format(argv[0]), file=stderr)
    input_data = stdin.read()
    ic = loads(input_data, parse_float=mpfr)['IC']
    print(input_data, file=stderr)
    bh = BhSymp(ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'], ic['cross'])
    step = ic['step']
    bh.solve(Symplectic(bh, step, ic['integrator'], ic['scheme']).method, step, ic['start'], ic['end'], ic['plotratio'])
else:
    print(__name__ + " module loaded", file=stderr)
