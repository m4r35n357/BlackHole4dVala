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
    def __init__(self, a, μ2, e, lz, cc, r0, θ0, xh):
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
        self.r = Dual.get(r0, variable=True)
        self.θ = Dual.get((make_mpfr(90) - θ0) * acos(make_mpfr(-1)) / make_mpfr(180), variable=True)
        self.φ = make_mpfr(0)
        self.cross = xh
        self.refresh()
        self.ur = - sqrt(self.R.val if self.R.val >= D0 else - self.R.val)
        self.uθ = - sqrt(self.Θ.val if self.Θ.val >= D0 else - self.Θ.val)

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
        self.ut = P_Δ * self.ra2.val - T.val * self.a
        self.uφ = P_Δ * self.a - T.val / self.sin2θ.val

    def p4_error(self, ut, ur, uθ, uφ):
        return (self.μ2 + self.sin2θ.val / self.Σ * (self.a * ut - self.ra2.val * uφ)**2 + self.Σ / self.Δ.val * ur**2
                + self.Σ * uθ**2 - self.Δ.val / self.Σ * (ut - self.a * self.sin2θ.val * uφ)**2)

    def q_update(self, c):
        self.t += c * self.ut
        self.r.val += c * self.ur
        self.θ.val += c * self.uθ
        self.φ += c * self.uφ
        self.refresh()

    def p_update(self, d):
        self.ur += 0.5 * d * self.R.der
        self.uθ += 0.5 * d * self.Θ.der

    def solve(self, method, h, start, end, tr):
        mino = τ = 0.0
        i = 0
        while (τ < end) and (self.cross or self.Δ.val > D0):
            if τ >= start and i % tr == 0:
                self.plot(mino, τ, self.ut / self.Σ, self.ur / self.Σ, self.uθ / self.Σ, self.uφ / self.Σ)
            method()
            i += 1
            mino = h * i
            τ += h * self.Σ
        self.plot(mino, τ, self.ut / self.Σ, self.ur / self.Σ, self.uθ / self.Σ, self.uφ / self.Σ)

    def plot(self, mino, τ, ut, ur, uθ, uφ):
        print(f'{{"mino":{mino:.9e},"tau":{τ:.9e},"v4e":{self.p4_error(ut, ur, uθ, uφ):.9e},'
              f'"ER":{ur**2 - self.R.val / self.Σ**2:.9e},"ETh":{uθ**2 - self.Θ.val / self.Σ**2:.9e},'
              f'"t":{self.t:.9e},"r":{self.r.val:.9e},"th":{self.θ.val:.9e},"ph":{self.φ:.9e}}}')


if __name__ == "__main__":
    #  Example: ./Bh3d.py initial-conditions.json  | ./filegraphics-pi.py initial-conditions.json
    print("Simulator: {}".format(argv[0]), file=stderr)
    input_data = open(argv[1]).read() if len(argv) == 2 else stdin.read()
    ic = loads(input_data, parse_float=mpfr)['IC']
    print(input_data, file=stderr)
    bh = BhSymp(ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'], ic['cross'])
    step = ic['step']
    bh.solve(Symplectic(bh, step, ic['integrator'], ic['scheme']).method, step, ic['start'], ic['end'], ic['plotratio'])
else:
    print(__name__ + " module loaded", file=stderr)
