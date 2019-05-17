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
get_context().precision = 236  # Set this BEFORE importing or defining any Taylor Series / Dual Number stuff!

#  ./Bh.py <initial-conditions.json | ./filegraphics-pi.py initial-conditions.json &
#  ./Bh.py <initial-conditions.json | ./plotErrors.py initial-conditions.json tau 1 &

# noinspection PyArgumentList
D05 = mpfr('0.5')

# noinspection PyArgumentList
def make_mpfr(x):
    return mpfr(str(x)) if isinstance(x, float) else mpfr(x) if isinstance(x, (int, str)) else x


class Dual:

    def __init__(self, value, derivative):
        self.val = value
        self.der = derivative

    @classmethod
    def get(cls, value, variable=False):
        return cls(make_mpfr(value), make_mpfr(1) if variable else make_mpfr(0))

    def __str__(self):
        return "{:+.6e} {:+.6e}".format(self.val, self.der)

    def __abs__(self):
        return self.sqr.sqrt

    def __pos__(self):
        return Dual(self.val, self.der)

    def __neg__(self):
        return Dual(- self.val, - self.der)

    def __add__(self, other):
        if isinstance(other, Dual):
            return Dual(self.val + other.val, self.der + other.der)
        else:
            return Dual(self.val + other, self.der)

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other):
        if isinstance(other, Dual):
            return Dual(self.val - other.val, self.der - other.der)
        else:
            return Dual(self.val - other, self.der)

    def __rsub__(self, other):
        return Dual(- self.val + other, - self.der)

    def __mul__(self, other):
        if isinstance(other, Dual):
            return Dual(self.val * other.val, self.der * other.val + self.val * other.der)
        else:
            return Dual(self.val * other, self.der * other)

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if isinstance(other, Dual):
            return Dual(self.val / other.val, (self.der * other.val - self.val * other.der) / other.val**2)
        else:
            return Dual(self.val / other, self.der / other)

    def __rtruediv__(self, other):
        return Dual(other / self.val, - other * self.der / self.val**2)

    def __pow__(self, a):
        return Dual(self.val**a, a * self.val**(a - 1) * self.der)

    @property
    def var(self):
        return Dual(self.val, make_mpfr(1))

    @property
    def sqr(self):
        return Dual(self.val**2, 2 * self.der * self.val)

    @property
    def sqrt(self):
        sqrt_val = sqrt(self.val)
        return Dual(sqrt_val, self.der / (2 * sqrt_val))

    @property
    def sin(self):
        return Dual(sin(self.val), self.der * cos(self.val))

    @property
    def cos(self):
        return Dual(cos(self.val), - self.der * sin(self.val))


class Kerr(object):

    def __init__(self, m, a, q, μ2, e, lz, cc, r0, θ0, ε):
        self.rs = 2 * m
        self.a = a
        self.q = q
        self.μ2 = μ2
        self.qt = make_mpfr(0)
        self.qr = Dual.get(r0)
        self.qθ = Dual.get((make_mpfr(90.0) - θ0) * acos(make_mpfr(-1)) / make_mpfr(180.0))
        self.qφ = make_mpfr(0)
        self.pt = Dual.get(- e)
        self.pr = Dual.get(make_mpfr(0))
        self.pθ = (cc - self.qθ.cos.sqr * (a**2 * (μ2 - e**2) + (lz / self.qθ.sin).sqr)).sqrt
        self.pφ = Dual.get(lz)
        self.qr_prev = self.qr.val + 0.001
        self.qθ_prev = self.qθ.val + 0.001
        self.pr_prev = self.pr.val + 0.001
        self.pθ_prev = self.pθ.val + 0.001
        self.ε = ε
        self.h0 = self.h(self.qr, self.qθ, self.pt, self.pr, self.pθ, self.pφ).val

    def h(self, qr, qθ, pt, pr, pθ, pφ):  # MTW p.900 equation 33.35
        Δ = qr.sqr - self.rs * qr + self.a**2 + self.q
        return D05 * (- ((qr.sqr + self.a**2) * pt + self.a * pφ).sqr / Δ + Δ * pr.sqr + pθ.sqr
                      + ((pφ + self.a * qθ.sin.sqr * pt) / qθ.sin).sqr) / (qr.sqr + self.a**2 * qθ.cos.sqr)

    def qt_update(self, δτ):
        return self.qt + D05 * δτ * self.h(self.qr, self.qθ, self.pt.var, self.pr, self.pθ, self.pφ).der

    def qr_implicit(self, δτ, r):
        return self.qr.val - r + D05 * δτ * self.h(Dual.get(r), self.qθ, self.pt, self.pr.var, self.pθ, self.pφ).der

    def qθ_implicit(self, δτ, θ):
        return self.qθ.val - θ + D05 * δτ * self.h(self.qr, Dual.get(θ), self.pt, self.pr, self.pθ.var, self.pφ).der

    def qφ_update(self, δτ):
        return self.qφ + D05 * δτ * self.h(self.qr, self.qθ, self.pt, self.pr, self.pθ, self.pφ.var).der

    def q_update_1(self, δτ):
        qr = secant(self.qr_implicit, self.qr.val, self.qr_prev, δτ, self.ε)
        qθ = secant(self.qθ_implicit, self.qθ.val, self.qθ_prev, δτ, self.ε)
        self.qr_prev = self.qr.val
        self.qθ_prev = self.qθ.val
        self.qr = Dual.get(qr)
        self.qθ = Dual.get(qθ)
        self.qt = self.qt_update(δτ)
        self.qφ = self.qφ_update(δτ)

    def q_update_2(self, δτ):
        qr = self.qr.val + D05 * δτ * self.h(self.qr, self.qθ, self.pt, self.pr.var, self.pθ, self.pφ).der
        qθ = self.qθ.val + D05 * δτ * self.h(self.qr, self.qθ, self.pt, self.pr, self.pθ.var, self.pφ).der
        self.qr = Dual.get(qr)
        self.qθ = Dual.get(qθ)
        self.qt = self.qt_update(δτ)
        self.qφ = self.qφ_update(δτ)

    def pr_implicit(self, δτ, pr):
        qr = self.qr.var
        return self.pr.val - pr - D05 * δτ * (self.h(qr, self.qθ, self.pt, self.pr, self.pθ, self.pφ).der
                                            + self.h(qr, self.qθ, self.pt, Dual.get(pr), self.pθ, self.pφ).der)

    def pθ_implicit(self, δτ, pθ):
        qθ = self.qθ.var
        return self.pθ.val - pθ - D05 * δτ * (self.h(self.qr, qθ, self.pt, self.pr, self.pθ, self.pφ).der
                                            + self.h(self.qr, qθ, self.pt, self.pr, Dual.get(pθ), self.pφ).der)

    def p_update(self, δτ):
        pr = secant(self.pr_implicit, self.pr.val, self.pr_prev, δτ, self.ε)
        pθ = secant(self.pθ_implicit, self.pθ.val, self.pθ_prev, δτ, self.ε)
        self.pr_prev = self.pr.val
        self.pθ_prev = self.pθ.val
        self.pr = Dual.get(pr)
        self.pθ = Dual.get(pθ)

    def stormer_verlet(self, δτ):
        self.q_update_1(δτ)
        self.p_update(δτ)
        self.q_update_2(δτ)

    def solve(self, δτ, start, end, tr):
        τ = make_mpfr(0.0)
        i = 0
        while τ < end:
            if τ >= start and i % tr == 0:
                self.plot(τ)
            self.stormer_verlet(δτ)
            i += 1
            τ = δτ * i
        self.plot(τ)

    def plot(self, τ):
        h = self.h(self.qr, self.qθ, self.pt, self.pr, self.pθ, self.pφ).val
        print(f'{{"tau":{τ:.9e},"v4e":{h - self.h0:.9e},"H":{h:.9e},"E":{- self.pt.val:.9e},"L":{self.pφ.val:.9e},'
              f'"Q":{(self.pθ.sqr + self.qθ.cos.sqr * (self.a**2 * (self.μ2 - self.pt.sqr) + (self.pφ / self.qθ.sin).sqr)).val:.9e},'
              f'"t":{self.qt:.9e},"r":{self.qr.val:.9e},"th":{self.qθ.val:.9e},"ph":{self.qφ:.9e}}}')


def secant(f, a, b, h, ε, limit=101):
    f_a = f(h, a)
    f_b = f(h, b)
    count = δx = c = f_c = 1
    while abs(f_c) > ε or abs(δx) > ε:
        if count == limit:
            raise RuntimeError("{}\n Giving up after {} iterations, value {}, previous {}".format(f, count - 1, a, b))
        c = (b * f_a - a * f_b) / (f_a - f_b)
        f_c = f(h, c)
        b = a
        f_b = f_a
        a = c
        f_a = f_c
        δx = b - a
        count += 1
    return c

if __name__ == "__main__":
    print("Simulator: {}".format(argv[0]), file=stderr)
    input_data = open(argv[1]).read() if len(argv) == 2 else stdin.read()
    ic = loads(input_data, parse_float=mpfr)['IC']
    print(input_data, file=stderr)
    step = ic['step']
    bh = Kerr(ic['M'], ic['a'], ic['q'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'], ic['tol'])
    bh.solve(step, ic['start'], ic['end'], ic['plotratio'])
else:
    print(__name__ + " module loaded", file=stderr)
