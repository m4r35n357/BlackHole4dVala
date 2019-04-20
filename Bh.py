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
    def from_number(cls, value, variable=False):
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

    def __init__(self, m, a, q, μ2, e, lz, cc, r_0, θ_0, tol):
        self.rs = 2 * m
        self.a = a
        self.q = q
        self.μ2 = μ2
        self.t = make_mpfr(0)
        self.r = Dual.from_number(r_0)
        self.θ = Dual.from_number((make_mpfr(90.0) - θ_0) * acos(make_mpfr(-1)) / make_mpfr(180.0))
        self.φ = make_mpfr(0)
        self.p_t = Dual.from_number(- e)
        self.p_r = Dual.from_number(make_mpfr(0))
        self.p_θ = (cc - self.θ.cos.sqr * (a**2 * (μ2 - e**2) + (lz / self.θ.sin).sqr)).sqrt
        self.p_φ = Dual.from_number(lz)
        self.r_prev = self.r.val + 0.001
        self.θ_prev = self.θ.val + 0.001
        self.p_r_prev = self.p_r.val + 0.001
        self.p_θ_prev = self.p_θ.val + 0.001
        self.ε = tol
        self.h0 = self.h(self.r, self.θ, self.p_t, self.p_r, self.p_θ, self.p_φ).val

    def h(self, r, θ, p_t, p_r, p_θ, p_φ):  # MTW p.900 equation 33.35
        Δ = r.sqr - self.rs * r + self.a**2 + self.q
        return D05 * (- ((r.sqr + self.a**2) * p_t + self.a * p_φ).sqr / Δ + Δ * p_r.sqr + p_θ.sqr
                      + ((p_φ + self.a * θ.sin.sqr * p_t) / θ.sin).sqr) / (r.sqr + self.a**2 * θ.cos.sqr)

    def q_t_update(self, h):
        return self.t + D05 * h * self.h(self.r, self.θ, self.p_t.var, self.p_r, self.p_θ, self.p_φ).der

    def q_r_implicit(self, h, r):
        return self.r.val - r + D05 * h * self.h(Dual.from_number(r), self.θ, self.p_t, self.p_r.var, self.p_θ, self.p_φ).der

    def q_θ_implicit(self, h, θ):
        return self.θ.val - θ + D05 * h * self.h(self.r, Dual.from_number(θ), self.p_t, self.p_r, self.p_θ.var, self.p_φ).der

    def q_φ_update(self, h):
        return self.φ + D05 * h * self.h(self.r, self.θ, self.p_t, self.p_r, self.p_θ, self.p_φ.var).der

    def q_update_1(self, h):
        r = secant(self.q_r_implicit, self.r.val, self.r_prev, h, self.ε)
        θ = secant(self.q_θ_implicit, self.θ.val, self.θ_prev, h, self.ε)
        self.r_prev = self.r.val
        self.θ_prev = self.θ.val
        self.t = self.q_t_update(h)
        self.r = Dual.from_number(r)
        self.θ = Dual.from_number(θ)
        self.φ = self.q_φ_update(h)

    def q_update_2(self, h):
        r = self.r.val + D05 * h * self.h(self.r, self.θ, self.p_t, self.p_r.var, self.p_θ, self.p_φ).der
        θ = self.θ.val + D05 * h * self.h(self.r, self.θ, self.p_t, self.p_r, self.p_θ.var, self.p_φ).der
        self.t = self.q_t_update(h)
        self.r = Dual.from_number(r)
        self.θ = Dual.from_number(θ)
        self.φ = self.q_φ_update(h)

    def p_r_implicit(self, h, p_r):
        r_var = self.r.var
        return self.p_r.val - p_r - D05 * h * (self.h(r_var, self.θ, self.p_t, self.p_r, self.p_θ, self.p_φ).der
                                             + self.h(r_var, self.θ, self.p_t, Dual.from_number(p_r), self.p_θ, self.p_φ).der)

    def p_θ_implicit(self, h, p_θ):
        θ_var = self.θ.var
        return self.p_θ.val - p_θ - D05 * h * (self.h(self.r, θ_var, self.p_t, self.p_r, self.p_θ, self.p_φ).der
                                             + self.h(self.r, θ_var, self.p_t, self.p_r, Dual.from_number(p_θ), self.p_φ).der)

    def p_update(self, h):
        pr = secant(self.p_r_implicit, self.p_r.val, self.p_r_prev, h, self.ε)
        pθ = secant(self.p_θ_implicit, self.p_θ.val, self.p_θ_prev, h, self.ε)
        self.p_r_prev = self.p_r.val
        self.p_θ_prev = self.p_θ.val
        self.p_r = Dual.from_number(pr)
        self.p_θ = Dual.from_number(pθ)

    def stormer_verlet(self, h):
        self.q_update_1(h)
        self.p_update(h)
        self.q_update_2(h)

    def solve(self, h, start, end, tr):
        τ = 0.0
        i = 0
        while τ < end:
            if τ >= start and i % tr == 0:
                self.plot(τ)
            self.stormer_verlet(h)
            i += 1
            τ = h * i
        self.plot(τ)

    def plot(self, τ):
        h = self.h(self.r, self.θ, self.p_t, self.p_r, self.p_θ, self.p_φ).val
        print('{{"tau":{:.9e},"v4e":{:.9e},"H":{:.9e},"E":{:.9e},"L":{:.9e},"Q":{:.9e},"t":{:.9e},"r":{:.9e},"th":{:.9e},"ph":{:.9e}}}'.format(
            τ, h - self.h0, h, - self.p_t.val, self.p_φ.val,
            (self.p_θ.sqr + self.θ.cos.sqr * (self.a**2 * (self.μ2 - self.p_t.sqr) + (self.p_φ / self.θ.sin).sqr)).val,
            self.t, self.r.val, self.θ.val, self.φ))


def secant(f, a, b, h, ε, limit=101):
    f_a = f(h, a)
    f_b = f(h, b)
    counter = delta = c = f_c = 1
    while abs(f_c) > ε or abs(delta) > ε:
        if counter == limit:
            raise RuntimeError("{}\n Giving up after {} iterations, value {}, previous {}".format(f, counter - 1, a, b))
        c = (b * f_a - a * f_b) / (f_a - f_b)
        f_c = f(h, c)
        b = a
        f_b = f_a
        a = c
        f_a = f_c
        delta = b - a
        counter += 1
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
