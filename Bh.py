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

    def __init__(self, m, a, q, mu2, e, lz, cc, r0, th0, tol):
        self.rs = 2 * m
        self.a = a
        self.q = q
        self.mu2 = mu2
        self.t = make_mpfr(0)
        self.r = Dual.from_number(r0)
        self.th = Dual.from_number(th0 + acos(make_mpfr(0)))
        self.phi = make_mpfr(0)
        self.p_t = Dual.from_number(- e)
        self.p_r = Dual.from_number(make_mpfr(0))
        self.p_th = (cc - self.th.cos.sqr * (a**2 * (mu2 - e**2) + (lz / self.th.sin).sqr)).sqrt
        self.p_phi = Dual.from_number(lz)
        self.r_prev = self.r.val + 0.001
        self.th_prev = self.th.val + 0.001
        self.p_r_prev = self.p_r.val + 0.001
        self.p_th_prev = self.p_th.val + 0.001
        self.tol = tol
        self.h0 = self.h(self.r, self.th, self.p_t, self.p_r, self.p_th, self.p_phi).val

    def h(self, r, th, p_t, p_r, p_th, p_phi):  # MTW p.900 equation 33.35
        delta = r.sqr - self.rs * r + self.a**2 + self.q
        return D05 * (- ((r.sqr + self.a**2) * p_t + self.a * p_phi).sqr / delta + delta * p_r.sqr + p_th.sqr
                      + ((p_phi + self.a * th.sin.sqr * p_t) / th.sin).sqr) / (r.sqr + self.a**2 * th.cos.sqr)

    def q_t_update(self, h):
        return self.t + D05 * h * self.h(self.r, self.th, self.p_t.var, self.p_r, self.p_th, self.p_phi).der

    def q_r_implicit(self, h, r):
        return self.r.val - r + D05 * h * self.h(Dual.from_number(r), self.th, self.p_t, self.p_r.var, self.p_th, self.p_phi).der

    def q_th_implicit(self, h, th):
        return self.th.val - th + D05 * h * self.h(self.r, Dual.from_number(th), self.p_t, self.p_r, self.p_th.var, self.p_phi).der

    def q_phi_update(self, h):
        return self.phi + D05 * h * self.h(self.r, self.th, self.p_t, self.p_r, self.p_th, self.p_phi.var).der

    def q_update_1(self, h):
        r = secant(self.q_r_implicit, self.r.val, self.r_prev, h, self.tol)
        th = secant(self.q_th_implicit, self.th.val, self.th_prev, h, self.tol)
        self.r_prev = self.r.val
        self.th_prev = self.th.val
        self.t = self.q_t_update(h)
        self.r = Dual.from_number(r)
        self.th = Dual.from_number(th)
        self.phi = self.q_phi_update(h)

    def q_update_2(self, h):
        r = self.r.val + D05 * h * self.h(self.r, self.th, self.p_t, self.p_r.var, self.p_th, self.p_phi).der
        th = self.th.val + D05 * h * self.h(self.r, self.th, self.p_t, self.p_r, self.p_th.var, self.p_phi).der
        self.t = self.q_t_update(h)
        self.r = Dual.from_number(r)
        self.th = Dual.from_number(th)
        self.phi = self.q_phi_update(h)

    def p_r_implicit(self, h, p_r):
        r_var = self.r.var
        return self.p_r.val - p_r - D05 * h * (self.h(r_var, self.th, self.p_t, self.p_r, self.p_th, self.p_phi).der
                                             + self.h(r_var, self.th, self.p_t, Dual.from_number(p_r), self.p_th, self.p_phi).der)

    def p_th_implicit(self, h, p_th):
        th_var = self.th.var
        return self.p_th.val - p_th - D05 * h * (self.h(self.r, th_var, self.p_t, self.p_r, self.p_th, self.p_phi).der
                                               + self.h(self.r, th_var, self.p_t, self.p_r, Dual.from_number(p_th), self.p_phi).der)

    def p_update(self, h):
        pr = secant(self.p_r_implicit, self.p_r.val, self.p_r_prev, h, self.tol)
        pth = secant(self.p_th_implicit, self.p_th.val, self.p_th_prev, h, self.tol)
        self.p_r_prev = self.p_r.val
        self.p_th_prev = self.p_th.val
        self.p_r = Dual.from_number(pr)
        self.p_th = Dual.from_number(pth)

    def stormer_verlet(self, h):
        self.q_update_1(h)
        self.p_update(h)
        self.q_update_2(h)

    def solve(self, h, start, end, tr):
        tau = 0.0
        i = 0
        while tau < end:
            if tau >= start and i % tr == 0:
                self.plot(tau)
            self.stormer_verlet(h)
            i += 1
            tau = h * i
        self.plot(tau)

    def plot(self, tau):
        h = self.h(self.r, self.th, self.p_t, self.p_r, self.p_th, self.p_phi).val
        print('{{"tau":{:.9e},"v4e":{:.9e},"H":{:.9e},"E":{:.9e},"L":{:.9e},"Q":{:.9e},"t":{:.9e},"r":{:.9e},"th":{:.9e},"ph":{:.9e}}}'.format(
            tau, h - self.h0, h, - self.p_t.val, self.p_phi.val,
            (self.p_th.sqr + self.th.cos.sqr * (self.a**2 * (self.mu2 - self.p_t.sqr) + (self.p_phi / self.th.sin).sqr)).val,
            self.t, self.r.val, self.th.val, self.phi))


def secant(f, a, b, h, tol, limit=101):
    f_a = f(h, a)
    f_b = f(h, b)
    counter = c = f_c = 1
    while abs(f_c) > tol:
        if counter == limit:
            raise RuntimeError("{}\n Giving up after {} iterations, value {}, previous {}".format(f, counter - 1, a, b))
        c = (b * f_a - a * f_b) / (f_a - f_b)
        f_c = f(h, c)
        b = a
        f_b = f_a
        a = c
        f_a = f_c
        counter += 1
    return c

if __name__ == "__main__":
    print("Simulator: {}".format(argv[0]), file=stderr)
    input_data = stdin.read()
    ic = loads(input_data, parse_float=mpfr)['IC']
    print(input_data, file=stderr)
    step = ic['step']
    bh = Kerr(ic['M'], ic['a'], ic['q'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'], ic['tol'])
    bh.solve(step, ic['start'], ic['end'], ic['plotratio'])
else:
    print(__name__ + " module loaded", file=stderr)
