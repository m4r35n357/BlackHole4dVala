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

#  ./DoublePendulum.py <initial-conditions.double-pendulum.json | ../../c/ODE-Playground/plotPi2d.py
#  ../../c/ODE-Playground/plotXY.py </tmp/data-Python 1 4 5 &

from json import loads
from sys import stdin, stderr, argv
from gmpy2 import get_context, mpfr, sin, cos, log10
get_context().precision = 236  # Set this BEFORE importing any Taylor Series stuff!
from dual import Dual

class DoublePendulum(object):
    def __init__(self, g, l1, m1, l2, m2, th1_0, pth1_0, th2_0, pth2_0, tol):
        self.g = g
        self.l1 = l1
        self.m1 = m1
        self.l2 = l2
        self.m2 = m2
        self.th1 = Dual.from_number(th1_0)
        self.pth1 = Dual.from_number(pth1_0)
        self.th2 = Dual.from_number(th2_0)
        self.pth2 = Dual.from_number(pth2_0)
        self.h0 = self.h(self.th1, self.pth1, self.th2, self.pth2).val
        self.tol = tol

    def h(self, th1, pth1, th2, pth2):
        return (self.l2**2 * self.m2 * pth1.sqr + self.l1**2 * (self.m1 + self.m2) * pth2.sqr
                - 2 * self.m2 * self.l1 * self.l2 * pth1 * pth2 * (th1 - th2).cos) \
               / (2 * self.l1**2 * self.l2**2 * self.m2 * (self.m1 + self.m2 * (th1 - th2).sin.sqr)) \
               - (self.m1 + self.m2) * self.g * self.l1 * th1.cos - self.m2 * self.g * self.l2 * th2.cos

    def qth1_update(self, c, th1):
        return self.th1.val - th1 + 0.5 * c * self.h(Dual.from_number(th1), self.pth1.var, self.th2, self.pth2).der

    def qth2_update(self, c, th2):
        return self.th2.val - th2 + 0.5 * c * self.h(self.th1, self.pth1, Dual.from_number(th2), self.pth2.var).der

    def q_update_1(self, c):
        th1 = secant(self.qth1_update, self.th1.val, c, self.tol)
        th2 = secant(self.qth2_update, self.th2.val, c, self.tol)
        self.th1 = Dual.from_number(th1)
        self.th2 = Dual.from_number(th2)

    def q_update_2(self, c):
        th1 = self.th1.val + 0.5 * c * self.h(self.th1, self.pth1.var, self.th2, self.pth2).der
        th2 = self.th2.val + 0.5 * c * self.h(self.th1, self.pth1, self.th2, self.pth2.var).der
        self.th1 = Dual.from_number(th1)
        self.th2 = Dual.from_number(th2)

    def pth1_update(self, d, pth1):
        th1_var = self.th1.var
        return self.pth1.val - pth1 - 0.5 * d * (self.h(th1_var, self.pth1, self.th2, self.pth2).der
                                               + self.h(th1_var, Dual.from_number(pth1), self.th2, self.pth2).der)

    def pth2_update(self, d, pth2):
        th2_var = self.th2.var
        return self.pth2.val - pth2 - 0.5 * d * (self.h(self.th1, self.pth1, th2_var, self.pth2).der
                                               + self.h(self.th1, self.pth1, th2_var, Dual.from_number(pth2)).der)

    def p_update(self, d):
        pth1 = secant(self.pth1_update, self.pth1.val, d, self.tol)
        pth2 = secant(self.pth2_update, self.pth2.val, d, self.tol)
        self.pth1 = Dual.from_number(pth1)
        self.pth2 = Dual.from_number(pth2)

    def stormer_verlet(self, h):
        self.q_update_1(h)
        self.p_update(h)
        self.q_update_2(h)

    def yoshida_4(self, h):
        z1 = 1.0 / (2.0 - 2.0**(1.0 / 3.0))
        self.stormer_verlet(h * z1)
        self.stormer_verlet(h * (1.0 - 2.0 * z1))
        self.stormer_verlet(h * z1)

    def suzuki_4(self, h):
        z1 = 1.0 / (4.0 - 4.0**(1.0 / 3.0))
        self.stormer_verlet(h * z1)
        self.stormer_verlet(h * z1)
        self.stormer_verlet(h * (1.0 - 4.0 * z1))
        self.stormer_verlet(h * z1)
        self.stormer_verlet(h * z1)

    def rk4(self, h):
        k1 = [0.0, 0.0, 0.0, 0.0]
        k2 = [0.0, 0.0, 0.0, 0.0]
        k3 = [0.0, 0.0, 0.0, 0.0]
        k4 = [0.0, 0.0, 0.0, 0.0]

        k1[0] = self.h(self.th1.var, self.pth1, self.th2, self.pth2).der
        k2[0] = self.h(self.pth1, self.pth1.var, self.th2, self.pth2).der
        k3[0] = self.h(self.th1, self.pth1, self.th2.var, self.pth2).der
        k4[0] = self.h(self.th1, self.pth1, self.th2, self.pth2.var).der

        k1[1] = self.h(
            (self.th1 + 0.5 * k2[0]).var, self.pth1 + 0.5 * k1[0], self.th2 + 0.5 * k4[0], self.pth2 + 0.5 * k3[0]).der
        k2[1] = self.h(
            self.th1 + 0.5 * k2[0], (self.pth1 + 0.5 * k1[0]).var, self.th2 + 0.5 * k4[0], self.pth2 + 0.5 * k3[0]).der
        k3[1] = self.h(
            self.th1 + 0.5 * k2[0], self.pth1 + 0.5 * k1[0], (self.th2 + 0.5 * k4[0]).var, self.pth2 + 0.5 * k3[0]).der
        k4[1] = self.h(
            self.th1 + 0.5 * k2[0], self.pth1 + 0.5 * k1[0], self.th2 + 0.5 * k4[0], (self.pth2 + 0.5 * k3[0]).var).der

        k1[2] = self.h(
            (self.th1 + 0.5 * k2[1]).var, self.pth1 + 0.5 * k1[1], self.th2 + 0.5 * k4[1], self.pth2 + 0.5 * k3[1]).der
        k2[2] = self.h(
            self.th1 + 0.5 * k2[1], (self.pth1 + 0.5 * k1[1]).var,  self.th2 + 0.5 * k4[1], self.pth2 + 0.5 * k3[1]).der
        k3[2] = self.h(
            self.th1 + 0.5 * k2[1], self.pth1 + 0.5 * k1[1], (self.th2 + 0.5 * k4[1]).var, self.pth2 + 0.5 * k3[1]).der
        k4[2] = self.h(
            self.th1 + 0.5 * k2[1], self.pth1 + 0.5 * k1[1], self.th2 + 0.5 * k4[1], (self.pth2 + 0.5 * k3[1]).var).der

        k1[3] = self.h(
            (self.th1 + k2[2]).var, self.pth1 + k1[2], self.th2 + k4[2], self.pth2 + k3[2]).der
        k2[3] = self.h(
            self.th1 + k2[2], (self.pth1 + k1[2]).var, self.th2 + k4[2], self.pth2 + k3[2]).der
        k3[3] = self.h(
            self.th1 + k2[2], self.pth1 + k1[2], (self.th2 + k4[2]).var, self.pth2 + k3[2]).der
        k4[3] = self.h(
            self.th1 + k2[2], self.pth1 + k1[2], self.th2 + k4[2], (self.pth2 + k3[2]).var).der

        self.th1 = Dual.from_number(self.th1.val + h * (k2[0] + 2 * (k2[1] + k2[2]) + k2[3]) / 6)
        self.pth1 = Dual.from_number(self.pth1.val - h * (k1[0] + 2 * (k1[1] + k1[2]) + k1[3]) / 6)
        self.th2 = Dual.from_number(self.th2.val + h * (k4[0] + 2 * (k4[1] + k4[2]) + k4[3]) / 6)
        self.pth2 = Dual.from_number(self.pth2.val - h * (k3[0] + 2 * (k3[1] + k3[2]) + k3[3]) / 6)

    def euler(self, h):
        self.th1 = Dual.from_number(self.th1.val + h * self.h(self.pth1, self.pth1.var, self.th2, self.pth2).der)
        self.pth1 = Dual.from_number(self.pth1.val - h * self.h(self.th1.var, self.pth1, self.th2, self.pth2).der)
        self.th2 = Dual.from_number(self.th2.val + h * self.h(self.th1, self.pth1, self.th2, self.pth2.var).der)
        self.pth2 = Dual.from_number(self.pth2.val - h * self.h(self.th1, self.pth1, self.th2.var, self.pth2).der)

    def solve(self, h, start, end, tr):
        t = 0.0
        i = 0
        while t < end:
            if t >= start and i % tr == 0:
                self.plot(t)
            # method()
            # self.euler(h)
            # self.rk4(h)
            self.stormer_verlet(h)
            i += 1
            t = h * i
        self.plot(t)

    def plot(self, time):
        x1 = self.l1 * sin(self.th1.val)
        y1 = - self.l1 * cos(self.th1.val)
        x2 = x1 + self.l2 * sin(self.th2.val)
        y2 = y1 - self.l2 * cos(self.th2.val)
        error = abs(self.h(self.th1, self.pth1, self.th2, self.pth2).val - self.h0)
        print("{:+.9e} {:+.9e} {:+.9e} {:+.9e} {:+.5e} {:+.9e}".format(
            x1, y1, x2, y2, time, 10 * log10(error if error > 1.0e-18 else 1.0e-18)))


def false_position(func, parameter, cd, tol, spread=0.001):
    a = parameter * (1 - spread) + (-1 if parameter * -1 > 0 else 1) * spread
    b = parameter * (1 + spread) + (1 if parameter * -1 < 0 else -1) * spread
    f_a = func(cd, a)
    f_b = func(cd, b)
    counter = c = f_c = 1
    while abs(f_c) > tol:
        c = a - f_a * (b - a) / (f_b - f_a)
        f_c = func(cd, c)
        if f_a * f_c < 0.0:
            b = c
            f_b = f_c
        else:
            a = c
            f_a = f_c
        counter += 1
        if counter == 1000:
            raise RuntimeError("Giving Up!")
    return c


def secant(func, parameter, cd, tol, spread=0.001):
    a = parameter * (1 - spread) + (-1 if parameter * -1 > 0 else 1) * spread
    b = parameter * (1 + spread) + (1 if parameter * -1 < 0 else -1) * spread
    f_a = func(cd, a)
    f_b = func(cd, b)
    counter = c = f_c = 1
    while abs(f_c) > tol:
        c = (b * f_a - a * f_b) / (f_a - f_b)
        f_c = func(cd, c)
        b = a
        f_b = f_a
        a = c
        f_a = f_c
        counter += 1
        if counter == 1000:
            raise RuntimeError("Giving Up!")
    return c


if __name__ == "__main__":
    print("Simulator: {}".format(argv[0]), file=stderr)
    input_data = stdin.read()
    ic = loads(input_data, parse_float=mpfr)['IC']
    print(input_data, file=stderr)
    dp = DoublePendulum(ic['g'], ic['l1'], ic['m1'], ic['l2'], ic['m2'], ic['th1'], ic['pth1'], ic['th2'], ic['pth2'], ic['tol'])
    step = ic['step']
    dp.solve(step, ic['start'], ic['end'], ic['plotratio'])
else:
    print(__name__ + " module loaded", file=stderr)
