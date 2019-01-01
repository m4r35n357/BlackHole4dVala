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
from gmpy2 import get_context, mpfr, sqrt, acos
get_context().precision = 113  # Set this BEFORE importing any Taylor Series stuff!
from Symplectic import Symplectic, D1, D0


class Newton(object):
    def __init__(self, l_fac, r0):
        self.PI_2 = acos(D0)
        self.ph = self.phDot = D0
        self.r = r0
        self.L = l_fac * sqrt(r0)
        self.L2 = self.L**2
        self.rDot = - sqrt(r0 - self.L2) / self.r
        self.H0 = self.h()

    def h(self):
        return 0.5 * (self.rDot**2 + self.L2 / (self.r**2)) - D1 / self.r

    def q_update(self, c):
        self.r += c * self.rDot
        self.phDot = self.L / self.r**2
        self.ph += c * self.phDot

    def p_update(self, d):
        self.rDot -= d * (D1 - self.L2 / self.r) / self.r**2

    def solve(self, method, h, start, end, tr):
        t = 0.0
        i = plot_count = 0
        while t < end:
            if t >= start and i % tr == 0:
                self.plot(t)
                plot_count += 1
            method()
            i += 1
            t = h * i
        self.plot(t)
        return i, plot_count

    def plot(self, time):
        print('{{"tau":{:.9e},"v4e":{:.9e},"t":{:.9e},"r":{:.9e},"th":{:.9e},"ph":{:.9e},"tP":{:.9e},"rP":{:.9e},"thP":{:.9e},"phP":{:.9e}}}'.format(
            time, self.h() - self.H0, time, self.r, self.PI_2, self.ph, D1, self.rDot, D0, self.phDot))


if __name__ == "__main__":
    print("Simulator: {}".format(argv[0]), file=stderr)
    input_data = stdin.read()
    ic = loads(input_data, parse_float=mpfr)['IC']
    print(input_data, file=stderr)
    bh = Newton(ic['Lfac'], ic['r0'])
    step = ic['step']
    bh.solve(Symplectic(bh, step, ic['integrator'], ic['scheme']).method, step, ic['start'], ic['end'], ic['plotratio'])
else:
    print(__name__ + " module loaded", file=stderr)
