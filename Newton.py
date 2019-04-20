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
from Symplectic import Symplectic
from dual import Dual, make_mpfr


class Newton(object):
    def __init__(self, g, m, l_fac, r0):
        zero = make_mpfr(0)
        self.π_2 = acos(zero)
        self.g = Dual.from_number(g)
        self.m = Dual.from_number(m)
        self.q_φ = Dual.from_number(zero)
        self.p_φ = Dual.from_number(l_fac * m * sqrt(r0))
        self.q_r = Dual.from_number(r0)
        self.p_r = Dual.from_number(zero)
        self.h0 = self.h(self.q_r, self.p_r, self.p_φ).val

    def h(self, q_r, p_r, p_φ):  # ph absent from Hamiltonian
        return (p_r**2 + p_φ**2 / q_r**2) / (2 * self.m) - self.g * self.m / q_r

    def q_update(self, c):
        q_r = c * self.h(self.q_r, self.p_r.var, self.p_φ).der
        q_φ = c * self.h(self.q_r, self.p_r, self.p_φ.var).der
        self.q_r = Dual.from_number(self.q_r.val + q_r)  # only update after all coordinates done!
        self.q_φ = Dual.from_number(self.q_φ.val + q_φ)

    def p_update(self, d):  # no self.p_ph update because ph absent from Hamiltonian
        self.p_r = Dual.from_number(self.p_r.val - d * self.h(self.q_r.var, self.p_r, self.p_φ).der)

    def solve(self, method, h, start, end, tr):
        t = 0.0
        i = 0
        while t < end:
            if t >= start and i % tr == 0:
                self.plot(t)
            method()
            i += 1
            t = h * i
        self.plot(t)

    def plot(self, time):
        print('{{"tau":{:.9e},"v4e":{:.9e},"t":{:.9e},"r":{:.9e},"th":{:.9e},"ph":{:.9e}}}'.format(
            time, self.h(self.q_r, self.p_r, self.p_φ).val - self.h0, time, self.q_r.val, self.π_2, self.q_φ.val))


if __name__ == "__main__":
    # ./Newton.py <initial-conditions.newton.json | ./filegraphics-pi.py initial-conditions.newton.json
    # ./Newton.py <initial-conditions.newton.json | ./plotErrors.py initial-conditions.newton.json t 1
    print("Simulator: {}".format(argv[0]), file=stderr)
    input_data = open(argv[1]).read() if len(argv) == 2 else stdin.read()
    ic = loads(input_data, parse_float=mpfr)['IC']
    print(input_data, file=stderr)
    bh = Newton(ic['g'], ic['m'], ic['Lfac'], ic['r0'])
    step = ic['step']
    integrator = Symplectic(bh, step, ic['integrator'], ic['scheme'])
    bh.solve(integrator.method, step, ic['start'], ic['end'], ic['plotratio'])
else:
    print(__name__ + " module loaded", file=stderr)
