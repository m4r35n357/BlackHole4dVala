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
from sys import stderr


class Symplectic(object):
    def __init__(self, model, h, order, debug=False):
        if debug:
            self.base2 = self.coefficients
        else:
            self.base2 = self.stormer_verlet
        self.total = 0.0
        self.total_abs = 0.0
        self.count = 0
        self.model = model
        self.h = h
        self.x1 = 1.0 / (4.0 - 4.0**(1.0 / 7.0))
        self.y1 = 1.0 / (4.0 - 4.0**(1.0 / 5.0))
        self.z1 = 1.0 / (4.0 - 4.0**(1.0 / 3.0))
        self.x3 = 1.0 - 4.0 * self.x1
        self.y3 = 1.0 - 4.0 * self.y1
        self.z3 = 1.0 - 4.0 * self.z1
        if order == 'b1':
            print >> stderr, "Python First order Symplectic Integrator"
            self.method = self.first_order
        elif order == 'b2':
            print >> stderr, "Python Second order Symplectic Integrator"
            self.method = self.second_order
        elif order == 'sb4':
            print >> stderr, "Python Fourth order Symplectic Integrator (using Suzuki Composition)"
            self.method = self.fourth_order
        elif order == 'sb6':
            print >> stderr, "Python Sixth order Symplectic Integrator (using Suzuki Composition)"
            self.method = self.sixth_order
        elif order == 'sb8':
            print >> stderr, "Python Eightth order Symplectic Integrator (using Suzuki Composition)"
            self.method = self.eightth_order
        else:
            raise Exception('>>> Integrator must be b1, b2, sb4, sb6 or sb8, was "{}" <<<'.format(order))
        if debug:
            print self.total

    def first_order(self):
        self.model.qUpdate(self.h)
        self.model.pUpdate(self.h)

    def coefficients(self, s):
        self.count += 1
        self.total += s
        self.total_abs += abs(s)
        print s, self.total, self.total_abs / self.count

    def stormer_verlet(self, s):
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


print >> stderr, __name__ + " module loaded"
