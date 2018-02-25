"""
Copyright (c) 2014-2018, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
from sys import stderr
from numpy import longfloat

D0 = longfloat(0.0)
D1 = longfloat(1.0)
D2 = longfloat(2.0)
D3 = longfloat(3.0)
D4 = longfloat(4.0)
D5 = longfloat(5.0)
D7 = longfloat(7.0)
D9 = longfloat(9.0)
D05 = longfloat(0.5)

class Symplectic(object):
    def __init__(self, model, h, order, scheme):
        self.model = model
        self.h = h
        if scheme == 'yoshida':
            print >> stderr, "Yoshida composition"
            self.scheme = self.yoshida
            self.scheme_root = D2
        elif scheme == 'suzuki':
            print >> stderr, "Suzuki composition"
            self.scheme = self.suzuki
            self.scheme_root = D4
        else:
            raise Exception('>>> Composition scheme must be yoshida or suzuki, was "{found}" <<<'.format(found=scheme))
        if order == 'b1':
            print >> stderr, "1st order (Euler-Cromer)"
            self.method = self.euler_cromer
        elif order == 'fr4':
            print >> stderr, "4th order (Forest-Ruth)"
            theta = D1 / (D2 - D2 ** (D1 / D3))
            self.c_d = [D05 * h * theta, h * theta, D05 * h * (D1 - theta), h * (D1 - D2 * theta)]
            self.method = self.forest_ruth
        elif order == 'b2':
            print >> stderr, "2nd order (Stormer-Verlet))"
            self.method = self.second_order
        elif order == 'b4':
            print >> stderr, "4th order (Composed)"
            self.method = self.fourth_order
        elif order == 'b6':
            print >> stderr, "6th order (Composed)"
            self.method = self.sixth_order
        elif order == 'b8':
            print >> stderr, "8th order (Composed)"
            self.method = self.eightth_order
        elif order == 'b10':
            print >> stderr, "10th order (Composed)"
            self.method = self.tenth_order
        else:
            raise Exception('>>> Integrator must be b1, b2, b4, b6, b8, b10 or fr4, was "{found}" <<<'.format(found=order))
        self.w_outer = D1 / (self.scheme_root - self.scheme_root**(D1 / D9))
        self.x_outer = D1 / (self.scheme_root - self.scheme_root**(D1 / D7))
        self.y_outer = D1 / (self.scheme_root - self.scheme_root**(D1 / D5))
        self.z_outer = D1 / (self.scheme_root - self.scheme_root**(D1 / D3))
        self.w_central = D1 - self.scheme_root * self.w_outer
        self.x_central = D1 - self.scheme_root * self.x_outer
        self.y_central = D1 - self.scheme_root * self.y_outer
        self.z_central = D1 - self.scheme_root * self.z_outer

    def euler_cromer(self):
        self.model.q_update(self.h)
        self.model.p_update(self.h)

    def forest_ruth(self):
        self.model.q_update(self.c_d[0])
        self.model.p_update(self.c_d[1])
        self.model.q_update(self.c_d[2])
        self.model.p_update(self.c_d[3])
        self.model.q_update(self.c_d[2])
        self.model.p_update(self.c_d[1])
        self.model.q_update(self.c_d[0])

    @staticmethod
    def yoshida(base_method, s, plus, minus):
        base_method(s * plus)
        base_method(s * minus)
        base_method(s * plus)

    @staticmethod
    def suzuki(base_method, s, plus, minus):
        base_method(s * plus)
        base_method(s * plus)
        base_method(s * minus)
        base_method(s * plus)
        base_method(s * plus)

    def stormer_verlet(self, s):
        self.model.q_update(self.h * s * D05)
        self.model.p_update(self.h * s)
        self.model.q_update(self.h * s * D05)

    def second_order(self):
        self.stormer_verlet(D1)

    def base4(self, s):
        self.scheme(self.stormer_verlet, s, self.z_outer, self.z_central)

    def fourth_order(self):
        # noinspection PyTypeChecker
        self.base4(D1)

    def base6(self, s):
        self.scheme(self.base4, s, self.y_outer, self.y_central)

    def sixth_order(self):
        # noinspection PyTypeChecker
        self.base6(D1)

    def base8(self, s):
        self.scheme(self.base6, s, self.x_outer, self.x_central)

    def eightth_order(self):
        # noinspection PyTypeChecker
        self.base8(D1)

    def tenth_order(self):
        # noinspection PyTypeChecker
        self.scheme(self.base8, D1, self.w_outer, self.w_central)


print >> stderr, __name__ + " module loaded"
