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
from gmpy2 import mpfr

D0 = mpfr("0.0")
D05 = mpfr("0.5")
D1 = mpfr("1.0")
D2 = mpfr("2.0")
D3 = mpfr("3.0")
D4 = mpfr("4.0")
D5 = mpfr("5.0")
D7 = mpfr("7.0")
D9 = mpfr("9.0")


class Symplectic(object):

    def __init__(self, model, h, order, scheme):
        self.model = model
        self.h = h
        if scheme == 'yoshida':
            print("Yoshida composition", file=stderr)
            self.scheme = self.yoshida
            self.scheme_root = D2
        elif scheme == 'suzuki':
            print("Suzuki composition", file=stderr)
            self.scheme = self.suzuki
            self.scheme_root = D4
        else:
            raise Exception('>>> Composition scheme must be yoshida or suzuki, was "{found}" <<<'.format(found=scheme))
        if order == 'b1':
            print("1st order (Euler-Cromer)", file=stderr)
            self.method = self.euler_cromer
        elif order == 'b2':
            print("2nd order (Stormer-Verlet))", file=stderr)
            self.method = self.second_order
        elif order == 'b4':
            print("4th order (Composed)", file=stderr)
            self.method = self.fourth_order
        elif order == 'f4':
            print("4th order (Forest-Ruth)", file=stderr)
            self.scheme = self.yoshida
            self.scheme_root = D2
            self.method = self.fourth_order_forest_ruth
        elif order == 'b6':
            print("6th order (Composed)", file=stderr)
            self.method = self.sixth_order
        elif order == 'f6':
            print("6th order (Forest-Ruth) (Composed)", file=stderr)
            self.scheme = self.yoshida
            self.scheme_root = D2
            self.method = self.sixth_order_forest_ruth
        elif order == 'b8':
            print("8th order (Composed)", file=stderr)
            self.method = self.eightth_order
        elif order == 'f8':
            print("8th order (Forest-Ruth) (Composed)", file=stderr)
            self.scheme = self.yoshida
            self.scheme_root = D2
            self.method = self.eightth_order_forest_ruth
        elif order == 'b10':
            print("10th order (Composed)", file=stderr)
            self.method = self.tenth_order
        elif order == 'f10':
            print("10th order (Forest-Ruth) (Composed)", file=stderr)
            self.scheme = self.yoshida
            self.scheme_root = D2
            self.method = self.tenth_order_forest_ruth
        elif order == 's4':
            print("4th order (Smith)", file=stderr)
            self.scheme = self.suzuki
            self.scheme_root = D4
            self.method = self.fourth_order_smith
        elif order == 's6':
            print("6th order (Smith)", file=stderr)
            self.scheme = self.suzuki
            self.scheme_root = D4
            self.method = self.sixth_order_smith
        elif order == 's8':
            print("8th order (Smith)", file=stderr)
            self.scheme = self.suzuki
            self.scheme_root = D4
            self.method = self.eightth_order_smith
        elif order == 's10':
            print("10th order (Smith) (Composed)", file=stderr)
            self.scheme = self.suzuki
            self.scheme_root = D4
            self.method = self.tenth_order_smith
        else:
            raise Exception(
                '>>> Integrator must be b1, b2, [bfs]4, [bfs]6, [bfs]8, or [bfs]10, was "{found}" <<<'.format(
                    found=order))
        self.w1 = D1 / (self.scheme_root - self.scheme_root ** (D1 / D9))
        self.x1 = D1 / (self.scheme_root - self.scheme_root ** (D1 / D7))
        self.y1 = D1 / (self.scheme_root - self.scheme_root ** (D1 / D5))
        self.z1 = D1 / (self.scheme_root - self.scheme_root ** (D1 / D3))
        self.w0 = D1 - self.scheme_root * self.w1
        self.x0 = D1 - self.scheme_root * self.x1
        self.y0 = D1 - self.scheme_root * self.y1
        self.z0 = D1 - self.scheme_root * self.z1
        self.cd_sv = [
            D05 * h, h
        ]
        self.cd_f4 = [
            D05 * h * self.z1, h * self.z1, D05 * h * (self.z0 + self.z1), h * self.z0
        ]
        self.cd_s4 = [
            D05 * h * self.z1, h * self.z1, h * self.z1, h * self.z1, D05 * h * (self.z1 + self.z0), h * self.z0
        ]
        self.cd_s6 = [
            D05 * h * self.z1 * self.y1,
            h * self.z1 * self.y1, h * self.z1 * self.y1, h * self.z1 * self.y1,
            D05 * h * (self.z1 + self.z0) * self.y1, h * self.z0 * self.y1, D05 * h * (self.z0 + self.z1) * self.y1,
            h * self.z1 * self.y1, h * self.z1 * self.y1, h * self.z1 * self.y1,
            h * self.z1 * self.y1,
            h * self.z1 * self.y1, h * self.z1 * self.y1, h * self.z1 * self.y1,
            D05 * h * (self.z1 + self.z0) * self.y1, h * self.z0 * self.y1, D05 * h * (self.z0 + self.z1) * self.y1,
            h * self.z1 * self.y1, h * self.z1 * self.y1, h * self.z1 * self.y1,
            D05 * h * self.z1 * (self.y1 + self.y0),
            h * self.z1 * self.y0, h * self.z1 * self.y0, h * self.z1 * self.y0,
            D05 * h * (self.z1 + self.z0) * self.y0, h * self.z0 * self.y0
        ]
        self.cd_s8 = [
            D05 * h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            D05 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, D05 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            D05 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, D05 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            D05 * h * self.z1 * (self.y1 + self.y0) * self.x1,
            h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1,
            D05 * h * (self.z1 + self.z0) * self.y0 * self.x1, h * self.z0 * self.y0 * self.x1, D05 * h * (self.z0 + self.z1) * self.y0 * self.x1,
            h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1,
            D05 * h * self.z1 * (self.y1 + self.y0) * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            D05 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, D05 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            D05 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, D05 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            D05 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, D05 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            D05 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, D05 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            D05 * h * self.z1 * (self.y1 + self.y0) * self.x1,
            h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1,
            D05 * h * (self.z1 + self.z0) * self.y0 * self.x1, h * self.z0 * self.y0 * self.x1, D05 * h * (self.z0 + self.z1) * self.y0 * self.x1,
            h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1, h * self.z1 * self.y0 * self.x1,
            D05 * h * self.z1 * (self.y1 + self.y0) * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            D05 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, D05 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            D05 * h * (self.z1 + self.z0) * self.y1 * self.x1, h * self.z0 * self.y1 * self.x1, D05 * h * (self.z0 + self.z1) * self.y1 * self.x1,
            h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1, h * self.z1 * self.y1 * self.x1,
            D05 * h * self.z1 * self.y1 * (self.x1 + self.x0),
            h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0,
            D05 * h * (self.z1 + self.z0) * self.y1 * self.x0, h * self.z0 * self.y1 * self.x0, D05 * h * (self.z0 + self.z1) * self.y1 * self.x0,
            h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0,
            h * self.z1 * self.y1 * self.x0,
            h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0,
            D05 * h * (self.z1 + self.z0) * self.y1 * self.x0, h * self.z0 * self.y1 * self.x0, D05 * h * (self.z0 + self.z1) * self.y1 * self.x0,
            h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0, h * self.z1 * self.y1 * self.x0,
            D05 * h * self.z1 * (self.y1 + self.y0) * self.x0,
            h * self.z1 * self.y0 * self.x0, h * self.z1 * self.y0 * self.x0, h * self.z1 * self.y0 * self.x0,
            D05 * h * (self.z1 + self.z0) * self.y0 * self.x0, h * self.z0 * self.y0 * self.x0
        ]

    def euler_cromer(self):
        self.model.q_update(self.h)
        self.model.p_update(self.h)

    def stormer_verlet(self, s):
        self.model.q_update(s * self.cd_sv[0])
        self.model.p_update(s * self.cd_sv[1])
        self.model.q_update(s * self.cd_sv[0])

    def second_order(self):
        self.stormer_verlet(D1)

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

    def base4(self, s):
        self.scheme(self.stormer_verlet, s, self.z1, self.z0)

    def base6(self, s):
        self.scheme(self.base4, s, self.y1, self.y0)

    def base8(self, s):
        self.scheme(self.base6, s, self.x1, self.x0)

    def fourth_order(self):
        self.base4(D1)

    def sixth_order(self):
        self.base6(D1)

    def eightth_order(self):
        self.base8(D1)

    def tenth_order(self):
        self.scheme(self.base8, D1, self.w1, self.w0)

    def forest_ruth_4(self, s):
        self.model.q_update(s * self.cd_f4[0])
        self.model.p_update(s * self.cd_f4[1])
        self.model.q_update(s * self.cd_f4[2])
        self.model.p_update(s * self.cd_f4[3])
        self.model.q_update(s * self.cd_f4[2])
        self.model.p_update(s * self.cd_f4[1])
        self.model.q_update(s * self.cd_f4[0])

    def fourth_order_forest_ruth(self):
        self.forest_ruth_4(D1)

    def base6_forest_ruth(self, s):
        self.scheme(self.forest_ruth_4, s, self.y1, self.y0)

    def sixth_order_forest_ruth(self):
        self.base6_forest_ruth(D1)

    def base8_forest_ruth(self, s):
        self.scheme(self.base6_forest_ruth, s, self.x1, self.x0)

    def eightth_order_forest_ruth(self):
        self.base8_forest_ruth(D1)

    def tenth_order_forest_ruth(self):
        self.scheme(self.base8_forest_ruth, D1, self.w1, self.w0)

    def smith(self, s):
        size = len(self.coefficients)
        for i in range(size):
            (self.model.q_update if i % 2 == 0 else self.model.p_update)(s * self.coefficients[i])
        for i in range(size - 2, -1, -1):
            (self.model.q_update if i % 2 == 0 else self.model.p_update)(s * self.coefficients[i])

    def fourth_order_smith(self):
        self.coefficients = self.cd_s4
        self.smith(D1)

    def sixth_order_smith(self):
        self.coefficients = self.cd_s6
        self.smith(D1)

    def eightth_order_smith(self):
        self.coefficients = self.cd_s8
        self.smith(D1)

    def tenth_order_smith(self):
        self.coefficients = self.cd_s8
        self.scheme(self.smith, D1, self.w1, self.w0)


print(__name__ + " module loaded", file=stderr)
