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
        elif order == 'b2':
            print >> stderr, "2nd order (Stormer-Verlet))"
            self.method = self.second_order
        elif order == 'b4':
            print >> stderr, "4th order (Composed)"
            self.method = self.fourth_order
        elif order == 'f4':
            print >> stderr, "4th order (Forest-Ruth)"
            self.scheme = self.yoshida
            self.scheme_root = D2
            self.method = self.fourth_order_forest_ruth
        elif order == 's4':
            print >> stderr, "4th order (Smith)"
            self.scheme = self.suzuki
            self.scheme_root = D4
            self.method = self.fourth_order_smith
        elif order == 'b6':
            print >> stderr, "6th order (Composed)"
            self.method = self.sixth_order
        elif order == 'f6':
            print >> stderr, "6th order (Forest-Ruth) (Composed)"
            self.scheme = self.yoshida
            self.scheme_root = D2
            self.method = self.sixth_order_forest_ruth
        elif order == 's6':
            print >> stderr, "6th order (Smith) (Composed)"
            self.scheme = self.suzuki
            self.scheme_root = D4
            self.method = self.sixth_order_smith
        elif order == 'b8':
            print >> stderr, "8th order (Composed)"
            self.method = self.eightth_order
        elif order == 'f8':
            print >> stderr, "8th order (Forest-Ruth) (Composed)"
            self.scheme = self.yoshida
            self.scheme_root = D2
            self.method = self.eightth_order_forest_ruth
        elif order == 's8':
            print >> stderr, "8th order (Smith) (Composed)"
            self.scheme = self.suzuki
            self.scheme_root = D4
            self.method = self.eightth_order_smith
        elif order == 'b10':
            print >> stderr, "10th order (Composed)"
            self.method = self.tenth_order
        elif order == 'f10':
            print >> stderr, "10th order (Forest-Ruth) (Composed)"
            self.scheme = self.yoshida
            self.scheme_root = D2
            self.method = self.tenth_order_forest_ruth
        elif order == 's10':
            print >> stderr, "10th order (Smith) (Composed)"
            self.scheme = self.suzuki
            self.scheme_root = D4
            self.method = self.tenth_order_smith
        else:
            raise Exception(
                '>>> Integrator must be [bfs]1, [bfs]2, [bfs]4, [bfs]6, [bfs]8, or [bfs]10, was "{found}" <<<'.format(
                    found=order))
        self.w1 = D1 / (self.scheme_root - self.scheme_root ** (D1 / D9))
        self.x1 = D1 / (self.scheme_root - self.scheme_root ** (D1 / D7))
        self.y1 = D1 / (self.scheme_root - self.scheme_root ** (D1 / D5))
        self.z1 = D1 / (self.scheme_root - self.scheme_root ** (D1 / D3))
        self.w0 = D1 - self.scheme_root * self.w1
        self.x0 = D1 - self.scheme_root * self.x1
        self.y0 = D1 - self.scheme_root * self.y1
        self.z0 = D1 - self.scheme_root * self.z1
        self.z0y0 = self.z0 * self.y0
        self.z0y1 = self.z0 * self.y1
        self.z1y0 = self.z1 * self.y0
        self.z1y1 = self.z1 * self.y1
        self.cd_sv = [D05 * h, h]
        self.cd_fr4 = [D05 * h * self.z1, h * self.z1, D05 * h * (self.z0 + self.z1), h * self.z0]
        self.cd_s4 = [D05 * h * self.z1, h * self.z1, h * self.z1, h * self.z1, D05 * h * (self.z0 + self.z1),
                      h * self.z0]
        self.cd_s6 = [D05 * h * self.z1y1, h * self.z1y1, h * self.z1y1, h * self.z1y1,
                      D05 * h * (self.z1y1 + self.z0y1), h * self.z0y1,
                      D05 * h * (self.z0y1 + self.z1y1), h * self.z1y1, h * self.z1y1, h * self.z1y1,
                      h * self.z1y1, h * self.z1y1, h * self.z1y1, h * self.z1y1,
                      D05 * h * (self.z1y1 + self.z0y1), h * self.z0y1,
                      D05 * h * (self.z0y1 + self.z1y1), h * self.z1y1, h * self.z1y1, h * self.z1y1,
                      D05 * h * (self.z1y1 + self.z1y0), h * self.z1y0, h * self.z1y0, h * self.z1y0,
                      D05 * h * (self.z1y0 + self.z0y0), h * self.z0y0]

    def euler_cromer(self):
        self.model.q_update(self.h)
        self.model.p_update(self.h)

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
        self.model.q_update(s * self.cd_sv[0])
        self.model.p_update(s * self.cd_sv[1])
        self.model.q_update(s * self.cd_sv[0])

    def forest_ruth_4(self, s):
        self.model.q_update(s * self.cd_fr4[0])
        self.model.p_update(s * self.cd_fr4[1])
        self.model.q_update(s * self.cd_fr4[2])
        self.model.p_update(s * self.cd_fr4[3])
        self.model.q_update(s * self.cd_fr4[2])
        self.model.p_update(s * self.cd_fr4[1])
        self.model.q_update(s * self.cd_fr4[0])

    def smith_4(self, s):
        self.model.q_update(s * self.cd_s4[0])
        self.model.p_update(s * self.cd_s4[1])
        self.model.q_update(s * self.cd_s4[2])
        self.model.p_update(s * self.cd_s4[3])
        self.model.q_update(s * self.cd_s4[4])
        self.model.p_update(s * self.cd_s4[5])
        self.model.q_update(s * self.cd_s4[4])
        self.model.p_update(s * self.cd_s4[3])
        self.model.q_update(s * self.cd_s4[2])
        self.model.p_update(s * self.cd_s4[1])
        self.model.q_update(s * self.cd_s4[0])

    def smith_6(self, s):
        self.model.q_update(s * self.cd_s6[0])
        self.model.p_update(s * self.cd_s6[1])
        self.model.q_update(s * self.cd_s6[2])
        self.model.p_update(s * self.cd_s6[3])
        self.model.q_update(s * self.cd_s6[4])
        self.model.p_update(s * self.cd_s6[5])
        self.model.q_update(s * self.cd_s6[6])
        self.model.p_update(s * self.cd_s6[7])
        self.model.q_update(s * self.cd_s6[8])
        self.model.p_update(s * self.cd_s6[9])
        self.model.q_update(s * self.cd_s6[10])
        self.model.p_update(s * self.cd_s6[11])
        self.model.q_update(s * self.cd_s6[12])
        self.model.p_update(s * self.cd_s6[13])
        self.model.q_update(s * self.cd_s6[14])
        self.model.p_update(s * self.cd_s6[15])
        self.model.q_update(s * self.cd_s6[16])
        self.model.p_update(s * self.cd_s6[17])
        self.model.q_update(s * self.cd_s6[18])
        self.model.p_update(s * self.cd_s6[19])
        self.model.q_update(s * self.cd_s6[20])
        self.model.p_update(s * self.cd_s6[21])
        self.model.q_update(s * self.cd_s6[22])
        self.model.p_update(s * self.cd_s6[23])
        self.model.q_update(s * self.cd_s6[24])
        self.model.p_update(s * self.cd_s6[25])
        self.model.q_update(s * self.cd_s6[24])
        self.model.p_update(s * self.cd_s6[23])
        self.model.q_update(s * self.cd_s6[22])
        self.model.p_update(s * self.cd_s6[21])
        self.model.q_update(s * self.cd_s6[20])
        self.model.p_update(s * self.cd_s6[19])
        self.model.q_update(s * self.cd_s6[18])
        self.model.p_update(s * self.cd_s6[17])
        self.model.q_update(s * self.cd_s6[16])
        self.model.p_update(s * self.cd_s6[15])
        self.model.q_update(s * self.cd_s6[14])
        self.model.p_update(s * self.cd_s6[13])
        self.model.q_update(s * self.cd_s6[12])
        self.model.p_update(s * self.cd_s6[11])
        self.model.q_update(s * self.cd_s6[10])
        self.model.p_update(s * self.cd_s6[9])
        self.model.q_update(s * self.cd_s6[8])
        self.model.p_update(s * self.cd_s6[7])
        self.model.q_update(s * self.cd_s6[6])
        self.model.p_update(s * self.cd_s6[5])
        self.model.q_update(s * self.cd_s6[4])
        self.model.p_update(s * self.cd_s6[3])
        self.model.q_update(s * self.cd_s6[2])
        self.model.p_update(s * self.cd_s6[1])
        self.model.q_update(s * self.cd_s6[0])

    def second_order(self):
        # noinspection PyTypeChecker
        self.stormer_verlet(D1)

    def base4(self, s):
        self.scheme(self.stormer_verlet, s, self.z1, self.z0)

    def fourth_order(self):
        # noinspection PyTypeChecker
        self.base4(D1)

    def base6(self, s):
        self.scheme(self.base4, s, self.y1, self.y0)

    def sixth_order(self):
        # noinspection PyTypeChecker
        self.base6(D1)

    def base8(self, s):
        self.scheme(self.base6, s, self.x1, self.x0)

    def eightth_order(self):
        # noinspection PyTypeChecker
        self.base8(D1)

    def tenth_order(self):
        # noinspection PyTypeChecker
        self.scheme(self.base8, D1, self.w1, self.w0)

    def fourth_order_forest_ruth(self):
        # noinspection PyTypeChecker
        self.forest_ruth_4(D1)

    def base6_forest_ruth(self, s):
        self.scheme(self.forest_ruth_4, s, self.y1, self.y0)

    def sixth_order_forest_ruth(self):
        # noinspection PyTypeChecker
        self.base6_forest_ruth(D1)

    def base8_forest_ruth(self, s):
        self.scheme(self.base6_forest_ruth, s, self.x1, self.x0)

    def eightth_order_forest_ruth(self):
        # noinspection PyTypeChecker
        self.base8_forest_ruth(D1)

    def tenth_order_forest_ruth(self):
        # noinspection PyTypeChecker
        self.scheme(self.base8_forest_ruth, D1, self.w1, self.w0)

    def fourth_order_smith(self):
        # noinspection PyTypeChecker
        self.smith_4(D1)

    def base6_smith(self, s):
        self.scheme(self.smith_4, s, self.y1, self.y0)

    def sixth_order_smith(self):
        # noinspection PyTypeChecker
        # self.base6_smith(D1)
        self.smith_6(D1)

    def base8_smith(self, s):
        # self.scheme(self.base6_smith, s, self.x1, self.x0)
        self.scheme(self.smith_6, s, self.x1, self.x0)

    def eightth_order_smith(self):
        # noinspection PyTypeChecker
        self.base8_smith(D1)

    def tenth_order_smith(self):
        # noinspection PyTypeChecker
        self.scheme(self.base8_smith, D1, self.w1, self.w0)


print >> stderr, __name__ + " module loaded"
