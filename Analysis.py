#!/usr/bin/env python
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
from numpy import longfloat, pi, sqrt, cos, sin
from sys import stdin, stderr, argv

from Symplectic import Symplectic, D1, D2, D0, D3


class Analysis(object):
    def __init__(self):
        self.c = 0.0
        self.d = 0.0


    def q_update(self, c):
        self.c += c
        print '{{"mino":{:.9e},"tau":{:.9e},"v4e":{:.9e},"ER":{:.9e},"ETh":{:.9e},"c":{:.9e},"d":{:.9e},"th":{:.9e},"ph":{:.9e},"tP":{:.9e},"rP":{:.9e},"thP":{:.9e},"phP":{:.9e}}}'.format(
            0.0, 0.0, 0.0, 0.0, 0.0, self.c, self.d, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    def p_update(self, d):
        self.d += d
        print '{{"mino":{:.9e},"tau":{:.9e},"v4e":{:.9e},"ER":{:.9e},"ETh":{:.9e},"c":{:.9e},"d":{:.9e},"th":{:.9e},"ph":{:.9e},"tP":{:.9e},"rP":{:.9e},"thP":{:.9e},"phP":{:.9e}}}'.format(
            0.0, 0.0, 0.0, 0.0, 0.0, self.c, self.d, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

    def solve(self, method, h, start, end, tr):
        tau = 0.0
        while tau < end:
            print '{{"mino":{:.9e},"tau":{:.9e},"v4e":{:.9e},"ER":{:.9e},"ETh":{:.9e},"c":{:.9e},"d":{:.9e},"th":{:.9e},"ph":{:.9e},"tP":{:.9e},"rP":{:.9e},"thP":{:.9e},"phP":{:.9e}}}'.format(
                0.0, 0.0, 0.0, 0.0, 0.0, self.c, self.d, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
            method()
            tau += 1.0


if __name__ == "__main__":
    print >> stderr, "Simulator: {}".format(argv[0])
    input_data = stdin.read()
    ic = loads(input_data, parse_float=longfloat)['IC']
    print >> stderr, input_data
    a = Analysis()
    step = ic['step']
    a.solve(Symplectic(a, step, ic['integrator'], ic['scheme']).method, step, ic['start'], ic['end'], ic['plotratio'])
else:
    print >> stderr, __name__ + " module loaded"
