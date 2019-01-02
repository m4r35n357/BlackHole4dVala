#!/usr/bin/env python3

#  Example: ./Analysis.py <analysis.json | ./plotXY.py 1 c d

"""
Copyright (c) 2014-2018, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

from gmpy2 import get_context, mpfr
get_context().precision = 113  # Set this BEFORE importing any Taylor Series stuff!
from json import loads
from sys import stdin, stderr, argv
from Symplectic import Symplectic


class Analysis(object):
    def __init__(self):
        self.c = 0.0
        self.d = 0.0
        self.count = 0

    def q_update(self, c):
        self.c += c
        self.count += 1
        self.output()

    def p_update(self, d):
        self.d += d
        self.count += 1
        self.output()

    def solve(self, method):
        self.output()
        method()

    def output(self):
        print('{{"count":{:d},"c":{:.9e},"d":{:.9e}}}'.format(self.count, self.c, self.d))

if __name__ == "__main__":
    print("Simulator: {}".format(argv[0]), file=stderr)
    input_data = stdin.read()
    ic = loads(input_data, parse_float=mpfr)['IC']
    print(input_data, file=stderr)
    a = Analysis()
    a.solve(Symplectic(a, 1.0, ic['integrator'], ic['scheme']).method)
else:
    print(__name__ + " module loaded", file=stderr)
