#!/usr/bin/env pypy
'''
Copyright (c) 2014, 2015, 2016, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''
from sys import stderr
from math import log10

class Integrator(object):
    def __init__(self, model, order):
        self.cbrt2 = 2.0**(1.0 / 3.0)
        self.f2 = 1.0 / (2.0 - self.cbrt2)
        self.coefficients = [0.5 * self.f2, self.f2, 0.5 * (1.0 - self.cbrt2) * self.f2, - self.cbrt2 * self.f2]
        if order == 'sb2':  # Second order base
            self.base = self.base2;
            self.w = [1.0]
        elif order == 'sc4':  # Fourth order, composed from Second order
            self.base = self.base2;
            self.w = [self.coefficients[1], self.coefficients[3], self.coefficients[1]]
        elif order == 'sb4':  # Fourth order base
            self.base = self.base4;
            self.w = [1.0]
        elif order == 'sc6':  # Sixth order, composed from Fourth order
            self.base = self.base4;
            fthrt2 = 2.0**(1.0 / 5.0)
            self.w = [1.0 / (2.0 - fthrt2), - fthrt2 / (2.0 - fthrt2), 1.0 / (2.0 - fthrt2)]
        else:  # Wrong value for integrator order
            raise Exception('>>> ERROR! Integrator order must be sb2, sc4, sb4, or sc6 <<<')
        self.wRange = range(len(self.w))
        self.model = model

    def base2 (self, w):  # Compose higher orders from this second-order symplectic base (d2 = 0.0)
        self.model.pUp(w * 0.5)  # c1 = 0.5
        self.model.qUp(w)        # d1 = 1.0
        self.model.pUp(w * 0.5)  # c2 = 0.5

    def base4 (self, w):  # Compose higher orders from this fourth-order symplectic base (d4 = 0.0)
        self.model.pUp(w * self.coefficients[0])  # w * c1
        self.model.qUp(w * self.coefficients[1])  # w * d1
        self.model.pUp(w * self.coefficients[2])  # w * c2
        self.model.qUp(w * self.coefficients[3])  # w * d2
        self.model.pUp(w * self.coefficients[2])  # w * c3
        self.model.qUp(w * self.coefficients[1])  # w * d3
        self.model.pUp(w * self.coefficients[0])  # w * c4

    def compose (self):
        for i in self.wRange:  # Composition happens in this loop
            self.base(self.w[i])

def logError (e):
    return 10.0 * log10(e if e > 1.0e-18 else 1.0e-18)
'''
if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"
'''

