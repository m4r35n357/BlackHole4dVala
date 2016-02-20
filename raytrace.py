#!/usr/bin/env python
'''
Copyright (c) 2014, 2015, 2016, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''
from sys import argv, stdin, stdout, stderr
from math import sqrt, sin, cos, pi
from json import loads

def main():
    if len(argv) < 1:
        raise Exception('>>> ERROR! Please supply a viewing angle in degrees <<<')
    angle = long(argv[1]) * pi / 180.0
    c = cos(angle)
    s = sin(angle)
    line = stdin.readline()
    while line:  # build raw data arrays line by line
        p = loads(line)
        print >> stdout, str(p['tau']) + ' 2 ' \
            + str(p['r']) + ' ' + str(cos(p['th'])) + ' ' + str(p['t']) + ' ' + str(p['ph']) + ' ' \
            + str(p['rP']) + ' ' + str(- sin(p['th']) * p['thP']) + ' ' + str(p['tP']) + ' ' + str(p['phP']) + ' ' \
            + str(-c) + ' ' + str(0) + ' ' + str(0) + ' ' + str(s) + ' ' \
            + str(s) + ' ' + str(0) + ' ' + str(0) + ' ' + str(c) + ' ' \
            + str(0) + ' ' + str(1) + ' ' + str(0) + ' ' + str(0)
        line = stdin.readline()

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"
