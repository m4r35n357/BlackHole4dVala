#!/usr/bin/env python
"""
Copyright (c) 2014, 2015, 2016, 2017, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
from sys import argv, stdin, stderr
from matplotlib import pyplot
from json import loads

def main():
    print "X-Y Plotter: {}".format(argv)
    if len(argv) < 3:
        raise Exception('>>> ERROR! Please supply a plotting interval and two quantities to plot <<<')
    interval = int(argv[1])
    coordinate1 = argv[2]
    coordinate2 = argv[3]
    line = stdin.readline()
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    ax1.set_xlabel(coordinate1, color='b')
    ax1.set_ylabel(coordinate2, color='b')
    # ax1.set_xlim(-1.0, 2.0)
    # ax1.set_ylim(-1.0, 2.0)
    n = 0
    x = []
    y = []
    while line:
        p = loads(line)
        if n % interval == 0:
            x.append(p[coordinate1])
            y.append(p[coordinate2])
        line = stdin.readline()
        n += 1
    ax1.plot(x, y, 'bo:', markersize=3)
    pyplot.show()

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"
