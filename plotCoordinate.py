#!/usr/bin/env python
"""
Copyright (c) 2014, 2015, 2016, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""
from sys import argv, stderr
from matplotlib import pyplot
from json import loads

def main():
    if len(argv) < 3:
        raise Exception('>>> ERROR! Please supply a data file name, a plotting interval, and a coordinate to plot <<<')
    dataFile = open(argv[1], 'r')
    interval = int(argv[2])
    coordinate = argv[3]
    line = dataFile.readline()
    ax1 = pyplot.figure().add_subplot(111)
    ax1.set_xlabel('tau', color='k')
    ax1.set_ylabel(coordinate, color='b')
    ax2 = ax1.twinx()
    ax2.set_ylabel(coordinate + 'Dot', color='r')
    n = 0
    while line:
        p = loads(line)
        if n % interval == 0:
            ax1.plot(p['tau'], p[coordinate], 'b.', markersize=2)
            ax2.plot(p['tau'], p[coordinate + 'P'], 'r.', markersize=2)
        line = dataFile.readline()
        n += 1
#        pyplot.legend(['t', 'r', 'th', 'ph'], loc='best')
    pyplot.show()

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"
