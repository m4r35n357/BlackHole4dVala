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
from sys import argv, stdin
from matplotlib import pyplot
from json import loads

def main():
    line = stdin.readline()
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, which='major', color='0.25', linestyle='-')
    ax1.set_xlabel('r, theta', color='0.20')
    ax1.set_ylabel('R(r)', color='b')
    ax1.set_ylim(-20, 20)
    ax2 = ax1.twinx()
    ax2.set_ylabel('THETA(theta)', color='r')
    ax2.set_ylim(-20, 20)
    n = 0
    while line:
        p = loads(line)
        ax1.plot(p['x'], p['R'], 'b.', markersize=2)
        ax2.plot(p['x'], p['THETA'], 'r.', markersize=2)
        line = stdin.readline()
        n += 1
    try:
        pyplot.show()
    except AttributeError as e:
        print('ATTRIBUTE ERROR: ' + str(argv[0]) + ': ' + str(e))
        exit(-1)

if __name__ == "__main__":
    main()
else:
    print >> sys.stderr, __name__ + " module loaded"
