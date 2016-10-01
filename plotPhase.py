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
from sys import argv, stdin, stdout, stderr
from matplotlib import pyplot
from json import loads
from array import array
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy as np

def main():
    if len(argv) < 3:
        raise Exception('>>> ERROR! Please supply a time variable, the number of points to plot, and a coordinate to plot <<<')
    timeCoordinate = str(argv[1])
    nPoints = int(argv[2])
    coordinate = argv[3]
    tau = array('d')
    c = array('d')
    cDot = array('d')
    tauMax = 0.0
    line = stdin.readline()
    while line:  # build raw data arrays
        p = loads(line)
        tauValue = p[timeCoordinate]
        tauMax = tauValue if tauValue > tauMax else tauMax
        tau.append(tauValue)
        c.append(p[coordinate])
        cDot.append(p[coordinate + 'P'])
        line = stdin.readline()
    try:  # interpolate here
        cI = InterpolatedUnivariateSpline(tau, c, k = 1)
        cDotI = InterpolatedUnivariateSpline(tau, cDot, k = 1)
    except ValueError as e:
        print('DATA ERROR: ' + str(argv[0]) + ': ' + str(e))
        exit(-2)
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.grid(b=True, which='major', color='0.25', linestyle='-')
    ax1.set_ylabel(coordinate + 'Dot', color='b')
    ax1.set_xlabel(coordinate, color='b')
    tauI = np.linspace(0, tauMax, num = nPoints)
    for i in range(len(tauI)):
        ax1.plot(cI(tauI[i]), cDotI(tauI[i]), color='magenta', marker='.', markersize=1)
    try:
        pyplot.show()
    except AttributeError as e:
        print('ATTRIBUTE ERROR: ' + str(argv[0]) + ':' + str(coordinate) + ': ' + str(e))
        exit(-3)

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"

