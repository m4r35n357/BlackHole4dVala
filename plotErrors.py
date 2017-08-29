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
from json import loads
from math import log10
from matplotlib import pyplot
from matplotlib.ticker import MultipleLocator
from sys import argv, stdin, stderr


def log_error(e):
    error = e if e >= 0.0 else -e
    return 10.0 * log10(error) if error > 1.0e-21 else -210.0


def main():
    print "Error Plotter: {}".format(argv)
    if len(argv) < 3:
        raise Exception('>>> ERROR! Please supply a time variable name and a plotting interval <<<')
    timeCoordinate = str(argv[1])
    interval = int(argv[2])
    ax1 = pyplot.figure().add_subplot(111)
    pyplot.minorticks_on()
    majorLocator = MultipleLocator(30)
    minorLocator = MultipleLocator(10)
    pyplot.grid(b=True, which='major', color='0.25', linestyle='-')
    pyplot.grid(b=True, which='minor', color='0.25', linestyle=':')
    ax1.yaxis.set_major_locator(majorLocator)
    ax1.yaxis.set_minor_locator(minorLocator)
    ax1.set_xlabel('Time: ' + timeCoordinate, color='0.20')
    ax1.set_ylabel('4-Velocity Norm Error', color='#006000')
    ax1.set_ylim(-210.0, 0.0)
    ax2 = ax1.twinx()
    ax2.yaxis.set_major_locator(majorLocator)
    ax2.yaxis.set_minor_locator(minorLocator)
    ax2.set_ylabel('Radial (blue) and Latitudinal (red) Errors, dB', color='0.25')
    ax2.set_ylim(-210.0, 0.0)
    pyplot.axhspan(-30.0, 0.0, facecolor='red', alpha=0.3)
    pyplot.axhspan(-60.0, -30.0, facecolor='orange', alpha=0.3)
    pyplot.axhspan(-90.0, -60.0, facecolor='yellow', alpha=0.3)
    pyplot.axhspan(-120.0, -90.0, facecolor='cyan', alpha=0.3)
    pyplot.axhspan(-210.0, -120.0, facecolor='green', alpha=0.3)
    count = 0
    eCum = ePk = 0.0
    line = stdin.readline()
    while line:  # build raw data arrays
        p = loads(line)
        e = p['v4e']
        e = e if e >= 0.0 else -e
        count += 1
        if count % interval == 0:
            timeValue = p[timeCoordinate]
            if 'ER' in p:
                ax2.plot(timeValue, log_error(p['ER']), color='blue', linestyle='-', marker='.', markersize=1)
            if 'ETh' in p:
                ax2.plot(timeValue, log_error(p['ETh']), color='red', linestyle='-', marker='.', markersize=1)
            ax1.plot(timeValue, log_error(eCum / count), color='black', linestyle='-', marker='.', markersize=1, zorder=10)
            ax1.plot(timeValue, log_error(e), color='#000f00', linestyle='-', marker='.', markersize=2)
        line = stdin.readline()
        eCum += e
        ePk = ePk if ePk > e else e
    print argv[0] + " - Average: " + str(log_error(eCum / count)) + ", Peak: " + str(log_error(ePk))
    try:
        pyplot.show()
    except AttributeError as e:
        print('ATTRIBUTE ERROR: ' + str(argv[0]) + ':' + str(e))
        exit(-3)

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"
