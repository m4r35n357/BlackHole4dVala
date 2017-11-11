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

import os


def log_error(e):
    """
    convert error to a pseudo db Value
    :param e: the numerical value of he error
    :return: the dB value, clamped to a minimum
    """
    error = e if e >= 0.0 else -e
    return 10.0 * log10(error) if error > 1.0e-21 else -210.0


def main():
    """
    Main method
    """
    print "Error Plotter: {}".format(argv)
    if len(argv) < 6:
        raise Exception(
            '>>> ERROR! Please supply an integrator type, number of composition stages, a time step, a time variable name and a plotting interval <<<')
    executable = os.environ['exe']
    integrator_type = str(argv[1])
    composition_stages = str(argv[2])
    time_step = str(argv[3])
    time_coordinate = str(argv[4])
    interval = int(argv[5])
    left = pyplot.figure().add_subplot(111)
    pyplot.minorticks_on()
    major_locator = MultipleLocator(30)
    minor_locator = MultipleLocator(10)
    pyplot.grid(b=True, color='0.25', linestyle='-')
    pyplot.grid(b=True, which='minor', color='0.25', linestyle=':')
    left.yaxis.set_major_locator(major_locator)
    left.yaxis.set_minor_locator(minor_locator)
    left.set_xlabel('Time: ' + time_coordinate, color='0.20')
    left.set_ylabel('4-Velocity Norm Error', color='#006000')
    left.set_ylim(-210.0, 0.0)
    right = left.twinx()
    right.yaxis.set_major_locator(major_locator)
    right.yaxis.set_minor_locator(minor_locator)
    right.set_ylabel('Radial (blue) and Latitudinal (red) Errors, dB', color='0.25')
    right.set_ylim(-210.0, 0.0)
    pyplot.axhspan(-30.0, 0.0, facecolor='red', alpha=0.3)
    pyplot.axhspan(-60.0, -30.0, facecolor='orange', alpha=0.3)
    pyplot.axhspan(-90.0, -60.0, facecolor='yellow', alpha=0.3)
    pyplot.axhspan(-120.0, -90.0, facecolor='cyan', alpha=0.3)
    pyplot.axhspan(-210.0, -120.0, facecolor='green', alpha=0.3)
    count = 0
    e_cum = e_pk = 0.0
    line = stdin.readline()
    while line:  # build raw data arrays
        p = loads(line)
        e = p['v4e']
        e = e if e >= 0.0 else -e
        count += 1
        if count % interval == 0:
            time_value = p[time_coordinate]
            if 'ER' in p:
                right.plot(time_value, log_error(p['ER']), color='blue', linestyle='-', marker='.', markersize=1)
            if 'ETh' in p:
                right.plot(time_value, log_error(p['ETh']), color='red', linestyle='-', marker='.', markersize=1)
            left.plot(time_value, log_error(e_cum / count), color='black', linestyle='-', marker='.', markersize=1,
                      zorder=10)
            left.plot(time_value, log_error(e), color='#000f00', linestyle='-', marker='.', markersize=2)
        line = stdin.readline()
        e_cum += e
        e_pk = e_pk if e_pk > e else e
    left.annotate(executable + " - " + integrator_type + " (" + composition_stages + " stages),  ts = " + time_step,
                  (0.0, 0.0), xytext=(0.15, 0.96), textcoords='figure fraction', color='0.20', )
    left.annotate("Peak: " + str(log_error(e_pk)) + ",   Average: " + str(log_error(e_cum / count)), (0.0, 0.0),
                  xytext=(0.25, 0.92), textcoords='figure fraction', color='0.20', )
    try:
        pyplot.show()
    except AttributeError as e:
        print('ATTRIBUTE ERROR: ' + str(argv[0]) + ':' + str(e))
        exit(-3)


if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"
