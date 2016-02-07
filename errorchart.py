#!/usr/bin/env python

from sys import argv, stdin
from math import fabs, log10
from matplotlib import pyplot
from matplotlib.ticker import MultipleLocator
from json import loads
from array import array
from scipy.interpolate import InterpolatedUnivariateSpline
import numpy as np

def main():
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
    ax1.set_ylim(-150.0, 0.0)
    ax2 = ax1.twinx()
    ax2.yaxis.set_major_locator(majorLocator)
    ax2.yaxis.set_minor_locator(minorLocator)
    ax2.set_ylabel('Radial (blue) and Latitudinal (red) Errors, dB', color='0.25')
    ax2.set_ylim(-150.0, 0.0)
    count = 0
    line = stdin.readline()
    while line:  # build raw data arrays
        p = loads(line)
        if (count % interval == 0):
            timeValue = p[timeCoordinate]
            ax1.plot(timeValue, float(p['v4e']), color='#006000', linestyle='-', marker='.', markersize=2, zorder=11)
            ax1.plot(timeValue, float(p['v4c']), color='#606060', linestyle='-', marker='.', markersize=2, zorder=10)
            ax2.plot(timeValue, float(p['ER']), color='blue', linestyle='-', marker=',', markersize=1, zorder=9)
            ax2.plot(timeValue, float(p['ETh']), color='red', linestyle='-', marker=',', markersize=1, zorder=8)
        line = stdin.readline()
        count += 1
    try:
        pyplot.show()
    except AttributeError as e:
        print('ATTRIBUTE ERROR: ' + str(argv[0]) + ':' + str(coordinate) + ': ' + str(e))
        exit(-3)

if __name__ == "__main__":
    main()
else:
    print >> sys.stderr, __name__ + " module loaded"
