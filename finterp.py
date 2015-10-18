#!/usr/bin/env python

from sys import argv, stdin, stdout
from math import fabs
from json import loads

def signChange (time, target, olderTime, newerTime):
    return (newerTime - target) > 0.0
#    return (newerTime - target) * (olderTime - target) < 0.0

def nearest (time, target, olderTime, newerTime, older, newer):
    if fabs(newerTime - target) <= fabs(olderTime - target):
        return newer
    else:
        return older

def main():
    if len(argv) > 3:
        raise Exception('>>> Please supply a time variable name (string) and a precision (float) <<<')
    timeVariable = argv[1]
    precision = float(argv[2])
    dataFile = stdin
    counter = 0
    previous = dataFile.readline()
    previousJson = loads(previous)
    previousTime = float(previousJson[timeVariable])
    latest = dataFile.readline()
    while latest:
        latestJson = loads(latest)
        latestTime = float(latestJson[timeVariable])
        target = counter * precision
        if signChange(timeVariable, target, previousTime, latestTime):
            print >> stdout, (previousTime)
            print >> stdout, (target)
            print >> stdout, (latestTime)
            print >> stdout, (nearest(timeVariable, target, previousTime, latestTime, previous, latest))
            counter += 1
        previous = latest
        previousJson = latestJson
        previousTime = latestTime
        latest = dataFile.readline()

if __name__ == "__main__":
    main()
else:
    print >> sys.stderr, __name__ + " module loaded"

