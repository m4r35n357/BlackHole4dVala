#!/usr/bin/env python

from sys import argv, stdin
from math import fabs
from json import loads

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
        if (latestTime - target) > 0.0:
            counter += 1
            if fabs(latestTime - target) <= fabs(previousTime - target):
                print latest,  # trailing comma suppresses newline
            else:
                print previous,
        previous = latest
        previousJson = latestJson
        previousTime = latestTime
        latest = dataFile.readline()

if __name__ == "__main__":
    main()
else:
    print >> sys.stderr, __name__ + " module loaded"

