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
from math import fabs
from json import loads

def main():
    print >> stderr, "Executable: {}".format(argv[0])
    if len(argv) != 3:
        raise Exception('>>> Please supply a time variable name (string) and a precision (float) <<<')
    timeVariable = argv[1]
    precision = float(argv[2])
    counter = 0
    previous = stdin.readline()
    previousJson = loads(previous)
    previousTime = float(previousJson[timeVariable])
    latest = stdin.readline()
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
        previousTime = latestTime
        latest = stdin.readline()

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"

