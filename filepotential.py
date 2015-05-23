#!/usr/bin/env python

from sys import argv
from matplotlib import pyplot
from json import loads

def main():
	if len(argv) < 2:
		raise Exception('>>> ERROR! Please supply a potential file name <<<')
	dataFile = open(argv[1], 'r')
	line = dataFile.readline()
	ax1 = pyplot.figure().add_subplot(111)
        pyplot.grid(b=True, which='major', color='k', linestyle='-')
	ax1.set_xlabel('r, theta')
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
		line = dataFile.readline()
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