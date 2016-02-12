#!/usr/bin/env pypy

from sys import stdin, stdout, stderr
from math import fabs, log10, sqrt
from json import loads
from array import array

class Particle(object):

	def __init__(self, qX, qY, qZ, pX, pY, pZ, mass):
		self.qX = qX
		self.qY = qY
		self.qZ = qZ
		self.pX = pX
		self.pY = pY
		self.pZ = pZ
		self.mass = mass

	def __str__(self):
		return "{\"qX\":%.6e,\"qY\":%.6e,\"qZ\":%.6e,\"pX\":%.6e,\"pY\":%.6e,\"pZ\":%.6e,\"mass\":%.3e}" % (self.qX, self.qY, self.qZ, self.pX, self.pY, self.pZ, self.mass)

class Symplectic(object):
	def __init__(self, g, runTime, timeStep, errorLimit, bodies, order):
		self.bodies = bodies
		self.pRange = range(len(bodies))
		self.g = g
		self.ts = timeStep
		self.eMax = errorLimit
		self.T = runTime
		if order == 2:  # Second order
			self.coeff = array('d', [1.0])
		elif order == 4:  # Fourth order
			cbrt2 = 2.0 ** (1.0 / 3.0)
			y = 1.0 / (2.0 - cbrt2)
			self.coeff = array('d', [y,- y * cbrt2])
		elif order == 6:  # Sixth order
			self.coeff = array('d', [0.78451361047755726381949763,
											0.23557321335935813368479318,
											-1.17767998417887100694641568,
											1.31518632068391121888424973])
		elif order == 8:  # Eighth order
			self.coeff = array('d', [0.74167036435061295344822780,
											-0.40910082580003159399730010,
											0.19075471029623837995387626,
											-0.57386247111608226665638773,
											0.29906418130365592384446354,
											0.33462491824529818378495798,
											0.31529309239676659663205666,
											-0.79688793935291635401978884])
		elif order == 10:  # Tenth order
			self.coeff = array('d', [0.09040619368607278492161150,
											0.53591815953030120213784983,
											0.35123257547493978187517736,
											-0.31116802097815835426086544,
											-0.52556314194263510431065549,
											0.14447909410225247647345695,
											0.02983588609748235818064083,
											0.17786179923739805133592238,
											0.09826906939341637652532377,
											0.46179986210411860873242126,
											-0.33377845599881851314531820,
											0.07095684836524793621031152,
											0.23666960070126868771909819,
											-0.49725977950660985445028388,
											-0.30399616617237257346546356,
											0.05246957188100069574521612,
											0.44373380805019087955111365])
		else:  # Wrong value for integrator order
			raise Exception('>>> ERROR! Integrator order must be 2, 4, 6, 8 or 10 <<<')
        	self.coefficientsUp = range(len(self.coeff) - 1)
        	self.coefficientsDown = range(len(self.coeff) - 1, -1, -1)

	@staticmethod
	def dist (xA, yA, zA, xB, yB, zB):  # Euclidean distance between point A and point B
		return sqrt((xB - xA)**2 + (yB - yA)**2 + (zB - zA)**2)

	def h (self):  # Conserved energy
		energy = 0.0
		for i in self.pRange:
			a = self.bodies[i]
			energy += 0.5 * (a.pX**2 + a.pY**2 + a.pZ**2) / a.mass
			for j in self.pRange:
				if i > j:
					b = self.bodies[j]
					energy -= self.g * a.mass * b.mass / self.dist(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ)
		return energy

	def updateQ (self, c):  # Update Positions
		for i in self.pRange:
			a = self.bodies[i]
			tmp = c / a.mass * self.ts
			a.qX += a.pX * tmp
			a.qY += a.pY * tmp
			a.qZ += a.pZ * tmp

	def updateP (self, c):  # Update Momenta
		for i in self.pRange:
			a = self.bodies[i]
			for j in self.pRange:
				if i > j:
					b = self.bodies[j]
					tmp = - c * self.g * a.mass * b.mass * self.ts / self.dist(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ)**3
					dPx = tmp * (b.qX - a.qX)
					dPy = tmp * (b.qY - a.qY)
					dPz = tmp * (b.qZ - a.qZ)
					a.pX -= dPx
					a.pY -= dPy
					a.pZ -= dPz
					b.pX += dPx
					b.pY += dPy
					b.pZ += dPz

	def solve (self):  # Generalized Symplectic Integrator
		def sympBase (y):  # Compose higher orders from this symmetrical second-order symplectic base
			self.updateQ(0.5 * y)
			self.updateP(y)
			self.updateQ(0.5 * y)		
		for i in self.coefficientsUp:  # Composition happens in these loops
			sympBase(self.coeff[i])
		for i in self.coefficientsDown:
			sympBase(self.coeff[i])

	def print_out (self, time, hNow, h0, hMin, hMax, dbValue):
		data = []
		for i in self.pRange:
			data.append(str(self.bodies[i]))
		print >> stdout, "[" + ','.join(data) + "]"  # Log data
		print >> stderr, '{"t":%.2f, "H":%.9e, "H0":%.9e, "H-":%.9e, "H+":%.9e, "ER":%.1f}' % (time, hNow, h0, hMin, hMax, dbValue)  # Log progress

def icJson ():
	ic = loads(stdin.read())
	bodies = []
	for a in ic['bodies']:
		if 'pX' in a and 'pY' in a and 'pZ' in a:  # momenta specified
			bodies.append(Particle(a['qX'], a['qY'], a['qZ'], a['pX'], a['pY'], a['pZ'], a['mass']))
		elif 'vX' in a and 'vY' in a and 'vZ' in a:  # velocities specified, convert to momenta
			mass = a['mass']
			bodies.append(Particle(a['qX'], a['qY'], a['qZ'], mass * a['vX'], mass * a['vY'], mass * a['vZ'], mass))
		else:
			raise Exception('>>> ERROR! Specify either momenta or velocites consistently <<<')
	return Symplectic(ic['g'], ic['simulationTime'], ic['timeStep'], ic['errorLimit'], bodies, ic['integratorOrder'])

def main ():  # Need to be inside a function to return . . .
	s = icJson()  # Create a symplectic integrator object from JSON input
	h0 = hMax = hMin = s.h()  # Set up error reporting
	s.print_out(0.0, h0, h0, h0, h0, -180.0)
        t = 0.0
	while True:
		s.solve()  # Perform one full integration step
		hNow = s.h()		
		tmp = fabs(hNow - h0)  # Protect logarithm against negative arguments
		dH = tmp if tmp > 0.0 else 1.0e-18  # Protect logarithm against small arguments
		if hNow < hMin:  # Low tide
			hMin = hNow
		elif hNow > hMax:  # High tide
			hMax = hNow
		dbValue = 10.0 * log10(fabs(dH / h0) + 1.0e-18)
		s.print_out(t, hNow, h0, hMin, hMax, dbValue)
		if fabs(t) > s.T or dbValue > s.eMax:
			return
		t += s.ts

if __name__ == "__main__":
	main()
else:
	print >> stderr, __name__ + " module loaded"

