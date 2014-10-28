#!/usr/bin/env pypy

from sys import stdin, stdout, stderr
from math import fabs, log10, sqrt, sin, cos, pi
from json import loads
from array import array

class BL(object):   # Boyer-Lindquist coordinates on the Kerr metric
    def __init__(self, mass, spin, pmass, energy, momentum, carter, r0, theta0, simtime, timestep, order):
    	self.m = 1.0
    	self.a = spin
        self.mu = pmass * pmass
    	self.E = energy
    	self.L = momentum
    	self.Q = carter
        self.t = 0.0
    	self.r = r0
    	self.theta = theta0
    	self.phi = 0.0
    	self.time = simtime
    	self.step = timestep
        self.mino = 0.0
        self.n = simtime / fabs(timestep)  # We can run backwards too!
        self.L_aE2 = (self.L - self.a * self.E)**2
        self.horizon = self.m * (1.0 + sqrt(1.0 - self.a**2))
        self.eCum = 0.0
	if order == 2:  # Second order
		self.coeff = array('d', [1.0])
	elif order == 4:  # Fourth order
		cbrt2 = 2.0 ** (1.0 / 3.0)
		y = 1.0 / (2.0 - cbrt2)
		self.coeff = array('d', [y, - y * cbrt2])
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

    def updatePotentials (self):  # Intermediate parameters
	self.delta = self.r**2 - 2.0 * self.r * self.m + self.a**2
	self.P1 = (self.r**2 + self.a**2) * self.E - self.a * self.L
	self.P2 = self.Q + self.L_aE2 + self.mu**2 * self.r**2
	self.R = self.P1**2 - self.delta * self.P2
#        self.R = self.R if self.R >= 0.0 else 0.0
	self.TH = self.a**2 * (self.mu**2 - self.E**2) + (self.L / sin(self.theta))**2
	self.THETA = self.Q - cos(self.theta)**2 * self.TH
#        self.THETA = self.THETA if self.THETA >= 0.0 else 0.0
	
    def errors (self):  # Error analysis
        e_r = abs(self.vR**2 - self.R) / 2.0
        e_th = abs(self.vTh**2 - self.THETA) / 2.0
        self.eR = 10.0 * log10(e_r + 1.0e-18)
        self.eTh = 10.0 * log10(e_th + 1.0e-18)
        self.e =  10.0 * log10(e_r + e_th + 1.0e-18)
        self.eCum += e_r + e_th	

    def update_t_phi (self):  # t and phi updates
        self.t += self.step * ((self.r**2 + self.a**2) * self.P1 / self.delta - self.a * (self.a * self.E * sin(self.theta)**2 - self.L))
        self.phi += self.step * (self.a * self.P1 / self.delta - self.a * self.E + self.L / sin(self.theta)**2)

    def qUpdate (self, c):  # r and theta updates
        self.r += c * self.step * self.vR
        self.theta += c * self.step * self.vTh
        self.updatePotentials()

    def qDotUpdate (self, c):  # Velocity updates
        self.vR += c * self.step * (2.0 * self.r * self.E * self.P1 - (self.r - self.m) * self.P2 - self.mu**2 * self.r * self.delta)
        self.vTh += c * self.step * (cos(self.theta) * sin(self.theta) * self.TH + self.L**2 * (cos(self.theta) / sin(self.theta))**3)

    def stormerVerlet (self, y):  # Compose higher orders from this symmetrical second-order symplectic base
	self.qUpdate(0.5 * y)
	self.qDotUpdate(y)
	self.qUpdate(0.5 * y)
			
    def solve (self):  # Generalized Symplectic Integrator
	for i in self.coefficientsUp:  # Composition happens in these loops
	    self.stormerVerlet(self.coeff[i])
	for i in self.coefficientsDown:
	    self.stormerVerlet(self.coeff[i])

def main ():  # Need to be inside a function to return . . .
    ic = loads(stdin.read())
    bl = BL(ic['M'], ic['a'], ic['mu'], ic['E'], ic['Lz'], ic['C'], ic['r'], ic['theta'], ic['time'], ic['step'], ic['integratorOrder'])
    bl.updatePotentials()
    bl.vR = -sqrt(bl.R if bl.R >= 0.0 else 0.0)
    bl.vTh = -sqrt(bl.THETA if bl.THETA >= 0.0 else 0.0)
    n = 1
    while n <= bl.n:
        bl.errors()
        ra = sqrt(bl.r**2 + bl.a**2)
	print >> stdout, '{"mino":%.9e, "tau":%.9e, "E":%.1f, "ER":%.1f, "ETh":%.1f, "EC":%.1f, "t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "R":%.9e, "THETA":%.9e, "x":%.9e, "y":%.9e, "z":%.9e}' % (bl.mino, bl.mino * (bl.r**2 + bl.a**2 * cos(bl.theta)**2), bl.e, bl.eR, bl.eTh, 10.0 * log10(bl.eCum + 1.0e-18), bl.t, bl.r, bl.theta, bl.phi, bl.R, bl.THETA, ra * sin(bl.theta) * cos(bl.phi), ra * sin(bl.theta) * sin(bl.phi), bl.r * cos(bl.theta))  # Log data
        bl.update_t_phi()  # Euler's method
        bl.solve()  # update r and theta with symplectic integrator
	if bl.eCum > 1.0e-0 or bl.r < bl.horizon:
	    break
        bl.mino += bl.step
	n += 1

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"

