#!/usr/bin/env pypy

from sys import stdin, stdout, stderr
from math import fabs, log10, sqrt, sin, cos
from json import loads
from array import array

class BL(object):

    def __init__(self, mass, spin, energy, momentum, carter, r0, theta0, simtime, timestep, order):
    	self.m = mass
    	self.a = spin
        self.mu = 1.0
    	self.E = energy
    	self.L = momentum
    	self.Q = carter
    	self.r = r0
    	self.theta = theta0
    	self.phi = 0.0
    	self.T = simtime
    	self.step = timestep
        self.horizon = self.m * (1.0 + sqrt(1.0 - self.a**2))
        self.tau = self.mino = self.eCum = 0.0
        self.t = 0.0
        self.n = simtime / fabs(timestep)  # We can run backwards too!
        self.eMax = -00.0
        self.L_aE = self.L - self.a * self.E
 	self.pR = -sqrt(self.R(self.r) if self.R(self.r) >= 0.0 else 0.0)
	self.pTh = -sqrt(self.THETA(self.theta) if self.THETA(self.theta) >= 0.0 else 0.0)
        self.nf = 1.0e-18
	if order == 2:  # Second order
		self.coefficients = array('d', [1.0])
	elif order == 4:  # Fourth order
		cbrt2 = 2.0 ** (1.0 / 3.0)
		y = 1.0 / (2.0 - cbrt2)
		self.coefficients = array('d', [y,- y * cbrt2])
	elif order == 6:  # Sixth order
		self.coefficients = array('d', [0.78451361047755726381949763,
						0.23557321335935813368479318,
						-1.17767998417887100694641568,
						1.31518632068391121888424973])
	elif order == 8:  # Eighth order
		self.coefficients = array('d', [0.74167036435061295344822780,
						-0.40910082580003159399730010,
						0.19075471029623837995387626,
						-0.57386247111608226665638773,
						0.29906418130365592384446354,
						0.33462491824529818378495798,
						0.31529309239676659663205666,
						-0.79688793935291635401978884])
	elif order == 10:  # Tenth order
		self.coefficients = array('d', [0.09040619368607278492161150,
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

    def clamp (self, potential):
        return potential if potential >= 0.0 else 0.0

# intermediates
    def delta (self, r):
    	return r**2 - 2.0 * r * self.m + self.a**2

    def P1 (self, r):
    	return (r**2 + self.a**2) * self.E - self.a * self.L

    def P2 (self, r):
    	return self.Q + self.L_aE**2 + self.mu**2 * r**2

    def R (self, r):
    	return self.P1(r)**2 - self.delta(r) * self.P2(r)

    def TH (self, theta):
    	return self.a**2 * (self.mu**2 - self.E**2) + self.L**2 / sin(theta)**2

    def THETA (self, theta):
    	return self.Q - cos(theta)**2 * self.TH(theta)

# hamiltonian
    def h (self, r, theta, pR, pTh):
        e_r = abs(pR**2 - (self.clamp(self.R(self.r)))) / 2.0
        e_th = abs(pTh**2 - (self.clamp(self.THETA(self.theta)))) / 2.0
        self.eCum += e_r + e_th
        self.eR = 10.0 * log10(e_r if e_r >= self.nf else self.nf)
        self.eTh = 10.0 * log10(e_th if e_th >= self.nf else self.nf)
        self.e =  10.0 * log10(e_r + e_th if e_r + e_th >= self.nf else self.nf)
        return self.e

# parameters
    def tDeriv (self, t, r, theta, phi, pR, pTh):
        return (r**2 + self.a**2) * self.P1(r) / self.delta(r) + self.a * self.L + self.a**2 * self.E * sin(theta)**2

    def rDeriv (self, t, r, theta, phi, pR, pTh):
        return pR

    def thetaDeriv (self, t, r, theta, phi, pR, pTh):
        return pTh

    def phiDeriv (self, t, r, theta, phi, pR, pTh):
        return self.a * self.P1(r) / self.delta(r) - self.a * self.E + self.L / sin(theta)**2

# derivatives
    def rDotDeriv (self, t, r, theta, phi, pR, pTh):
        return 2.0 * r * self.E * self.P1(r) - (r - self.m) * self.P2(r) - self.mu**2 * r * self.delta(r)

    def thDotDeriv (self, t, r, theta, phi, pR, pTh):
        return cos(theta) * sin(theta) * self.TH(theta) + self.L**2 * cos(theta)**3 / sin(theta)**3

# Integrators
    def sympBase (self, y, t, r, theta, phi, pR, pTh):  # Compose higher orders from this symmetrical second-order symplectic base
	halfY = 0.5 * y
	self.updateQ(halfY, t, r, theta, phi, pR, pTh)
	self.updateP(y, t, r, theta, phi, pR, pTh)
	self.updateQ(halfY, t, r, theta, phi, pR, pTh)
			
    def solve (self, t, r, theta, phi, pR, pTh):  # Generalized Symplectic Integrator
	tmp = len(self.coefficients) - 1
	for i in range(tmp):  # Composition happens in these loops
		self.sympBase(self.coefficients[i], t, r, theta, phi, pR, pTh)
	for i in range(tmp, -1, -1):
		self.sympBase(self.coefficients[i], t, r, theta, phi, pR, pTh)

    def updateQ (self, c, t, r, theta, phi, pR, pTh):
        hstep = c * self.step
        self.t += hstep * self.tDeriv (t, r, theta, phi, pR, pTh)
        self.r += hstep * self.rDeriv (t, r, theta, phi, pR, pTh)
        self.theta += hstep * self.thetaDeriv (t, r, theta, phi, pR, pTh)
        self.phi += hstep * self.phiDeriv (t, r, theta, phi, pR, pTh)

    def updateP (self, c, t, r, theta, phi, pR, pTh):
        hstep = c * self.step
        self.pR += hstep * self.rDotDeriv (t, r, theta, phi, pR, pTh)
        self.pTh += hstep * self.thDotDeriv (t, r, theta, phi, pR, pTh)

    def euler (self, t, r, theta, phi, pR, pTh):
        self.t += self.step * self.tDeriv (t, r, theta, phi, pR, pTh)
        self.r += self.step * self.rDeriv (t, r, theta, phi, pR, pTh)
        self.theta += self.step * self.thetaDeriv (t, r, theta, phi, pR, pTh)
        self.phi += self.step * self.phiDeriv (t, r, theta, phi, pR, pTh)
        self.pR += self.step * self.rDotDeriv (t, r, theta, phi, pR, pTh)
        self.pTh += self.step * self.thDotDeriv (t, r, theta, phi, pR, pTh)

    def rk4 (self, t, r, theta, phi, pR, pTh):
        hstep = 0.5 * self.step
        a = array('d', [self.tDeriv (t, r, theta, phi, pR, pTh),
                 self.rDeriv (t, r, theta, phi, pR, pTh),
                 self.thetaDeriv (t, r, theta, phi, pR, pTh),
                 self.phiDeriv (t, r, theta, phi, pR, pTh),
                 self.rDotDeriv (t, r, theta, phi, pR, pTh),
                 self.thDotDeriv (t, r, theta, phi, pR, pTh)])
        b = array('d', [self.tDeriv (t + hstep * a[0], r + hstep * a[1], theta + hstep * a[2], phi + hstep * a[3], pR + hstep * a[4], pTh + hstep * a[5]),
                 self.rDeriv (t + hstep * a[0], r + hstep * a[1], theta + hstep * a[2], phi + hstep * a[3], pR + hstep * a[4], pTh + hstep * a[5]),
                 self.thetaDeriv (t + hstep * a[0], r + hstep * a[1], theta + hstep * a[2], phi + hstep * a[3], pR + hstep * a[4], pTh + hstep * a[5]),
                 self.phiDeriv (t + hstep * a[0], r + hstep * a[1], theta + hstep * a[2], phi + hstep * a[3], pR + hstep * a[4], pTh + hstep * a[5]),
                 self.rDotDeriv (t + hstep * a[0], r + hstep * a[1], theta + hstep * a[2], phi + hstep * a[3], pR + hstep * a[4], pTh + hstep * a[5]),
                 self.thDotDeriv (t + hstep * a[0], r + hstep * a[1], theta + hstep * a[2], phi + hstep * a[3], pR + hstep * a[4], pTh + hstep * a[5])])
        c = array('d', [self.tDeriv (t + hstep * b[0], r + hstep * b[1], theta + hstep * b[2], phi + hstep * b[3], pR + hstep * b[4], pTh + hstep * b[5]),
                 self.rDeriv (t + hstep * b[0], r + hstep * b[1], theta + hstep * b[2], phi + hstep * b[3], pR + hstep * b[4], pTh + hstep * b[5]),
                 self.thetaDeriv (t + hstep * b[0], r + hstep * b[1], theta + hstep * b[2], phi + hstep * b[3], pR + hstep * b[4], pTh + hstep * b[5]),
                 self.phiDeriv (t + hstep * b[0], r + hstep * b[1], theta + hstep * b[2], phi + hstep * b[3], pR + hstep * b[4], pTh + hstep * b[5]),
                 self.rDotDeriv (t + hstep * b[0], r + hstep * b[1], theta + hstep * b[2], phi + hstep * b[3], pR + hstep * b[4], pTh + hstep * b[5]),
                 self.thDotDeriv (t + hstep * b[0], r + hstep * b[1], theta + hstep * b[2], phi + hstep * b[3], pR + hstep * b[4], pTh + hstep * b[5])])
        d = array('d', [self.tDeriv (t + self.step * c[0], r + self.step * c[1], theta + self.step * c[2], phi + self.step * c[3], pR + self.step * c[4], pTh + self.step * c[5]),
                 self.rDeriv (t + self.step * c[0], r + self.step * c[1], theta + self.step * c[2], phi + self.step * c[3], pR + self.step * c[4], pTh + self.step * c[5]),
                 self.thetaDeriv (t + self.step * c[0], r + self.step * c[1], theta + self.step * c[2], phi + self.step * c[3], pR + self.step * c[4], pTh + self.step * c[5]),
                 self.phiDeriv (t + self.step * c[0], r + self.step * c[1], theta + self.step * c[2], phi + self.step * c[3], pR + self.step * c[4], pTh + self.step * c[5]),
                 self.rDotDeriv (t + self.step * c[0], r + self.step * c[1], theta + self.step * c[2], phi + self.step * c[3], pR + self.step * c[4], pTh + self.step * c[5]),
                 self.thDotDeriv (t + self.step * c[0], r + self.step * c[1], theta + self.step * c[2], phi + self.step * c[3], pR + self.step * c[4], pTh + self.step * c[5])])
        self.t += self.step * (a[0] + 2.0 * b[0] + 2.0 * c[0] + d[0]) / 6.0
        self.r += self.step * (a[1] + 2.0 * b[1] + 2.0 * c[1] + d[1]) / 6.0
        self.theta += self.step * (a[2] + 2.0 * b[2] + 2.0 * c[2] + d[2]) / 6.0
        self.phi += self.step * (a[3] + 2.0 * b[3] + 2.0 * c[3] + d[3]) / 6.0
        self.pR += self.step * (a[4] + 2.0 * b[4] + 2.0 * c[4] + d[4]) / 6.0
        self.pTh += self.step * (a[5] + 2.0 * b[5] + 2.0 * c[5] + d[5]) / 6.0

# parse input
def icJson ():
	ic = loads(stdin.read())
        return BL(ic['M'], ic['a'], ic['E'], ic['Lz'], ic['C'], ic['r'], ic['theta'], ic['time'], ic['step'], ic['integratorOrder'])

def main ():  # Need to be inside a function to return . . .
    bh = icJson()
    h0 = hMax = hMin = bh.h(bh.r,  bh.theta, bh.pR, bh.pTh)  # Set up error reporting
    while True:
        ra = sqrt(bh.r**2 + bh.a**2)
        x = ra * sin(bh.theta) * cos(bh.phi)
        y = ra * sin(bh.theta) * sin(bh.phi)
        z = bh.r * cos(bh.theta)
	hNow = bh.h(bh.r, bh.theta, bh.pR, bh.pTh)		
	dbValue = hNow
	print >> stdout, '{"mino":%.9e, "tau":%.9e, "E":%.1f, "ER":%.1f, "ETh":%.1f, "EC":%.1f, "t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "x":%.9e, "y":%.9e, "z":%.9e}' % (bh.mino, bh.tau, dbValue, bh.eR, bh.eTh, 10.0 * log10(abs(bh.eCum) + 1.0e-18), bh.t, bh.r, bh.theta, bh.phi, x, y, z)  # Log data
        bh.euler(bh.t, bh.r, bh.theta, bh.phi, bh.pR, bh.pTh)
#        bh.rk4(bh.t, bh.r, bh.theta, bh.phi, bh.pR, bh.pTh)
#        bh.solve(bh.t, bh.r, bh.theta, bh.phi, bh.pR, bh.pTh)
	if abs(bh.mino) > bh.T or bh.eCum > 1.0e-0 or bh.r < bh.horizon:
	    break
        bh.mino += bh.step
        bh.tau += bh.step * (bh.r**2 + bh.a**2 * cos(bh.theta)**2)

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"

