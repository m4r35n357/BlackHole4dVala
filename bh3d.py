#!/usr/bin/env pypy

from sys import stdin, stdout, stderr
from math import fabs, log10, sqrt, sin, cos
from json import loads
from array import array

class BL(object):

    def __init__(self, mass, spin, energy, momentum, carter, r0, theta0, phi0, simtime, timestep):
    	self.m = mass
    	self.a = spin
        self.mu = 1.0
    	self.E = energy
    	self.L = momentum
    	self.Q = carter
    	self.r = r0
    	self.theta = theta0
    	self.phi = phi0
    	self.time = simtime
    	self.step = timestep
        self.tau = 0.0
        self.t = 0.0
        self.n = simtime / fabs(timestep)  # We can run backwards too!
        self.eMax = -30.0
        self.L_aE = self.L - self.a * self.E
        self.pR = 0.0
        self.pTh = sqrt(self.Q)

# intermediates
    def sigma (self, r, theta):
#    	return r**2 + self.a**2 * cos(theta)**2
    	return 1.0

    def delta (self, r):
    	return r**2 - 2.0 * r * self.m + self.a**2

    def deltaDash (self, r):
    	return 2.0 * (r - self.m)

    def P1 (self, r):
    	return (r**2 + self.a**2) * self.E - self.a * self.L

    def P2 (self, r):
    	return self.Q + self.L_aE**2 + self.mu**2 * r**2

    def R (self, r):
    	return self.P1(r)**2 - self.delta(r) * self.P2(r)

    def RDash (self, r):
    	return - self.deltaDash(r) * self.P2(r) + 4.0 * r * self. E * self.P1(r) - 2.0 * r * self.delta(r)

    def T1 (self, theta):
    	return self.L**2 / sin(theta)**2 + self.a**2 * (self.mu**2 - self.E**2) 

    def THETA (self, theta):
    	return self.Q - cos(theta)**2 * self.T1(theta)

# hamiltonian
    def h (self, r, theta):
        return 0.5 * (self.delta(r) * self.pR**2 + self.pTh**2 - self.R(r) / self.delta(r) - self.THETA(theta)) / self.sigma(r, theta) - 0.5

# parameters
    def tDeriv (self, t, r, theta, phi, pR, pTh):
        return - 0.5 * ((r**2 + self.a**2) * self.P1(r) / self.delta(r) + self.a * self.L_aE + cos(theta)**2 * self.a**2 * self.E) / self.sigma(r, theta)

    def rDeriv (self, t, r, theta, phi, pR, pTh):
        return pR * self.delta(r) / self.sigma(r, theta) 

    def thetaDeriv (self, t, r, theta, phi, pR, pTh):
        return pTh / self.sigma(r, theta) 

    def phiDeriv (self, t, r, theta, phi, pR, pTh):
        return - 0.5 * (self.a * self.P1(r) / self.delta(r) - 2.0 * self.a * self.L_aE + self.L * cos(theta)**2 / sin(theta)**2) / self.sigma(r, theta)

# derivatives
    def rDotDeriv (self, t, r, theta, phi, pR, pTh):
        return 0.5 * (self.deltaDash(r) * self.R(r) / self.delta(r)**2 - self.RDash(r) / self.delta(r) - self.deltaDash(r) * pR**2) / self.sigma(r, theta)

    def thDotDeriv (self, t, r, theta, phi, pR, pTh):
        return (cos(theta) * sin(theta) * self.T1(theta) + self.L**2 * cos(theta)**3 / sin(theta)**3) / self.sigma(r, theta) 

# Integrators
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
        b = array('d', [self.tDeriv (t + hstep * a[0], r + 0.5 * a[1], theta + hstep * a[2], phi + hstep * a[3], pR + hstep * a[4], pTh + hstep * a[5]),
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
        return BL(ic['M'], ic['a'], ic['E'], ic['Lz'], ic['C'], ic['r'], ic['theta'], ic['phi'], ic['time'], ic['step'])

def main ():  # Need to be inside a function to return . . .
    bh = icJson()
    h0 = hMax = hMin = bh.h(bh.r,  bh.theta)  # Set up error reporting
    n = 1
    while n <= bh.n:
        ra = sqrt(bh.r**2 + bh.a**2)
        x = ra * sin(bh.theta) * cos(bh.phi)
        y = ra * sin(bh.theta) * sin(bh.phi)
        z = bh.r * cos(bh.theta)
	hNow = bh.h(bh.r, bh.theta)		
	tmp = fabs(hNow - h0)  # Protect logarithm against negative arguments
	dH = tmp if tmp > 0.0 else 1.0e-18  # Protect logarithm against small arguments
	if hNow < hMin:  # Low tide
		hMin = hNow
	elif hNow > hMax:  # High tide
		hMax = hNow
	dbValue = 10.0 * log10(dH)
	print >> stdout, '{"tau":%.9e, "H":%.9e, "H0":%.9e, "H-":%.9e, "H+":%.9e, "ER":%.1f, "t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "x":%.9e, "y":%.9e, "z":%.9e}' % (-bh.tau, hNow, h0, hMin, hMax, dbValue, bh.t, bh.r, bh.theta, bh.phi, x, y, z)  # Log data
	print >> stderr, '{"tau":%.9f, "H":%.9e, "H0":%.9e, "H-":%.9e, "H+":%.9e, "ER":%.1f}' % (bh.tau, hNow, h0, hMin, hMax, dbValue)  # Log progress
        bh.euler(bh.t, bh.r, bh.theta, bh.phi, bh.pR, bh.pTh)
#        bh.rk4(bh.t, bh.r, bh.theta, bh.phi, bh.pR, bh.pTh)
	if dbValue > bh.eMax:
		return
	n += 1
        bh.tau += bh.step

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"

