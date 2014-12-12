#!/usr/bin/env pypy
'''
Copyright (c) 2014, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
'''
from sys import stdin, stdout, stderr
from math import fabs, log10, sqrt, sin, cos, pi
from json import loads
from array import array

class BL(object):   # Boyer-Lindquist coordinates on the Kerr metric
    def __init__(self, mass, spin, pmass, energy, momentum, carter, r0, theta0, simtime, timestep, order):
    	self.m = 1.0
    	self.a = spin
    	self.a2 = self.a**2
        self.mu2 = pmass**2
    	self.E = energy
        self.E2 = self.E**2
        self.aE = self.a * self.E
        self.a2E = self.a2 * self.E
        self.a2mu2_E2 = self.a2 * (self.mu2 - self.E2)
    	self.L = momentum
        self.L2 = self.L**2
        self.aL = self.a * self.L
        self.L_aE2 = (self.L - self.aE)**2
    	self.Q = carter
    	self.r = r0
    	self.th = theta0
    	self.time = simtime
    	self.h = timestep
        self.T = abs(simtime)
        self.t = self.ph = self.mino = self.tau = self.eCum = 0.0
        self.nf = 1.0e-18
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

    def clamp (self, potential):
        return potential if potential >= 0.0 else 0.0

    def updatePotentials (self):  # Intermediate parameters
        self.sth = sin(self.th)
        self.cth = cos(self.th)
        self.sth2 = self.sth**2
	self.delta = self.r**2 - 2.0 * self.r * self.m + self.a2
	self.P1 = (self.r**2 + self.a2) * self.E - self.aL
	self.P2 = self.Q + self.L_aE2 + self.mu2 * self.r**2
	self.R = self.P1**2 - self.delta * self.P2
	self.TH = self.a2mu2_E2 + (self.L / self.sth)**2
	self.THETA = self.Q - self.cth**2 * self.TH
	
    def errors (self):  # Error analysis
        e_r = abs(self.vR**2 - self.clamp(self.R)) / 2.0
        e_th = abs(self.vTh**2 - self.clamp(self.THETA)) / 2.0
        self.eR = 10.0 * log10(e_r if e_r >= self.nf else self.nf)
        self.eTh = 10.0 * log10(e_th if e_th >= self.nf else self.nf)
        self.e = 10.0 * log10(e_r + e_th if e_r + e_th >= self.nf else self.nf)
        self.eCum += e_r + e_th	

    def update_t_phi_Dot (self):  # t and phi updates
        self.tDot = (self.r**2 + self.a2) * self.P1 / self.delta + self.aL - self.a2E * self.sth2
        self.phDot = self.a * self.P1 / self.delta - self.aE + self.L / self.sth2

    def update_t_phi (self):  # t and phi updates
        self.update_t_phi_Dot()
        self.t += self.h * self.tDot
        self.ph += self.h * self.phDot

    def qUp (self, c):  # r and theta updates
        self.r += c * self.h * self.vR
        self.th += c * self.h * self.vTh
        self.updatePotentials()

    def qDotUp (self, c):  # Velocity updates
        self.vR += c * self.h * (2.0 * self.r * self.E * self.P1 - (self.r - self.m) * self.P2 - self.mu2 * self.r * self.delta)
        self.vTh += c * self.h * (self.cth * self.sth * self.TH + self.L2 * (self.cth / self.sth)**3)

    def solve (self):  # Generalized Symplectic Integrator
        def sv (y):  # Compose higher orders from this symmetrical second-order symplectic base
	    self.qUp(0.5 * y)
	    self.qDotUp(y)
	    self.qUp(0.5 * y)
	for i in self.coefficientsUp:  # Composition happens in these loops
	    sv(self.coeff[i])
	for i in self.coefficientsDown:
	    sv(self.coeff[i])

def main ():  # Need to be inside a function to return . . .
    ic = loads(stdin.read())
    bl = BL(ic['M'], ic['a'], ic['mu'], ic['E'], ic['Lz'], ic['C'], ic['r'], ic['theta'], ic['time'], ic['step'], ic['integratorOrder'])
    bl.updatePotentials()
    bl.vR = - sqrt(bl.clamp(bl.R))
    bl.vTh = - sqrt(bl.clamp(bl.THETA))
    bl.update_t_phi_Dot()
    while True:
        bl.errors()
        ra = sqrt(bl.r**2 + bl.a2)
	print >> stdout, '{"mino":%.9e, "tau":%.9e, "E":%.1f, "ER":%.1f, "ETh":%.1f, "EC":%.1f, "t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "tDot":%.9e, "rDot":%.9e, "thDot":%.9e, "phDot":%.9e, "x":%.9e, "y":%.9e, "z":%.9e}' % (bl.mino, bl.tau, bl.e, bl.eR, bl.eTh, 10.0 * log10(bl.eCum if bl.eCum >= bl.nf else bl.nf), bl.t, bl.r, bl.th, bl.ph, bl.tDot, bl.vR, bl.vTh, bl.phDot, ra * bl.sth * cos(bl.ph), ra * bl.sth * sin(bl.ph), bl.r * bl.cth)  # Log data
        bl.update_t_phi()  # Euler's method
        bl.solve()  # update r and theta with symplectic integrator
	if abs(bl.mino) > bl.T or bl.eCum > 1.0e-0:
	    break
        bl.mino += bl.h
        bl.tau += bl.h * (bl.r**2 + bl.a2 * cos(bl.th)**2)

if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"

