/*
Copyright (c) 2014, 2015, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

using GLib.Math;

namespace Kerr {

    public class Particle : GLib.Object {

        public double qX;
        public double qY;
        public double qZ;
        public double pX;
        public double pY;
        public double pZ;
        public double mass;

        public Particle (double qX, double qY, double qZ, double pX, double pY, double pZ, double mass) {
            this.qX = qX;
            this.qY = qY;
            this.qZ = qZ;
            this.pX = pX;
            this.pY = pY;
            this.pZ = pZ;
            this.mass = mass;
        }

        public void output (double mino, double tau) {
            stdout.printf("{\"qX\":%.6e,\"qY\":%.6e,\"qZ\":%.6e,\"pX\":%.6e,\"pY\":%.6e,\"pZ\":%.6e,\"mass\":%.3e}\n", qX, qY, qZ, pX, pY, pZ, mass);
        }
    }

    public class System : GLib.Object, IModel {

        public double[] bodies;
        public double pRange;
        public double g;
        public double pX;
        public double pY;
        public double pZ;
        public double mass;
        public double h { get; set; }
        public Integrator integrator;
        public double T;

        public System (double[] bodies, double g, double timeStep, double errorLimit, double runTime, string type) {
			this.bodies = bodies;
			this.pRange = bodies.length;
			this.g = g;
            this.h = timeStep;
			this.T = runTime;
            this.integrator = Integrator.getIntegrator(this, type);
        }

        public void pUp (double c) {
        }

        public void qUp (double d) {
        }

        public static System icJson () {
            var ic = getJson();
            foreach (var node in ic.get_array_member("bodies").get_elements()) {
               
            }
/*
		    return new System(ic.get_array_member("bodies"), ic.get_member("a").get_double(), ic.get_member("mu").get_double(), ic.get_member("E").get_double(), ic.get_member("Lz").get_double(), ic.get_member("C").get_double(), ic.get_member("r").get_double(), ic.get_member("theta").get_double(), ic.get_member("start").get_double(), ic.get_member("duration").get_double(), ic.get_member("step").get_double(), ic.get_member("integrator").get_string());
*/
            return new System({ 0.0 }, 0.0, 0.0, 0.0, 0.0, "");
        }

        public void output (double mino, double tau) {
/*
            stdout.printf("{\"mino\":%.9e, \"tau\":%.9e, \"v4e\":%.1f, \"v4c\":%.1f, \"ER\":%.1f, \"ETh\":%.1f, \"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, \"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e}\n", mino, tau, v4e, v4c, eR, eTh, t, r, th, ph, tP / S, rP / S, thP / S, phP / S);
*/
        }

		static int main (string[] args) {
		    System s = icJson();
/*
			var mino = 0.0;
			var tau = 0.0;
			while (! (fabs(mino) > g.endtime)) {
				g.errors();
				if (fabs(mino) > g.starttime) {
					g.output(mino, tau);
		        }
                g.integrator.compose();
				mino += g.h;
				tau += g.h * g.S;
		    }
*/
        	return 0; 
    	}
	}
}
/*
class Particle(object):

	def __init__(self, qX, qY, qZ, pX, pY, pZ, mass):
		this.qX = qX
		this.qY = qY
		this.qZ = qZ
		this.pX = pX
		this.pY = pY
		this.pZ = pZ
		this.mass = mass

	def __str__(self):
		return "{\"qX\":%.6e,\"qY\":%.6e,\"qZ\":%.6e,\"pX\":%.6e,\"pY\":%.6e,\"pZ\":%.6e,\"mass\":%.3e}" % (this.qX, this.qY, this.qZ, this.pX, this.pY, this.pZ, this.mass)

class Symplectic(object):
	def __init__(self, g, runTime, timeStep, errorLimit, bodies, order):
		this.bodies = bodies
		this.pRange = range(len(bodies))
		this.g = g
		this.ts = timeStep
		this.eMax = errorLimit
		this.T = runTime
		if order == 2:  # Second order
			this.coeff = array('d', [1.0])
		elif order == 4:  # Fourth order
			cbrt2 = 2.0 ** (1.0 / 3.0)
			y = 1.0 / (2.0 - cbrt2)
			this.coeff = array('d', [y,- y * cbrt2])
		elif order == 6:  # Sixth order
			this.coeff = array('d', [0.78451361047755726381949763,
											0.23557321335935813368479318,
											-1.17767998417887100694641568,
											1.31518632068391121888424973])
		elif order == 8:  # Eighth order
			this.coeff = array('d', [0.74167036435061295344822780,
											-0.40910082580003159399730010,
											0.19075471029623837995387626,
											-0.57386247111608226665638773,
											0.29906418130365592384446354,
											0.33462491824529818378495798,
											0.31529309239676659663205666,
											-0.79688793935291635401978884])
		elif order == 10:  # Tenth order
			this.coeff = array('d', [0.09040619368607278492161150,
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
        	this.coefficientsUp = range(len(this.coeff) - 1)
        	this.coefficientsDown = range(len(this.coeff) - 1, -1, -1)

	@staticmethod
	def dist (xA, yA, zA, xB, yB, zB):  # Euclidean distance between point A and point B
		return sqrt((xB - xA)**2 + (yB - yA)**2 + (zB - zA)**2)

	def h (self):  # Conserved energy
		energy = 0.0
		for i in this.pRange:
			a = this.bodies[i]
			energy += 0.5 * (a.pX**2 + a.pY**2 + a.pZ**2) / a.mass
			for j in this.pRange:
				if i > j:
					b = this.bodies[j]
					energy -= this.g * a.mass * b.mass / this.dist(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ)
		return energy

	def updateQ (self, c):  # Update Positions
		for i in this.pRange:
			a = this.bodies[i]
			tmp = c / a.mass * this.ts
			a.qX += a.pX * tmp
			a.qY += a.pY * tmp
			a.qZ += a.pZ * tmp

	def updateP (self, c):  # Update Momenta
		for i in this.pRange:
			a = this.bodies[i]
			for j in this.pRange:
				if i > j:
					b = this.bodies[j]
					tmp = - c * this.g * a.mass * b.mass * this.ts / this.dist(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ)**3
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
			this.updateQ(0.5 * y)
			this.updateP(y)
			this.updateQ(0.5 * y)		
		for i in this.coefficientsUp:  # Composition happens in these loops
			sympBase(this.coeff[i])
		for i in this.coefficientsDown:
			sympBase(this.coeff[i])

	def print_out (self, time, hNow, h0, hMin, hMax, dbValue):
		data = []
		for i in this.pRange:
			data.append(str(this.bodies[i]))
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
*/
