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

    public class Orbit : GLib.Object, IModel {

        private double E;
        private double L;
        private double r;
        private double th;
        public double starttime;
        public double endtime;
        public double h { get; set; }
        private double ph;
        private double L2;
        private double rP;
        private double phP;
		private double R;
		private double H0;
		private double eR;
		private double eTh;
        private ISymplectic integrator;

        public Orbit (double spin, double pMass2, double energy, double momentum, double carter, double r0, double thetaMin, 
                         double starttime, double duration, double timestep, string type) {
            this.E = energy;
            this.L = momentum;
            this.r = r0;
            this.th = thetaMin;
            this.starttime = starttime;
            this.endtime = starttime + duration;
            this.h = timestep;
            this.integrator = Integrator.getIntegrator(this, type);
		    this.L2 = L * L;
			refresh();
            this.H0 = H(rP, R);
			//this.rP = - sqrt(fabs(E - 1.0 - R));
        }

		private double logError (double e) {
		    return 10.0 * log10(e > 1.0e-18 ? e : 1.0e-18);
		}

		private double H (double xDot, double X) {
		    return 0.5 * (xDot * xDot + X);
		}

		public void errors () {
		    eR = logError(H(rP, R));
		}

		private void refresh () {
            R = - 1.0 / r + L2 / (2.0 * r * r);
		    phP = L / (r * r);
//            stderr.printf("{\"L\":%.9e}\n", L);
//            stderr.printf("{\"r\":%.9e}\n", r);
//            stderr.printf("{\"phP\":%.9e}\n", phP);
        }

        public void pUp (double c) {
		    rP += c * (L2 / (r * r * r) - 1.0 / (r * r));
        }

        public void qUp (double d) {
		    r += d * rP;
		    ph += d * L / (r * r);
		    refresh();
        }

        public void evolve () {
            integrator.compose();
        }

        public void solve () {
			var mino = 0.0;
			while (! (mino > endtime)) {
				errors();
				if (mino > starttime) {
					output(mino, 0.0);
			    }
		        evolve();
				mino += h;
			}
        }

        public static Orbit fromJson () {
            var ic = getJson();
		    return new Orbit(ic.get_double_member("a"),
                                ic.get_double_member("mu"),
                                ic.get_double_member("E"),
                                ic.get_double_member("Lz"),
                                ic.get_double_member("C"),
                                ic.get_double_member("r"),
                                ic.get_double_member("theta"),
                                ic.get_double_member("start"),
                                ic.get_double_member("duration"),
                                ic.get_double_member("step"),
                                ic.get_string_member("integrator"));
        }

        public void output (double mino, double tau) {
            //stderr.printf("{\"phP\":%.9e}\n", phP);
            stdout.printf("{\"mino\":%.9e, \"tau\":%.9e, \"v4e\":%.1f, \"v4c\":%.1f, \"ER\":%.1f, \"ETh\":%.1f, ", mino, mino, H(rP, R), 0.0, eR, eTh);
            stdout.printf("\"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, ", mino, r, th, ph);
            stdout.printf("\"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e}\n", 0.0, rP, 0.0, phP);
        }
	}

	static int main (string[] args) {
	    Orbit.fromJson().solve();
    	return 0; 
	}
}

