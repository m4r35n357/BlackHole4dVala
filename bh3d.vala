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

    public class Geodesic : GLib.Object, IModel {

        public double a;
        public double mu2;
        public double E;
        public double L;
        public double Q;
        public double r;
        public double th;
        public double starttime;
        public double endtime;
        public double h { get; set; }
        public double t;
        public double ph;
        public double eCum;
        public int count;
        public double a2;
        public double aE;
        public double a2E;
        public double L2;
        public double aL;
        public double a2xE2_mu2;
        public double[] cr;
        public double tP;
        public double rP;
        public double thP;
        public double phP;
        public double sth;
        public double cth;
        public double sth2;
        public double ra2;
		public double D;
		public double S;
		public double R;
		public double TH;
		public double THETA;
		public double eR;
		public double eTh;
		public double v4Cum;
		public double v4c;
		public double v4e;
        public ISymplectic integrator;

        public Geodesic (double bhMass, double spin, double pMass2, double energy, double momentum, double carter, double r0, double thetaMin, 
                         double starttime, double duration, double timestep, string type) {
            this.a = spin;
            this.mu2 = pMass2;
            this.E = energy;
            this.L = momentum;
            this.Q = carter;
            this.r = r0;
            this.th = thetaMin;
            this.starttime = starttime;
            this.endtime = starttime + fabs(duration);
            this.h = timestep;
            this.integrator = Integrator.getIntegrator(this, type);
			this.a2 = a * a;
		    this.aE = a * E;
		    this.a2E = a2 * E;
		    this.L2 = L * L;
		    this.aL = a * L;
		    var E2_mu2 = E * E - mu2;
		    this.cr = { E2_mu2, 2.0 * mu2, a2 * E2_mu2 - L2 - Q, 2.0 * ((aE - L) * (aE - L) + Q), - a2 * Q };
		    this.a2xE2_mu2 = - a2 * E2_mu2;
			refresh();
			this.rP = sqrt(fabs(R));
			this.thP = sqrt(fabs(THETA));
        }

		private double logError (double e) {
		    return 10.0 * log10(e > 1.0e-18 ? e : 1.0e-18);
		}

		private double modH (double xDot, double X) {
		    return 0.5 * fabs(xDot * xDot - X);
		}

		private double v4Error () {
		    var tmp1 = a * tP - ra2 * phP;
		    var tmp2 = tP - a * sth2 * phP;
		    return fabs(mu2 + sth2 / S * tmp1 * tmp1 + S / D * rP * rP + S * thP * thP - D / S * tmp2 * tmp2);
		}

		public void errors () {
		    eR = logError(modH(rP, R));
		    eTh = logError(modH(thP, THETA));
		    var error = v4Error();
		    v4e = logError(error);
		    v4Cum += error;
		    v4c = logError(v4Cum / count);
			count++;
		}

		public void refresh () {
            var r2 = r * r;
		    sth = sin(th);
		    cth = cos(th);
		    sth2 = sth * sth;
		    var cth2 = 1.0 - sth2;
		    ra2 = r2 + a2;
			D = ra2 - 2.0 * r;
			S = r2 + a2 * cth2;
		    R = (((cr[0] * r + cr[1]) * r + cr[2]) * r + cr[3]) * r + cr[4];
			TH = a2xE2_mu2 + L2 / sth2;
			THETA = Q - cth2 * TH;
			var P_D = (ra2 * E - aL) / D;
		    tP = ra2 * P_D + aL - a2E * sth2;
		    phP = a * P_D - aE + L / sth2;
        }

        public void pUp (double c) {
		    rP += c * (((4.0 * cr[0] * r + 3.0 * cr[1]) * r + 2.0 * cr[2]) * r + cr[3]) * 0.5;
            var cot = cth / sth;
		    thP += c * (cth * sth * TH + L2 * cot * cot * cot);
        }

        public void qUp (double d) {
		    t += d * tP;
		    r += d * rP;
		    th += d * thP;
		    ph += d * phP;
		    refresh();
        }

        public static Geodesic icJson () {
            var ic = getJson();
		    return new Geodesic(ic.get_double_member("M"),
                                ic.get_double_member("a"),
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
            stdout.printf("{\"mino\":%.9e, \"tau\":%.9e, \"v4e\":%.1f, \"v4c\":%.1f, \"ER\":%.1f, \"ETh\":%.1f ", mino, tau, v4e, v4c, eR, eTh);
            stdout.printf("{\"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, ", t, r, th, ph);
            stdout.printf("{\"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e}\n", tP / S, rP / S, thP / S, phP / S);
        }

		static int main (string[] args) {
		    Geodesic g = icJson();
			var mino = 0.0;
			var tau = 0.0;
			while (! (fabs(mino) > fabs(g.endtime))) {
				g.errors();
				if (fabs(mino) > fabs(g.starttime)) {
					g.output(mino, tau);
		        }
                g.integrator.compose();
				mino += g.h;
				tau += g.h * g.S;
		    }
        	return 0; 
    	}
	}
}

