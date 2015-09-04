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
using Json;

namespace BH.BL {

    public string read_stdin () {
        var input = new StringBuilder ();
        var buffer = new char[1024];
        while (!stdin.eof ()) {
            string read_chunk = stdin.gets (buffer);
            if (read_chunk != null) {
                input.append (read_chunk);
            }
        }
        return input.str;
    }

    public class Geodesic : GLib.Object {

        public double M;
        public double a;
        public double mu2;
        public double E;
        public double L;
        public double Q;
        public double r;
        public double th;
        public double starttime;
        public double endtime;
        public double h;
        public string order;
        public double cbrt2;
        public double f2;
        public double[] coefficients;
        public double[] w;
        public int wRange;
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

        public Geodesic (double bhMass, double spin, double pMass2, double energy, double momentum, double carter, double r0, double thetaMin, double starttime, double duration, double timestep, string order) {
            this.M = bhMass;
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
            this.order = order;
            this.cbrt2 = pow(2.0, (1.0 / 3.0));
            this.f2 = 1.0 / (2.0 - cbrt2);
            this.coefficients = { 0.5 * f2, f2, 0.5 * (1.0 - cbrt2) * f2, - cbrt2 * f2 };
            switch (order) {
                case "sb2":
                    w = { 1.0 };
                    break;
                case "sc4":
                    w = { coefficients[1], coefficients[3], coefficients[1] };
                    break;
                case "sb4":
                    w = { 1.0 };
                    break;
                case "sc6":
                    var fthrt2 = pow(2.0, (1.0 / 5.0));
                    w = { 1.0 / (2.0 - fthrt2), - fthrt2 / (2.0 - fthrt2), 1.0 / (2.0 - fthrt2) };
                    break;
                default:
				    stderr.printf ("Unknown integrator type: %s\n", order);
                    break;
            }
            this.wRange = w.length;
            this.t = 0.0;
            this.ph = 0.0;
            this.eCum = 0.0;
		    this.count = 0;
			this.a2 = a * a;
		    this.aE = a * E;
		    this.a2E = a2 * E;
		    this.L2 = L * L;
		    this.aL = a * L;
		    var E2_mu2 = E * E - mu2;
		    this.cr = { E2_mu2, 2.0 * mu2, a2 * E2_mu2 - L2 - Q, 2.0 * ((aE - L) * (aE - L) + Q), - a2 * Q };
		    this.a2xE2_mu2 = - a2 * E2_mu2;
        }

		private double logError (double e) {
		    return 10.0 * log10(e > 1.0e-18 ? e : 1.0e-18);
		}

		private double modH (double xDot, double X) {
		    return 0.5 * fabs(xDot * xDot - X);
		}

		private double v4Error (double tP, double rP, double thP, double phP) {
		    double tmp1 = a * tP - ra2 * phP;
		    double tmp2 = tP - a * sth2 * phP;
		    return fabs(mu2 + sth2 / S * tmp1 * tmp1 + S / D * rP * rP + S * thP * thP - D / S * tmp2 * tmp2);
		}

		public void errors (double R, double THETA, double tP, double rP, double thP, double phP) {
		    eR = logError(modH(rP, R));
		    eTh = logError(modH(thP, THETA));
		    double error = v4Error(tP / S, rP / S, thP / S, phP / S);
		    v4Cum += error;
		    v4c = logError(v4Cum / count);
		    v4e = logError(error);
		}

		public void refresh (double r, double th) {
		    sth = sin(th);
		    cth = cos(th);
		    sth2 = sth * sth;
		    double cth2 = 1.0 - sth2;
		    ra2 = r * r + a2;
			D = ra2 - 2.0 * r;
			S = r * r + a2 * cth2;
		    R = (((cr[0] * r + cr[1]) * r + cr[2]) * r + cr[3]) * r + cr[4];
			TH = a2xE2_mu2 + L2 / sth2;
			THETA = Q - cth2 * TH;
			double P_D = (ra2 * E - aL) / D;
		    tP = ra2 * P_D + aL - a2E * sth2;
		    phP = a * P_D - aE + L / sth2;
        }

        private void pUp (double c) {
		    rP += c * (((4.0 * cr[0] * r + 3.0 * cr[1]) * r + 2.0 * cr[2]) * r + cr[3]) * 0.5;
            double cot = cth / sth;
		    thP += c * (cth * sth * TH + L2 * cot * cot * cot);
        }

        private void qUp (double d) {
		    t += d * tP;
		    r += d * rP;
		    th += d * thP;
		    ph += d * phP;
		    refresh(r, th);
        }

        public void base2 (double w) {
		    pUp(w * 0.5);
		    qUp(w);
		    pUp(w * 0.5);
        }

        public void base4 (double w) {
		    pUp(w * coefficients[0]);
		    qUp(w * coefficients[1]);
		    pUp(w * coefficients[2]);
		    qUp(w * coefficients[3]);
		    pUp(w * coefficients[2]);
		    qUp(w * coefficients[1]);
		    pUp(w * coefficients[0]);
        }
    }

    static int main (string[] args) {
        unowned Json.Object obj;
        Json.Parser parser = new Json.Parser ();
        try {
            parser.load_from_data(read_stdin());
            obj = parser.get_root().get_object ();
        } catch (Error e) {
            stdout.printf ("Unable to parse the data file: %s\n", e.message);
            return -1;
        }
        Geodesic g = new Geodesic(obj.get_member("M").get_double(), obj.get_member("a").get_double(), obj.get_member("mu").get_double(), obj.get_member("E").get_double(), obj.get_member("Lz").get_double(), obj.get_member("C").get_double(), obj.get_member("r").get_double(), obj.get_member("theta").get_double(), obj.get_member("start").get_double(), obj.get_member("duration").get_double(), obj.get_member("step").get_double(), obj.get_member("integrator").get_string());
		g.refresh(g.r, g.th);
		g.rP = sqrt(fabs(g.R));
		g.thP = sqrt(fabs(g.THETA));
		double mino = 0.0;
		double tau = 0.0;
    	while (! (fabs(mino) > g.endtime)) {
		    g.count += 1;
		    g.errors(g.R, g.THETA, g.tP, g.rP, g.thP, g.phP);
		    if (fabs(mino) > g.starttime) {
			    stdout.printf("{\"mino\":%.9e, \"tau\":%.9e, \"v4e\":%.1f, \"v4c\":%.1f, \"ER\":%.1f, \"ETh\":%.1f, \"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, \"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e}\n", mino, tau, g.v4e, g.v4c, g.eR, g.eTh, g.t, g.r, g.th, g.ph, g.tP / g.S, g.rP / g.S, g.thP / g.S, g.phP / g.S);
            }
		    for (int i = 0; i < g.wRange; i++) {
		        g.base2(g.w[i] * g.h);
            }
		    mino += g.h;
		    tau += g.h * g.S;
        }
        return 0; 
    }
}

