/*
Copyright (c) 2014, 2015, 2016, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
using GLib.Math;

namespace Simulations {

    public class KerrGeodesic : GLib.Object, IModel {
        // Constants from IC file
        private double a;
        private double mu2;
        private double E;
        private double L;
        private double Q;
        private double starttime;
        private double endtime;
        private ISymplectic integrator;
        // Derived Constants
        private double a2;
        private double horizon;
        private double aE;
        private double a2E;
        private double L2;
        private double aL;
        private double a2xE2_mu2;
        private double[] cr;
        // Variables
        private double h;
        private double sth;
        private double cth;
        private double sth2;
        private double ra2;
        private double D;
        private double S;
        private double R;
        private double TH;
        private double THETA;
        private double eR;
        private double eTh;
        private double v4Cum;
        private double v4c;
        private double v4e;
        // Boyer-Lindquist Coordinates
        private double t;
        private double r;
        private double th;
        private double ph;
        private double tDot;
        private double rP;
        private double thP;
        private double phDot;

        /**
         * Private constructor, use the static factory
         */
        private KerrGeodesic (double spin, double pMass2, double energy, double momentum, double carter, double r0, double thetaMin,
                         double starttime, double duration, double timestep, string type) {
            // Constants
            this.a = spin;
            this.mu2 = pMass2;
            this.E = energy;
            this.L = momentum;
            this.Q = carter;
            this.a2 = a * a;
            this.horizon = 1.0 + sqrt(1.0 - a2);
            this.aE = a * E;
            this.a2E = a2 * E;
            this.L2 = L * L;
            this.aL = a * L;
            var E2_mu2 = E * E - mu2;
            this.a2xE2_mu2 = - a2 * E2_mu2;
            this.cr = { E2_mu2, 2.0 * mu2, - a2xE2_mu2 - L2 - Q, 2.0 * ((aE - L) * (aE - L) + Q), - a2 * Q };  // see Wilkins
            this.starttime = starttime;
            this.endtime = starttime + duration;
            this.h = timestep;
            this.integrator = Integrator.getIntegrator(this, type);
            // Coordinates
            this.r = r0;
            this.th = thetaMin;
            refresh();
            this.rP = sqrt(fabs(R));
            this.thP = sqrt(fabs(THETA));
        }

        /**
         * Static factory
         */
        public static KerrGeodesic fromJson () {
            var ic = getJson();
            return new KerrGeodesic(ic.get_double_member("a"),
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

        private double modH (double xDot, double X) {
            return 0.5 * fabs(xDot * xDot - X);
        }

        private double v4Error (double tDot, double rDot, double thDot, double phDot) {
            var tmp1 = a * tDot - ra2 * phDot;
            var tmp2 = tDot - a * sth2 * phDot;
            return fabs(mu2 + sth2 / S * tmp1 * tmp1 + S / D * rDot * rDot + S * thDot * thDot - D / S * tmp2 * tmp2);
        }

        private void errors (int count) {
            eR = logError(modH(rP, R));
            eTh = logError(modH(thP, THETA));
            var error = v4Error(tDot / S, rP / S, thP / S, phDot / S);
            v4e = logError(error);
            v4Cum += error;
            v4c = logError(v4Cum / (count + 1));
        }

        private void refresh () {  // update intermediate variables
            var r2 = r * r;
            sth = sin(th);
            cth = cos(th);
            sth2 = sth * sth;
            var cth2 = 1.0 - sth2;
            ra2 = r2 + a2;
            D = ra2 - 2.0 * r;
            S = r2 + a2 * cth2;
            R = (((cr[0] * r + cr[1]) * r + cr[2]) * r + cr[3]) * r + cr[4];  // see Wilkins
            TH = a2xE2_mu2 + L2 / sth2;
            THETA = Q - cth2 * TH;  // see Wilkins
            var P_D = (ra2 * E - aL) / D;
            tDot = ra2 * P_D + aL - a2E * sth2;  // MTW eq.33.32d
            phDot = a * P_D - aE + L / sth2;  // MTW eq.33.32c
        }

        public void qUp (double d) {  // dx/dTau = dH/dxP
            t += d * h * tDot;
            r += d * h * rP;
            th += d * h * thP;
            ph += d * h * phDot;
            refresh();
        }

        public void pUp (double c) {  // dxP/dTau = - dH/dx
            rP += c * h * (((4.0 * cr[0] * r + 3.0 * cr[1]) * r + 2.0 * cr[2]) * r + cr[3]) * 0.5;  // dR/dr
            var cot = cth / sth;
            thP += c * h * (cth * sth * TH + L2 * cot * cot * cot);  // dTheta/dtheta see Maxima file maths.wxm, "My Equations (Mino Time)"
        }

        /**
         * Sole user method
         */
        public void solve () {
            var mino = 0.0;
            var tau = 0.0;
            var count = 0;
            while (mino <= endtime) {
                errors(count);
                if (mino >= starttime) {
                    output(mino, tau);
                }
                integrator.compose();
                mino += h;
                tau += h * S;
                count++;
            }
        }

        public void output (double mino, double tau) {
            stdout.printf("{");
            stdout.printf("\"mino\":%.9e, \"tau\":%.9e, ", mino, tau);                                                      // time variables
            stdout.printf("\"v4e\":%.1f, \"v4c\":%.1f, \"ER\":%.1f, \"ETh\":%.1f, ", v4e, v4c, eR, eTh);                    // errors
            stdout.printf("\"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, ", t, r, th, ph);                              // coordinates
            stdout.printf("\"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e", tDot / S, rP / S, thP / S, phDot / S);    // coordinate derivatives
            stdout.printf("}\n");
        }
    }
}

