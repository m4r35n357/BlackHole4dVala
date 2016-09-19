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
        private double l_3;
        private double a;
        private double mu2;
        private double E;
        private double L;
        private double start;
        private double end;
        private double h;
        private double tr;
        private ISymplectic integrator;
        // Derived Constants
        private double a2;
        private double a2l_3;
        private double a2mu2;
        private double X2;
        private double horizon;
        private double aE;
        private double aL;
        private double K;
        // Variables
        private double sth;
        private double cth;
        private double sth2;
        private double cth2;
        private double r2;
        private double ra2;
        private double D_r;
        private double D_th;
        private double S;
        private double P;
        private double T;
        private double R;
        private double TH;
        private double eR;
        private double eTh;
        private double v4Cum;
        private double v4c;
        private double v4e;
        // Boyer-Lindquist Coordinates
        private double t = 0.0;
        private double r;
        private double th;
        private double ph = 0.0;
        private double tDot;
        private double rP;
        private double thP;
        private double phDot;

        /**
         * Private constructor, use the static factory
         */
        private KerrGeodesic (double lambda, double a, double mu2, double E, double L, double Q, double r0, double th0, double tau0, double deltaTau, double tStep, int64 tRatio, string type) {
            // Constants
            this.l_3 = lambda / 3.0;
            this.a = a;
            this.mu2 = mu2;
            this.E = E;
            this.L = L;
            this.a2 = a * a;
            this.a2l_3 = a2 * l_3;
            this.a2mu2 = a2 * mu2;
            this.horizon = 1.0 + sqrt(1.0 - a2);
            this.aE = a * E;
            this.aL = a * L;
            this.X2 = (1.0 + a2l_3) * (1.0 + a2l_3);
            this.K = Q + X2 * (L - aE) * (L - aE);
            this.start = tau0;
            this.end = tau0 + deltaTau;
            this.h = tStep;
            this.tr = tRatio;
            this.integrator = Integrator.getIntegrator(this, type);
            // Coordinates
            this.r = r0;
            this.th = (90.0 - th0) * PI / 180.0;
            refresh(r, th);
            this.rP = - sqrt(fabs(R));
            this.thP = - sqrt(fabs(TH));
        }

        /**
         * Static factory
         */
        public static KerrGeodesic fromJson () {
            var ic = getJson().get_object_member("IC");
            return new KerrGeodesic(ic.get_double_member("lambda"), ic.get_double_member("a"),
                                    ic.get_double_member("mu"), ic.get_double_member("E"), ic.get_double_member("L"), ic.get_double_member("Q"),
                                    ic.get_double_member("r0"), ic.get_double_member("th0"),
                                    ic.get_double_member("start"), ic.get_double_member("duration"), ic.get_double_member("step"),
                                    ic.get_int_member("plotratio"), ic.get_string_member("integrator"));
        }

        private double modH (double xDot, double X) {
            return 0.5 * fabs(xDot * xDot - X);
        }

        private double v4Error (double tDot, double rDot, double thDot, double phDot) {
            var U1 = a * tDot - ra2 * phDot;
            var U4 = tDot - a * sth2 * phDot;
            var SX2 = S * X2;
            return fabs(mu2 + sth2 * D_th / SX2 * U1 * U1 + S / D_r * rDot * rDot + S / D_th * thDot * thDot - D_r / SX2 * U4 * U4);
        }

        private void errors (int count) {
            eR = logError(modH(rP, R));
            eTh = logError(modH(thP, TH));
            var error = v4Error(tDot / S, rP / S, thP / S, phDot / S);
            v4e = logError(error);
            v4Cum += error;
            v4c = logError(v4Cum / count);
        }

        private void refresh (double radius, double theta) {  // update intermediate variables, see Maxima file maths.wxm, "Kerr-deSitter"
            // R potential
            r2 = radius * radius;
            ra2 = r2 + a2;
            P = ra2 * E - aL;
            D_r = (1.0 - l_3 * r2) * ra2 - 2.0 * radius;
            R = X2 * P * P - D_r * (mu2 * r2 + K);
            // THETA potential
            sth = sin(theta);
            cth = cos(theta);
            sth2 = sth * sth;
            cth2 = 1.0 - sth2;
            T = aE * sth2 - L;
            D_th = 1.0 + a2l_3 * cth2;
            TH = D_th * (K - a2mu2 * cth2) - X2 / sth2 * T * T;
            // Equations of motion
            var P_Dr = P / D_r;
            var T_Dth = T / D_th;
            S = r2 + a2 * cth2;
            tDot = (P_Dr * ra2 - T_Dth * a) * X2;
            phDot = (P_Dr * a - T_Dth / sth2) * X2;
        }

        public void qUp (double d) {  // dx/dTau = dH/dxP
            t += d * h * tDot;
            r += d * h * rP;
            th += d * h * thP;
            ph += d * h * phDot;
            refresh(r, th);
        }

        public void pUp (double c) {  // dxP/dTau = - dH/dx
            // dR/dr, see Maxima file maths.wxm, "Kerr-deSitter"
            rP += c * h * (2.0 * r * E * P * X2 - (r * (1.0 - l_3 * r2) - 1.0 - l_3 * r * ra2) * (K + mu2 * r2) - mu2 * r * D_r);
            // dTheta/dtheta
            thP += c * h * (cth * sth * a2 * (mu2 * D_th - l_3 * (K - a2mu2 * cth2)) + cth * X2 * T / sth * (T / sth2 - 2.0 * a * E ));
        }

        /**
         * Sole user method
         */
        public void solve () {
            var mino = 0.0;
            var tau = 0.0;
            var count = 1;
            var tmp = 0;
            while ((tau <= end) && (r >= horizon) && (D_r >= 0.0)) {
                if ((tau >= start) && (tau >= tmp * tr * h)) {
                    errors(count);
                    output(mino, tau);
                    tmp += 1;
                }
                integrator.compose();
                mino += h;
                tau += h * S;
                count += 1;
            }
            output(mino, tau);
        }

        public void output (double mino, double tau) {
            stdout.printf("{\"mino\":%.9e, \"tau\":%.9e, ", mino, tau);                                                      // time variables
            stdout.printf("\"v4e\":%.1f, \"v4c\":%.1f, \"ER\":%.1f, \"ETh\":%.1f, ", v4e, v4c, eR, eTh);                    // errors
            stdout.printf("\"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, ", t, r, th, ph);                              // coordinates
            stdout.printf("\"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e}\n", tDot / S, rP / S, thP / S, phDot / S);    // coordinate derivatives
        }
    }

    public static int main (string[] args) {
        KerrGeodesic.fromJson().solve();
        return 0;
    }
}

