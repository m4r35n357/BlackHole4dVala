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

    public class KerrGeodesic : KerrBase, IModel {
        /**
         * All fields are private
         */
        // Constants from IC file
        private ISymplectic integrator;
        // Variables
        private double sth;
        private double cth;
        private double cth2;
        private double r2;
        private double P;
        private double T;
        private double eR;
        private double eTh;
        private double v4Cum;
        private double v4c;
        private double v4e;

        /**
         * Private constructor, use the static factory
         */
        private KerrGeodesic (double lambda, double a, double mu2, double E, double L, double Q, double r0, double th0, double tau0, double deltaTau, double tStep, int64 tRatio, string type) {
            base(lambda, a, mu2, E, L, Q, r0, th0, tau0, deltaTau, tStep, tRatio);
            this.integrator = Integrator.getIntegrator(this, type);
            refresh(r, th);
            this.Ur = - sqrt(fabs(R));
            this.Uth = - sqrt(fabs(TH));
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
            eR = logError(modH(Ur, R));
            eTh = logError(modH(Uth, TH));
            var error = v4Error(Ut / S, Ur / S, Uth / S, Uph / S);
            v4e = logError(error);
            v4Cum += error;
            v4c = logError(v4Cum / count);
        }

        /**i
         * Calculate potentials & coordinate velocites
         */
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
            Ut = (P_Dr * ra2 - T_Dth * a) * X2;
            Uph = (P_Dr * a - T_Dth / sth2) * X2;
        }

        public void qUp (double d) {  // dx/dTau = dH/dxP
            t += d * h * Ut;
            r += d * h * Ur;
            th += d * h * Uth;
            ph += d * h * Uph;
            refresh(r, th);
        }

        public void pUp (double c) {  // dxP/dTau = - dH/dx (dR/dr & dTheta/dtheta), see Maxima file maths.wxm, "Kerr-deSitter"
            Ur += c * h * (2.0 * r * E * P * X2 - (r * (1.0 - l_3 * r2) - 1.0 - l_3 * r * ra2) * (K + mu2 * r2) - mu2 * r * D_r);
            Uth += c * h * (cth * sth * a2 * (mu2 * D_th - l_3 * (K - a2mu2 * cth2)) + cth * X2 * T / sth * (T / sth2 - 2.0 * a * E ));
        }

        /**
         * Externally visible method, sets up and controls the simulation
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

        /**
         * Write the simulated data to STDOUT
         */
        public void output (double mino, double tau) {
            stdout.printf("{\"mino\":%.9e, \"tau\":%.9e, ", mino, tau);                                                      // time variables
            stdout.printf("\"v4e\":%.1f, \"v4c\":%.1f, \"ER\":%.1f, \"ETh\":%.1f, ", v4e, v4c, eR, eTh);                    // errors
            stdout.printf("\"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, ", t, r, th, ph);                              // coordinates
            stdout.printf("\"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e}\n", Ut / S, Ur / S, Uth / S, Uph / S);    // coordinate derivatives
        }

        /**
         * Static factory from STDIN in JSON format
         */
        public static KerrGeodesic fromJson () {
            var ic = getJson().get_object_member("IC");
            return new KerrGeodesic(ic.get_double_member("lambda"), ic.get_double_member("a"),
                                    ic.get_double_member("mu"), ic.get_double_member("E"), ic.get_double_member("L"), ic.get_double_member("Q"),
                                    ic.get_double_member("r0"), ic.get_double_member("th0"),
                                    ic.get_double_member("start"), ic.get_double_member("duration"), ic.get_double_member("step"),
                                    ic.get_int_member("plotratio"), ic.get_string_member("integrator"));
        }
    }

    public static int main (string[] args) {
        KerrGeodesic.fromJson().solve();
        return 0;
    }
}

