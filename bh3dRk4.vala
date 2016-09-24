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

    public class BL : KerrBase {
        /**
         * All fields are private
         */
        // variables
        private double[] kt = { 0.0, 0.0, 0.0, 0.0 };
        private double[] kr = { 0.0, 0.0, 0.0, 0.0 };
        private double[] kth = { 0.0, 0.0, 0.0, 0.0 };
        private double[] kph = { 0.0, 0.0, 0.0, 0.0 };
        private double sgnR = -1.0;
        private double sgnTH = -1.0;

        /**
         * Private constructor, use the static factory
         */
        private BL (double lambda, double a, double mu2, double E, double L, double Q, double r0, double th0,
                    double tau0, double deltaTau, double tStep, int64 tRatio) {
            base(lambda, a, mu2, E, L, Q, r0, th0, tau0, deltaTau, tStep, tRatio);
            f(r, th, 0);
        }

        /**
         * Calculate potentials & coordinate velocities, and populate RK4 arrays for each stage
         */
        private void f (double radius, double theta, int stage) {  // see maths.wxm in Maxima
            // R potential
            var r2 = radius * radius;
            ra2 = r2 + a2;
            var P = ra2 * E - aL;
            D_r = (1.0 - l_3 * r2) * ra2 - 2.0 * radius;
            R = X2 * P * P - D_r * (mu2 * r2 + K);
            // THETA potential
            var sth = sin(theta);
            sth2 = sth * sth;
            var cth2 = 1.0 - sth2;
            var T = aE * sth2 - L;
            D_th = 1.0 + a2l_3 * cth2;
            TH = D_th * (K - a2mu2 * cth2) - X2 / sth2 * T * T;
            // Equations of motion
            var P_Dr = P / D_r;
            var T_Dth = T / D_th;
            S = r2 + a2 * cth2;
            var X2_S = X2 / S;
            Ut = (P_Dr * ra2 - T_Dth * a) * X2_S;
            Ur = sqrt(R > 0.0 ? R : -R) / S;
            Uth = sqrt(TH > 0.0 ? TH : -TH) / S;
            Uph = (P_Dr * a - T_Dth / sth2) * X2_S;
            // RK4 arrays
            kt[stage] = h * Ut;
            kr[stage] = h * Ur;
            kth[stage] = h * Uth;
            kph[stage] = h * Uph;
        }

        /**
         * Sum the RK4 terms for a single coordinate
         */
        private double sumK (double[] kx) {
            return (kx[0] + 2.0 * (kx[1] + kx[2]) + kx[3]) / 6.0;
        }

        /**
         * The RK4 algorithm including turning point handling
         */
        private void rk4Step () {
            sgnR = R > 0.0 ? sgnR : -sgnR;
            sgnTH = TH > 0.0 ? sgnTH : -sgnTH;
            f(r + 0.5 * kr[0], th + 0.5 * kth[0], 1);
            f(r + 0.5 * kr[1], th + 0.5 * kth[1], 2);
            f(r + kr[2], th + kth[2], 3);
            t += sumK(kt);
            r += sumK(kr) * sgnR;
            th += sumK(kth) * sgnTH;
            ph += sumK(kph);
            f(r, th, 0);
        }

        /**
         * Externally visible method, sets up and controls the simulation
         */
        public void solve () {
            int64 count = 0;
            while ((tau <= end) && (r >= horizon) && (D_r >= 0.0)) {
                if ((tau >= start) && (count % tr == 0)) {
                    output();
                }
                rk4Step();
                count += 1;
                tau += h;
            }
            output();
        }

        /**
         * Write the simulated data to STDOUT
         */
        private void output () {
            var U1 = a * Ut - ra2 * Uph;
            var U4 = Ut - a * sth2 * Uph;
            var SX2 = S * X2;
            var e = mu2 + sth2 * D_th / SX2 * U1 * U1 + S / D_r * Ur * Ur + S / D_th * Uth * Uth - D_r / SX2 * U4 * U4;
            stdout.printf("{\"tau\":%.9e, \"v4e\":%.1f, \"v4c\":%.1f, \"ER\":%.1f, \"ETh\":%.1f, ",
                            tau, logError(e = e > 0.0 ? e : -e), -180.0, -180.0, -180.0);
            stdout.printf("\"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, \"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e}\n",
                            t, r, th, ph, Ut, Ur, Uth, Uph);
        }

        /**
         * Static factory from STDIN in JSON format
         */
        public static BL fromJson () {
            var o = getJson().get_object_member("IC");
            return new BL(o.get_double_member("lambda"), o.get_double_member("a"), o.get_double_member("mu"), o.get_double_member("E"),
                          o.get_double_member("L"), o.get_double_member("Q"), o.get_double_member("r0"), o.get_double_member("th0"),
                          o.get_double_member("start"), o.get_double_member("duration"), o.get_double_member("step"), o.get_int_member("plotratio"));
        }
    }

    public static int main (string[] args) {
        BL.fromJson().solve();
        return 0;
    }
}
