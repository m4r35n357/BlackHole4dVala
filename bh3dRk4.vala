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
using Json;

namespace Sim {

    public class BL : GLib.Object {

        /**
         * All fields are private
         */
        private double lambda;
        private double a;
        private double mu2;
        private double E;
        private double L;
        private double K;
        private double start;
        private double end;
        private double ts;
        private int64 tr;
        private double a2;
        private double horizon;
        private double aE;
        private double aL;
        private double sth2;
        private double ra2;
        private double chi2;
        private double Delta_r;
        private double Delta_th;
        private double S;
        private double R;
        private double THETA;
        private double t = 0.0;
        private double r;
        private double th;
        private double ph = 0.0;
        private double tDot;
        private double rDot;
        private double thDot;
        private double phDot;
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
            this.lambda = lambda;
            this.a = a;
            this.mu2 = mu2;
            this.E = E;
            this.L = L;
            this.a2 = a * a;
            this.horizon = 1.0 + sqrt(1.0 - a2);
            this.aE = a * E;
            this.aL = a * L;
            this.chi2 = (1.0 + a2 * lambda / 3.0) * (1.0 + a2 * lambda / 3.0);
            this.K = Q + chi2 * (L - a * E) * (L - a * E);
            this.start = tau0;
            this.end = tau0 + deltaTau;
            this.ts = tStep;
            this.tr = tRatio;
            this.r = r0;
            this.th = th0;
            f(r, th, 0);
        }

        /**
         * Calculate potentials & coordinate velocites, and populate RK4 arrays for each stage
         */
        private void f (double radius, double theta, int stage) {  // update intermediate variables, see Wilkins 1972 eqns. 1 to 5
            var r2 = radius * radius;
            sth2 = sin(theta) * sin(theta);
            var cth2 = 1.0 - sth2;
            ra2 = r2 + a2;
            S = r2 + a2 * cth2;
            var P1 = ra2 * E - aL;
            Delta_r = (1.0 - lambda / 3 * r2) * ra2 - 2 * r;
            R = chi2 * P1 * P1 - Delta_r * (mu2 * r2 + K);
            var T1 = aE * sth2 - L;
            Delta_th = 1.0 + lambda / 3 * a2 * cth2;
            THETA = Delta_th * (K - mu2 * a2 * cth2) - chi2 / sth2 * T1 * T1;
            //THETA = Q + (L - a * E) * (L - a * E) - (a * E * sth2 - L) * (a * E * sth2 - L) / sth2 - mu2 * a2 * cth2;
            tDot = (P1 * ra2 / Delta_r - T1 * a / Delta_th) * chi2 / S;
            rDot = sqrt(R > 0.0 ? R : -R) / S;
            thDot = sqrt(THETA > 0.0 ? THETA : -THETA) / S;
            phDot = (P1 * a / Delta_r - T1 / (Delta_th * sth2)) * chi2 / S;

            kt[stage] = ts * tDot;
            kr[stage] = ts * rDot;
            kth[stage] = ts * thDot;
            kph[stage] = ts * phDot;
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
            sgnTH = THETA > 0.0 ? sgnTH : -sgnTH;
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
         * Write the simulated data to STDOUT
         */
        private void output (double tau) {
            var tmp1 = a * tDot - ra2 * phDot;
            var tmp2 = tDot - a * sth2 * phDot;
            var e = mu2 + sth2 * Delta_th / (S * chi2) * tmp1 * tmp1 + S / Delta_r * rDot * rDot
                        + S / Delta_th * thDot * thDot - Delta_r / (S * chi2) * tmp2 * tmp2;
            e = e > 0.0 ? e : -e;
            stdout.printf("{\"tau\":%.9e, \"v4e\":%.1f, \"v4c\":%.1f, \"ER\":%.1f, \"ETh\":%.1f, ",
                            tau, 10.0 * log10(e > 1.0e-18 ? e : 1.0e-18), -180.0, -180.0, -180.0);
            stdout.printf("\"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, ", t, r, th, ph);
            stdout.printf("\"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e}\n", tDot, rDot, thDot, phDot);
        }

        /**
         * Externally visible method, sets up and controls the simulation
         */
        public void solve () {
            int64 count = 0;
            var tau = 0.0;
            while ((tau <= end) && (r >= horizon)) {
                if ((tau >= start) && (count % tr == 0)) {
                    output(tau);
                }
                rk4Step();
                count += 1;
                tau += ts;
            }
            output(tau);
        }

        /**
         * Static factory from STDIN in JSON format
         */
        public static BL fromJson () {
            var input = new StringBuilder();
            var buffer = new char[1024];
            while (!stdin.eof()) {
                var chunk = stdin.gets(buffer);
                if (chunk != null) {
                    input.append(chunk);
                }
            }
            unowned Json.Object o;
            var p = new Parser();
            try {
                p.load_from_data(input.str);
                o = p.get_root().get_object();
            } catch (Error e) {
                stderr.printf("Unable to parse the input data: %s\n", e.message);
                return_if_reached();
            }
            return new BL(o.get_double_member("lambda"), o.get_double_member("a"), o.get_double_member("mu"),
                          o.get_double_member("E"), o.get_double_member("L"), o.get_double_member("Q"),
                          o.get_double_member("r0"), (90.0 - o.get_double_member("th0")) * PI / 180.0,
                          o.get_double_member("start"), o.get_double_member("duration"), o.get_double_member("step"), o.get_int_member("plotratio"));
        }
    }

    public static int main (string[] args) {
        BL.fromJson().solve();
        return 0;
    }
}
