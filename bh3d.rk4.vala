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
        private double a;
        private double mu2;
        private double E;
        private double L;
        private double Q;
        private double start;
        private double end;
        private double ts;
        private int64 tr;
        private double a2;
        private double horizon;
        private double aE;
        private double a2E;
        private double L2;
        private double aL;
        private double L_aE2;
        private double a2xE2_mu2;
        private double sth2;
        private double ra2;
        private double D;
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

        private BL (double a, double mu2, double E, double L, double Q, double r0, double th0, double tau0, double deltaTau, double tStep, int64 tRatio) {
            this.a = a;
            this.mu2 = mu2;
            this.E = E;
            this.L = L;
            this.Q = Q;
            this.a2 = a * a;
            this.horizon = 1.0 + sqrt(1.0 - a2);
            this.aE = a * E;
            this.a2E = a2 * E;
            this.L2 = L * L;
            this.aL = a * L;
            this.L_aE2 = (L - aE) * (L - aE);
            this.a2xE2_mu2 = - a2 * (E * E - mu2);
            this.start = tau0;
            this.end = tau0 + deltaTau;
            this.ts = tStep;
            this.tr = tRatio;
            this.r = r0;
            this.th = th0;
            derivatives(r, th, 0);
        }

        private void derivatives (double radius, double theta, int stage) {  // update intermediate variables, see Wilkins
            var r2 = radius * radius;
            sth2 = sin(theta) * sin(theta);
            var cth2 = 1.0 - sth2;
            ra2 = r2 + a2;
            D = ra2 - 2.0 * radius;
            S = r2 + a2 * cth2;
            R = (ra2 * E - aL) * (ra2 * E - aL) - D * (Q + L_aE2 + mu2 * r2);
            THETA = Q - cth2 * (a2xE2_mu2 + L2 / sth2);
            var P_D = (ra2 * E - aL) / D;
            tDot = (ra2 * P_D + aL - a2E * sth2) / S;
            rDot = sqrt(R > 0.0 ? R : -R) / S;
            thDot = sqrt(THETA > 0.0 ? THETA : -THETA) / S;
            phDot = (a * P_D - aE + L / sth2) / S;
            kt[stage] = tDot;
            kr[stage] = rDot;
            kth[stage] = thDot;
            kph[stage] = phDot;
        }

        private double sumK (double[] kx) {
            return (kx[0] + 3.0 * (kx[1] + kx[2]) + kx[3]) * ts / 8.0;
        }

        private void rk4Step () {
            sgnR = R > 0.0 ? sgnR : -sgnR;
            sgnTH = THETA > 0.0 ? sgnTH : -sgnTH;
            derivatives(r + 1.0 / 3.0 * sgnR * kr[0] * ts, th + 1.0 / 3.0 * sgnTH * kth[0] * ts, 1);
            derivatives(r + 2.0 / 3.0 * sgnR * kr[1] * ts, th + 2.0 / 3.0 * sgnTH * kth[1] * ts, 2);
            derivatives(r + sgnR * kr[2] * ts, th + sgnTH * kth[2] * ts, 3);
            t += sumK(kt);
            r += sumK(kr) * sgnR;
            th += sumK(kth) * sgnTH;
            ph += sumK(kph);
            derivatives(r, th, 0);
        }

        private void output (double tau) {
            var tmp1 = a * tDot - ra2 * phDot;
            var tmp2 = tDot - a * sth2 * phDot;
            var e = mu2 + sth2 / S * tmp1 * tmp1 + S / D * rDot * rDot + S * thDot * thDot - D / S * tmp2 * tmp2;
            e = e > 0.0 ? e : -e;
            stdout.printf("{\"tau\":%.9e, \"v4e\":%.1f, \"v4c\":%.1f, \"ER\":%.1f, \"ETh\":%.1f, ", tau, 10.0 * log10(e > 1.0e-18 ? e : 1.0e-18), -180.0, -180.0, -180.0);
            stdout.printf("\"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, ", t, r, th, ph);
            stdout.printf("\"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e}\n", tDot, rDot, thDot, phDot);
        }

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
            return new BL(o.get_double_member("a"), o.get_double_member("mu"), o.get_double_member("E"), o.get_double_member("Lz"), o.get_double_member("C"),
                          o.get_double_member("r"), o.get_double_member("theta"),
                          o.get_double_member("start"), o.get_double_member("duration"), o.get_double_member("step"), o.get_int_member("plotratio"));
        }
    }

    public static int main (string[] args) {
        BL.fromJson().solve();
        return 0;
    }
}
