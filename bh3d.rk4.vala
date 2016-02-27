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

namespace Simulations {

    public class KerrGeodesicRk4 : GLib.Object {
        private double a;
        private double mu2;
        private double E;
        private double L;
        private double Q;
        private double starttime;
        private double endtime;
        private double a2;
        private double aE;
        private double a2E;
        private double L2;
        private double aL;
        private double L_aE2;
        private double a2xE2_mu2;
        private double h;
        private double sth2;
        private double ra2;
        private double D;
        private double S;
        private double R;
        private double THETA;
        private double v4Cum;
        private double v4c;
        private double v4e;
        private double t = 0.0;
        private double r;
        private double th;
        private double ph = 0.0;
        private double tDot;
        private double rP;
        private double thP;
        private double phDot;
        private double[] kt = { 0.0, 0.0, 0.0, 0.0 };
        private double[] kr = { 0.0, 0.0, 0.0, 0.0 };
        private double[] kth = { 0.0, 0.0, 0.0, 0.0 };
        private double[] kph = { 0.0, 0.0, 0.0, 0.0 };
        private double sgnR = 1.0;
        private double sgnTHETA = 1.0;

        private KerrGeodesicRk4 (double spin, double pMass2, double energy, double momentum, double carter, double r0, double thetaMin,
                         double starttime, double duration, double timestep) {
            this.a = spin;
            this.mu2 = pMass2;
            this.E = energy;
            this.L = momentum;
            this.Q = carter;
            this.a2 = a * a;
            this.aE = a * E;
            this.a2E = a2 * E;
            this.L2 = L * L;
            this.aL = a * L;
            this.L_aE2 = (L - aE) * (L - aE);
            this.a2xE2_mu2 = - a2 * (E * E - mu2);
            this.starttime = starttime;
            this.endtime = starttime + duration;
            this.h = timestep;
            this.r = r0;
            this.th = thetaMin;
            refresh(r, th);
        }

        public static KerrGeodesicRk4 fromJson () {
            var ic = getJson();
            return new KerrGeodesicRk4(ic.get_double_member("a"),
                                    ic.get_double_member("mu"),
                                    ic.get_double_member("E"),
                                    ic.get_double_member("Lz"),
                                    ic.get_double_member("C"),
                                    ic.get_double_member("r"),
                                    ic.get_double_member("theta"),
                                    ic.get_double_member("start"),
                                    ic.get_double_member("duration"),
                                    ic.get_double_member("step"));
        }

        private double v4Error (double tDot, double rDot, double thDot, double phDot) {
            var tmp1 = a * tDot - ra2 * phDot;
            var tmp2 = tDot - a * sth2 * phDot;
            return fabs(mu2 + sth2 / S * tmp1 * tmp1 + S / D * rDot * rDot + S * thDot * thDot - D / S * tmp2 * tmp2);
        }

        private void errors (int count) {
            var error = v4Error(tDot / S, rP / S, thP / S, phDot / S);
            v4e = logError(error);
            v4Cum += error;
            v4c = logError(v4Cum / count);
        }

        private void refresh (double r, double th) {  // update intermediate variables
            var r2 = r * r;
            sth2 = sin(th) * sin(th);
            var cth2 = 1.0 - sth2;
            ra2 = r2 + a2;
            D = ra2 - 2.0 * r;
            S = r2 + a2 * cth2;
            var P = ra2 * E - aL;
            R = P * P - D * (Q + L_aE2 + mu2 * r2);  // see Wilkins
            THETA = Q - cth2 * (a2xE2_mu2 + L2 / sth2);  // see Wilkins
            var P_D = (ra2 * E - aL) / D;
            tDot = ra2 * P_D + aL - a2E * sth2;  // MTW eq.33.32d
            rP = sqrt(fabs(R));
            thP = sqrt(fabs(THETA));
            phDot = a * P_D - aE + L / sth2;  // MTW eq.33.32c
        }

        private void k (int i) {
            kt[i] = h * tDot;
            kr[i] = h * rP;
            kth[i] = h * thP;
            kph[i] = h * phDot;
        }

        private double rk4update (double[] k) {
            return (k[0] + 3.0 * (k[1] + k[2]) + k[3]) / 8.0;
        }

        private void rk4 () {
            sgnR = R > 0.0 ? sgnR : -sgnR;
            sgnTHETA = THETA > 0.0 ? sgnTHETA : -sgnTHETA;
            k(0);
            refresh(r + 1.0 / 3.0 * sgnR * kr[0], th + 1.0 / 3.0 * sgnTHETA * kth[0]);
            k(1);
            refresh(r + 2.0 / 3.0 * sgnR * kr[1], th + 2.0 / 3.0 * sgnTHETA * kth[1]);
            k(2);
            refresh(r + sgnR * kr[2], th + sgnTHETA * kth[2]);
            k(3);
            t += rk4update(kt);
            r += sgnR * rk4update(kr);
            th += sgnTHETA * rk4update(kth);
            ph += rk4update(kph);
            refresh(r, th);
        }

        public void solve () {
            var mino = 0.0;
            var tau = 0.0;
            var count = 1;
            while (mino <= endtime) {
                errors(count);
                if (mino >= starttime) {
                    output(mino, tau);
                }
                rk4();
                mino += h;
                tau += h * S;
                count++;
            }
        }

        public void output (double mino, double tau) {
            stdout.printf("{");
            stdout.printf("\"mino\":%.9e, \"tau\":%.9e, ", mino, tau);                                                      // time variables
            stdout.printf("\"v4e\":%.1f, \"v4c\":%.1f, \"ER\":%.1f, \"ETh\":%.1f, ", v4e, v4c, -180.0, -180.0);             // errors
            stdout.printf("\"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, ", t, r, th, ph);                              // coordinates
            stdout.printf("\"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e", tDot / S, rP / S, thP / S, phDot / S);    // coordinate derivatives
            stdout.printf("}\n");
        }
    }

    private static Json.Object getJson () {
        var input = new StringBuilder();
        var buffer = new char[1024];
        while (!stdin.eof()) {
            var chunk = stdin.gets(buffer);
            if (chunk != null) {
                input.append(chunk);
            }
        }
        unowned Json.Object obj;
        var p = new Parser();
        try {
            p.load_from_data(input.str);
            obj = p.get_root().get_object();
        } catch (Error e) {
            stderr.printf("Unable to parse the data file: %s\n", e.message);
            return new Json.Object();
        }
        return obj;
    }

    private static double logError (double e) {
        return 10.0 * log10(e > 1.0e-18 ? e : 1.0e-18);
    }

    public static int main (string[] args) {
        KerrGeodesicRk4.fromJson().solve();
        return 0;
    }
}

