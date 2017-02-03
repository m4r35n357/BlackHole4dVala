/*
Copyright (c) 2014, 2015, 2016, 2017 Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
using GLib.Math;

namespace Simulations {

    protected abstract class KdSBase : GLib.Object {
        /**
         * All fields are protected
         */
        // Constants from IC file
        protected double l_3;
        protected double a;
        protected double a2;
        protected double a2l_3;
        protected double mu2;
        protected double a2mu2;
        protected double X2;
        protected double E;
        protected double L;
        protected double K;
        protected double aE;
        protected double aL;
        // Variables
        protected double r2;
        protected double ra2;
        protected double sth2;
        protected double cth2;
        protected double D_r;
        protected double D_th;
        protected double S;
        protected double R;
        protected double P;
        protected double T;
        protected double TH;
        // State
        protected double t = 0.0;
        protected double r;
        protected double th;
        protected double ph = 0.0;
        protected double Ut;
        protected double Ur;
        protected double Uth;
        protected double Uph;

        /**
         * Protected constructor, use the static factory
         */
        protected KdSBase (double lambda, double a, double mu2, double E, double L, double Q, double r0, double th0) {
            stderr.printf("Kerr-deSitter Geodesic\n");
            this.l_3 = lambda / 3.0;
            this.a = a;
            this.mu2 = mu2;
            this.E = E;
            this.L = L;
            this.a2 = a * a;
            this.a2l_3 = a2 * l_3;
            this.a2mu2 = a2 * mu2;
            this.aE = a * E;
            this.aL = a * L;
            this.X2 = (1.0 + a2l_3) * (1.0 + a2l_3);
            this.K = Q + X2 * (L - aE) * (L - aE);
            this.r = r0;
            this.th = (90.0 - th0) * PI / 180.0;
        }

        /**
         * Calculate potentials & coordinate velocites
         */
        protected void refresh (double radius, double theta) {  // update intermediate variables, see Maxima file maths.wxm, "Kerr-deSitter"
            // R potential
            r2 = radius * radius;
            ra2 = r2 + a2;
            P = ra2 * E - aL;
            D_r = (1.0 - l_3 * r2) * ra2 - 2.0 * radius;
            R = X2 * P * P - D_r * (mu2 * r2 + K);
            // THETA potential
            sth2 = sin(theta) * sin(theta);
            cth2 = 1.0 - sth2;
            T = aE * sth2 - L;
            D_th = 1.0 + a2l_3 * cth2;
            TH = D_th * (K - a2mu2 * cth2) - X2 * T * T / sth2;
            // Equations of motion (Mino time, partial)
            var P_Dr = P / D_r;
            var T_Dth = T / D_th;
            S = r2 + a2 * cth2;
            Ut = (P_Dr * ra2 - T_Dth * a) * X2;
            Uph = (P_Dr * a - T_Dth / sth2) * X2;
        }

        /**
         * Calculate 4 velocity norm error
         */
        protected double v4Error (double tDot, double rDot, double thDot, double phDot) {
            var U1 = a * tDot - ra2 * phDot;
            var U4 = tDot - a * sth2 * phDot;
            var SX2 = S * X2;
            return fabs(mu2 + sth2 * D_th / SX2 * U1 * U1 + S / D_r * rDot * rDot + S / D_th * thDot * thDot - D_r / SX2 * U4 * U4);
        }

        /**
         *  Static factory
         */
        public static ISolver newInstance(double lambda, double a, double mu2, double E, double L, double Q, double r0, double th0, string t) {
            if (("rk4" == t) || ("rk438" == t)) {
                return new BhRk4(lambda, a, mu2, E, L, Q, r0, th0, t);
            } else if (("sb1" == t) || ("sb2" == t) || ("sb4" == t) || ("sc4" == t) || ("sc6" == t) || ("sh6" == t) || ("sh8" == t) || ("sh10" == t)) {
                return new BhSymp(lambda, a, mu2, E, L, Q, r0, th0, t);
            } else {
                stderr.printf("Bad integrator; should be [ rk4 | rk438 | sb2 | sb4 | sc4 | sc6 | sh6 | sh8 | sh10 ], found {%s}\n", t);
                assert_not_reached();
            }
        }
    }
}
