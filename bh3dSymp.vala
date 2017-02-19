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

    public class BhSymp : IModel, ISolver, GLib.Object {

        private double l_3;  // Constants
        private double a;
        private double a2;
        private double a2l_3;
        private double mu2;
        private double a2mu2;
        private double X2;
        private double E;
        private double L;
        private double K;
        private double aE;
        private double aL;
        private double r2;  // Variables
        private double ra2;
        private double sth;
        private double sth2;
        private double cth2;
        private double D_r;
        private double D_th;
        private double S;
        private double R;
        private double P;
        private double T;
        private double TH;
        private double t = 0.0;  // State
        private double r;
        private double th;
        private double ph = 0.0;
        private double Ut;
        private double Ur;
        private double Uth;
        private double Uph;

        public BhSymp (double lambda, double a, double mu2, double E, double L, double Q, double r0, double th0) {
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
        private void refresh () {  // updates intermediate variables, see Maxima file maths.wxm, "Kerr-deSitter"
            r2 = r * r;  // R potential
            ra2 = r2 + a2;
            P = ra2 * E - aL;
            D_r = (1.0 - l_3 * r2) * ra2 - 2.0 * r;
            R = X2 * P * P - D_r * (mu2 * r2 + K);
            sth = sin(th);  // THETA potential
            sth2 = sth * sth;
            cth2 = 1.0 - sth2;
            T = aE * sth2 - L;
            D_th = 1.0 + a2l_3 * cth2;
            TH = D_th * (K - a2mu2 * cth2) - X2 * T * T / sth2;
            var P_Dr = P / D_r;  // Equations of motion (Mino time, t & ph)
            var T_Dth = T / D_th;
            Ut = (P_Dr * ra2 - T_Dth * a) * X2;
            Uph = (P_Dr * a - T_Dth / sth2) * X2;
            S = r2 + a2 * cth2;  // Sigma
        }

        public void qUpdate (double c) {  // dq/dTau = dT/dp
            t += c * Ut;
            r += c * Ur;
            th += c * Uth;
            ph += c * Uph;
            refresh();
        }

        public void pUpdate (double d) {  // dp/dTau = - dV/dq (for dR/dr & dTheta/dtheta, see Maxima file maths.wxm, "Kerr-deSitter")
            Ur += d * (r * (2.0 * E * P * X2 - mu2 * D_r) - (r * (1.0 - l_3 * (r2 + ra2)) - 1.0) * (K + mu2 * r2));
            Uth += d * cos(th) * (sth * a2 * (mu2 * D_th - l_3 * (K - a2mu2 * cth2)) + X2 * T / sth * (T / sth2 - 2.0 * aE));
        }

        /**
         * Externally visible method, sets up and controls the simulation
         */
        public int[] solve (ISymplectic integrator, double h, double start, double end, int64 tr) {
            var mino = 0.0;
            var tau = 0.0;
            var i = 0;
            var plotCount = 0;
            refresh();
            this.Ur = - sqrt(R >= 0.0 ? R : -R);
            this.Uth = - sqrt(TH >= 0.0 ? TH : -TH);
            while (tau < end) {
                if ((tau >= start) && (i % tr == 0)) {
                    plot(mino, tau, Ut / S, Ur / S, Uth / S, Uph / S);
                    plotCount += 1;
                }
                integrator.step(this);
                i += 1;
                mino = h * i;
                tau += h * S;
            }
            plot(mino, tau, Ut / S, Ur / S, Uth / S, Uph / S);
            return { i, plotCount };  // for testing
        }

        /**
         * 4 velocity norm error (Hamiltonian - mass)
         */
        private double v4Error (double Ut, double Ur, double Uth, double Uph) {
            var U1 = a * Ut - ra2 * Uph;
            var U4 = Ut - a * sth2 * Uph;
            var SX2 = S * X2;
            return mu2 + sth2 * D_th / SX2 * U1 * U1 + S / D_r * Ur * Ur + S / D_th * Uth * Uth - D_r / SX2 * U4 * U4;
        }

        /**
         * Write the simulated data to STDOUT
         */
        private void plot (double mino, double tau, double Ut, double Ur, double Uth, double Uph) {
            var S2 = S * S;
            var eR = 0.5 * (Ur * Ur - R / S2);
            var eTh = 0.5 * (Uth * Uth - TH / S2);
            var v4e = v4Error(Ut, Ur, Uth, Uph);
            stdout.printf("{\"mino\":%.9e,\"tau\":%.9e,\"v4e\":%.9e,\"ER\":%.9e,\"ETh\":%.9e,\"t\":%.9e,\"r\":%.9e,\"th\":%.9e,\"ph\":%.9e,\"tP\":%.9e,\"rP\":%.9e,\"thP\":%.9e,\"phP\":%.9e}\n", mino,tau, v4e,eR,eTh, t,r,th,ph, Ut/S,Ur/S,Uth/S,Uph/S);
        }
    }
}
