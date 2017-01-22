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

    protected class BhSymp : KdSBase, IModel, ISolver {
        /**
         * All fields are private
         */
        // Constants from IC file
        private ISymplectic integrator;

        /**
         * Private constructor, use a static factory
         */
        protected BhSymp (double lambda, double a, double mu2, double E, double L, double Q, double r0, double th0,
                      double tau0, double tauN, double tStep, int64 tRatio, string type) {
            base(lambda, a, mu2, E, L, Q, r0, th0, tau0, tauN, tStep, tRatio);
            this.integrator = Symplectic.getIntegrator(this, type);
            refresh(r, th);
            this.Ur = - sqrt(fabs(R));
            this.Uth = - sqrt(fabs(TH));
        }

        public void qUp (double c) {  // dx/dTau = dH/dxP
            t += c * h * Ut;
            r += c * h * Ur;
            th += c * h * Uth;
            ph += c * h * Uph;
            refresh(r, th);
        }

        public void pUp (double d) {  // dxP/dTau = - dH/dx (dR/dr & dTheta/dtheta), see Maxima file maths.wxm, "Kerr-deSitter"
            Ur += d * h * (r * (2.0 * E * P * X2 - mu2 * D_r) - (r * (1.0 - l_3 * r2) - l_3 * r * ra2 - 1.0) * (K + mu2 * r2));
            var sth = sin(th);
            Uth += d * h * cos(th) * (sth * a2 * (mu2 * D_th - l_3 * (K - a2mu2 * cth2)) + X2 * T / sth * (T / sth2 - 2.0 * aE));
        }

        /**
         * Externally visible method, sets up and controls the simulation
         */
        public int[] solve () {
            var mino = 0.0;
            var tau = 0.0;
            var iterationCount = 0;
            var plotCount = 0;
            while (tau < end) {
                if ((tau >= start) && (iterationCount % tr == 0)) {
                    plot(mino, tau);
                    plotCount += 1;
                }
                integrator.compose();
                iterationCount += 1;
                mino = iterationCount * h;
                tau += h * S;
            }
            plot(mino, tau);
            return { iterationCount, plotCount };  // for testing
        }

        private double modH (double xDot, double X) {
            return 0.5 * fabs(xDot * xDot - X);
        }

        /**
         * Write the simulated data to STDOUT
         */
        private void plot (double mino, double tau) {
            var eR = logError(modH(Ur, R));
            var eTh = logError(modH(Uth, TH));
            var v4e = logError(v4Error(Ut / S, Ur / S, Uth / S, Uph / S));
            stdout.printf("{\"mino\":%.9e, \"tau\":%.9e,", mino, tau);                                                    // time variables
            stdout.printf(" \"v4e\":%.1f, \"ER\":%.1f, \"ETh\":%.1f, \"D_r\":%.9e,", v4e, eR, eTh, D_r);                  // errors
            stdout.printf(" \"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e,", t, r, th, ph);                            // coordinates
            stdout.printf(" \"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e}\n", Ut / S, Ur / S, Uth / S, Uph / S);  // coordinate derivatives
        }
    }
}
