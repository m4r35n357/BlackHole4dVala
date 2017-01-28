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

    protected class BhRk4 : KdSBase, ISolver {
        /**
         * For handling different RK4 implementations
         */
        delegate double UpdaterType(double[] x);
        delegate void EvaluatorType();

        /**
         * All fields are private
         */
        // Constants derived from IC file
        private EvaluatorType evaluator;
        private UpdaterType updater;
        // Variables
        private double[] kt = { 0.0, 0.0, 0.0, 0.0 };
        private double[] kr = { 0.0, 0.0, 0.0, 0.0 };
        private double[] kth = { 0.0, 0.0, 0.0, 0.0 };
        private double[] kph = { 0.0, 0.0, 0.0, 0.0 };
        private int sgnR = -1;
        private int sgnTH = -1;

        /**
         * Protected constructor, use a static factory
         */
        protected BhRk4 (double lambda, double a, double mu2, double E, double L, double Q, double r0, double th0,
                     double tau0, double tauN, double tStep, int64 tRatio, string type) {
            base(lambda, a, mu2, E, L, Q, r0, th0, tau0, tauN, tStep, tRatio);
            if ("rk4" == type) {
                this.evaluator = evaluator4;
                this.updater = updater4;
                stderr.printf("RK4 Integrator\n");
            } else if ("rk438" == type){
                this.evaluator = evaluator438;
                this.updater = updater438;
                stderr.printf("RK4 (3/8) Integrator\n");
            } else {
                stderr.printf("Integrator not recognized: %s\n", type);
                assert_not_reached();
            }
            f(r, th, 0);

        }

        /**
         * Calculate potentials & coordinate velocities, and populate RK4 arrays for each stage
         */
        private void f (double radius, double theta, int stage) {  // see maths.wxm in Maxima
            refresh(radius, theta);
            // Equations of motion; the base class gives expressions for Mino time, so divide by Sigma
            Ut /= S;
            Ur = sqrt(R > 0.0 ? R : -R) / S;
            Uth = sqrt(TH > 0.0 ? TH : -TH) / S;
            Uph /= S;
            // RK4 arrays
            kt[stage] = h * Ut;
            kr[stage] = h * Ur;
            kth[stage] = h * Uth;
            kph[stage] = h * Uph;
        }

        /**
         * updater - Sum the RK4 terms for a single coordinate
         */
        private double updater4 (double[] kx) {
            return (kx[0] + 2.0 * (kx[1] + kx[2]) + kx[3]) / 6.0;
        }

        /**
         * updater - Sum the RK4 terms for a single coordinate (3/8 method)
         */
        private double updater438 (double[] kx) {
            return (kx[0] + 3.0 * (kx[1] + kx[2]) + kx[3]) / 8.0;
        }

        /**
         * evaluator - (most of) the RK4 algorithm
         */
        private void evaluator4 () {
            f(r + 0.5 * kr[0], th + 0.5 * kth[0], 1);
            f(r + 0.5 * kr[1], th + 0.5 * kth[1], 2);
            f(r + kr[2], th + kth[2], 3);
        }

        /**
         * evaluator - (most of) the RK4 algorithm (3/8 method)
         */
        private void evaluator438 () {
            f(r + kr[0] / 3.0, th + kth[0] / 3.0, 1);
            f(r - kr[0] / 3.0 + kr[1], th - kth[0] / 3.0 + kth[1], 2);
            f(r + kr[0] - kr[1] + kr[2], th + kth[0] - kth[1] + kth[2], 3);
        }

        /**
         * RK4 common code including turning point handling
         */
        private void iterate () {
            sgnR = R > 0.0 ? sgnR : -sgnR;  // establish direction of R motion
            sgnTH = TH > 0.0 ? sgnTH : -sgnTH;  // establish direction of TH motion
            evaluator();
            t += updater(kt);
            r += updater(kr) * sgnR;
            th += updater(kth) * sgnTH;
            ph += updater(kph);
            f(r, th, 0);
        }

        /**
         * Externally visible method, sets up and controls the simulation
         */
        public int[] solve () {
            var tau = 0.0;
            var iterationCount = 0;
            var plotCount = 0;
            while (tau < end) {
                if ((tau >= start) && (iterationCount % tr == 0)) {
                    plot(tau);
                    plotCount += 1;
                }
                iterate();
                iterationCount += 1;
                tau = iterationCount * h;
            }
            plot(tau);
            return { iterationCount, plotCount };  // for testing
        }

        /**
         * Write the simulated data to STDOUT
         */
        private void plot (double tau) {
            stdout.printf("{\"tau\":%.9e, \"v4e\":%.1f, \"D_r\":%.9e, \"R\":%.9e, \"TH\":%.9e,",
                            tau, logError(v4Error(Ut, Ur, Uth, Uph)), D_r, R, TH);
            stdout.printf(" \"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, \"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e}\n",
                            t, r, th, ph, Ut, Ur, Uth, Uph);
        }
    }
}
