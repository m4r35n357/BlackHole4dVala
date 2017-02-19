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

    public class Newton : IModel, ISolver, GLib.Object {

        private double L;
        private double r;
        private double ph = 0.0;
        private double rDot;
        private double phDot;
        private double L2;
        private double H0;

        public Newton (double lFac, double r0) {
            stderr.printf("Newtonian Orbit\n");
            this.r = r0;
            this.L = sqrt(r0);
            this.L2 = L * L;
            var E0 = 0.5 * L2 / (r * r) - 1.0 / r;
            this.L = lFac * L;
            this.L2 = L * L;
            var V = 0.5 * L2 / (r * r) - 1.0 / r;
            this.rDot = - sqrt(2.0 * (E0 > V ? E0 - V : V - E0));
            this.H0 = H();
        }

       /**
         * Total (kinetic + potential) energy of the system, the Hamiltonian
         */
        private double H () {
            return 0.5 * (rDot * rDot + L2 / (r * r)) - 1.0 / r;
        }

        public void qUpdate (double d) {
            r += d * rDot;
            phDot = L / (r * r);
            ph += d * phDot;
        }

        public void pUpdate (double c) {
            rDot -= c * (1.0 / (r * r) - L2 / (r * r * r));
        }

        /**
         * Externally visible method, sets up and controls the simulation
         */
        public int[] solve (ISymplectic integrator, double h, double start, double end, int64 tr) {
            var i = 0;
            var t = 0.0;
            while (t < end) {
                if ((t > start) && (i % tr == 0)) {
                    output(t);
                }
                integrator.step(this);
                i += 1;
                t = i * h;
            }
            output(t);
            return { i };
        }

        /**
         * Write the simulated data to STDOUT
         */
        private void output (double time) {
            stdout.printf("{\"tau\":%.9e,\"v4e\":%.9e,", time, H() - H0);
            stdout.printf("\"t\":%.9e,\"r\":%.9e,\"th\":%.9e,\"ph\":%.9e,", time, r, PI_2, ph);
            stdout.printf("\"tP\":%.9e,\"rP\":%.9e,\"thP\":%.9e,\"phP\":%.9e}\n", 1.0, rDot, 0.0, phDot);
        }
    }
}

