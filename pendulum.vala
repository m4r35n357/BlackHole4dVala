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

    /**
     * Hamiltonian model for a simple pendulum
     */
    public class Pendulum : IModel, ISolver, GLib.Object {

        private double H0;
        private double l;
        private double l2;
        private double g;
        private double m;
        private double theta;
        private double Utheta = 0.0;

        /**
         * Public constructor
         *
         * @param length of the pendulum
         * @param theta0 initial angle
         */
        public Pendulum (double g, double m, double length, double theta0) {
            stderr.printf("Simple Pendulum\n");
            this.g = g;
            this.m = m;
            this.l = length;
            this.l2 = length * length;
            this.theta = theta0;
            this.H0 = H();
        }

       /**
        *
        * Total (kinetic + potential) energy of the system, T + V
        *
        * @return the Hamiltonian
        */
        private double H () {
            return 0.5 * m * l2 * Utheta * Utheta - m * g * l * cos(theta);
        }

        /**
         * {@inheritDoc}
         * @see IModel.qUpdate
         */
        public void qUpdate (double c) {
            theta += c * m * l2 * Utheta;
        }

        /**
         * {@inheritDoc}
         * @see IModel.pUpdate
         */
        public void pUpdate (double d) {
            Utheta -= d * m * g * l * sin(theta);
        }

        /**
         * {@inheritDoc}
         * @see ISolver.solve
         */
        public int64[] solve (Integrators.ISymplectic integrator, double h, double start, double end, int64 tr) {
            int64 i = 0;
            double t = 0.0;
            while (t < end) {
                if ((t > start) && (i % tr == 0)) {
                    output(t);
                }
                integrator.step();
                i += 1;
                t = i * h;
            }
            output(t);
            return { i };
        }

        /**
         * Write the simulated data to STDOUT as a JSON string
         *
         * @param time absolute Newtonian time
         */
        private void output (double time) {
            stdout.printf("{\"tau\":%.9e,\"v4e\":%.9e,", time, H() - H0);
            stdout.printf("\"t\":%.9e,\"th\":%.9e,", time, theta);
            stdout.printf("\"thP\":%.9e}\n", Utheta);
        }
    }
}

