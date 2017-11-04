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
     * Description of a single Newtonian body using its coordinates, momenta and mass
     */
    public class Body : GLib.Object {

        public double qX;
        public double qY;
        public double qZ;
        public double pX;
        public double pY;
        public double pZ;
        public double mass;

        /**
         * Public constructor
         *
         * @param qX X coordinate
         * @param qY Y coordinate
         * @param qZ Z coordinate
         * @param pX X momentum
         * @param pY Y momentum
         * @param pZ Z momentum
         * @param mass particle mass
         */
        public Body (double qX, double qY, double qZ, double pX, double pY, double pZ, double mass) {
            this.qX = qX;
            this.qY = qY;
            this.qZ = qZ;
            this.pX = pX;
            this.pY = pY;
            this.pZ = pZ;
            this.mass = mass;
        }

        public string toString () {
            return("{\"qX\":%.9e,\"qY\":%.9e,\"qZ\":%.9e,\"pX\":%.9e,\"pY\":%.9e,\"pZ\":%.9e,\"mass\":%.9e}".printf(qX, qY, qZ, pX, pY, pZ, mass));
        }
    }

    /**
     * Hamiltonian model for Newtonian N body orbits
     */
    public class NBody : IModel, ISolver, GLib.Object {

        private Body[] bodies;
        private int n;
        private double m = 0.0;
        private double g;
        private double errorLimit;

        /**
         * Public constructor
         *
         * @param bodies array of {@link Body} instances
         * @param g Newton's gravitational constant
         * @param errorLimit termination condition
         */
        public NBody (Body[] bodies, double g, double errorLimit) {
            stderr.printf("Newtonian N-Body Simulation\n");
            this.bodies = bodies;
            this.n = bodies.length;
            for (var i = 0; i < n; i++) {
                this.m += bodies[i].mass;
            }
            this.g = g;
            this.errorLimit = errorLimit;
        }

        /**
         * Euclidean distance between two bodies, i and j
         */
        private double separation (Body i, Body j) {
            var sX = i.qX - j.qX;
            var sY = i.qY - j.qY;
            var sZ = i.qZ - j.qZ;
            return sqrt(sX * sX + sY * sY + sZ * sZ);
        }

        /**
         * Compensate for drift in centre of mass
         */
        private void trackCentreOfMass () {
            var mX = 0.0;
            var mY = 0.0;
            var mZ = 0.0;
            for (var i = 0; i < n; i++) {
                var a = bodies[i];
                mX += a.qX * a.mass;
                mY += a.qY * a.mass;
                mZ += a.qZ * a.mass;
            }
            for (var i = 0; i < n; i++) {
                var a = bodies[i];
                a.qX -= mX / m;
                a.qY -= mY / m;
                a.qZ -= mZ / m;
            }
        }

        /**
         *
         * Total (kinetic + potential) energy of the system, T + V
         *
         * @return the Hamiltonian
         */
        private double h () {
            double energy = 0.0;
            for (var i = 0; i < n; i++) {
                var a = bodies[i];
                energy += 0.5 * (a.pX * a.pX + a.pY * a.pY + a.pZ * a.pZ) / a.mass;
                for (var j = i + 1; j < n; j++) {
                    var b = bodies[j];
                    energy -= g * a.mass * b.mass / separation(a, b);
                }
            }
            return energy;
        }

        /**
         * {@inheritDoc}
         * @see IModel.qUpdate
         */
        public void qUpdate (double d) {
            for (var i = 0; i < n; i++) {
                var a = bodies[i];
                var tmp = d / a.mass;
                a.qX += a.pX * tmp;
                a.qY += a.pY * tmp;
                a.qZ += a.pZ * tmp;
            }
        }

        /**
         * {@inheritDoc}
         * @see IModel.pUpdate
         */
        public void pUpdate (double c) {
            for (var i = 0; i < n; i++) {
                var a = bodies[i];
                for (var j = i + 1; j < n; j++) {
                    var b = bodies[j];
                    var s = separation(a, b);
                    var tmp = - c * g * a.mass * b.mass / (s * s * s);
                    var dPx = (b.qX - a.qX) * tmp;
                    var dPy = (b.qY - a.qY) * tmp;
                    var dPz = (b.qZ - a.qZ) * tmp;
                    a.pX -= dPx;
                    a.pY -= dPy;
                    a.pZ -= dPz;
                    b.pX += dPx;
                    b.pY += dPy;
                    b.pZ += dPz;
                }
            }
        }

        /**
         * {@inheritDoc}
         * @see ISolver.solve
         */
        public int64[] solve (Integrators.ISymplectic integrator, double step, double start, double end, int64 tr) {
            var h0 = h();
            int64 i = 0;
            double t = 0.0;
            double error = 0.0;
            while ((t < end) && (error < errorLimit)) {
                var hNow = h();
                error = hNow / h0 - 1.0;
                if ((t > start) && (i % tr == 0)) {
                    output(t, hNow, h0, error);
                }
                integrator.step();
                i += 1;
                t = i * step;
            }
            return { i };
        }

        /**
         * Write the simulated data to STDOUT as a JSON string
         *
         * @param time absolute Newtonian time
         * @param hNow current Hamiltonian
         * @param h0 original Hamiltonian
         * @param error based on the difference
         */
        private void output (double time, double hNow, double h0, double error) {
            trackCentreOfMass();
            string[] data = {};
            foreach (var body in bodies) {
                data += body.toString();
            }
            stdout.printf("[".concat(string.joinv(",", data),"]\n"));
            stderr.printf("{\"t\":%.2f,\"H\":%.9e,\"H0\":%.9e,\"ER\":%.9e}\n", time, hNow, h0, error);
        }
    }
}
