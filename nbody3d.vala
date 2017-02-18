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

    public class Body : GLib.Object {

        public double qX;
        public double qY;
        public double qZ;
        public double pX;
        public double pY;
        public double pZ;
        public double mass;

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

    public class NBody : IModel, ISolver, GLib.Object {

        private Body[] bodies;
        private int np;
        private double g;
        private double errorLimit;

        public NBody (Body[] bodies, double g, double errorLimit) {
            stderr.printf("Newtonian N-Body Simulation\n");
            this.bodies = bodies;
            this.np = bodies.length;
            this.g = g;
            this.errorLimit = errorLimit;
        }

        /**
         * Euclidean distance between points A and B
         */
        private double distance (double xA, double yA, double zA, double xB, double yB, double zB) {
            return sqrt((xB - xA) * (xB - xA) + (yB - yA) * (yB - yA) + (zB - zA) * (zB - zA));
        }

       /**
         * Total (kinetic + potential) energy of the system, the Hamiltonian
         */
        private double h () {
            var energy = 0.0;
            for (var i = 0; i < np; i++) {
                var a = bodies[i];
                energy += 0.5 * (a.pX * a.pX + a.pY * a.pY + a.pZ * a.pZ) / a.mass;
                for (var j = i + 1; j < np; j++) {
                    var b = bodies[j];
                    energy -= g * a.mass * b.mass / distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ);
                }
            }
            return energy;
        }

        public void qUp (double d) {
            for (var i = 0; i < np; i++) {
                var a = bodies[i];
                var tmp = d / a.mass;
                a.qX += a.pX * tmp;
                a.qY += a.pY * tmp;
                a.qZ += a.pZ * tmp;
            }
        }

        public void pUp (double c) {
            for (var i = 0; i < np; i++) {
                var a = bodies[i];
                for (var j = i + 1; j < np; j++) {
                    var b = bodies[j];
                    var d = distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ);
                    var tmp = - c * g * a.mass * b.mass / (d * d * d);
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
         * Externally visible method, sets up and controls the simulation
         */
        public int[] solve (ISymplectic integrator, double ts, double start, double end, int64 tr) {
            var h0 = h();
            var i = 0;
            var t = 0.0;
            var error = 0.0;
            while ((t < end) && (error < errorLimit)) {
                var hNow = h();
                error = hNow / h0 - 1.0;
                if ((t > start) && (i % tr == 0)) {
                    output(t, hNow, h0, error);
                }
                integrator.integrate(this);
                i += 1;
                t = i * ts;
            }
            return { i };
        }

        /**
         * Write the simulated data to STDOUT
         */
        private void output (double time, double hNow, double h0, double error) {
            string[] data = {};
            foreach (var body in bodies) {
                data += body.toString();
            }
            stdout.printf("[".concat(string.joinv(",", data),"]\n"));
            stderr.printf("{\"t\":%.2f,\"H\":%.9e,\"H0\":%.9e,\"ER\":%.9e}\n", time, hNow, h0, error);
        }
    }
}
