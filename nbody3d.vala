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

namespace Simulations {

    public class Particle : GLib.Object {

        public double qX;
        public double qY;
        public double qZ;
        public double pX;
        public double pY;
        public double pZ;
        public double mass;

        public Particle (double qX, double qY, double qZ, double pX, double pY, double pZ, double mass) {
            this.qX = qX;
            this.qY = qY;
            this.qZ = qZ;
            this.pX = pX;
            this.pY = pY;
            this.pZ = pZ;
            this.mass = mass;
        }

        public string toString () {
            return("{\"qX\":%.6e,\"qY\":%.6e,\"qZ\":%.6e,\"pX\":%.6e,\"pY\":%.6e,\"pZ\":%.6e,\"mass\":%.3e}".printf(qX, qY, qZ, pX, pY, pZ, mass));
        }
    }

    public class NBody : GLib.Object, IModel {

        private Particle[] bodies;
        private int np;
        private double g;
        private double ts;
        private double errorLimit;
        private double simulationTime;
        private ISymplectic integrator;

        /**
         * Private constructor, use the static factory
         */
        private NBody (Particle[] bodies, double g, double timeStep, double errorLimit, double simulationTime, string type) {
            this.bodies = bodies;
            this.np = bodies.length;
            this.g = g;
            this.ts = timeStep;
            this.errorLimit = errorLimit;
            this.simulationTime = simulationTime;
            this.integrator = Integrator.getIntegrator(this, type);
        }

        /**
         * Static factory
         */
        public static NBody fromJson () {
            Particle[] bodies = {};
            var ic = getJson();
            foreach (var node in ic.get_array_member("bodies").get_elements()) {
                var body = node.get_object();
                if (body.has_member("pX") && body.has_member("pY") && body.has_member("pZ")) {
                    bodies += new Particle(body.get_double_member("qX"),
                                           body.get_double_member("qY"),
                                           body.get_double_member("qZ"),
                                           body.get_double_member("pX"),
                                           body.get_double_member("pY"),
                                           body.get_double_member("pZ"),
                                           body.get_double_member("mass"));
                } else if (body.has_member("vX") && body.has_member("vY") && body.has_member("vZ")) {
                    var mass = body.get_double_member("mass");
                    bodies += new Particle(body.get_double_member("qX"),
                                           body.get_double_member("qY"),
                                           body.get_double_member("qZ"),
                                           body.get_double_member("vX") * mass,
                                           body.get_double_member("vY") * mass,
                                           body.get_double_member("vZ") * mass,
                                           mass);
                } else {
                    stderr.printf("Mixed use of momenta and velocity\n");
                }
            }
            return new NBody(bodies,
                              ic.get_double_member("g"),
                              ic.get_double_member("timeStep"),
                              ic.get_double_member("errorLimit"),
                              ic.get_double_member("simulationTime"),
                              ic.get_string_member("integratorOrder"));
        }

        /**
         * Euclidean distance between points A and B
         */
        double distance (double xA, double yA, double zA, double xB, double yB, double zB) {
            return sqrt(pow(xB - xA, 2) + pow(yB - yA, 2) + pow(zB - zA, 2));
        }

       /**
         * Total (kinetic + potential) energy of the system
         */
        public double h () {
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

        /**
         * Position update implements dH/dp, which in this case is a function of p only
         */
        public void qUp (double d) {
            for (var i = 0; i < np; i++) {
                var a = bodies[i];
                var tmp = d * ts / a.mass;
                a.qX += a.pX * tmp;
                a.qY += a.pY * tmp;
                a.qZ += a.pZ * tmp;
            }
        }

        /**
         * Momentum update implements -dH/dq, which in this case is a function of q only
         */
        public void pUp (double c) {
            for (var i = 0; i < np; i++) {
                var a = bodies[i];
                for (var j = i + 1; j < np; j++) {
                    var b = bodies[j];
                    var tmp = - c * ts * g * a.mass * b.mass / pow(distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ), 3);
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
         * Sole user method
         */
        public void solve () {
            var h0 = h();
            var t = 0.0;
            while (true) {
                var hNow = h();
                var dbValue = logError(fabs(hNow / h0 - 1.0));
                output(t, hNow, h0, dbValue);
                if (fabs(t) > simulationTime || dbValue > errorLimit) {
                    return;
                }
                integrator.compose();
                t += ts;
            }
        }

        public void output (double time, double hNow, double h0, double dbValue) {
            string[] data = {};
            foreach (var particle in bodies) {
                data += particle.toString();
            }
            stdout.printf("[".concat(string.joinv(",", data), "]\n"));
            stderr.printf("{\"t\":%.2f, \"H\":%.9e, \"H0\":%.9e, \"ER\":%.1f}\n", time, hNow, h0, dbValue);
        }
    }

    public static int main (string[] args) {
        NBody.fromJson().solve();
        return 0;
    }
}
