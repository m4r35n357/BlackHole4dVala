/*
Copyright (c) 2014, 2015, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

using GLib.Math;
using Gee;

namespace Kerr {

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

        public void output () {
            stdout.printf("{\"qX\":%.6e,\"qY\":%.6e,\"qZ\":%.6e,\"pX\":%.6e,\"pY\":%.6e,\"pZ\":%.6e,\"mass\":%.3e}\n", qX, qY, qZ, pX, pY, pZ, mass);
        }
    }

    public class Simulation : GLib.Object, IModel {

        public ArrayList<Particle> bodies;
        private int np;
        public double iterations;
        public double g;
        public double timeStep;
        public double errorLimit;
        public double h { get; set; }
        public ISymplectic integrator;
        public double T;

        public Simulation (ArrayList<Particle> bodies, double g, double timeStep, double errorLimit, double simulationTime, string type) {
            this.bodies = bodies;
            this.np = bodies.size;
            this.g = g;
            this.h = timeStep;
            this.iterations = simulationTime / timeStep;
            this.integrator = Integrator.getIntegrator(this, type);
        }

        /**
         * Separation between points A and B
         * @param xA x coordinate of point A
         * @param yA y coordinate of point A
         * @param zA z coordinate of point A
         * @param xB x coordinate of point B
         * @param yB y coordinate of point B
         * @param zB z coordinate of point B
         * @return the Euclidean distance between points A and B
         */
        double distance (double xA, double yA, double zA, double xB, double yB, double zB) {
            return Math.sqrt(Math.pow(xB - xA, 2) + Math.pow(yB - yA, 2) + Math.pow(zB - zA, 2));
        }

        /**
         * Total (kinetic + potential) energy of the system
         * @return the total energy
         */
        public double getH () {
            double energy = 0.0;
            for (int i = 0; i < np; i++) {
                var a = bodies[i];
                energy += (0.5 * (a.pX * a.pX + a.pY * a.pY + a.pZ * a.pZ)) / a.mass;
                for (int j = 0; j < np; j++) {
                    if (i > j) {
                        var b = bodies[j];
                        energy -= (g * a.mass * b.mass) / (distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ));
                    }
                }
            }
            return energy;
        }

        /**
         * Position update implements dH/dp, which in this case is a function of p only
         * @param c composition coefficient
         */
        public void qUp (double c) {
            for (int i = 0; i < np; i++) {
                var a = bodies[i];
                double tmp = c * timeStep / a.mass;
                a.qX += a.pX * tmp;
                a.qY += a.pY * tmp;
                a.qZ += a.pZ * tmp;
            }
        }

        /**
         * Momentum update implements -dH/dq, which in this case is a function of q only
         * @param c composition coefficient
         */
        public void pUp (double c) {
            for (int i = 0; i < np; i++) {
                var a = bodies[i];
                for (int j = 0; j < np; j++) {
                    if (i > j) {
                        var b = bodies[j];
                        double tmp = - c * g * a.mass * b.mass * timeStep / Math.pow(distance(a.qX, a.qY, a.qZ, b.qX, b.qY, b.qZ), 3);
                        double dPx = (b.qX - a.qX) * tmp;
                        double dPy = (b.qY - a.qY) * tmp;
                        double dPz = (b.qZ - a.qZ) * tmp;
                        a.pX -= dPx;
                        a.pY -= dPy;
                        a.pZ -= dPz;
                        b.pX += dPx;
                        b.pY += dPy;
                        b.pZ += dPz;
                    }
                }
            }
        }

        public void evolve () {
            integrator.compose();
        }

        public void particlesJson () {
            stdout.printf("[");
            foreach (var particle in bodies) {
                particle.output();
            }
            stdout.printf("]");
        }

        /**
         * Sole user method
         */
        public void solve (string[] args) {
            double h0 = getH();
            double hMin = h0;
            double hMax = h0;
            particlesJson();
            long n = 1;
            output(n, 0.0, h0, h0, h0, h0, -999.9);
            while (n <= iterations) {
                evolve();
                double hNow = getH();
                double tmp = fabs(hNow - h0);
                double dH = tmp > 0.0 ? tmp : 1.0e-18;
                if (hNow < hMin) {
                    hMin = hNow;
                } else if (hNow > hMax) {
                    hMax = hNow;
                }
                double dbValue = 10.0 * Math.log10(dH);
                particlesJson();
                output(n, timeStep, hNow, h0, hMin, hMax, dbValue);
                if (dbValue > errorLimit) {
                    return;
                }
                n += 1;
            }
        }

        public static Simulation fromJson () {
            var bodies = new ArrayList<Particle> ();
            var ic = getJson();
            foreach (var body in ic.get_array_member("bodies").get_elements()) {
                bodies.add(new Particle(ic.get_double_member("qX"), ic.get_double_member("qY"), ic.get_double_member("qZ"), ic.get_double_member("pX"), ic.get_double_member("pY"), ic.get_double_member("pZ"), ic.get_double_member("mass")));
            }
            return new Simulation(bodies, ic.get_double_member("g"), ic.get_double_member("simulationTime"), ic.get_double_member("timeStep"), ic.get_double_member("errorLimit"), ic.get_string_member("integrator"));
        }

        public void output (long n, double timeStep, double hNow, double h0, double hMin, double hMax, double dbValue) {
            stderr.printf("{\"t\":%.2f, \"H\":%.9e, \"H0\":%.9e, \"H-\":%.9e, \"H+\":%.9e, \"ER\":%.1f}\n", n * timeStep, hNow, h0, hMin, hMax, dbValue);
        }

    }

    static int main (string[] args) {
        Simulation.fromJson().solve(args);
        return 0;
    }
}
