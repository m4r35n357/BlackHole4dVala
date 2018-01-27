/*
Copyright (c) 2014-2018, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
using GLib.Math;

namespace Models {

    /**
     * Hamiltonian model for a Harmonic Oscillator
     */
    public class Oscillator : IModel, ISolver, GLib.Object {

        private double H0;
        private double k;
        private double m;
        private double x;
        private double Ux = 0.0;

        /**
         * Public constructor
         *
         * @param k the ""spring constant"
         * @param x0 the initial displacement
         */
        public Oscillator (double m, double k, double x0) {
            stderr.printf("Harmonic oscillator\n");
            this.m = m;
            this.k = k;
            this.x = x0;
            H0 = H();
        }

       /**
        *
        * Total (kinetic + potential) energy of the system, T + V
        *
        * @return the Hamiltonian
        */
        private double H () {
            return 0.5 * (m * Ux * Ux + k * x * x);
        }

        /**
         * {@inheritDoc}
         * @see IModel.qUpdate
         */
        public void qUpdate (double c) {
            x += c * m * Ux;
        }

        /**
         * {@inheritDoc}
         * @see IModel.pUpdate
         */
        public void pUpdate (double d) {
            Ux -= d * k * x;
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
            stdout.printf("{\"tau\":%.9e,\"v4e\":%.9e,\"t\":%.9e,\"x\":%.9e,\"xP\":%.9e}\n", time, H() - H0, time, x, Ux);
        }
    }

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

    /**
     * Hamiltonian model for Newtonian central body orbit
     */
    public class Newton : IModel, ISolver, GLib.Object {

        private double L;
        private double r;
        private double ph = 0.0;
        private double rDot;
        private double phDot;
        private double L2;
        private double H0;

        /**
         * Public constructor
         *
         * @param lFac fraction of azimuthal angular momentum requred for circular orbit at a radius of r0
         * @param r0 initial r coordinate
         */
        public Newton (double lFac, double r0) {
            stderr.printf("Newtonian Orbit\n");
            this.r = r0;
            this.L = lFac * sqrt(r0);
            this.L2 = L * L;
            this.rDot = - sqrt(r0 - L2) / r;
            this.H0 = H();
        }

       /**
        *
        * Total (kinetic + potential) energy of the system, T + V
        *
        * @return the Hamiltonian
        */
        private double H () {
            return 0.5 * (rDot * rDot + L2 / (r * r)) - 1.0 / r;
        }

        /**
         * {@inheritDoc}
         * @see IModel.qUpdate
         */
        public void qUpdate (double d) {
            r += d * rDot;
            phDot = L / (r * r);
            ph += d * phDot;
        }

        /**
         * {@inheritDoc}
         * @see IModel.pUpdate
         */
        public void pUpdate (double c) {
            rDot -= c * (1.0 - L2 / r) / (r * r);
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
            stdout.printf("\"t\":%.9e,\"r\":%.9e,\"th\":%.9e,\"ph\":%.9e,", time, r, PI_2, ph);
            stdout.printf("\"tP\":%.9e,\"rP\":%.9e,\"thP\":%.9e,\"phP\":%.9e}\n", 1.0, rDot, 0.0, phDot);
        }
    }

    /**
     * Hamiltonian model for Henon-Heiles Potential
     */
    public class HenonHeiles : IModel, ISolver, GLib.Object {

        private double H0;
        private double lambda = 1.0;
        private double x;
        private double y;
        private double Ux = 0.0;
        private double Uy = 0.0;

        /**
         * Public constructor
         *
         * @param x0 the initial displacement in x
         * @param y0 the initial displacement in y
         */
        public HenonHeiles (double x0, double y0) {
            stderr.printf("Henon-Heiles Potential\n");
            this.x = x0;
            this.y = y0;
            H0 = H();
        }

       /**
        *
        * Total (kinetic + potential) energy of the system, T + V
        *
        * @return the Hamiltonian
        */
        private double H () {
            return 0.5 * (Ux * Ux + Uy * Uy + x * x + y * y) + lambda * (x * x * y - y * y * y / 3.0);
        }

        /**
         * {@inheritDoc}
         * @see IModel.qUpdate
         */
        public void qUpdate (double c) {
            x += c * Ux;
            y += c * Uy;
        }

        /**
         * {@inheritDoc}
         * @see IModel.pUpdate
         */
        public void pUpdate (double d) {
            Ux -= d * (x + 2.0 * lambda * x * y);
            Uy -= d * (y + lambda * (x * x - y * y));
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
            stdout.printf("{\"tau\":%.9e,\"v4e\":%.9e,\"t\":%.9e,\"x\":%.9e,\"xP\":%.9e,\"y\":%.9e,\"yP\":%.9e}\n",
                            time, H() - H0, time, x, Ux, y, Uy);
        }
    }

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

    /**
     * A model for generating geodesics in the Kerr-deSitter spacetime
     */
    public class Bh3d : IModel, ISolver, GLib.Object {

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
        private double two_EX2;
        private double two_aE;
        private double aL;
        private double r2;  // Variables
        private double ra2;
        private double sth;
        private double cth;
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
        private bool cross;

        /**
         * Public constructor
         *
         * @param lambda the cosmological constant
         * @param spin black hole spin parameter = J / M
         * @param pMass2 squared particle mass (0 or 1) = mu2
         * @param energy constant of motion = E / mu
         * @param angMom constant of motion (azimuthal angular momentum) = L / (mu M)
         * @param Q constant of motion (Carter's constant) = C / (mu2 M2)
         * @param r0 initial r coordinate = r / M
         * @param th0 initial theta coordinate
         */
        public Bh3d (double lambda, double spin, double pMass2, double energy, double angMom, double Q, double r0, double th0, bool xh) {
            stderr.printf("Kerr-deSitter Geodesic\n");
            l_3 = lambda / 3.0;
            a = spin;
            mu2 = pMass2;
            E = energy;
            L = angMom;
            a2 = a * a;
            a2l_3 = a2 * l_3;
            a2mu2 = a2 * mu2;
            aE = a * E;
            aL = a * L;
            X2 = (1.0 + a2l_3) * (1.0 + a2l_3);
            two_EX2 = 2.0 * E * X2;
            two_aE = 2.0 * aE;
            K = Q + X2 * (L - aE) * (L - aE);
            r = r0;
            th = (90.0 - th0) * PI / 180.0;
            cross = xh;
            refresh();
            Ur = - sqrt(R >= 0.0 ? R : -R);
            Uth = - sqrt(TH >= 0.0 ? TH : -TH);
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
            sincos(th, out sth, out cth);  // THETA potential
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

        /**
         * {@inheritDoc}
         * @see IModel.qUpdate
         */
        public void qUpdate (double c) {  // dq/dTau = dT/dp
            t += c * Ut;
            r += c * Ur;
            th += c * Uth;
            ph += c * Uph;
            refresh();
        }

        /**
         * {@inheritDoc}
         * @see IModel.pUpdate
         */
        public void pUpdate (double d) {  // dp/dTau = - dV/dq (for dR/dr & dTheta/dtheta, see Maxima file maths.wxm, "Kerr-deSitter")
            Ur += d * (r * (two_EX2 * P - mu2 * D_r) - (r * (1.0 - l_3 * (r2 + ra2)) - 1.0) * (K + mu2 * r2));
            Uth += d * cth * (sth * a2 * (mu2 * D_th - l_3 * (K - a2mu2 * cth2)) + X2 * T / sth * (T / sth2 - two_aE));
        }

        /**
         * {@inheritDoc}
         * @see ISolver.solve
         */
        public int64[] solve (Integrators.ISymplectic integrator, double h, double start, double end, int64 tr) {
            double mino = 0.0;
            double tau = 0.0;
            int64 i = 0;
            int64 plotCount = 0;
            while ((tau < end) && (cross || D_r > 0.0)) {
                if ((tau >= start) && (i % tr == 0)) {
                    plot(mino, tau, Ut / S, Ur / S, Uth / S, Uph / S);
                    plotCount += 1;
                }
                integrator.step();
                i += 1;
                mino = h * i;
                tau += h * S;
            }
            plot(mino, tau, Ut / S, Ur / S, Uth / S, Uph / S);
            return { i, plotCount };  // for testing
        }

        /**
         * 4 velocity norm error (Hamiltonian - mass)
         *
         * @param Ut the time component of the four-velocity
         * @param Ur the radial component of the four-velocity
         * @param Uth the elevation component of the four-velocity
         * @param Uph the azimuth component of the four-velocity
         *
         * @return difference between the mass and the hamiltonian
         */
        private double v4Error (double Ut, double Ur, double Uth, double Uph) {
            var U1 = a * Ut - ra2 * Uph;
            var U4 = Ut - a * sth2 * Uph;
            var SX2 = S * X2;
            return mu2 + sth2 * D_th / SX2 * U1 * U1 + S / D_r * Ur * Ur + S / D_th * Uth * Uth - D_r / SX2 * U4 * U4;
        }

        /**
         * Write the simulated data to STDOUT as a JSON string
         *
         * @param mino the particle's mino time
         * @param tau the particle's proper time
         * @param Ut the time component of the four-velocity
         * @param Ur the radial component of the four-velocity
         * @param Uth the elevation component of the four-velocity
         * @param Uph the azimuth component of the four-velocity
         */
        private void plot (double mino, double tau, double Ut, double Ur, double Uth, double Uph) {
            var S2 = S * S;
            stdout.printf("{\"mino\":%.9e,\"tau\":%.9e,\"v4e\":%.9e,\"ER\":%.9e,\"ETh\":%.9e,\"t\":%.9e,\"r\":%.9e,\"th\":%.9e,\"ph\":%.9e,\"tP\":%.9e,\"rP\":%.9e,\"thP\":%.9e,\"phP\":%.9e}\n", mino, tau, v4Error(Ut,Ur,Uth,Uph), Ur*Ur-R/S2, Uth*Uth-TH/S2, t,r,th,ph, Ut,Ur,Uth,Uph);
        }
    }
}
