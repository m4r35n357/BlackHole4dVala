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
using Json;
using Gsl;

namespace Sim {

    public class IcGenerator : GLib.Object {

        /**
         * Fixed parameters and constraints
         */
        private struct Params {
            public double mu2;
            public double rMin;
            public double rMax;
            public double elevation;
            public double a;
            public double Lfac;
        }

        /**
         * Array indices for variables
         */
        private enum X {
            E = 0, L = 1, Q = 2;
        }

        /**
         * Array indices for objective functions
         */
        private enum F {
            R1 = 0, R2 = 1, TH = 2;
        }

        /**
         * R potential
         */
        private static double R (double r, double E, double L, double Q, Params* p) {
            var a = p->a;
            var mu2 = p->mu2;
            return (E * (r * r + a * a) - a * L) * (E * (r * r + a * a) - a * L)
                    - (r * r + a * a - 2.0 * r) * (mu2 * r * r + Q + (L - a * E) * (L - a * E));
        }

        /**
         * Derivative of R potential wrt R
         */
        private static double dRdr (double r, double E, double L, double Q, Params* p) {
            var a = p->a;
            var mu2 = p->mu2;
            return 4.0 * r * E * (E * (r * r + a * a) - a * L)
                    - (2.0 * r - 2.0) * (mu2 * r * r + Q + (L - a * E) * (L - a * E)) - 2.0 * mu2 * r * (r * r + a * a - 2.0 * r);
        }

        /**
         * THETA potential
         */
        private static double THETA (double theta, double E, double L, double Q, Params* p) {
            var a = p->a;
            var mu2 = p->mu2;
            return Q - cos(theta) * cos(theta) * (a * a * (mu2 - E * E) + L * L / (sin(theta) * sin(theta)));
        }

        /**
         * Spherical orbits are determined by R potential at constant radius, its derivative, and THETA potential
         */
        private static int sphericalOrbit (Vector x, Params* p, Vector f) {
            var E = x.get(X.E);
            var L = x.get(X.L);
            var Q = x.get(X.Q);

            f.set(F.R1, R(p->rMax, E, L, Q, p));
            f.set(F.R2, dRdr(p->rMax, E, L, Q, p));
            f.set(F.TH, THETA(p->elevation, E, L, Q, p));

            return Status.SUCCESS;
        }

        /**
         * Non-spherical orbits are determined by R potentials at maximum and minimum radii, and THETA potential
         */
        private static int nonSphericalOrbit (Vector x, Params* p, Vector f) {
            var E = x.get(X.E);
            var L = x.get(X.L);
            var Q = x.get(X.Q);

            f.set(F.R1, R(p->rMin, E, L, Q, p));
            f.set(F.R2, R(p->rMax, E, L, Q, p));
            f.set(F.TH, THETA(p->elevation, E, L, Q, p));

            return Status.SUCCESS;
        }

        /**
         * Write the initial conditions file to STDOUT and potential data to STDERR for plotting
         */
        private void printOutput (MultirootFsolver s, size_t iterations) {
            var p = (Params*) s.function.params;
            var E = s.x.get(X.E);
            var L = s.x.get(X.L) * p->Lfac;
            var Q = s.x.get(X.Q);
            stdout.printf("{\n");
            stdout.printf("  \"generator\" : \"icgenParticle %s\",\n", s.name());
            stdout.printf("  \"iterations\" : %zu,\n", iterations);
            stdout.printf("  \"residuals\" : \"R1: %.1e, R2: %.1e, TH: %.1e\",\n", s.f.get(F.R1), s.f.get(F.R2), s.f.get(F.TH));
            stdout.printf("  \"deltas\" : \"dE: %.1e, dL: %.1e, dQ: %.1e\",\n", s.dx.get(X.E), s.dx.get(X.L), s.dx.get(X.Q));
            if (p->a * L >= 0.0) {
                stdout.printf("  \"direction\" : \"PROGRADE\",\n");
            } else {
                stdout.printf("  \"direction\" : \"RETROGRADE\",\n");
            }
            stdout.printf("  \"M\" : %.1f,\n", 1.0);
            stdout.printf("  \"a\" : %.1f,\n", p->a);
            stdout.printf("  \"mu\" : %.1f,\n", p->mu2);
            stdout.printf("  \"E\" : %.17g,\n", E);
            stdout.printf("  \"L\" : %.17g,\n", L);
            stdout.printf("  \"Q\" : %.17g,\n", Q);
            stdout.printf("  \"r0\" : %.1f,\n", 0.5 * (p->rMin + p->rMax));
            stdout.printf("  \"th0\" : %.0f,\n", 0.0);
            stdout.printf("  \"start\" : %.1f,\n", 0.0);
            stdout.printf("  \"duration\" : %.1f,\n", 5000.0);
            stdout.printf("  \"step\" : %.3f,\n", 0.001);
            stdout.printf("  \"plotratio\" : %.1d,\n", 500);
            stdout.printf("  \"integrator\" : \"%s\"\n", "sc4");
            stdout.printf("}\n");
            for (var x = 1; x <= 1001; x++) {
                var xValue = 1.0 * x / 1001;
                stderr.printf("{ \"x\" : %.6f, \"R\" : %.6f, \"THETA\" : %.6f }\n",
                                xValue * p->rMax * 1.1, R(xValue * p->rMax * 1.1, E, L, Q, p), THETA(xValue * PI, E, L, Q, p));
            }
        }

        private Vector initializeVariables (Json.Object input) {
            var initialValues = new Vector(3);
            initialValues.set(X.E, input.has_member("E0") ? input.get_double_member("E0") : 1.0);
            initialValues.set(X.L, input.has_member("L0") ? input.get_double_member("L0") : 5.0);
            initialValues.set(X.Q, input.has_member("Q0") ? input.get_double_member("Q0") : 0.0);
            return initialValues;
        }

        private Params initializeParams (Json.Object input, double rMin, double rMax) {
            return Params() {
                mu2 = 1.0,
                rMin = rMin,
                rMax = rMax,
                elevation = (1.0 - (input.has_member("elevation") ? (90.0 - input.get_double_member("elevation")) / 180.0 : 0.5)) * PI,
                a = input.has_member("spin") ? input.get_double_member("spin") : 0.0,
                Lfac = input.has_member("Lfac") ? input.get_double_member("Lfac") : 1.0
            };
        }

        /**
         * Externally visible method, sets up and controls the solver
         */
        public void generate (Json.Object input) {
            var initialValues = initializeVariables(input);
            var nDim = initialValues.size;

            // choose a solver
            MultirootFsolver solver;
            switch (input.has_member("method") ? input.get_string_member("method") : "dnewton") {
                case "dnewton":
                    solver = new MultirootFsolver(MultirootFsolverTypes.dnewton, nDim);
                    break;
                case "broyden":
                    solver = new MultirootFsolver(MultirootFsolverTypes.broyden, nDim);
                    break;
                case "hybrid":
                    solver = new MultirootFsolver(MultirootFsolverTypes.hybrid, nDim);
                    break;
                case "hybrids":
                    solver = new MultirootFsolver(MultirootFsolverTypes.hybrids, nDim);
                    break;
                default:
                    stderr.printf("ERROR: Invalid solver name!");
                    return_if_reached();
            }

            // configure the solver
            MultirootFunction objectiveFunctionData;
            Params parameters;
            if (input.has_member("r") && ! input.has_member("rMin") && ! input.has_member("rMax")) {
                parameters = initializeParams(input, input.get_double_member("r"), input.get_double_member("r"));
                objectiveFunctionData = MultirootFunction() {
                    f = sphericalOrbit,
                    n = initialValues.size,
                    params = &parameters
                };
            } else if (! input.has_member("r") && input.has_member("rMin") && input.has_member("rMax")) {
                parameters = initializeParams(input, input.get_double_member("rMin"), input.get_double_member("rMax"));
                objectiveFunctionData = MultirootFunction() {
                    f = nonSphericalOrbit,
                    n = initialValues.size,
                    params = &parameters
                };
            } else {
                stderr.printf("ERROR: Invalid radius constraint: use either r alone or rMin & rMax together");
                return_if_reached();
            }
            solver.set(&objectiveFunctionData, initialValues);

            // run the solver
            var epsabs = input.has_member("epsabs") ? input.get_double_member("epsabs") : 1.0e-12;
            var epsrel = input.has_member("epsrel") ? input.get_double_member("epsrel") : 1.0e-12;
            var maxIterations = input.has_member("maxIterations") ? input.get_int_member("maxIterations") : 1000;
            var termination = input.has_member("termination") ? input.get_string_member("termination") : "deltas";
            bool continuing = true;
            var iterations = 0;
            do {
                iterations++;
                var solverStatus = solver.iterate();
                if (solverStatus == Status.ENOPROG || solverStatus == Status.EBADFUNC) {
                    break;
                }
                var residualStatus = MultirootTest.residual(solver.f, epsabs);
                var deltaStatus = MultirootTest.delta(solver.dx, solver.x, epsabs, epsrel);
                switch (termination) {
                    case "residuals":
                        continuing = residualStatus == Status.CONTINUE;
                        break;
                    case "deltas":
                        continuing = deltaStatus == Status.CONTINUE;
                        break;
                    case "either":
                        continuing = residualStatus == Status.CONTINUE && deltaStatus == Status.CONTINUE;
                        break;
                    case "both":
                        continuing = residualStatus == Status.CONTINUE || deltaStatus == Status.CONTINUE;
                        break;
                    default:
                        stderr.printf("ERROR: Invalid termination criteria!");
                        return_if_reached();
                }
            } while (continuing && iterations < maxIterations);

            // generate output
            printOutput(solver, iterations);
        }
    }

    public static void main (string[] args) {
        new IcGenerator().generate(Simulations.getJson());
    }
}
