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
using Gsl;
using Json;

/**
 * Generators of initial conditions for simulating geodesics in the Kerr spacetime, and also of potential data
 */
namespace Generators {

    /**
     * Internal interface for the parameter generators
     */
    public interface IGenerator : GLib.Object {
        /**
         * Sets up and runs the solver
         *
         * @param input the JSON IC generation data
         */
        public abstract void generateInitialConditions (Json.Object input);

        /**
         * Writes the potential data to STDOUT for plotting
         *
         * @param input the generated JSON simulation data
         */
        public abstract void printPotentials (Json.Object input);
    }

    /**
    * Turns fixed parameters and other constraints into initial conditions data
    * suitable as input to {@link Simulations.Bh3d}
    * (passing some items directly through to the output data),
    * alternatively takes the data and creates potential plots from it.
     */
    public class Particle : IGenerator, GLib.Object {

        /**
         * The fixed parameters
         */
        private struct Params {
            public double a;
            public double mu2;
            public double rMin;
            public double rMax;
            public double elevation;
            public bool cross;
            public double start;
            public double end;
            public double step;
            public int64 plotratio;
            public string integrator;
        }

        /**
         * Array indices for the variables
         */
        private enum X {
            E = 0, L = 1, Q = 2;
        }

        /**
         * Array indices for the objective functions
         */
        private enum F {
            R1 = 0, R2 = 1, TH = 2;
        }

        /**
         * R potential
         */
        private static double R (double r, double E, double L, double Q, Params* p) {
            return R_base(r, E, L, Q, p->a, p->mu2);
        }

        private static double R_base (double r, double E, double L, double Q, double a, double mu2) {
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
            return THETA_base(theta, E, L, Q, p->a, p->mu2);
        }

        private static double THETA_base (double theta, double E, double L, double Q, double a, double mu2) {
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
         * Write the initial conditions file to STDOUT
         */
        private void printInitialConditions (MultirootFsolver s, size_t iterations) {
            var p = (Params*) s.function.params;
            stdout.printf("{\n");
            stdout.printf("  \"Generator\" : {\n");
            stdout.printf("    \"name\" : \"icgenParticle %s\",\n", s.name());
            stdout.printf("    \"iterations\" : %zu,\n", iterations);
            stdout.printf("    \"residuals\" : \"R1: %.1e, R2: %.1e, TH: %.1e\",\n", s.f.get(F.R1), s.f.get(F.R2), s.f.get(F.TH));
            stdout.printf("    \"deltas\" : \"dE: %.1e, dL: %.1e, dQ: %.1e\",\n", s.dx.get(X.E), s.dx.get(X.L), s.dx.get(X.Q));
            if (p->a * s.x.get(X.L) >= 0.0) {
                stdout.printf("    \"direction\" : \"PROGRADE\"\n");
            } else {
                stdout.printf("    \"direction\" : \"RETROGRADE\"\n");
            }
            stdout.printf("  },\n");
            stdout.printf("  \"Simulator\" : \"KerrDeSitter\",\n");
            stdout.printf("  \"IC\" : {\n");
            stdout.printf("    \"M\" : %.17g,\n", 1.0);
            stdout.printf("    \"a\" : %.17g,\n", p->a);
            stdout.printf("    \"lambda\" : %.17g,\n", 0.0);
            stdout.printf("    \"mu\" : %.17g,\n", p->mu2);
            stdout.printf("    \"E\" : %.17g,\n", s.x.get(X.E));
            stdout.printf("    \"Efac\" : %.17g,\n", 1.0);
            stdout.printf("    \"L\" : %.17g,\n", s.x.get(X.L));
            stdout.printf("    \"Lfac\" : %.17g,\n", 1.0);
            stdout.printf("    \"Q\" : %.17g,\n", s.x.get(X.Q));
            stdout.printf("    \"Qfac\" : %.17g,\n", 1.0);
            stdout.printf("    \"r0\" : %.17g,\n", 0.5 * (p->rMin + p->rMax));
            stdout.printf("    \"th0\" : %.17g,\n", 0.0);
            stdout.printf("    \"cross\" : %s,\n", p->cross ? "true" : "false");
            stdout.printf("    \"start\" : %.17g,\n", p->start);
            stdout.printf("    \"end\" : %.17g,\n", p->end);
            stdout.printf("    \"step\" : %.17g,\n", p->step);
            stdout.printf("    \"plotratio\" : %.1d,\n", (int)p->plotratio);
            stdout.printf("    \"integrator\" : \"%s\"\n", p->integrator);
            stdout.printf("  }\n");
            stdout.printf("}\n");
        }

        /**
         * These are the quantities to vary: E, L & Q
         */
        private Vector initializeVariables (Json.Object input) {
            var initialValues = new Vector(3);
            initialValues.set(X.E, input.has_member("E0") ? input.get_double_member("E0") : 1.0);
            initialValues.set(X.L, input.has_member("L0") ? input.get_double_member("L0") : 5.0);
            initialValues.set(X.Q, input.has_member("Q0") ? input.get_double_member("Q0") : 0.0);
            return initialValues;
        }

        /**
         * These are the fixed quantities
         */
        private Params initializeParams (Json.Object input, double rMin, double rMax) {
            return Params() {
                mu2 = 1.0,
                rMin = rMin,
                rMax = rMax,
                elevation = (1.0 - (input.has_member("elevation") ? (90.0 - input.get_double_member("elevation")) / 180.0 : 0.5)) * PI,
                a = input.has_member("spin") ? input.get_double_member("spin") : 0.0,
                cross = input.has_member("cross") ? input.get_boolean_member("cross") : false,
                start = input.has_member("start") ? input.get_double_member("start") : 0.0,
                end = input.has_member("end") ? input.get_double_member("end") : 1000.0,
                step = input.has_member("step") ? input.get_double_member("step") : 0.01,
                plotratio = input.has_member("plotratio") ? input.get_int_member("plotratio") : 1,
                integrator = input.has_member("integrator") ? input.get_string_member("integrator") : "sb2"
            };
        }

        /**
         * {@inheritDoc}
         * @see IGenerator.generateInitialConditions
         */
        public void generateInitialConditions (Json.Object input) {
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
                    stderr.printf("\nERROR: Invalid solver name!\n");
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
                stderr.printf("\nERROR: Invalid radius constraint: use either r alone or rMin & rMax together\n");
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
                        stderr.printf("\nERROR: Invalid termination criteria!\n");
                        return_if_reached();
                }
            } while (continuing && iterations < maxIterations);

            // generate output
            printInitialConditions(solver, iterations);
        }

        /**
         * {@inheritDoc}
         * @see IGenerator.printPotentials
         */
        public void printPotentials (Json.Object input) {
            var a = input.get_double_member("a");
            var mu2 = input.get_double_member("mu");
            var E = input.get_double_member("E");
            var L = input.get_double_member("L");
            var Q = input.get_double_member("Q");
            var rMax = input.get_double_member("r0") * 2.0;
            for (var x = 1; x < 1000; x++) {
                var xValue = 1.0 * x / 1001;
                stdout.printf("{ \"x\" : %.6f, \"R\" : %.6f, \"y\" : %.6f, \"THETA\" : %.6f }\n",
                    xValue * rMax, R_base(xValue * rMax, E, L, Q, a, mu2), xValue * PI, THETA_base(xValue * PI, E, L, Q, a, mu2));
            }
        }
    }

    /**
     * Generates initial conditions and potentials for light
     */
    public class Light : IGenerator, GLib.Object {

        private static double R (double r, double a, double E, double L, double Q) {
            return (E * (r * r + a * a) - a * L) * (E * (r * r + a * a) - a * L) - (r * r + a * a - 2.0 * r) * (Q + (L - a * E) * (L - a * E));
        }

        private static double THETA (double theta, double a, double E, double L, double Q) {
            return Q - cos(theta) * cos(theta) * (L * L / (sin(theta) * sin(theta)) - a * a * E * E);
        }

        private double r1 (double a) {
            return 2.0 * (1.0 + cos(2.0 / 3.0 * acos(- fabs(a))));
        }

        private double r2 (double a) {
            return 2.0 * (1.0 + cos(2.0 / 3.0 * acos(fabs(a))));
        }

        private double L (double r, double a) {
            return - (r * r * r - 3.0 * r * r + a * a * r + a * a) / (a * (r - 1.0));
        }

        private double Q (double r, double a) {
            return - r * r * r * (r * r * r - 6.0 * r * r + 9.0 * r - 4.0 * a * a) / (a * a * (r - 1.0) * (r - 1.0));
        }

        /**
         * Write the initial conditions file to STDOUT and potential data to STDERR for plotting
         */
        private void printOutput (double r, double a, bool cross, double start, double end, double step, int64 plotratio, string integrator) {
            var E = 1.0;
            var L = L(r,a);
            var Q = Q(r,a);
            stdout.printf("{\n");
            stdout.printf("  \"Generator\" : {\n");
            stdout.printf("    \"name\" : \"icgenLight\",\n");
            if (a * L >= 0.0) {
                stdout.printf("    \"direction\" : \"PROGRADE\"\n");
            } else {
                stdout.printf("    \"direction\" : \"RETROGRADE\"\n");
            }
            stdout.printf("  },\n");
            stdout.printf("  \"Simulator\" : \"KerrDeSitter\",\n");
            stdout.printf("  \"IC\" : {\n");
            stdout.printf("    \"M\" : %.17g,\n", 1.0);
            stdout.printf("    \"a\" : %.17g,\n", a);
            stdout.printf("    \"lambda\" : %.17g,\n", 0.0);
            stdout.printf("    \"mu\" : %.17g,\n", 0.0);
            stdout.printf("    \"E\" : %.17g,\n", E);
            stdout.printf("    \"L\" : %.17g,\n", L);
            stdout.printf("    \"Q\" : %.17g,\n", Q);
            stdout.printf("    \"r0\" : %.17g,\n", r);
            stdout.printf("    \"r1\" : %.17g,\n", r1(a));
            stdout.printf("    \"r2\" : %.17g,\n", r2(a));
            stdout.printf("    \"th0\" : %.17g,\n", 0.0);
            stdout.printf("    \"cross\" : %s,\n", cross ? "true" : "false");
            stdout.printf("    \"start\" : %.17g,\n", start);
            stdout.printf("    \"end\" : %.17g,\n", end);
            stdout.printf("    \"step\" : %.17g,\n", step);
            stdout.printf("    \"plotratio\" : %.1d,\n", (int)plotratio);
            stdout.printf("    \"integrator\" : \"%s\"\n", integrator);
            stdout.printf("  }\n");
            stdout.printf("}\n");
        }

        /**
         * {@inheritDoc}
         * @see IGenerator.generateInitialConditions
         */
        public void generateInitialConditions (Json.Object input) {
            // generate output
            printOutput(input.has_member("r") ? input.get_double_member("r") : 3.0,
                        input.has_member("spin") ? input.get_double_member("spin") : 1.0,
                        input.has_member("cross") ? input.get_boolean_member("cross") : false,
                        input.has_member("start") ? input.get_double_member("start") : 0.0,
                        input.has_member("end") ? input.get_double_member("end") : 1000.0,
                        input.has_member("step") ? input.get_double_member("step") : 0.001,
                        input.has_member("plotratio") ? input.get_int_member("plotratio") : 50,
                        input.has_member("integrator") ? input.get_string_member("integrator") : "sb2");
        }

        /**
         * {@inheritDoc}
         * @see IGenerator.printPotentials
         */
        public void printPotentials (Json.Object input) {
            var a = input.get_double_member("a");
            var E = input.get_double_member("E");
            var L = input.get_double_member("L");
            var Q = input.get_double_member("Q");
            for (var x = 1; x <= 1001; x++) {
                var xValue = 1.0 * x / 1001;
                stdout.printf("{ \"x\" : %.6f, \"R\" : %.6f, \"y\" : %.6f, \"THETA\" : %.6f }\n",
                                xValue * r2(a) * 1.1, R(xValue * r2(a) * 1.1, a, E, L, Q), xValue * PI, THETA(xValue * PI, a, E, L, Q));
            }
        }
    }
}
