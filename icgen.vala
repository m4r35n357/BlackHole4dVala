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

        private struct IcGenParam {
            public double mu2;
            public double rMin;
            public double rMax;
            public double thMin;
            public double a;
            public double Efac;
            public double Lfac;
            public double Qfac;
        }

        private enum ConstantOfMotion {
            E, L, Q;
        }

        private enum Objective {
            R1, R2, THETA;
        }

        private static double R (double r, double E, double L, double Q, void* params) {
            var a = ((IcGenParam*) params) -> a;
            var mu2 = ((IcGenParam*) params) -> mu2;
            return (E * (r * r + a * a) - a * L) * (E * (r * r + a * a) - a * L)
                    - (r * r + a * a - 2.0 * r) * (mu2 * r * r + Q + (L - a * E) * (L - a * E));
        }

        private static double dRdr (double r, double E, double L, double Q, void* params) {
            var a = ((IcGenParam*) params) -> a;
            var mu2 = ((IcGenParam*) params) -> mu2;
            return 4.0 * r * E * (E * (r * r + a * a) - a * L)
                    - (2.0 * r - 2.0) * (mu2 * r * r + Q + (L - a * E) * (L - a * E)) - 2.0 * mu2 * r * (r * r + a * a - 2.0 * r);
        }

        private static double THETA (double theta, double E, double L, double Q, void* params) {
            var a = ((IcGenParam*) params) -> a;
            var mu2 = ((IcGenParam*) params) -> mu2;
            return Q - cos(theta) * cos(theta) * (a * a * (mu2 - E * E) + L * L / (sin(theta) * sin(theta)));
        }

        private static int sphericalOrbit (Vector x, void* params, Vector f) {
            var E = x.get(ConstantOfMotion.E);
            var L = x.get(ConstantOfMotion.L);
            var Q = x.get(ConstantOfMotion.Q);

            f.set(Objective.R1, R(((IcGenParam*) params) -> rMax, E, L, Q, params));
            f.set(Objective.R2, dRdr(((IcGenParam*) params) -> rMax, E, L, Q, params));
            f.set(Objective.THETA, THETA(((IcGenParam*) params) -> thMin, E, L, Q, params));

            return Status.SUCCESS;
        }

        private static int nonSphericalOrbit (Vector x, void* params, Vector f) {
            var E = x.get(ConstantOfMotion.E);
            var L = x.get(ConstantOfMotion.L);
            var Q = x.get(ConstantOfMotion.Q);

            f.set(Objective.R1, R(((IcGenParam*) params) -> rMin, E, L, Q, params));
            f.set(Objective.R2, R(((IcGenParam*) params) -> rMax, E, L, Q, params));
            f.set(Objective.THETA, THETA(((IcGenParam*) params) -> thMin, E, L, Q, params));

            return Status.SUCCESS;
        }

        private void print_potential (MultirootFsolver s, void* params) {
            var rMax = ((IcGenParam*) params) -> rMax;
            var E = s.x.get(ConstantOfMotion.E) * ((IcGenParam*) params) -> Efac;
            var L = s.x.get(ConstantOfMotion.L) * ((IcGenParam*) params) -> Lfac;
            var Q = s.x.get(ConstantOfMotion.Q) * ((IcGenParam*) params) -> Qfac;
            for (var x = 1; x <= 1001; x++) {
                var xValue = 1.0 * x / 1001;
                stderr.printf("{ \"x\" : %.6f, \"R\" : %.6f, \"THETA\" : %.6f }\n",
                                xValue * rMax * 1.1, R(xValue * rMax * 1.1, E, L, Q, params), THETA(xValue * PI, E, L, Q, params));
            }
        }

        private void print_inital_conditions (MultirootFsolver s, void* params, size_t iterations) {
            stdout.printf("{ \"solver\" : \"%s\",\n", s.name());
            stdout.printf("  \"iterations\" : %zu,\n", iterations);
            stdout.printf("  \"errors\" : \"%.3e %.3e %.3e\",\n", s.f.get(Objective.R1), s.f.get(Objective.R2), s.f.get(Objective.THETA));
            stdout.printf("  \"M\" : %.1f,\n", 1.0);
            stdout.printf("  \"a\" : %.1f,\n", ((IcGenParam*) params) -> a);
            stdout.printf("  \"mu\" : %.1f,\n", ((IcGenParam*) params) -> mu2);
            stdout.printf("  \"E\" : %.17g,\n", s.x.get(ConstantOfMotion.E) * ((IcGenParam*) params) -> Efac);
            stdout.printf("  \"Lz\" : %.17g,\n", s.x.get(ConstantOfMotion.L) * ((IcGenParam*) params) -> Lfac);
            stdout.printf("  \"C\" : %.17g,\n", s.x.get(ConstantOfMotion.Q) * ((IcGenParam*) params) -> Qfac);
            stdout.printf("  \"r\" : %.1f,\n", 0.5 * (((IcGenParam*) params) -> rMin + ((IcGenParam*) params) -> rMax));
            stdout.printf("  \"theta\" : %.9f,\n", 0.5 * PI);
            stdout.printf("  \"start\" : %.1f,\n", 0.0);
            stdout.printf("  \"duration\" : %.1f,\n", 5000.0);
            stdout.printf("  \"step\" : %.3f,\n", 0.001);
            stdout.printf("  \"plotratio\" : %.1d,\n", 500);
            stdout.printf("  \"integrator\" : \"%s\"\n", "sc4");
            stdout.printf("}\n");
        }

        public void generate (Json.Object input) {
            size_t nDim = 3;

            var initialValues = new Vector(nDim);
            initialValues.set(ConstantOfMotion.E, 1.0);
            initialValues.set(ConstantOfMotion.L, 5.0);
            initialValues.set(ConstantOfMotion.Q, 0.0);

            MultirootFsolver solver;
            switch (input.get_string_member("method")) {
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
                    return_if_reached();
            }

            IcGenParam parameters;
            MultirootFunction objectiveFunctionData;
            switch (input.get_size()) {
                case 8:
                    parameters = IcGenParam() {
                        mu2 = 1.0,
                        rMin = input.get_double_member("rMin"),
                        rMax = input.get_double_member("rMax"),
                        thMin = (1.0 - input.get_double_member("thMin")) * PI,
                        a = input.get_double_member("spin"),
                        Efac = input.get_double_member("Efac"),
                        Lfac = input.get_double_member("Lfac"),
                        Qfac = input.get_double_member("Qfac")
                    };
                    objectiveFunctionData = MultirootFunction() {
                        f = nonSphericalOrbit,
                        n = nDim,
                        params = &parameters
                    };
                    break;
                case 7:
                    parameters = IcGenParam() {
                        mu2 = 1.0,
                        rMin = input.get_double_member("r"),
                        rMax = input.get_double_member("r"),
                        thMin = (1.0 - input.get_double_member("thMin")) * PI,
                        a = input.get_double_member("spin"),
                        Efac = input.get_double_member("Efac"),
                        Lfac = input.get_double_member("Lfac"),
                        Qfac = input.get_double_member("Qfac")
                    };
                    objectiveFunctionData = MultirootFunction() {
                        f = sphericalOrbit,
                        n = nDim,
                        params = &parameters
                    };
                    break;
                default:
                    return_if_reached();
            }
            solver.set(&objectiveFunctionData, initialValues);

            int status = 0;
            size_t iterations = 0;
            do {
                iterations++;
                status = solver.iterate();
                if ((bool) status) {
                    break;
                }
                status = MultirootTest.residual(solver.f, 1.0e-12);
            } while (status == Status.CONTINUE && iterations < 1000);

            print_inital_conditions(solver, &parameters, iterations);
            print_potential(solver, &parameters);
        }

        public static Json.Object fromJson () {
            var input = new StringBuilder();
            var buffer = new char[1024];
            while (!stdin.eof()) {
                var chunk = stdin.gets(buffer);
                if (chunk != null) {
                    input.append(chunk);
                }
            }
            unowned Json.Object o;
            var p = new Parser();
            try {
                p.load_from_data(input.str);
                o = p.get_root().get_object();
            } catch (GLib.Error e) {
                stderr.printf("Unable to parse the input data: %s\n", e.message);
                return_if_reached();
            }
            return o;
        }
    }

    public static void main (string[] args) {
        new IcGenerator().generate(IcGenerator.fromJson());
    }
}
