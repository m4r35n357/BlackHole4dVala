/*
Copyright (c) 2014, 2015, 2016, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
using GLib;
using GLib.Math;
using Gsl;

namespace Sim {

    public class IcGenerator : GLib.Object {

        private struct IcGenParams {
            public double mu2;
            public double rMin;
            public double rMax;
            public double thMin;
            public double a;
        }

        private static double R (double r, double E, double L, double Q, void* params) {
            var a = ((IcGenParams*) params) -> a;
            var mu2 = ((IcGenParams*) params) -> mu2;
            return (E * (r * r + a * a) - a * L) * (E * (r * r + a * a) - a * L)
                    - (r * r + a * a - 2.0 * r) * (mu2 * r * r + Q + (L - a * E) * (L - a * E));
        }

        private static double dRdr (double r, double E, double L, double Q, void* params) {
            var a = ((IcGenParams*) params) -> a;
            var mu2 = ((IcGenParams*) params) -> mu2;
            return 4.0 * r * E * (E * (r * r + a * a) - a * L)
                    - (2.0 * r - 2.0) * (mu2 * r * r + Q + (L - a * E) * (L - a * E)) - 2.0 * mu2 * r * (r * r + a * a - 2.0 * r);
        }

        private static double THETA (double theta, double E, double L, double Q, void* params) {
            var a = ((IcGenParams*) params) -> a;
            var mu2 = ((IcGenParams*) params) -> mu2;
            return Q - cos(theta) * cos(theta) * (a * a * (mu2 - E * E) + L * L / (sin(theta) * sin(theta)));
        }

        private static int spherical (Vector x, void* params, Vector f) {
            var E = x.get(0);
            var L = x.get(1);
            var Q = x.get(2);

            f.set(0, R(((IcGenParams*) params) -> rMax, E, L, Q, params));
            f.set(1, dRdr(((IcGenParams*) params) -> rMax, E, L, Q, params));
            f.set(2, THETA(((IcGenParams*) params) -> thMin, E, L, Q, params));

            return Status.SUCCESS;
        }

        private static int nonSpherical (Vector x, void* params, Vector f) {
            var E = x.get(0);
            var L = x.get(1);
            var Q = x.get(2);

            f.set(0, R(((IcGenParams*) params) -> rMin, E, L, Q, params));
            f.set(1, R(((IcGenParams*) params) -> rMax, E, L, Q, params));
            f.set(2, THETA(((IcGenParams*) params) -> thMin, E, L, Q, params));

            return Status.SUCCESS;
        }

        private void print_potential (MultirootFsolver s, void* params) {
            var rMax = ((IcGenParams*) params) -> rMax;
            var E = s.x.get(0);
            var L = s.x.get(1);
            var Q = s.x.get(2);
            for (var x = 1; x <= 1001; x++) {
                var xValue = 1.0 * x / 1001;
                stderr.printf("{ \"x\" : %.6f, \"R\" : %.6f, \"THETA\" : %.6f }\n",
                                xValue * rMax * 1.1, R(xValue * rMax * 1.1, E, L, Q, params), THETA(xValue * PI, E, L, Q, params));
            }
        }

        private void print_inital_conditions (MultirootFsolver s, void* params, size_t iterations) {
            stdout.printf("{ \"solver\" : \"%s\",\n", s.name());
            stdout.printf("  \"iterations\" : %zu,\n", iterations);
            stdout.printf("  \"errors\" : \"%.3e %.3e %.3e\",\n", s.f.get(0), s.f.get(1), s.f.get(2));
            stdout.printf("  \"M\" : %.1f,\n", 1.0);
            stdout.printf("  \"a\" : %.1f,\n", ((IcGenParams*) params) -> a);
            stdout.printf("  \"mu\" : %.1f,\n", ((IcGenParams*) params) -> mu2);
            stdout.printf("  \"E\" : %.17g,\n", s.x.get(0));
            stdout.printf("  \"Lz\" : %.17g,\n", s.x.get(1));
            stdout.printf("  \"C\" : %.17g,\n", s.x.get(2));
            stdout.printf("  \"r\" : %.1f,\n", 0.5 * ((((IcGenParams*) params) -> rMin) + ((IcGenParams*) params) -> rMax));
            stdout.printf("  \"theta\" : %.9f,\n", 0.5 * PI);
            stdout.printf("  \"start\" : %.1f,\n", 0.0);
            stdout.printf("  \"duration\" : %.1f,\n", 5000.0);
            stdout.printf("  \"step\" : %.3f,\n", 0.001);
            stdout.printf("  \"plotratio\" : %.1d,\n", 500);
            stdout.printf("  \"integrator\" : \"%s\"\n", "sc4");
            stdout.printf("}\n");
        }

        public void generate (string[] args) {
            size_t nDim = 3;
            var x = new Vector(nDim);
            x.set(0, 1.0);  // E
            x.set(1, 5.0);  // L
            x.set(2, 0.0);  // Q

            MultirootFsolver solver;
            switch (args[1]) {
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
            IcGenParams params;
            MultirootFunction f;
            switch (args.length) {
                case 6:
                    params = { 1.0, double.parse(args[2]), double.parse(args[3]), (1.0 - double.parse(args[4])) * PI, double.parse(args[5]) };
                    f = MultirootFunction() { f = nonSpherical, n = nDim, params = &params };
                    break;
                case 5:
                    params = { 1.0, double.parse(args[2]), double.parse(args[2]), (1.0 - double.parse(args[3])) * PI, double.parse(args[4]) };
                    f = MultirootFunction() { f = spherical, n = nDim, params = &params };
                    break;
                default:
                    return_if_reached();
            }
            solver.set(&f, x);

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

            print_inital_conditions(solver, &params, iterations);
            print_potential(solver, &params);
        }
    }

    public static void main (string[] args) {
        new IcGenerator().generate(args);
    }
}
