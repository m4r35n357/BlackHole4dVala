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
using Json;

/**
 * Generators of initial conditions for simulating geodesics in the Kerr spacetime, and also of potential data
 */
namespace Generators {

    /**
     * Generates initial conditions and potentials for light
     */
    public class Light : Simulations.IGenerator, GLib.Object {

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
        private void printOutput (double r, double a, double start, double end, double step, int64 plotratio, string integrator) {
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
            stdout.printf("  \"IC\" : {\n");
            stdout.printf("    \"M\" : %.1f,\n", 1.0);
            stdout.printf("    \"a\" : %.1f,\n", a);
            stdout.printf("    \"lambda\" : %.17g,\n", 0.0);
            stdout.printf("    \"mu\" : %.1f,\n", 0.0);
            stdout.printf("    \"E\" : %.17g,\n", E);
            stdout.printf("    \"L\" : %.17g,\n", L);
            stdout.printf("    \"Q\" : %.17g,\n", Q);
            stdout.printf("    \"r0\" : %.17g,\n", r);
            stdout.printf("    \"r1\" : %.3f,\n", r1(a));
            stdout.printf("    \"r2\" : %.3f,\n", r2(a));
            stdout.printf("    \"th0\" : %.0f,\n", 0.0);
            stdout.printf("    \"start\" : %.1f,\n", start);
            stdout.printf("    \"end\" : %.1f,\n", end);
            stdout.printf("    \"step\" : %.3f,\n", step);
            stdout.printf("    \"plotratio\" : %.1d,\n", (int)plotratio);
            stdout.printf("    \"integrator\" : \"%s\"\n", integrator);
            stdout.printf("  }\n");
            stdout.printf("}\n");
        }

        /**
         * {@inheritDoc}
         * @see Simulations.IGenerator.generateInitialConditions
         */
        public void generateInitialConditions (Json.Object input) {
            // generate output
            printOutput(input.has_member("r") ? input.get_double_member("r") : 3.0,
                        input.has_member("spin") ? input.get_double_member("spin") : 1.0,
                        input.has_member("start") ? input.get_double_member("start") : 0.0,
                        input.has_member("end") ? input.get_double_member("end") : 1000.0,
                        input.has_member("step") ? input.get_double_member("step") : 0.001,
                        input.has_member("plotratio") ? input.get_int_member("plotratio") : 50,
                        input.has_member("integrator") ? input.get_string_member("integrator") : "sb2");
        }

        /**
         * {@inheritDoc}
         * @see Simulations.IGenerator.printPotentials
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
