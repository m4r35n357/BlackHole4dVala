/*
Copyright (c) 2014, 2015, 2016, 2017 Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
using Json;
using GLib.Math;

namespace Simulations {

    /**
     * Common entry point for all simulators based on symplectic integration, reads in JSON data from stdin and executes a program according to the following table:
     *
     * || ''program name'' || ''key JSON variable'' || ''class(es))'' || ''simulation'' ||
     * || Simulate || "Lfac" || {@link Newton} || Newtonian central body problem ||
     * || Simulate || "a" || {@link BhSymp} || Kerr-deSitter geodesics ||
     * || Simulate || "bodies" || {@link NBody} {@link Body} || Newtonian N body problem ||
     * || GenParticle || none || {@link Generators.Particle} || Generate initial conditions for Kerr-deSitter particle geodesics ||
     * || GenLight || none || {@link Generators.Light} || Generate initial conditions for Kerr-deSitter light geodesics ||
     *
     * @param args the program is hard linked, so args[0] is used as the program name
     */
    public static void main (string[] args) {
        var executable = args[0];
        stderr.printf("Executable: %s\n", executable);
        if ("Simulate" in executable) {
            var o = getJson().get_object_member("IC");
            var type = o.get_string_member("integrator");
            if (("sb1" == type) || ("sb2" == type) || ("sb4" == type) || ("sb6" == type) || ("sb8" == type)) {
                if (o.has_member("Lfac")) {
                    var model = new Newton(o.get_double_member("Lfac"), o.get_double_member("r0"));
                    var integrator = Simulations.getIntegrator(model, o.get_double_member("step"), type);
                    model.solve(integrator, o.get_double_member("step"),
                                            o.get_double_member("start"),
                                            o.get_double_member("end"),
                                            o.get_int_member("plotratio"));
                } else if (o.has_member("a")) {
                    var model = new BhSymp(o.get_double_member("lambda"),
                                           o.get_double_member("a"),
                                           o.get_double_member("mu"),
                                           o.get_double_member("E"),
                                           o.get_double_member("L"),
                                           o.get_double_member("Q"),
                                           o.get_double_member("r0"),
                                           o.get_double_member("th0"),
                                           o.get_boolean_member("cross"));
                    var integrator = Simulations.getIntegrator(model, o.get_double_member("step"), type);
                    model.solve(integrator, o.get_double_member("step"),
                                            o.get_double_member("start"),
                                            o.get_double_member("end"),
                                            o.get_int_member("plotratio"));
                } else if (o.has_member("bodies")) {
                    Body[] bodies = {};
                    foreach (var node in o.get_array_member("bodies").get_elements()) {
                        var body = node.get_object();
                        if (body.has_member("pX") && body.has_member("pY") && body.has_member("pZ")) {
                            bodies += new Body(body.get_double_member("qX"),
                                               body.get_double_member("qY"),
                                               body.get_double_member("qZ"),
                                               body.get_double_member("pX"),
                                               body.get_double_member("pY"),
                                               body.get_double_member("pZ"),
                                               body.get_double_member("mass"));
                        } else if (body.has_member("vX") && body.has_member("vY") && body.has_member("vZ")) {
                            var m = body.get_double_member("mass");
                            bodies += new Body(body.get_double_member("qX"),
                                               body.get_double_member("qY"),
                                               body.get_double_member("qZ"),
                                               body.get_double_member("vX") * m,
                                               body.get_double_member("vY") * m,
                                               body.get_double_member("vZ") * m,
                                               m);
                        } else {
                            stderr.printf("Mixed use of momenta and velocity\n");
                        }
                    }
                    var model = new NBody(bodies, o.get_double_member("g"), o.get_double_member("errorLimit"));
                    var integrator = Simulations.getIntegrator(model, o.get_double_member("step"), type);
                    model.solve(integrator, o.get_double_member("step"),
                                            o.get_double_member("start"),
                                            o.get_double_member("end"),
                                            o.get_int_member("plotratio"));
                }
            } else {
                stderr.printf("Bad integrator; should be [ sb1 | sb2 | sb4 | sb6 | sb8 ], found {%s}\n", type);
                assert_not_reached();
            }
        } else if ("GenParticle" in executable) {
            var json = getJson();
            if (json.has_member("IC")) {
                new Generators.Particle().printPotentials(json.get_object_member("IC"));
            } else {
                new Generators.Particle().generateInitialConditions(json);
            }
        } else if ("GenLight" in executable) {
            var json = getJson();
            if (json.has_member("IC")) {
                new Generators.Light().printPotentials(json.get_object_member("IC"));
            } else {
                new Generators.Light().generateInitialConditions(json);
            }
        } else {
            stderr.printf("Executable not recognized: %s\n", executable);
        }
    }
}
