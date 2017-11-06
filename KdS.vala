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

/**
 * Common entry point for all simulators based on symplectic integration, reads in JSON data from stdin and executes a program according to the following table:
 *
 * || ''program name'' || ''key JSON variable'' || ''class(es))'' || ''simulation'' ||
 * || Simulate || "Newton" || {@link Simulations.Newton} || Newtonian central body problem ||
 * || Simulate || "KerrDeSitter" || {@link Simulations.Bh3d} || Kerr-deSitter geodesics ||
 * || Simulate || "NBody" || {@link Simulations.NBody} {@link Simulations.Body} || Newtonian N body problem ||
 * || GenParticle || N/A || {@link Generators.Particle} || Generate initial conditions for Kerr-deSitter particle geodesics ||
 * || GenLight || N/A || {@link Generators.Light} || Generate initial conditions for Kerr-deSitter light geodesics ||
 *
 * @param args the program is hard linked, so args[0] is used as the program name
 */
public static void main (string[] args) {
    var executable = args[0];
    stderr.printf("Executable: %s\n", executable);
    if ("Simulate" in executable) {
        var content = JsonParser.getJson();
        var simulator = content.get_string_member("Simulator");
        var o = content.get_object_member("IC");
        var type = o.get_string_member("integrator");
        var step = o.get_double_member("step");
        var start = o.get_double_member("start");
        var end = o.get_double_member("end");
        var plotratio = o.get_int_member("plotratio");
        if (("b1" == type) || ("b2" == type) || ("sb4" == type) || ("sb6" == type) || ("sb8" == type) ||
            ("yb4" == type) || ("yb6" == type) || ("yb8" == type) || ("kl6" == type) || ("kl8" == type)) {
            switch (simulator) {
                case "Oscillator":
                    var model = new Simulations.Oscillator(o.get_double_member("m"), o.get_double_member("k"), o.get_double_member("x"));
                    model.solve(Integrators.getIntegrator(model, step, type), step, start, end, plotratio);
                    break;
                case "Pendulum":
                    var model = new Simulations.Pendulum(o.get_double_member("g"),
                                                         o.get_double_member("m"),
                                                         o.get_double_member("length"),
                                                         o.get_double_member("angle"));
                    model.solve(Integrators.getIntegrator(model, step, type), step, start, end, plotratio);
                    break;
                case "Newton":
                    var model = new Simulations.Newton(o.get_double_member("Lfac"), o.get_double_member("r0"));
                    model.solve(Integrators.getIntegrator(model, step, type), step, start, end, plotratio);
                    break;
                case "HenonHeiles":
                    var model = new Simulations.HenonHeiles(o.get_double_member("x0"), o.get_double_member("y0"));
                    model.solve(Integrators.getIntegrator(model, step, type), step, start, end, plotratio);
                    break;
                case "NBody":
                    Simulations.Body[] bodies = {};
                    foreach (var node in o.get_array_member("bodies").get_elements()) {
                        var body = node.get_object();
                        var m = body.get_double_member("mass");
                        if (body.has_member("pX") && body.has_member("pY") && body.has_member("pZ")) {
                            bodies += new Simulations.Body(body.get_double_member("qX"),
                                                           body.get_double_member("qY"),
                                                           body.get_double_member("qZ"),
                                                           body.get_double_member("pX"),
                                                           body.get_double_member("pY"),
                                                           body.get_double_member("pZ"),
                                                           m);
                        } else if (body.has_member("vX") && body.has_member("vY") && body.has_member("vZ")) {
                            bodies += new Simulations.Body(body.get_double_member("qX"),
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
                    var model = new Simulations.NBody(bodies, o.get_double_member("g"), o.get_double_member("errorLimit"));
                    model.solve(Integrators.getIntegrator(model, step, type), step, start, end, plotratio);
                    break;
                case "KerrDeSitter":
                    var model = new Simulations.Bh3d(o.get_double_member("lambda"),
                                           o.get_double_member("a"),
                                           o.get_double_member("mu"),
                                           o.has_member("Efac") ? o.get_double_member("Efac")*o.get_double_member("E") : o.get_double_member("E"),
                                           o.has_member("Lfac") ? o.get_double_member("Lfac")*o.get_double_member("L") : o.get_double_member("L"),
                                           o.has_member("Qfac") ? o.get_double_member("Qfac")*o.get_double_member("Q") : o.get_double_member("Q"),
                                           o.get_double_member("r0"),
                                           o.get_double_member("th0"),
                                           o.get_boolean_member("cross"));
                    model.solve(Integrators.getIntegrator(model, step, type), step, start, end, plotratio);
                    break;
                default:
                    stderr.printf("Bad simulator; should be [ Newton | KerrDeSitter | NBody ], found {%s}\n", simulator);
                    assert_not_reached();
             }
        } else {
            stderr.printf("Bad integrator; should be [ b1 | b2 | sb4 | sb6 | sb8 | yb4 | yb6 | yb8 | kl6 | kl8 ], found {%s}\n", type);
            assert_not_reached();
        }
    } else if ("GenParticle" in executable) {
        var json = JsonParser.getJson();
        if (json.has_member("IC")) {
            new Generators.Particle().printPotentials(json.get_object_member("IC"));
        } else {
            new Generators.Particle().generateInitialConditions(json);
        }
    } else if ("GenLight" in executable) {
        var json = JsonParser.getJson();
        if (json.has_member("IC")) {
            new Generators.Light().printPotentials(json.get_object_member("IC"));
        } else {
            new Generators.Light().generateInitialConditions(json);
        }
    } else {
        stderr.printf("Executable not recognized: %s\n", executable);
    }
}
