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
     * External user interface for the physical model
     */
    public interface ISolver : GLib.Object {
        /**
         * Externally visible method. It sets up, controls and terminates the simulation,
         * writing its data to stdout.
         *
         * Sole method called by main(), calls {@link Integrators.Symplectic.integrator} on the selected integrator once per iteration
         * @param integrator the selected implementation
         * @param h the time step
         * @param start start time
         * @param end finish time
         * @param tr plot ratio i.e. only plot every tr-th point
         * @return an array of iteration counters
         */
        public abstract int64[] solve (Integrators.Symplectic.Integrator integrator, double h, double start, double end, int64 tr);
    }

    /**
     * Internal interface for the symplectic integrator to access and update the model
     */
    public interface IModel : GLib.Object {
        /**
         * Hamiltonian equations of motion - coordinate updates (dT/dp), called by {@link Integrators.Symplectic.integrator}
         * @param d scaled time step
         */
        public abstract void qUpdate (double d);

        /**
         * Hamiltonian equations of motion - momentum updates (dV/dq), called by {@link Integrators.Symplectic.integrator}
         * @param c scaled time step
         */
        public abstract void pUpdate (double c);
    }
}

namespace Integrators {

    /**
     * Symplectic integrator class
     */
    public class Symplectic : GLib.Object {

        /**
         * Signature for integrator methods
         */
        public delegate void Integrator ();

        /**
         * The integrator method
         */
        public Integrator integrator;

        /**
         * The physical model
         */
        private Models.IModel model;

        /**
         * Simulation time step,
         */
        private double h;

        /**
         * Composition coefficients
         */
        private double[] cd;

        /**
         * This constructor produces instances from its label and scheme arguments according to the following tables:
         *
         * || ''label'' || ''Integrator Order'' ||  ''Description'' ||
         * || "b1" || {@link firstOrder} || 1st Order, Symplectic, NOT Reversible ||
         * || "b2" || {@link secondOrder} || 2nd Order, Symplectic, Reversible ||
         * || "b4" || {@link smith} || 4th Order, Symplectic, Reversible ||
         * || "b6" || {@link smith} || 6th Order, Symplectic, Reversible ||
         * || "b8" || {@link smith} || 8th Order, Symplectic, Reversible ||
         *
         * @param model the model
         * @param h the time step
         * @param label the integrator order
         * @param scheme the composition scheme
         */
        public Symplectic (Models.IModel model, double h, string label, string scheme) {
            this.model = model;
            this.h = h;
            var x1 = 1.0 / (4.0 - pow(4.0, (1.0 / 7.0)));
            var x0 = 1.0 - 4.0 * x1;
            var y1 = 1.0 / (4.0 - pow(4.0, (1.0 / 5.0)));
            var y0 = 1.0 - 4.0 * y1;
            var z1 = 1.0 / (4.0 - pow(4.0, (1.0 / 3.0)));
            var z0 = 1.0 - 4.0 * z1;
            switch (label) {
                case "b1":
                    stderr.printf("1st Order (Euler-Cromer)\n");
                    integrator = firstOrder;
                    break;
                case "b2":
                    stderr.printf("2nd Order (Stormer-Verlet)\n");
                    integrator = secondOrder;
                    break;
                case "b4":
                    stderr.printf("4th Order (Smith)\n");
                    integrator = smith;
                    cd = {
                      0.5 * h * z1,
                      h * z1, h * z1, h * z1,
                      0.5 * h * (z1 + z0), h * z0
                    };
                    break;
                case "b6":
                    stderr.printf("6th Order (Smith)\n");
                    integrator = smith;
                    cd = {
                      0.5 * h * z1 * y1,
                      h * z1 * y1, h * z1 * y1, h * z1 * y1,
                      0.5 * h * (z1 + z0) * y1, h * z0 * y1, 0.5 * h * (z0 + z1) * y1,
                      h * z1 * y1, h * z1 * y1, h * z1 * y1,
                      h * z1 * y1,
                      h * z1 * y1, h * z1 * y1, h * z1 * y1,
                      0.5 * h * (z1 + z0) * y1, h * z0 * y1, 0.5 * h * (z0 + z1) * y1,
                      h * z1 * y1, h * z1 * y1, h * z1 * y1,
                      0.5 * h * z1 * (y1 + y0),
                      h * z1 * y0, h * z1 * y0, h * z1 * y0,
                      0.5 * h * (z1 + z0) * y0, h * z0 * y0
                    };
                    break;
                case "b8":
                    stderr.printf("8th Order (Smith)\n");
                    integrator = smith;
                    cd = {
                      0.5 * h * z1 * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      h * z1 * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * z1 * (y1 + y0) * x1,
                      h * z1 * y0 * x1, h * z1 * y0 * x1, h * z1 * y0 * x1,
                      0.5 * h * (z1 + z0) * y0 * x1, h * z0 * y0 * x1, 0.5 * h * (z0 + z1) * y0 * x1,
                      h * z1 * y0 * x1, h * z1 * y0 * x1, h * z1 * y0 * x1,
                      0.5 * h * z1 * (y1 + y0) * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      h * z1 * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      h * z1 * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      h * z1 * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * z1 * (y1 + y0) * x1,
                      h * z1 * y0 * x1, h * z1 * y0 * x1, h * z1 * y0 * x1,
                      0.5 * h * (z1 + z0) * y0 * x1, h * z0 * y0 * x1, 0.5 * h * (z0 + z1) * y0 * x1,
                      h * z1 * y0 * x1, h * z1 * y0 * x1, h * z1 * y0 * x1,
                      0.5 * h * z1 * (y1 + y0) * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      h * z1 * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * z1 * y1 * (x1 + x0),
                      h * z1 * y1 * x0, h * z1 * y1 * x0, h * z1 * y1 * x0,
                      0.5 * h * (z1 + z0) * y1 * x0, h * z0 * y1 * x0, 0.5 * h * (z0 + z1) * y1 * x0,
                      h * z1 * y1 * x0, h * z1 * y1 * x0, h * z1 * y1 * x0,
                      h * z1 * y1 * x0,
                      h * z1 * y1 * x0, h * z1 * y1 * x0, h * z1 * y1 * x0,
                      0.5 * h * (z1 + z0) * y1 * x0, h * z0 * y1 * x0, 0.5 * h * (z0 + z1) * y1 * x0,
                      h * z1 * y1 * x0, h * z1 * y1 * x0, h * z1 * y1 * x0,
                      0.5 * h * z1 * (y1 + y0) * x0,
                      h * z1 * y0 * x0, h * z1 * y0 * x0, h * z1 * y0 * x0,
                      0.5 * h * (z1 + z0) * y0 * x0, h * z0 * y0 * x0
                    };
                    break;
                default:
                    stderr.printf("Integrator not recognized: %s\n", label);
                    assert_not_reached();
            }
        }

        /**
         * Euler-Cromer 1st order integrator.
         *
         * Performs the following calls on {@link Models.IModel} per iteration:
         *
         * {{{
         * qUpdate(h)
         * pUpdate(h)
         * }}}
         *
         * where h is the time step.
         */
        private void firstOrder () {
            model.qUpdate(h);
            model.pUpdate(h);
        }

        /**
         * Stormer-Verlet 2nd order integrator.
         *
         * Performs the following calls on {@link Models.IModel} per iteration:
         *
         * {{{
         * qUpdate(h/2)
         * pUpdate(h)
         * qUpdate(h/2)
         * }}}
         *
         * where h is the time step.
         */
        private void secondOrder () {
            model.qUpdate(h * 0.5);
            model.pUpdate(h);
            model.qUpdate(h * 0.5);
        }

        /**
         * Direct integrator.
         *
         * Performs the following calls on {@link Models.IModel} per iteration:
         *
         * {{{
         * qUpdate(cd[0])
         * pUpdate(cd[1])
         * ...
         * pUpdate(cd[size - 1])
         * ...
         * pUpdate(cd[1])
         * qUpdate(cd[0])
         * }}}
         */
        private void smith () {
            var size = cd.length;
            for (int i = 0; i < size; i++) {
                if (i % 2 == 0) {
                    model.qUpdate(cd[i]);
                } else {
                    model.pUpdate(cd[i]);
                }
            }
            for (int i = size - 2; i > -1; i--) {
                if (i % 2 == 0) {
                    model.qUpdate(cd[i]);
                } else {
                    model.pUpdate(cd[i]);
                }
            }
        }
    }
}
