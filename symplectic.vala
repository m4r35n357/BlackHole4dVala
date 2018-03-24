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
        private double x1;
        private double x0;

        /**
         * Base method coefficients
         */
        private double[] cd_s4;
        private double[] cd_s6;

        /**
         * This constructor produces instances from its label and scheme arguments according to the following tables:
         *
         * || ''label'' || ''Integrator Order'' ||  ''Description'' ||
         * || "b1" || {@link firstOrder} || 1st Order, Symplectic, NOT Reversible ||
         * || "b2" || {@link secondOrder} || 2nd Order, Symplectic, Reversible ||
         * || "b4" || {@link fourthOrder} || 4th Order, Symplectic, Reversible ||
         * || "b6" || {@link sixthOrder} || 6th Order, Symplectic, Reversible ||
         * || "b8" || {@link eightthOrder} || 8th Order, Symplectic, Reversible ||
         *
         * @param model the model
         * @param h the time step
         * @param label the integrator order
         * @param scheme the composition scheme
         */
        public Symplectic (Models.IModel model, double h, string label, string scheme) {
            this.model = model;
            this.h = h;
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
                    integrator = fourthOrder;
                    break;
                case "b6":
                    stderr.printf("6th Order (Smith)\n");
                    integrator = sixthOrder;
                    break;
                case "b8":
                    stderr.printf("8th Order (Suzuki Composition)\n");
                    integrator = eightthOrder;
                    break;
                default:
                    stderr.printf("Integrator not recognized: %s\n", label);
                    assert_not_reached();
            }
            x1 = 1.0 / (4.0 - pow(4.0, (1.0 / 7.0)));
            x0 = 1.0 - 4.0 * x1;
            var y1 = 1.0 / (4.0 - pow(4.0, (1.0 / 5.0)));
            var y0 = 1.0 - 4.0 * y1;
            var z1 = 1.0 / (4.0 - pow(4.0, (1.0 / 3.0)));
            var z0 = 1.0 - 4.0 * z1;
            cd_s4 = { 0.5 * h * z1,
                      h * z1, h * z1, h * z1,
                      0.5 * h * (z1 + z0), h * z0 };
            cd_s6 = { 0.5 * h * z1 * y1,
                      h * z1 * y1, h * z1 * y1, h * z1 * y1,
                      0.5 * h * (z1 + z0) * y1, h * z0 * y1, 0.5 * h * (z0 + z1) * y1,
                      h * z1 * y1, h * z1 * y1, h * z1 * y1,
                      h * z1 * y1,
                      h * z1 * y1, h * z1 * y1, h * z1 * y1,
                      0.5 * h * (z1 + z0) * y1, h * z0 * y1, 0.5 * h * (z0 + z1) * y1,
                      h * z1 * y1, h * z1 * y1, h * z1 * y1,
                      0.5 * h * z1 * (y1 + y0),
                      h * z1 * y0, h * z1 * y0, h * z1 * y0,
                      0.5 * h * (z1 + z0) * y0, h * z0 * y0 };
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
         * Direct fourth order integrator (my own - Suzuki composition of Stormer-Verlet).
         *
         * Performs the following calls on {@link Models.IModel} per iteration:
         *
         * {{{
         * qUpdate(h * x1 / 2)
         * pUpdate(h * x1)
         * qUpdate(h * x1)
         * pUpdate(h * x1)
         * qUpdate(h * (x1 + x0) / 2 )
         * pUpdate(h * x0)
         * qUpdate(h * (x1 + x0) / 2 )
         * pUpdate(h * x1)
         * qUpdate(h * x1)
         * pUpdate(h * x1)
         * qUpdate(h * x1 / 2)
         * }}}
         */
        private void fourthOrder () {
            model.qUpdate(cd_s4[0]);
            model.pUpdate(cd_s4[1]);
            model.qUpdate(cd_s4[2]);
            model.pUpdate(cd_s4[3]);
            model.qUpdate(cd_s4[4]);
            model.pUpdate(cd_s4[5]);
            model.qUpdate(cd_s4[4]);
            model.pUpdate(cd_s4[3]);
            model.qUpdate(cd_s4[2]);
            model.pUpdate(cd_s4[1]);
            model.qUpdate(cd_s4[0]);
        }

        /**
         * Direct sixth order base method (my own - two Suzuki compositions of Stormer-Verlet).
         *
         * @param s the current multipler
         */
        private void base6 (double s) {
            model.qUpdate(s * cd_s6[0]);
            model.pUpdate(s * cd_s6[1]);
            model.qUpdate(s * cd_s6[2]);
            model.pUpdate(s * cd_s6[3]);
            model.qUpdate(s * cd_s6[4]);
            model.pUpdate(s * cd_s6[5]);
            model.qUpdate(s * cd_s6[6]);
            model.pUpdate(s * cd_s6[7]);
            model.qUpdate(s * cd_s6[8]);
            model.pUpdate(s * cd_s6[9]);
            model.qUpdate(s * cd_s6[10]);
            model.pUpdate(s * cd_s6[11]);
            model.qUpdate(s * cd_s6[12]);
            model.pUpdate(s * cd_s6[13]);
            model.qUpdate(s * cd_s6[14]);
            model.pUpdate(s * cd_s6[15]);
            model.qUpdate(s * cd_s6[16]);
            model.pUpdate(s * cd_s6[17]);
            model.qUpdate(s * cd_s6[18]);
            model.pUpdate(s * cd_s6[19]);
            model.qUpdate(s * cd_s6[20]);
            model.pUpdate(s * cd_s6[21]);
            model.qUpdate(s * cd_s6[22]);
            model.pUpdate(s * cd_s6[23]);
            model.qUpdate(s * cd_s6[24]);
            model.pUpdate(s * cd_s6[25]);
            model.qUpdate(s * cd_s6[24]);
            model.pUpdate(s * cd_s6[23]);
            model.qUpdate(s * cd_s6[22]);
            model.pUpdate(s * cd_s6[21]);
            model.qUpdate(s * cd_s6[20]);
            model.pUpdate(s * cd_s6[19]);
            model.qUpdate(s * cd_s6[18]);
            model.pUpdate(s * cd_s6[17]);
            model.qUpdate(s * cd_s6[16]);
            model.pUpdate(s * cd_s6[15]);
            model.qUpdate(s * cd_s6[14]);
            model.pUpdate(s * cd_s6[13]);
            model.qUpdate(s * cd_s6[12]);
            model.pUpdate(s * cd_s6[11]);
            model.qUpdate(s * cd_s6[10]);
            model.pUpdate(s * cd_s6[9]);
            model.qUpdate(s * cd_s6[8]);
            model.pUpdate(s * cd_s6[7]);
            model.qUpdate(s * cd_s6[6]);
            model.pUpdate(s * cd_s6[5]);
            model.qUpdate(s * cd_s6[4]);
            model.pUpdate(s * cd_s6[3]);
            model.qUpdate(s * cd_s6[2]);
            model.pUpdate(s * cd_s6[1]);
            model.qUpdate(s * cd_s6[0]);
        }

        /**
         * 6th order integration step.  Delegates to a 6th order base.
         *
         * Calls {@link base6} with s = 1.
         */
        private void sixthOrder () {
            base6(1.0);
        }

        /**
         * 8th order integration step using Suzuki composition from a 6th order base.
         *
         * Performs the following calls to {@link base6} per iteration:
         * {{{
         * base6(1 / (4 - 4^(1/7)))
         * base6(1 / (4 - 4^(1/7)))
         * base6(1 - 4 * 1 / (4 - 4^(1/7)))
         * base6(1 / (4 - 4^(1/7)))
         * base6(1 / (4 - 4^(1/7)))
         * }}}
         */
        private void eightthOrder () {
            base6(x1);
            base6(x1);
            base6(x0);
            base6(x1);
            base6(x1);
        }
    }
}
