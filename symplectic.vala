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
         * Signature for composable base methods
         */
        private delegate void BaseMethod (double s);

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
        private double yOuter;
        private double xOuter;
        private double yCentral;
        private double xCentral;

        /**
         * Base method coefficients
         */
        private double[] c_d;

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
         * || ''scheme'' || ''Composition Method'' ||  ''Description'' ||
         * || "yoshida" || {@link composeYoshida} || Yoshida Composition ||
         * || "suzuki" || {@link composeSuzuki} || Suzuki Composition ||
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
                    stderr.printf("6th Order (Suzuki Composition)\n");
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
            var zOuter = 1.0 / (4.0 - pow(4.0, (1.0 / 3.0)));
            yOuter = 1.0 / (4.0 - pow(4.0, (1.0 / 5.0)));
            xOuter = 1.0 / (4.0 - pow(4.0, (1.0 / 7.0)));
            var zCentral = 1.0 - 4.0 * zOuter;
            yCentral = 1.0 - 4.0 * yOuter;
            xCentral = 1.0 - 4.0 * xOuter;
            c_d = { 0.5 * h * zOuter, h * zOuter, h * zOuter, h * zOuter, 0.5 * h * (zOuter + zCentral), h * zCentral };
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
         * Suzuki composition
         *
         * Performs the following calls to {@link BaseMethod} per iteration:
         *
         * {{{
         * method(s * outer)
         * method(s * outer)
         * method(s * central)
         * method(s * outer)
         * method(s * outer)
         * }}}
         *
         * where outer = 1 / (4 - 4^(1/X)), central = 1 - 4 * outer, and X = 3, 5 or 7
         * @param method the method being composed to a higher order
         * @param s the current multipler
         * @param outer size of the forward steps
         * @param central size of the backward step
         */
        private void composeSuzuki (BaseMethod method, double s, double outer, double central) {
            method(s * outer);
            method(s * outer);
            method(s * central);
            method(s * outer);
            method(s * outer);
        }

        /**
         * Composition from 2nd order to 4th order.
         *
         * @param s the current multipler
         */
        private void base4 (double s) {
            model.qUpdate(s * c_d[0]);
            model.pUpdate(s * c_d[1]);
            model.qUpdate(s * c_d[2]);
            model.pUpdate(s * c_d[3]);
            model.qUpdate(s * c_d[4]);
            model.pUpdate(s * c_d[5]);
            model.qUpdate(s * c_d[4]);
            model.pUpdate(s * c_d[3]);
            model.qUpdate(s * c_d[2]);
            model.pUpdate(s * c_d[1]);
            model.qUpdate(s * c_d[0]);
        }

        /**
         * 4th order integration step.
         *
         * Calls {@link base4} with s = 1.
         */
        private void fourthOrder () {
            base4(1.0);
        }

        /**
         * Composition from 4th order to 6th order.
         *
         * @param s the current multipler
         */
        private void base6 (double s) {
            composeSuzuki(base4, s, yOuter, yCentral);
        }

        /**
         * 6th order integration step.
         *
         * Calls {@link base6} with s = 1.
         */
        private void sixthOrder () {
            base6(1.0);
        }

        /**
         * 8th order integration step.
         *
         * Calls {@link base8} with s = 1.
         */
        private void eightthOrder () {
            composeSuzuki(base6, 1.0, xOuter, xCentral);
        }
    }
}
