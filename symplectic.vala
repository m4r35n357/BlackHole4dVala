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
     * Symplectic integrator abstract superclass, leaves integration method selection to subclasses
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
        private double zOuter = 1.0 / (4.0 - pow(4.0, (1.0 / 3.0)));
        private double yOuter = 1.0 / (4.0 - pow(4.0, (1.0 / 5.0)));
        private double xOuter = 1.0 / (4.0 - pow(4.0, (1.0 / 7.0)));
        private double zCentral;
        private double yCentral;
        private double xCentral;

        /**
         * This constructor produces instances from its label argument according to the following table:
         *
         * || ''label'' || ''Subclass'' ||  ''Description'' ||
         * || "b1" || {@link firstOrder} || 1st Order, Symplectic, NOT Reversible ||
         * || "b2" || {@link secondOrder} || 2nd Order, Symplectic, Reversible ||
         * || "b4" || {@link fourthOrder} || 4th Order, Symplectic, Reversible ||
         * || "b6" || {@link sixthOrder} || 6th Order, Symplectic, Reversible ||
         * || "b8" || {@link eightthOrder} || 8th Order, Symplectic, Reversible ||
         * @param model the model
         * @param h the time step
         * @param label the integrator order
         */
        public Symplectic (Models.IModel model, double h, string label) {
            this.model = model;
            this.h = h;
            zCentral = 1.0 - 4.0 * zOuter;
            yCentral = 1.0 - 4.0 * yOuter;
            xCentral = 1.0 - 4.0 * xOuter;
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
                    stderr.printf("4th Order (Suzuki composition)\n");
                    integrator = fourthOrder;
                    break;
                case "b6":
                    stderr.printf("6th Order (Suzuki composition)\n");
                    integrator = sixthOrder;
                    break;
                case "b8":
                    stderr.printf("8th Order (Suzuki composition)\n");
                    integrator = eightthOrder;
                    break;
                default:
                    stderr.printf("Integrator not recognized: %s\n", label);
                    break;
            }
        }

        /**
         * Euler-Cromer 1st order integrator.
         *
         * Performs the following calls on {@link Models.IModel} per iteration:
         *
         * {{{
         * qUpdate(s * h)
         * pUpdate(s * h)
         * }}}
         *
         * where h is the time step.
         */
        private void firstOrder () {
            model.qUpdate(h);
            model.pUpdate(h);
        }

        /**
         * Suzuki composition
         * @param baseMethod the method being composed into a higher order
         * @param s the current multipler
         * @param outer size of the forward steps
         * @param central size of the backward step
         */
        private void composeSuzuki (BaseMethod baseMethod, double s, double outer, double central) {
            baseMethod(s * outer);
            baseMethod(s * outer);
            baseMethod(s * central);
            baseMethod(s * outer);
            baseMethod(s * outer);
        }

        /**
         * Stormer-Verlet base integrator.
         *
         * Performs the following calls on {@link Models.IModel} per iteration:
         *
         * {{{
         * qUpdate(s * h/2)
         * pUpdate(s * h)
         * qUpdate(s * h/2)
         * }}}
         *
         * where h is the time step.
         */
        private void base2 (double s) {
            model.qUpdate(h * s * 0.5);
            model.pUpdate(h * s);
            model.qUpdate(h * s * 0.5);
        }

        /**
         * 2nd order integration step.
         *
         * Calls {@link base2} with s = 1.
         */
        private void secondOrder () {
            base2(1.0);
        }

        /**
         * Suzuki composition from 2nd order to 4th order.
         *
         * Performs the following calls to {@link base2} per iteration:
         *
         * {{{
         * base2(s * zOuter)
         * base2(s * zOuter)
         * base2(s * zCentral)
         * base2(s * zOuter)
         * base2(s * zOuter)
         * }}}
         *
         * where zOuter = 1 / (4 - 4**(1/3)), and zCentral = 1 - 4 * zOuter
         * @param s the current multipler
         */
        private void base4 (double s) {
            composeSuzuki(base2, s, zOuter, zCentral);
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
         * Suzuki composition from 4th order to 6th order.
         *
         * Performs the following calls to {@link base4} per iteration:
         *
         * {{{
         * base4(s * yOuter)
         * base4(s * yOuter)
         * base4(s * yCentral)
         * base4(s * yOuter)
         * base4(s * yOuter)
         * }}}
         *
         * where yOuter = 1 / (4 - 4**(1/5)), and yCentral = 1 - 4 * yOuter
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
         * Suzuki composition from 6th order to 8th order.
         *
         * Performs the following calls to {@link base6} per iteration:
         *
         * {{{
         * base6(s * xOuter)
         * base6(s * xOuter)
         * base6(s * xCentral)
         * base6(s * xOuter)
         * base6(s * xOuter)
         * }}}
         *
         * where xOuter = 1 / (4 - 4**(1/7)), and xCentral = 1 - 4 * xOuter
         * @param s the current multipler
         */
        private void base8 (double s) {
            composeSuzuki(base6, s, xOuter, xCentral);
        }

        /**
         * 8th order integration step.
         *
         * Calls {@link base8} with s = 1.
         */
        private void eightthOrder () {
            base8(1.0);
        }
    }
}
