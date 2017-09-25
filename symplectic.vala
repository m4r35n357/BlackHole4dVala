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

namespace Simulations {

    /**
     * Symplectic integrator abstract superclass, leaves integration details to subclasses
     */
    protected abstract class Symplectic : ISymplectic, GLib.Object {

        /**
         * The physical model, defined at subclass construction
         */
        protected IModel model;

        /**
         * Simulation time step, defined at subclass construction
         */
        protected double h;

        /**
         * Suzuki composition coefficients
         */
        protected double x1;
        protected double y1;
        protected double z1;
        protected double x3;
        protected double y3;
        protected double z3;

        /**
         * This is a protected constructor; use the static factory {@link getIntegrator} to obtain a concrete subclass.
         *
         * @param h the time step
         */
        protected Symplectic (IModel model, double h) {
            this.model = model;
            this.h = h;
            x1 = 1.0 / (4.0 - pow(4.0, (1.0 / 7.0)));
            y1 = 1.0 / (4.0 - pow(4.0, (1.0 / 5.0)));
            z1 = 1.0 / (4.0 - pow(4.0, (1.0 / 3.0)));
            x3 = 1.0 - 4.0 * x1;
            y3 = 1.0 - 4.0 * y1;
            z3 = 1.0 - 4.0 * z1;
        }

        /**
         * Stormer-Verlet base integrator.  Performs the following calls on {@link IModel} per iteration:
         *
         * {{{
         * qUpdate(h/2)
         * pUpdate(h)
         * qUpdate(h/2)
         * }}}
         *
         * where h is the time step.
         */
        protected void base2 (double s) {
            model.qUpdate(h * s * 0.5);
            model.pUpdate(h * s);
            model.qUpdate(h * s * 0.5);
        }

        /**
         * Suzuki composition from 2nd order to 4th order
         */
        protected void base4 (double s) {
            base2(s * z1);
            base2(s * z1);
            base2(s * z3);
            base2(s * z1);
            base2(s * z1);
        }

        /**
         * Suzuki composition from 4th order to 6th order
         */
        protected void base6 (double s) {
            base4(s * y1);
            base4(s * y1);
            base4(s * y3);
            base4(s * y1);
            base4(s * y1);
        }

        /**
         * Suzuki composition from 6th order to 8th order
         */
        protected void base8 (double s) {
            base6(s * x1);
            base6(s * x1);
            base6(s * x3);
            base6(s * x1);
            base6(s * x1);
        }
        /**
         * Subclasses should perform one integration step by executing alternating {@link IModel.qUpdate} and {@link IModel.pUpdate}
         * methods.
         *
         * @see ISymplectic.step
         */
        protected abstract void step ();
    }

    /**
     * First-order symplectic integrator concrete subclass.  This integrator is NOT time symmetrical.
     */
    public class FirstOrder : Symplectic {

        /**
         * {@inheritDoc}
         * @see Symplectic.Symplectic
         */
        public FirstOrder (IModel model, double h) {
            base(model, h);
        }

        /**
         * 1st order integration step.  Performs the following calls on {@link IModel} per iteration:
         *
         * {{{
         * qUpdate(h)
         * pUpdate(h)
         * }}}
         *
         * where h is the time step.
         *
         * @see Symplectic.step
         */
        public override void step () {
            model.qUpdate(h);
            model.pUpdate(h);
        }
    }

    /**
     * Second-order symplectic integrator concrete subclass.  This integrator is time symmetrical.
     */
    public class SecondOrder : Symplectic {

        /**
         * {@inheritDoc}
         * @see Symplectic.Symplectic
         */
        public SecondOrder (IModel model, double h) {
            base(model, h);
        }

        /**
         * 2nd order integration step
         *
         * @see Symplectic.step
         */
        public override void step () {
            base2(1.0);
        }
    }

    /**
     * Fourth-order symplectic integrator concrete subclass.  This integrator is time symmetrical.
     */
    public class FourthOrder : Symplectic {

        /**
         * {@inheritDoc}
         * @see Symplectic.Symplectic
         */
        public FourthOrder (IModel model, double h) {
            base(model, h);
        }

        /**
         * 4th order integration step
         *
         * @see Symplectic.step
         */
        public override void step () {
            base4(1.0);
        }
    }

    /**
     * Sixth-order symplectic integrator concrete subclass.  This integrator is time symmetrical.
     */
    public class SixthOrder : Symplectic {

        /**
         * {@inheritDoc}
         * @see Symplectic.Symplectic
         */
        public SixthOrder (IModel model, double h) {
            base(model, h);
        }

        /**
         * 6th order integration step
         *
         * @see Symplectic.step
         */
        public override void step () {
            base6(1.0);
        }
    }

    /**
     * Eightth-order symplectic integrator concrete subclass.  This integrator is time symmetrical.
     */
    public class EightthOrder : Symplectic {

        /**
         * {@inheritDoc}
         * @see Symplectic.Symplectic
         */
        public EightthOrder (IModel model, double h) {
            base(model, h);
        }

        /**
         * 8th order integration step
         *
         * @see Symplectic.step
         */
        public override void step () {
            base8(1.0);
        }
    }

    /**
     * Static factory for producing subclass instances from its type argument according to the following table:
     *
     * || ''type'' || ''Subclass'' ||  ''Description'' ||
     * || "sb1" || {@link FirstOrder} || 1st Order, Symplectic ||
     * || "sb2" || {@link SecondOrder} || 2nd Order, Symplectic, Reversible ||
     * || "sb4" || {@link FourthOrder} || 4th Order, Symplectic, Reversible ||
     * || "sb6" || {@link SixthOrder} || 6th Order, Symplectic, Reversible ||
     * || "sb8" || {@link EightthOrder} || 8th Order, Symplectic, Reversible ||
     *
     * @param h the time step
     * @param type the selected implementation
     *
     * @return concrete implementation of a symplectic integrator
     */
    public static ISymplectic getIntegrator (IModel model, double h, string type) {
        switch (type) {
            case "sb1":
                stderr.printf("1st Order Symplectic Integrator\n");
                return new FirstOrder(model, h);
            case "sb2":
                stderr.printf("2nd Order Symplectic Integrator\n");
                return new SecondOrder(model, h);
            case "sb4":
                stderr.printf("4th Order Symplectic Integrator\n");
                return new FourthOrder(model, h);
            case "sb6":
                stderr.printf("6th Order Symplectic Integrator\n");
                return new SixthOrder(model, h);
            case "sb8":
                stderr.printf("8th Order Symplectic Integrator\n");
                return new EightthOrder(model, h);
        }
        stderr.printf("Integrator not recognized: %s\n", type);
        assert_not_reached();
    }
}
