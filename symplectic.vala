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
         * Simulation time step, defined at subclass construction
         */
        private double h;

        /**
         * Suzuki composition coefficients
         */
        protected double[] x;
        protected double[] y;
        protected double[] z;

        /**
         * This is a protected constructor; use the static factory {@link getIntegrator} to obtain a concrete subclass.
         *
         * @param h the time step
         */
        protected Symplectic (double h) {
            this.h = h;
            var rt3 = pow(4.0, (1.0 / 3.0));
            var rt5 = pow(4.0, (1.0 / 5.0));
            var rt7 = pow(4.0, (1.0 / 7.0));
            this.x = { 1.0 / (4.0 - rt7), - rt7 / (4.0 - rt7) };
            this.y = { 1.0 / (4.0 - rt5), - rt5 / (4.0 - rt5) };
            this.z = { 1.0 / (4.0 - rt3), - rt3 / (4.0 - rt3) };
        }

        protected void base2 (IModel model, double s) {
            model.qUpdate(h * s * 0.5);
            model.pUpdate(h * s);
            model.qUpdate(h * s * 0.5);
        }

        protected void base4 (IModel model, double s) {
            this.base2(model, s * z[0]);
            this.base2(model, s * z[0]);
            this.base2(model, s * z[1]);
            this.base2(model, s * z[0]);
            this.base2(model, s * z[0]);
        }

        protected void base6 (IModel model, double s) {
            this.base4(model, s * y[0]);
            this.base4(model, s * y[0]);
            this.base4(model, s * y[1]);
            this.base4(model, s * y[0]);
            this.base4(model, s * y[0]);
        }

        /**
         * Subclasses should perform one integration step by executing alternating {@link IModel.qUpdate} and {@link IModel.pUpdate}
         * methods.
         *
         * @see ISymplectic.step
         */
        protected abstract void step (IModel model);

        /**
         * {@inheritDoc}
         * @see ISymplectic.getH
         */
        public double getH () {
            return h;
        }
    }

    /**
     * First-order symplectic integrator concrete subclass.  This integrator is NOT time symmetrical.
     */
    public class FirstOrder : Symplectic {

        /**
         * {@inheritDoc}
         * @see Symplectic.Symplectic
         */
        public FirstOrder (double h) {
            base(h);
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
        public override void step (IModel model) {
            model.qUpdate(1.0);
            model.pUpdate(1.0);
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
        public SecondOrder (double h) {
            base(h);
        }

        /**
         * 2nd order integration step  Performs the following calls on {@link IModel} per iteration:
         *
         * {{{
         * qUpdate(h/2)
         * pUpdate(h)
         * qUpdate(h/2)
         * }}}
         *
         * where h is the time step.
         *
         * @see Symplectic.step
         */
        public override void step (IModel model) {
            base2(model, 1.0);
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
        public FourthOrder (double h) {
            base(h);
        }

        /**
         * 4th order integration step
         *
         * @see Symplectic.step
         */
        public override void step (IModel model) {
            base4(model, 1.0);
        }
    }

    public class SixthOrder : Symplectic {

        /**
         * {@inheritDoc}
         * @see Symplectic.Symplectic
         */
        public SixthOrder (double h) {
            base(h);
        }

        /**
         * 6th order integration step
         *
         * @see Symplectic.step
         */
        public override void step (IModel model) {
            base6(model, 1.0);
        }
    }

    public class EightthOrder : Symplectic {

        /**
         * {@inheritDoc}
         * @see Symplectic.Symplectic
         */
        public EightthOrder (double h) {
            base(h);
        }

        /**
         * 6th order integration step
         *
         * @see Symplectic.step
         */
        public override void step (IModel model) {
            this.base6(model, x[0]);
            this.base6(model, x[0]);
            this.base6(model, x[1]);
            this.base6(model, x[0]);
            this.base6(model, x[0]);
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
    public static ISymplectic getIntegrator (double h, string type) {
        switch (type) {
            case "sb1":
                stderr.printf("1st Order Symplectic Integrator\n");
                return new FirstOrder(h);
            case "sb2":
                stderr.printf("2nd Order Symplectic Integrator\n");
                return new SecondOrder(h);
            case "sb4":
                stderr.printf("4th Order Symplectic Integrator\n");
                return new FourthOrder(h);
            case "sb6":
                stderr.printf("6th Order Symplectic Integrator\n");
                return new SixthOrder(h);
            case "sb8":
                stderr.printf("8th Order Symplectic Integrator\n");
                return new EightthOrder(h);
        }
        stderr.printf("Integrator not recognized: %s\n", type);
        assert_not_reached();
    }
}
