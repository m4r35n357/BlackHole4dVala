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
     * External user interface for the physical model
     */
    public interface ISolver : GLib.Object {
        /**
         * Externally visible method. It sets up, controls and terminates the simulation,
         * writing its data to stdout.
         *
         * Sole method called by main(), calls {@link ISymplectic.step} on the selected integrator once per iteration
         *
         * @param integrator the selected implementation
         * @param h the time step
         * @param start start time
         * @param end finish time
         * @param tr plot ratio i.e. only plot every tr-th point
         *
         * @return an array of iteration counters
         */
        public abstract int64[] solve (ISymplectic integrator, double h, double start, double end, int64 tr);
    }

    /**
     * Internal interface for the symplectic integrator to access and update the model
     */
    public interface IModel : GLib.Object {
        /**
         * Hamiltonian equations of motion - coordinate updates (dT/dp), called by {@link ISymplectic.step}
         *
         * @param d scaled time step
         */
        public abstract void qUpdate (double d);

        /**
         * Hamiltonian equations of motion - momentum updates (dV/dq), called by {@link ISymplectic.step}
         *
         * @param c scaled time step
         */
        public abstract void pUpdate (double c);
    }

    /**
     * Internal interface for a model to drive the symplectic integrator
     */
    public interface ISymplectic : GLib.Object {
        /**
         * Should be called by {@link ISolver.solve} as needed, in turn calls {@link IModel.qUpdate} and {@link IModel.pUpdate}
         */
        public abstract void step ();
    }

    /**
     * Symplectic integrator abstract superclass, leaves integration method selection to subclasses
     */
    protected abstract class SymplecticBase : ISymplectic, GLib.Object {

        /**
         * The physical model, defined at subclass construction
         */
        protected IModel model;

        /**
         * Simulation time step, defined at subclass construction
         */
        protected double h;

        /**
         * This is a protected constructor; use the static factory {@link getIntegrator} to obtain a concrete subclass.
         *
         * @param h the time step
         */
        protected SymplecticBase (IModel model, double h) {
            this.model = model;
            this.h = h;
        }

        /**
         * Stormer-Verlet base integrator.  Performs the following calls on {@link IModel} per iteration:
         *
         * {{{
         * qUpdate(s * h/2)
         * pUpdate(s * h)
         * qUpdate(s * h/2)
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
    public class FirstOrder : SymplecticBase {

        /**
         * {@inheritDoc}
         * @see SymplecticBase.SymplecticBase
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
         * @see SymplecticBase.step
         */
        public override void step () {
            model.qUpdate(h);
            model.pUpdate(h);
        }
    }

    /**
     * Second-order symplectic integrator concrete subclass.  This integrator is time symmetrical.
     */
    public class SecondOrder : SymplecticBase {

        /**
         * {@inheritDoc}
         * @see SymplecticBase.SymplecticBase
         */
        public SecondOrder (IModel model, double h) {
            base(model, h);
        }

        /**
         * 2nd order integration step.  Calls {@link SymplecticBase.base2} with s = 1.
         *
         * @see SymplecticBase.step
         */
        public override void step () {
            base2(1.0);
        }
    }

    /**
     * Suzuki composed integrator abstract superclass, leaves integration method selection to subclasses
     */
    protected abstract class Suzuki : SymplecticBase {

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
        protected Suzuki (IModel model, double h) {
            base(model, h);
            x1 = 1.0 / (4.0 - pow(4.0, (1.0 / 7.0)));
            y1 = 1.0 / (4.0 - pow(4.0, (1.0 / 5.0)));
            z1 = 1.0 / (4.0 - pow(4.0, (1.0 / 3.0)));
            x3 = 1.0 - 4.0 * x1;
            y3 = 1.0 - 4.0 * y1;
            z3 = 1.0 - 4.0 * z1;
        }

        /**
         * Suzuki composition from 2nd order to 4th order.  Performs the following calls to {@link SymplecticBase.base2} per iteration:
         *
         * {{{
         * base2(s * z1)
         * base2(s * z1)
         * base2(s * z3)
         * base2(s * z1)
         * base2(s * z1)
         * }}}
         *
         * where z1 = 1 / (4 - 4**(1/3)), and z3 = 1 - 4 * z1
         */
        protected void base4 (double s) {
            base2(s * z1);
            base2(s * z1);
            base2(s * z3);
            base2(s * z1);
            base2(s * z1);
        }

        /**
         * Suzuki composition from 4th order to 6th order.  Performs the following calls per iteration:
         *
         * {{{
         * base4(s * y1)
         * base4(s * y1)
         * base4(s * y3)
         * base4(s * y1)
         * base4(s * y1)
         * }}}
         *
         * where y1 = 1 / (4 - 4**(1/5)), and y3 = 1 - 4 * y1
         */
        protected void base6 (double s) {
            base4(s * y1);
            base4(s * y1);
            base4(s * y3);
            base4(s * y1);
            base4(s * y1);
        }

        /**
         * Suzuki composition from 6th order to 8th order.  Performs the following calls per iteration:
         *
         * {{{
         * base6(s * x1)
         * base6(s * x1)
         * base6(s * x3)
         * base6(s * x1)
         * base6(s * x1)
         * }}}
         *
         * where x1 = 1 / (4 - 4**(1/7)), and x3 = 1 - 4 * x1
         */
        protected void base8 (double s) {
            base6(s * x1);
            base6(s * x1);
            base6(s * x3);
            base6(s * x1);
            base6(s * x1);
        }
    }

    /**
     * Yoshida composed integrator abstract superclass, leaves integration method selection to subclasses
     */
    protected abstract class Yoshida : SymplecticBase {

        /**
         * Yoshida composition coefficients
         */
        protected double x1;
        protected double y1;
        protected double z1;
        protected double x2;
        protected double y2;
        protected double z2;

        /**
         * This is a protected constructor; use the static factory {@link getIntegrator} to obtain a concrete subclass.
         *
         * @param h the time step
         */
        protected Yoshida (IModel model, double h) {
            base(model, h);
            x1 = 1.0 / (2.0 - pow(2.0, (1.0 / 7.0)));
            y1 = 1.0 / (2.0 - pow(2.0, (1.0 / 5.0)));
            z1 = 1.0 / (2.0 - pow(2.0, (1.0 / 3.0)));
            x2 = 1.0 - 2.0 * x1;
            y2 = 1.0 - 2.0 * y1;
            z2 = 1.0 - 2.0 * z1;
        }

        /**
         * Yoshida composition from 2nd order to 4th order.  Performs the following calls to {@link SymplecticBase.base2} per iteration:
         *
         * {{{
         * base2(s * z1)
         * base2(s * z2)
         * base2(s * z1)
         * }}}
         *
         * where z1 = 1 / (2 - 2**(1/3)), and z2 = 1 - 2 * z1
         */
        protected void base4 (double s) {
            base2(s * z1);
            base2(s * z2);
            base2(s * z1);
        }

        /**
         * Yoshida composition from 4th order to 6th order.  Performs the following calls per iteration:
         *
         * {{{
         * base4(s * y1)
         * base4(s * y2)
         * base4(s * y1)
         * }}}
         *
         * where y1 = 1 / (2 - 2**(1/5)), and y2 = 1 - 2 * y1
         */
        protected void base6 (double s) {
            base4(s * y1);
            base4(s * y2);
            base4(s * y1);
        }

        /**
         * Yoshida composition from 6th order to 8th order.  Performs the following calls per iteration:
         *
         * {{{
         * base6(s * x1)
         * base6(s * x2)
         * base6(s * x1)
         * }}}
         *
         * where x1 = 1 / (2 - 2**(1/7)), and x2 = 1 - 2 * x1
         */
        protected void base8 (double s) {
            base6(s * x1);
            base6(s * x2);
            base6(s * x1);
        }
    }

    /**
     * Generic composed integrator abstract superclass, leaves integration method selection to subclasses
     */
    protected abstract class GenericComposer : SymplecticBase {
        /**
         * Generic composition coefficients
         */
        protected double[] delta;

        /**
         * Size of coefficients array
         */
        protected int n;

        /**
         * This is a protected constructor; use the static factory {@link getIntegrator} to obtain a concrete subclass.
         *
         * @param h the time step
         */
        protected GenericComposer (IModel model, double h, double[] coefficients) {
            base(model, h);
            n = 2 * coefficients.length - 1;
            delta = new double[n];
            var m = coefficients.length - 1;
            for (var i = 0; i < m; i++) {
                delta[i] = coefficients[i];
                delta[n - 1 - i] = coefficients[i];
            }
            delta[m] = coefficients[m];
        }

        /**
         * Generic integration step.  Calls {@link SymplecticBase.base2} with delta array.
         *
         * @see SymplecticBase.step
         */
        public override void step () {
            for (var i = 0; i < n; i++) {
                base2(delta[i]);
            }
        }
    }

    /**
     * Fourth-order symplectic integrator concrete subclass using Suzuki composition.  This integrator is time symmetrical.
     */
    public class FourthOrderSuzuki : Suzuki {

        /**
         * {@inheritDoc}
         * @see Suzuki.Suzuki
         */
        public FourthOrderSuzuki (IModel model, double h) {
            base(model, h);
        }

        /**
         * 4th order integration step.  Calls {@link Suzuki.base4} with s = 1.
         *
         * @see SymplecticBase.step
         */
        public override void step () {
            base4(1.0);
        }
    }

    /**
     * Fourth-order symplectic integrator concrete subclass using Yoshida composition.  This integrator is time symmetrical.
     */
    public class FourthOrderYoshida : Yoshida {

        /**
         * {@inheritDoc}
         * @see Yoshida.Yoshida
         */
        public FourthOrderYoshida (IModel model, double h) {
            base(model, h);
        }

        /**
         * 4th order integration step.  Calls {@link Yoshida.base4} with s = 1.
         *
         * @see SymplecticBase.step
         */
        public override void step () {
            base4(1.0);
        }
    }

    /**
     * Sixth-order symplectic integrator concrete subclass using Suzuki composition.  This integrator is time symmetrical.
     */
    public class SixthOrderSuzuki : Suzuki {

        /**
         * {@inheritDoc}
         * @see Suzuki.Suzuki
         */
        public SixthOrderSuzuki (IModel model, double h) {
            base(model, h);
        }

        /**
         * 6th order integration step.  Calls {@link Suzuki.base6} with s = 1.
         *
         * @see SymplecticBase.step
         */
        public override void step () {
            base6(1.0);
        }
    }

    /**
     * Sixth-order symplectic integrator concrete subclass using Yoshida composition.  This integrator is time symmetrical.
     */
    public class SixthOrderYoshida : Yoshida {

        /**
         * {@inheritDoc}
         * @see Yoshida.Yoshida
         */
        public SixthOrderYoshida (IModel model, double h) {
            base(model, h);
        }

        /**
         * 6th order integration step.  Calls {@link Yoshida.base6} with s = 1.
         *
         * @see SymplecticBase.step
         */
        public override void step () {
            base6(1.0);
        }
    }

    /**
     * Sixth-order symplectic integrator concrete subclass using Kahan-Li composition s9odr6b.  This integrator is time symmetrical.
     */
    public class SixthOrderKahanLi9o6b : GenericComposer {

        /**
         * {@inheritDoc}
         * @see GenericComposer.GenericComposer
         */
        public SixthOrderKahanLi9o6b (IModel model, double h) {
            base(model, h, {
                            0.39103020330868478817,
                            0.33403728961113601749,
                            -0.70622728118756134346,
                            0.081877549648059445768,
                            0.79856447723936218406
                           });
        }
    }

    /**
     * Eightth-order symplectic integrator concrete subclass using Suzuki composition.  This integrator is time symmetrical.
     */
    public class EightthOrderSuzuki : Suzuki {

        /**
         * {@inheritDoc}
         * @see Suzuki.Suzuki
         */
        public EightthOrderSuzuki (IModel model, double h) {
            base(model, h);
        }

        /**
         * 8th order integration step.  Calls {@link Suzuki.base8} with s = 1.
         *
         * @see SymplecticBase.step
         */
        public override void step () {
            base8(1.0);
        }
    }

    /**
     * Eightth-order symplectic integrator concrete subclass using Yoshida composition.  This integrator is time symmetrical.
     */
    public class EightthOrderYoshida : Yoshida {

        /**
         * {@inheritDoc}
         * @see Yoshida.Yoshida
         */
        public EightthOrderYoshida (IModel model, double h) {
            base(model, h);
        }

        /**
         * 8th order integration step.  Calls {@link Yoshida.base8} with s = 1.
         *
         * @see SymplecticBase.step
         */
        public override void step () {
            base8(1.0);
        }
    }

    /**
     * Eightth-order symplectic integrator concrete subclass using Kahan-Li composition s17odr8b.  This integrator is time symmetrical.
     */
    public class EightthOrderKahanLi17o8b : GenericComposer {

        /**
         * {@inheritDoc}
         * @see GenericComposer.GenericComposer
         */
        public EightthOrderKahanLi17o8b (IModel model, double h) {
            base(model, h, {
                            0.12713692773487857916,
                            0.56170253798880269972,
                            -0.38253471994883018888,
                            0.16007605629464743119,
                            -0.40181637432680696673,
                            0.18736671654227849724,
                            0.26070870920779240570,
                            0.29039738812516162389,
                            -0.60607448323584816258
                           });
        }
    }

    /**
     * Static factory for producing subclass instances from its type argument according to the following table:
     *
     * || ''type'' || ''Subclass'' ||  ''Description'' ||
     * || "b1" || {@link FirstOrder} || 1st Order, Symplectic ||
     * || "b2" || {@link SecondOrder} || 2nd Order, Symplectic, Reversible ||
     * || "sb4" || {@link FourthOrderSuzuki} || 4th Order, Symplectic, Reversible ||
     * || "sb6" || {@link SixthOrderSuzuki} || 6th Order, Symplectic, Reversible ||
     * || "sb8" || {@link EightthOrderSuzuki} || 8th Order, Symplectic, Reversible ||
     * || "yb4" || {@link FourthOrderYoshida} || 4th Order, Symplectic, Reversible ||
     * || "yb6" || {@link SixthOrderYoshida} || 6th Order, Symplectic, Reversible ||
     * || "yb8" || {@link EightthOrderYoshida} || 8th Order, Symplectic, Reversible ||
     * || "kl6" || {@link SixthOrderKahanLi9o6b} || 6th Order, Symplectic, Reversible ||
     * || "kl8" || {@link EightthOrderKahanLi17o8b} || 8th Order, Symplectic, Reversible ||
     *
     * @param h the time step
     * @param type the selected implementation
     *
     * @return concrete implementation of a symplectic integrator
     */
    public static ISymplectic getIntegrator (IModel model, double h, string type) {
        switch (type) {
            case "b1":
                stderr.printf("1st Order Symplectic Integrator\n");
                return new FirstOrder(model, h);
            case "b2":
                stderr.printf("2nd Order Symplectic Integrator\n");
                return new SecondOrder(model, h);
            case "sb4":
                stderr.printf("4th Order Symplectic Integrator (using Suzuki Composition)\n");
                return new FourthOrderSuzuki(model, h);
            case "sb6":
                stderr.printf("6th Order Symplectic Integrator (using Suzuki Composition)\n");
                return new SixthOrderSuzuki(model, h);
            case "sb8":
                stderr.printf("8th Order Symplectic Integrator (using Suzuki Composition)\n");
                return new EightthOrderSuzuki(model, h);
            case "yb4":
                stderr.printf("4th Order Symplectic Integrator (using Yoshida Composition)\n");
                return new FourthOrderYoshida(model, h);
            case "yb6":
                stderr.printf("6th Order Symplectic Integrator (using Yoshida Composition)\n");
                return new SixthOrderYoshida(model, h);
            case "yb8":
                stderr.printf("8th Order Symplectic Integrator (using Yoshida Composition)\n");
                return new EightthOrderYoshida(model, h);
            case "kl6":
                stderr.printf("6th Order Symplectic Integrator (using Kahan-Li s9odr6b Composition)\n");
                return new SixthOrderKahanLi9o6b(model, h);
            case "kl8":
                stderr.printf("8th Order Symplectic Integrator (using Kahan-Li s17odr8b Composition)\n");
                return new EightthOrderKahanLi17o8b(model, h);
        }
        stderr.printf("Integrator not recognized: %s\n", type);
        assert_not_reached();
    }
}
