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
         * Sole method called by main(), calls {@link Integrators.ISymplectic.step} on the selected integrator once per iteration
         *
         * @param integrator the selected implementation
         * @param h the time step
         * @param start start time
         * @param end finish time
         * @param tr plot ratio i.e. only plot every tr-th point
         *
         * @return an array of iteration counters
         */
        public abstract int64[] solve (Integrators.ISymplectic integrator, double h, double start, double end, int64 tr);
    }

    /**
     * Internal interface for the symplectic integrator to access and update the model
     */
    public interface IModel : GLib.Object {
        /**
         * Hamiltonian equations of motion - coordinate updates (dT/dp), called by {@link Integrators.ISymplectic.step}
         *
         * @param d scaled time step
         */
        public abstract void qUpdate (double d);

        /**
         * Hamiltonian equations of motion - momentum updates (dV/dq), called by {@link Integrators.ISymplectic.step}
         *
         * @param c scaled time step
         */
        public abstract void pUpdate (double c);
    }
}

namespace Integrators {

    /**
     * Internal interface for a model to drive the symplectic integrator
     */
    public interface ISymplectic : GLib.Object {
        /**
         * Should be called by {@link Models.ISolver.solve} as needed, in turn calls {@link Models.IModel.qUpdate} and {@link Models.IModel.pUpdate}
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
        protected Models.IModel model;

        /**
         * Simulation time step, defined at subclass construction
         */
        protected double h;

        /**
         * This is a protected constructor; use the static factory {@link getIntegrator} to obtain a concrete subclass.
         *
         * @param h the time step
         */
        protected SymplecticBase (Models.IModel model, double h) {
            this.model = model;
            this.h = h;
        }

        /**
         * Subclasses should perform one integration step by executing alternating {@link Models.IModel.qUpdate} and {@link Models.IModel.pUpdate}
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
        public FirstOrder (Models.IModel model, double h) {
            base(model, h);
        }

        /**
         * 1st order integration step.  Performs the following calls on {@link Models.IModel} per iteration:
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
     * Second-order symplectic integrator concrete subclass.
     */
    public class SecondOrder : SymplecticBase {

        /**
         * {@inheritDoc}
         * @see SymplecticBase.SymplecticBase
         */
        public SecondOrder (Models.IModel model, double h) {
            base(model, h);
        }

        /**
         * Stormer-Verlet base integrator.  Performs the following calls on {@link Models.IModel} per iteration:
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
         * 2nd order integration step.  Calls {@link base2} with s = 1.
         *
         * @see SymplecticBase.step
         */
        public override void step () {
            base2(1.0);
        }
    }

    /**
     * Abstract base for compositions that are defined by explicit formulas (Yoshida, Suzuki and others)).
     */
    protected abstract class ExplicitComposed : SecondOrder {

        /**
         * Root of exponential formulas
         */
        protected int64 root;

        /**
         * For looping over the outer stages
         */
        protected int64 nOuter;

        /**
         * {@inheritDoc}
         * @see SecondOrder.SecondOrder
         */
        public ExplicitComposed (Models.IModel model, double h, int64 stages) {
            base(model, h);
            root = stages - 1;
            nOuter = root / 2;
        }
    }

    /**
     * Fourth-order symplectic integrator concrete subclass using various compositions.
     */
    public class FourthOrder : ExplicitComposed {

        /**
         * Composition coefficients
         */
        protected double zOuter;
        protected double zCentral;

        /**
         * {@inheritDoc}
         * @see SymplecticBase.SymplecticBase
         */
        public FourthOrder (Models.IModel model, double h, int64 stages) {
            base(model, h, stages);
            zOuter = 1.0 / (root - pow(root, (1.0 / 3.0)));
            zCentral = 1.0 - root * zOuter;
        }

        /**
         * General composition from 2nd order to 4th order.  Performs the following calls to {@link SecondOrder.base2} per iteration:
         *
         * {{{
         * base2(s * zOuter)
         * ...
         * base2(s * zCentral)
         * ...
         * base2(s * zOuter)
         * }}}
         *
         * where zOuter = 1 / (4 - 4**(1/3)), and zCentral = 1 - 4 * zOuter
         */
        protected void base4 (double s) {
            for (var _ = 0; _ < nOuter; _++) {
                base2(s * zOuter);
            }
            base2(s * zCentral);
            for (var _ = 0; _ < nOuter; _++) {
                base2(s * zOuter);
            }
        }

        /**
         * 4th order integration step.  Calls {@link base4} with s = 1.
         *
         * @see SymplecticBase.step
         */
        public override void step () {
            base4(1.0);
        }
    }

    /**
     * Sixth-order symplectic integrator concrete subclass using various compositions.
     */
    public class SixthOrder : FourthOrder {

        /**
         * Composition coefficients
         */
        protected double yOuter;
        protected double yCentral;

        /**
         * {@inheritDoc}
         * @see FourthOrder.FourthOrder
         */
        public SixthOrder (Models.IModel model, double h, int64 stages) {
            base(model, h, stages);
            yOuter = 1.0 / (root - pow(root, (1.0 / 5.0)));
            yCentral = 1.0 - root * yOuter;
        }

        /**
         * General composition from 4th order to 6th order.  Performs the following calls per iteration:
         *
         * {{{
         * base4(s * yOuter)
         * ...
         * base4(s * yCentral)
         * ...
         * base4(s * yOuter)
         * }}}
         *
         * where yOuter = 1 / (4 - 4**(1/5)), and yCentral = 1 - 4 * yOuter
         */
        protected void base6 (double s) {
            for (var _ = 0; _ < nOuter; _++) {
                base4(s * yOuter);
            }
            base4(s * yCentral);
            for (var _ = 0; _ < nOuter; _++) {
                base4(s * yOuter);
            }
        }

        /**
         * 6th order integration step.  Calls {@link base6} with s = 1.
         *
         * @see SymplecticBase.step
         */
        public override void step () {
            base6(1.0);
        }
    }

    /**
     * Eightth-order symplectic integrator concrete subclass using various compositions.
     */
    public class EightthOrder : SixthOrder {

        /**
         * Composition coefficients
         */
        protected double xOuter;
        protected double xCentral;

        /**
         * {@inheritDoc}
         * @see SixthOrder.SixthOrder
         */
        public EightthOrder (Models.IModel model, double h, int64 stages) {
            base(model, h, stages);
            xOuter = 1.0 / (root - pow(root, (1.0 / 7.0)));
            xCentral = 1.0 - root * xOuter;
        }

        /**
         * General composition from 6th order to 8th order.  Performs the following calls per iteration:
         *
         * {{{
         * base6(s * xOuter)
         * ...
         * base6(s * xCentral)
         * ...
         * base6(s * xOuter)
         * }}}
         *
         * where xOuter = 1 / (4 - 4**(1/7)), and xCentral = 1 - 4 * xOuter
         */
        protected void base8 (double s) {
            for (var _ = 0; _ < nOuter; _++) {
                base6(s * xOuter);
            }
            base6(s * xCentral);
            for (var _ = 0; _ < nOuter; _++) {
                base6(s * xOuter);
            }
        }

        /**
         * 8th order integration step.  Calls {@link base8} with s = 1.
         *
         * @see SymplecticBase.step
         */
        public override void step () {
            base8(1.0);
        }
    }

    /**
     * Generic composed integrator abstract superclass, leaves integration method selection to subclasses
     */
    protected abstract class GenericComposer : SecondOrder {
        /**
         * Generic composition coefficients
         */
        protected double[] delta;

        /**
         * Size of coefficients array
         */
        protected int n;

        /**
         * Populates a symmetric delta array from the asymmetric coefficients passed in.
         *
         * This is a protected constructor; use the static factory {@link getIntegrator} to obtain a concrete subclass.
         *
         * @param h the time step
         */
        protected GenericComposer (Models.IModel model, double h, double[] coefficients) {
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
         * Generic integration step.  Calls {@link SecondOrder.base2} with delta array.
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
     * Sixth-order symplectic integrator concrete subclass using Kahan-Li composition s9odr6b.
     */
    public class SixthOrderKahanLi9o6b : GenericComposer {

        /**
         * {@inheritDoc}
         * @see GenericComposer.GenericComposer
         */
        public SixthOrderKahanLi9o6b (Models.IModel model, double h) {
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
     * Eightth-order symplectic integrator concrete subclass using Kahan-Li composition s17odr8b.
     */
    public class EightthOrderKahanLi17o8b : GenericComposer {

        /**
         * {@inheritDoc}
         * @see GenericComposer.GenericComposer
         */
        public EightthOrderKahanLi17o8b (Models.IModel model, double h) {
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
     * || "b4" || {@link FourthOrder} || 4th Order, Symplectic, Reversible ||
     * || "b6" || {@link SixthOrder} || 6th Order, Symplectic, Reversible ||
     * || "b8" || {@link EightthOrder} || 8th Order, Symplectic, Reversible ||
     * || "kl6" || {@link SixthOrderKahanLi9o6b} || 6th Order, Symplectic, Reversible ||
     * || "kl8" || {@link EightthOrderKahanLi17o8b} || 8th Order, Symplectic, Reversible ||
     *
     * @param h the time step
     * @param type the selected implementation
     *
     * @return concrete implementation of a symplectic integrator
     */
    public static ISymplectic getIntegrator (Models.IModel model, double h, string type, int64 stages) {
        if ((stages < 3) || (stages % 2 == 0)) {
            stderr.printf("'stages' should be odd and at least 3\n");
            assert_not_reached();
        }
        switch (type) {
            case "b1":
                stderr.printf("1st Order Symplectic Integrator\n");
                return new FirstOrder(model, h);
            case "b2":
                stderr.printf("2nd Order Symplectic Integrator\n");
                return new SecondOrder(model, h);
            case "b4":
                stderr.printf("4th Order Symplectic Integrator (using explicit composition)\n");
                return new FourthOrder(model, h, stages);
            case "b6":
                stderr.printf("6th Order Symplectic Integrator (using explicit composition)\n");
                return new SixthOrder(model, h, stages);
            case "b8":
                stderr.printf("8th Order Symplectic Integrator (using explicit composition)\n");
                return new EightthOrder(model, h, stages);
            case "kl6":
                stderr.printf("6th Order Symplectic Integrator (using Kahan-Li s9odr6b composition)\n");
                return new SixthOrderKahanLi9o6b(model, h);
            case "kl8":
                stderr.printf("8th Order Symplectic Integrator (using Kahan-Li s17odr8b composition)\n");
                return new EightthOrderKahanLi17o8b(model, h);
        }
        stderr.printf("Integrator not recognized: %s\n", type);
        assert_not_reached();
    }
}
