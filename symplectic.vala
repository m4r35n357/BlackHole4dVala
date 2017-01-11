/*
Copyright (c) 2014, 2015, 2016, Ian Smith (m4r35n357)
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
     * Interface for the physical model
     */
    public interface ISolver : GLib.Object {
        /**
         * Sole method called by main(), calls iterate() on RK4, and ISymplectic.compose() on the Symplectics
         */
        public abstract int[] solve ();
    }

    /**
     * Interface for the symplectic integrators (client)
     */
    public interface IModel : GLib.Object {
        /**
         * Coordinate updates, called by ISymplectic.compose()
         */
        public abstract void qUp (double d);

        /**
         * Momentum updates, called by ISymplectic.compose()
         */
        public abstract void pUp (double c);
    }

    /**
     * Interface for the integrators themselves
     */
    public interface ISymplectic : GLib.Object {
        /**
         * Should be called by IModel.solve() as needed, calls IModel.pUp() and IModel.qUp()
         */
        public abstract void compose ();
    }

    /**
     * Integrator superclass, controls composition but leaves integration details to subclass
     */
    protected abstract class Integrator : ISymplectic, GLib.Object {

        private double[] compositionWeights;
        private int wRange;
        protected IModel model;
        protected double[] baseMethodCoefficients;

        /**
         * Protected constructor, use the static factory
         */
        protected Integrator (IModel model, double[] compositionWeights) {
            this.compositionWeights = compositionWeights;
            this.wRange = compositionWeights.length;
            this.model = model;
        }

        /**
         * Static factory
         */
        public static ISymplectic getIntegrator (IModel model, string type) {
            switch (type) {
                case "sb2":  // second order, basic
                    return new Base2(model, { 1.0 });
                case "sb4":  // fourth order, basic
                    return new Base4(model, { 1.0 });
                case "sc4":  // fourth order, composed
                    var CUBE_ROOT_2 = pow(2.0, (1.0 / 3.0));
                    return new Base2(model, { 1.0 / (2.0 - CUBE_ROOT_2),
                                            - CUBE_ROOT_2 / (2.0 - CUBE_ROOT_2),
                                            1.0 / (2.0 - CUBE_ROOT_2) });
                case "sc6":  // sixth order, composed
                    var FIFTH_ROOT_2 = pow(2.0, (1.0 / 5.0));
                    return new Base4(model, { 1.0 / (2.0 - FIFTH_ROOT_2),
                                            - FIFTH_ROOT_2 / (2.0 - FIFTH_ROOT_2),
                                            1.0 / (2.0 - FIFTH_ROOT_2) });
            }
            assert_not_reached();
        }

        /**
         * Subclasses should perform one integration step with the current compositional weight
         */
        protected abstract void integrate (double compositionWeight);

        /**
         * Perform a composition of weighted integration steps
         */
        public void compose () {
            for (int i = 0; i < wRange; i++) {
                integrate(compositionWeights[i]);
            }
        }
    }

    /**
     * Second-order symplectic integrator subclass
     */
    public class Base2 : Integrator {

        /**
         * Protected constructor, use the static factory in superclass
         */
        protected Base2 (IModel model, double[] compositionWeights) {
            base(model, compositionWeights);
            this.baseMethodCoefficients = { 0.5, 1.0 };
        }

        /**
         * Weighted 2nd order integration step
         */
        protected override void integrate (double compositionWeight) {
            model.pUp(baseMethodCoefficients[0] * compositionWeight);
            model.qUp(baseMethodCoefficients[1] * compositionWeight);
            model.pUp(baseMethodCoefficients[0] * compositionWeight);
        }
    }

    /**
     * Fourth-order symplectic integrator subclass
     */
    public class Base4 : Integrator {

        /**
         * Protected constructor, use the static factory in superclass
         */
        protected Base4 (IModel model, double[] compositionWeights) {
            base(model, compositionWeights);
            var CUBE_ROOT_2 = pow(2.0, (1.0 / 3.0));
            this.baseMethodCoefficients = { 0.5 / (2.0 - CUBE_ROOT_2),
                                            1.0 / (2.0 - CUBE_ROOT_2),
                                            0.5 * (1.0 - CUBE_ROOT_2) / (2.0 - CUBE_ROOT_2),
                                            - CUBE_ROOT_2 / (2.0 - CUBE_ROOT_2) };
        }

        /**
         * Weighted 4th order integration step
         */
        protected override void integrate (double compositionWeight) {
            model.qUp(baseMethodCoefficients[0] * compositionWeight);
            model.pUp(baseMethodCoefficients[1] * compositionWeight);
            model.qUp(baseMethodCoefficients[2] * compositionWeight);
            model.pUp(baseMethodCoefficients[3] * compositionWeight);
            model.qUp(baseMethodCoefficients[2] * compositionWeight);
            model.pUp(baseMethodCoefficients[1] * compositionWeight);
            model.qUp(baseMethodCoefficients[0] * compositionWeight);
        }
    }
}
