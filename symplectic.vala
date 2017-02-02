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
     * Symplectic integrator abstract superclass, controls composition but leaves integration details to subclass
     */
    protected abstract class Symplectic : ISymplectic, GLib.Object {

        private double[] compositionWeights;
        protected IModel model;
        protected double[] coefficients;

        /**
         * Protected constructor, use the static factory
         */
        protected Symplectic (IModel model, double[] compositionWeights) {
            this.compositionWeights = compositionWeights;
            this.model = model;
        }

        /**
         * Static factory
         */
        public static ISymplectic getIntegrator (IModel model, string type) {
            switch (type) {
                case "sb2":
                    stderr.printf("2nd Order Symplectic Integrator\n");
                    return new Base2(model, { 1.0 });
                case "sb4":
                    stderr.printf("4th Order Symplectic Integrator\n");
                    return new Base4(model, { 1.0 });
                case "sc4":
                    stderr.printf("4th Order Symplectic Integrator (Composed from 2nd order)\n");
                    var CUBRT2 = pow(2.0, 1.0/3.0);
                    return new Base2(model, { 1.0/(2.0-CUBRT2), -CUBRT2/(2.0-CUBRT2), 1.0/(2.0-CUBRT2) });
                case "sc6":
                    stderr.printf("6th Order Symplectic Integrator (Composed from 4th order)\n");
                    var FTHRT2 = pow(2.0, 1.0/5.0);
                    return new Base4(model, { 1.0/(2.0-FTHRT2), -FTHRT2/(2.0-FTHRT2), 1.0/(2.0-FTHRT2) });
                case "sh6":
                    stderr.printf("6th Order Symplectic Integrator (Composed from 2nd order)\n");
                    return new Base2(model, {
                                                0.78451361047755726381949763,
                                                0.23557321335935813368479318,
                                                -1.17767998417887100694641568,
                                                1.31518632068391121888424973,
                                                -1.17767998417887100694641568,
                                                0.23557321335935813368479318,
                                                0.78451361047755726381949763
                                            });
                case "sh8":
                    stderr.printf("8th Order Symplectic Integrator (Composed from 2nd order)\n");
                    return new Base2(model, {
                                                0.74167036435061295344822780,
                                                -0.40910082580003159399730010,
                                                0.19075471029623837995387626,
                                                -0.57386247111608226665638773,
                                                0.29906418130365592384446354,
                                                0.33462491824529818378495798,
                                                0.31529309239676659663205666,
                                                -0.79688793935291635401978884,
                                                0.31529309239676659663205666,
                                                0.33462491824529818378495798,
                                                0.29906418130365592384446354,
                                                -0.57386247111608226665638773,
                                                0.19075471029623837995387626,
                                                -0.40910082580003159399730010,
                                                0.74167036435061295344822780
                                            });
                case "sh10":
                    stderr.printf("10th Order Symplectic Integrator (Composed from 2nd order)\n");
                    return new Base2(model, {
                                                0.09040619368607278492161150,
                                                0.53591815953030120213784983,
                                                0.35123257547493978187517736,
                                                -0.31116802097815835426086544,
                                                -0.52556314194263510431065549,
                                                0.14447909410225247647345695,
                                                0.02983588609748235818064083,
                                                0.17786179923739805133592238,
                                                0.09826906939341637652532377,
                                                0.46179986210411860873242126,
                                                -0.33377845599881851314531820,
                                                0.07095684836524793621031152,
                                                0.23666960070126868771909819,
                                                -0.49725977950660985445028388,
                                                -0.30399616617237257346546356,
                                                0.05246957188100069574521612,
                                                0.44373380805019087955111365,
                                                0.05246957188100069574521612,
                                                -0.30399616617237257346546356,
                                                -0.49725977950660985445028388,
                                                0.23666960070126868771909819,
                                                0.07095684836524793621031152,
                                                -0.33377845599881851314531820,
                                                0.46179986210411860873242126,
                                                0.09826906939341637652532377,
                                                0.17786179923739805133592238,
                                                0.02983588609748235818064083,
                                                0.14447909410225247647345695,
                                                -0.52556314194263510431065549,
                                                -0.31116802097815835426086544,
                                                0.35123257547493978187517736,
                                                0.53591815953030120213784983,
                                                0.09040619368607278492161150
                                            });
            }
            stderr.printf("Integrator not recognized: %s\n", type);
            assert_not_reached();
        }

        /**
         * Subclasses should perform one integration step with the current compositional weight
         */
        protected abstract void integrate (double compositionWeight, double h);

        /**
         * Perform a composition of weighted integration steps
         */
        public void compose (double h) {
            for (int i = 0; i < compositionWeights.length; i += 1) {
                integrate(compositionWeights[i], h);
            }
        }
    }

    /**
     * Second-order symplectic integrator (position Verlet) concrete subclass
     */
    public class Base2 : Symplectic {

        /**
         * Protected constructor, use the static factory in superclass
         */
        protected Base2 (IModel model, double[] compositionWeights) {
            base(model, compositionWeights);
            this.coefficients = { 0.5, 1.0 };
        }

        /**
         * Weighted 2nd order integration step
         */
        protected override void integrate (double weight, double h) {
            model.qUp(coefficients[0] * weight, h);
            model.pUp(coefficients[1] * weight, h);
            model.qUp(coefficients[0] * weight, h);
        }
    }

    /**
     * Fourth-order symplectic integrator (position Forest-Ruth) concrete subclass
     */
    public class Base4 : Symplectic {

        /**
         * Protected constructor, use the static factory in superclass
         */
        protected Base4 (IModel model, double[] compositionWeights) {
            base(model, compositionWeights);
            var CUBRT2 = pow(2.0, (1.0 / 3.0));
            this.coefficients = { 0.5/(2.0-CUBRT2), 1.0/(2.0-CUBRT2), 0.5*(1.0-CUBRT2)/(2.0-CUBRT2), -CUBRT2/(2.0-CUBRT2) };
        }

        /**
         * Weighted 4th order integration step
         */
        protected override void integrate (double weight, double h) {
            model.qUp(coefficients[0] * weight, h);
            model.pUp(coefficients[1] * weight, h);
            model.qUp(coefficients[2] * weight, h);
            model.pUp(coefficients[3] * weight, h);
            model.qUp(coefficients[2] * weight, h);
            model.pUp(coefficients[1] * weight, h);
            model.qUp(coefficients[0] * weight, h);
        }
    }
}
