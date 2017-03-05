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

        private double h;  // the time step
        protected double[] c_d;  // the update coefficients

        /**
         * Protected constructor, use the static factory
         * Subclass constructors should populate a c_d coefficient array, premultiplied by the step size
         */
        protected Symplectic (double h) {
            this.h = h;
        }

        /**
         * Static factory
         */
        public static ISymplectic getIntegrator (double h, string type) {
            switch (type) {
                case "sb1":
                    stderr.printf("1st Order Symplectic Integrator\n");
                    return new EulerCromer(h);
                case "sb2":
                    stderr.printf("2nd Order Symplectic Integrator\n");
                    return new StormerVerlet(h);
                case "sb4":
                    stderr.printf("4th Order Symplectic Integrator\n");
                    return new ForestRuth(h);
            }
            stderr.printf("Integrator not recognized: %s\n", type);
            assert_not_reached();
        }

        /**
         * Subclasses should perform one integration step by executing alternating qUpdate() and pUpdate()
         * methods with the c_d array elements as parameters
         */
        protected abstract void step (IModel model);
    }

    /**
     * First-order symplectic integrator concrete subclass
     */
    public class EulerCromer : Symplectic {

        /**
         * Protected constructor, use the static factory in superclass
         */
        protected EulerCromer (double h) {
            base(h);
            this.c_d = { h };
        }

        /**
         * 1st order integration step
         */
        protected override void step (IModel model) {
            model.qUpdate(c_d[0]);
            model.pUpdate(c_d[0]);
        }
    }

    /**
     * Second-order symplectic integrator concrete subclass
     */
    public class StormerVerlet : Symplectic {

        /**
         * Protected constructor, use the static factory in superclass
         */
        protected StormerVerlet (double h) {
            base(h);
            this.c_d = { 0.5 * h, h };
        }

        /**
         * 2nd order integration step
         */
        protected override void step (IModel model) {
            model.qUpdate(c_d[0]);
            model.pUpdate(c_d[1]);
            model.qUpdate(c_d[0]);
        }
    }

    /**
     * Fourth-order symplectic integrator concrete subclass
     */
    public class ForestRuth : Symplectic {

        /**
         * Protected constructor, use the static factory in superclass
         */
        protected ForestRuth (double h) {
            base(h);
            var theta = 1.0 / (2.0 - pow(2.0, (1.0 / 3.0)));
            this.c_d = { 0.5 * h * theta, h * theta, 0.5 * h * (1.0 - theta), h * (1.0 - 2.0 * theta) };
        }

        /**
         * 4th order integration step
         */
        protected override void step (IModel model) {
            model.qUpdate(c_d[0]);
            model.pUpdate(c_d[1]);
            model.qUpdate(c_d[2]);
            model.pUpdate(c_d[3]);
            model.qUpdate(c_d[2]);
            model.pUpdate(c_d[1]);
            model.qUpdate(c_d[0]);
        }
    }
}
