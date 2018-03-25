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
         * Base method coefficients
         */
        private double[] cd_s4;
        private double[] cd_s6;
        private double[] cd_s8;

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
            var x1 = 1.0 / (4.0 - pow(4.0, (1.0 / 7.0)));
            var x0 = 1.0 - 4.0 * x1;
            var y1 = 1.0 / (4.0 - pow(4.0, (1.0 / 5.0)));
            var y0 = 1.0 - 4.0 * y1;
            var z1 = 1.0 / (4.0 - pow(4.0, (1.0 / 3.0)));
            var z0 = 1.0 - 4.0 * z1;
            cd_s4 = {
                      0.5 * h * z1,
                      h * z1, h * z1, h * z1,
                      0.5 * h * (z1 + z0), h * z0
                    };
            cd_s6 = {
                      0.5 * h * z1 * y1,
                      h * z1 * y1, h * z1 * y1, h * z1 * y1,
                      0.5 * h * (z1 + z0) * y1, h * z0 * y1, 0.5 * h * (z0 + z1) * y1,
                      h * z1 * y1, h * z1 * y1, h * z1 * y1,
                      h * z1 * y1,
                      h * z1 * y1, h * z1 * y1, h * z1 * y1,
                      0.5 * h * (z1 + z0) * y1, h * z0 * y1, 0.5 * h * (z0 + z1) * y1,
                      h * z1 * y1, h * z1 * y1, h * z1 * y1,
                      0.5 * h * z1 * (y1 + y0),
                      h * z1 * y0, h * z1 * y0, h * z1 * y0,
                      0.5 * h * (z1 + z0) * y0, h * z0 * y0
                    };
            cd_s8 = {
                      0.5 * h * z1 * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      h * z1 * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * z1 * (y1 + y0) * x1,
                      h * z1 * y0 * x1, h * z1 * y0 * x1, h * z1 * y0 * x1,
                      0.5 * h * (z1 + z0) * y0 * x1, h * z0 * y0 * x1, 0.5 * h * (z0 + z1) * y0 * x1,
                      h * z1 * y0 * x1, h * z1 * y0 * x1, h * z1 * y0 * x1,
                      0.5 * h * z1 * (y1 + y0) * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      h * z1 * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      h * z1 * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      h * z1 * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * z1 * (y1 + y0) * x1,
                      h * z1 * y0 * x1, h * z1 * y0 * x1, h * z1 * y0 * x1,
                      0.5 * h * (z1 + z0) * y0 * x1, h * z0 * y0 * x1, 0.5 * h * (z0 + z1) * y0 * x1,
                      h * z1 * y0 * x1, h * z1 * y0 * x1, h * z1 * y0 * x1,
                      0.5 * h * z1 * (y1 + y0) * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      h * z1 * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * (z1 + z0) * y1 * x1, h * z0 * y1 * x1, 0.5 * h * (z0 + z1) * y1 * x1,
                      h * z1 * y1 * x1, h * z1 * y1 * x1, h * z1 * y1 * x1,
                      0.5 * h * z1 * y1 * (x1 + x0),
                      h * z1 * y1 * x0, h * z1 * y1 * x0, h * z1 * y1 * x0,
                      0.5 * h * (z1 + z0) * y1 * x0, h * z0 * y1 * x0, 0.5 * h * (z0 + z1) * y1 * x0,
                      h * z1 * y1 * x0, h * z1 * y1 * x0, h * z1 * y1 * x0,
                      h * z1 * y1 * x0,
                      h * z1 * y1 * x0, h * z1 * y1 * x0, h * z1 * y1 * x0,
                      0.5 * h * (z1 + z0) * y1 * x0, h * z0 * y1 * x0, 0.5 * h * (z0 + z1) * y1 * x0,
                      h * z1 * y1 * x0, h * z1 * y1 * x0, h * z1 * y1 * x0,
                      0.5 * h * z1 * (y1 + y0) * x0,
                      h * z1 * y0 * x0, h * z1 * y0 * x0, h * z1 * y0 * x0,
                      0.5 * h * (z1 + z0) * y0 * x0, h * z0 * y0 * x0
                    };
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
         * Direct fourth order integrator (my own - equivalent to a Suzuki composition of Stormer-Verlet).
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
         * Direct sixth order base method (my own - equivalent to two Suzuki compositions of Stormer-Verlet).
         */
        private void sixthOrder () {
            model.qUpdate(cd_s6[0]);
            model.pUpdate(cd_s6[1]);
            model.qUpdate(cd_s6[2]);
            model.pUpdate(cd_s6[3]);
            model.qUpdate(cd_s6[4]);
            model.pUpdate(cd_s6[5]);
            model.qUpdate(cd_s6[6]);
            model.pUpdate(cd_s6[7]);
            model.qUpdate(cd_s6[8]);
            model.pUpdate(cd_s6[9]);
            model.qUpdate(cd_s6[10]);
            model.pUpdate(cd_s6[11]);
            model.qUpdate(cd_s6[12]);
            model.pUpdate(cd_s6[13]);
            model.qUpdate(cd_s6[14]);
            model.pUpdate(cd_s6[15]);
            model.qUpdate(cd_s6[16]);
            model.pUpdate(cd_s6[17]);
            model.qUpdate(cd_s6[18]);
            model.pUpdate(cd_s6[19]);
            model.qUpdate(cd_s6[20]);
            model.pUpdate(cd_s6[21]);
            model.qUpdate(cd_s6[22]);
            model.pUpdate(cd_s6[23]);
            model.qUpdate(cd_s6[24]);
            model.pUpdate(cd_s6[25]);
            model.qUpdate(cd_s6[24]);
            model.pUpdate(cd_s6[23]);
            model.qUpdate(cd_s6[22]);
            model.pUpdate(cd_s6[21]);
            model.qUpdate(cd_s6[20]);
            model.pUpdate(cd_s6[19]);
            model.qUpdate(cd_s6[18]);
            model.pUpdate(cd_s6[17]);
            model.qUpdate(cd_s6[16]);
            model.pUpdate(cd_s6[15]);
            model.qUpdate(cd_s6[14]);
            model.pUpdate(cd_s6[13]);
            model.qUpdate(cd_s6[12]);
            model.pUpdate(cd_s6[11]);
            model.qUpdate(cd_s6[10]);
            model.pUpdate(cd_s6[9]);
            model.qUpdate(cd_s6[8]);
            model.pUpdate(cd_s6[7]);
            model.qUpdate(cd_s6[6]);
            model.pUpdate(cd_s6[5]);
            model.qUpdate(cd_s6[4]);
            model.pUpdate(cd_s6[3]);
            model.qUpdate(cd_s6[2]);
            model.pUpdate(cd_s6[1]);
            model.qUpdate(cd_s6[0]);
        }

        /**
         * Direct eigtth order base method (my own - equivalent to three Suzuki compositions of Stormer-Verlet).
         */
        private void eightthOrder () {
            model.qUpdate(cd_s8[0]);
            model.pUpdate(cd_s8[1]);
            model.qUpdate(cd_s8[2]);
            model.pUpdate(cd_s8[3]);
            model.qUpdate(cd_s8[4]);
            model.pUpdate(cd_s8[5]);
            model.qUpdate(cd_s8[6]);
            model.pUpdate(cd_s8[7]);
            model.qUpdate(cd_s8[8]);
            model.pUpdate(cd_s8[9]);
            model.qUpdate(cd_s8[10]);
            model.pUpdate(cd_s8[11]);
            model.qUpdate(cd_s8[12]);
            model.pUpdate(cd_s8[13]);
            model.qUpdate(cd_s8[14]);
            model.pUpdate(cd_s8[15]);
            model.qUpdate(cd_s8[16]);
            model.pUpdate(cd_s8[17]);
            model.qUpdate(cd_s8[18]);
            model.pUpdate(cd_s8[19]);
            model.qUpdate(cd_s8[20]);
            model.pUpdate(cd_s8[21]);
            model.qUpdate(cd_s8[22]);
            model.pUpdate(cd_s8[23]);
            model.qUpdate(cd_s8[24]);
            model.pUpdate(cd_s8[25]);
            model.qUpdate(cd_s8[26]);
            model.pUpdate(cd_s8[27]);
            model.qUpdate(cd_s8[28]);
            model.pUpdate(cd_s8[29]);
            model.qUpdate(cd_s8[30]);
            model.pUpdate(cd_s8[31]);
            model.qUpdate(cd_s8[32]);
            model.pUpdate(cd_s8[33]);
            model.qUpdate(cd_s8[34]);
            model.pUpdate(cd_s8[35]);
            model.qUpdate(cd_s8[36]);
            model.pUpdate(cd_s8[37]);
            model.qUpdate(cd_s8[38]);
            model.pUpdate(cd_s8[39]);
            model.qUpdate(cd_s8[40]);
            model.pUpdate(cd_s8[41]);
            model.qUpdate(cd_s8[42]);
            model.pUpdate(cd_s8[43]);
            model.qUpdate(cd_s8[44]);
            model.pUpdate(cd_s8[45]);
            model.qUpdate(cd_s8[46]);
            model.pUpdate(cd_s8[47]);
            model.qUpdate(cd_s8[48]);
            model.pUpdate(cd_s8[49]);
            model.qUpdate(cd_s8[50]);
            model.pUpdate(cd_s8[51]);
            model.qUpdate(cd_s8[52]);
            model.pUpdate(cd_s8[53]);
            model.qUpdate(cd_s8[54]);
            model.pUpdate(cd_s8[55]);
            model.qUpdate(cd_s8[56]);
            model.pUpdate(cd_s8[57]);
            model.qUpdate(cd_s8[58]);
            model.pUpdate(cd_s8[59]);
            model.qUpdate(cd_s8[60]);
            model.pUpdate(cd_s8[61]);
            model.qUpdate(cd_s8[62]);
            model.pUpdate(cd_s8[63]);
            model.qUpdate(cd_s8[64]);
            model.pUpdate(cd_s8[65]);
            model.qUpdate(cd_s8[66]);
            model.pUpdate(cd_s8[67]);
            model.qUpdate(cd_s8[68]);
            model.pUpdate(cd_s8[69]);
            model.qUpdate(cd_s8[70]);
            model.pUpdate(cd_s8[71]);
            model.qUpdate(cd_s8[72]);
            model.pUpdate(cd_s8[73]);
            model.qUpdate(cd_s8[74]);
            model.pUpdate(cd_s8[75]);
            model.qUpdate(cd_s8[76]);
            model.pUpdate(cd_s8[77]);
            model.qUpdate(cd_s8[78]);
            model.pUpdate(cd_s8[79]);
            model.qUpdate(cd_s8[80]);
            model.pUpdate(cd_s8[81]);
            model.qUpdate(cd_s8[82]);
            model.pUpdate(cd_s8[83]);
            model.qUpdate(cd_s8[84]);
            model.pUpdate(cd_s8[85]);
            model.qUpdate(cd_s8[86]);
            model.pUpdate(cd_s8[87]);
            model.qUpdate(cd_s8[88]);
            model.pUpdate(cd_s8[89]);
            model.qUpdate(cd_s8[90]);
            model.pUpdate(cd_s8[91]);
            model.qUpdate(cd_s8[92]);
            model.pUpdate(cd_s8[93]);
            model.qUpdate(cd_s8[94]);
            model.pUpdate(cd_s8[95]);
            model.qUpdate(cd_s8[96]);
            model.pUpdate(cd_s8[97]);
            model.qUpdate(cd_s8[98]);
            model.pUpdate(cd_s8[99]);
            model.qUpdate(cd_s8[100]);
            model.pUpdate(cd_s8[101]);
            model.qUpdate(cd_s8[102]);
            model.pUpdate(cd_s8[103]);
            model.qUpdate(cd_s8[104]);
            model.pUpdate(cd_s8[105]);
            model.qUpdate(cd_s8[106]);
            model.pUpdate(cd_s8[107]);
            model.qUpdate(cd_s8[108]);
            model.pUpdate(cd_s8[109]);
            model.qUpdate(cd_s8[110]);
            model.pUpdate(cd_s8[111]);
            model.qUpdate(cd_s8[112]);
            model.pUpdate(cd_s8[113]);
            model.qUpdate(cd_s8[114]);
            model.pUpdate(cd_s8[115]);
            model.qUpdate(cd_s8[116]);
            model.pUpdate(cd_s8[117]);
            model.qUpdate(cd_s8[118]);
            model.pUpdate(cd_s8[119]);
            model.qUpdate(cd_s8[120]);
            model.pUpdate(cd_s8[121]);
            model.qUpdate(cd_s8[122]);
            model.pUpdate(cd_s8[123]);
            model.qUpdate(cd_s8[124]);
            model.pUpdate(cd_s8[125]);
            model.qUpdate(cd_s8[124]);
            model.pUpdate(cd_s8[123]);
            model.qUpdate(cd_s8[122]);
            model.pUpdate(cd_s8[121]);
            model.qUpdate(cd_s8[120]);
            model.pUpdate(cd_s8[119]);
            model.qUpdate(cd_s8[118]);
            model.pUpdate(cd_s8[117]);
            model.qUpdate(cd_s8[116]);
            model.pUpdate(cd_s8[115]);
            model.qUpdate(cd_s8[114]);
            model.pUpdate(cd_s8[113]);
            model.qUpdate(cd_s8[112]);
            model.pUpdate(cd_s8[111]);
            model.qUpdate(cd_s8[110]);
            model.pUpdate(cd_s8[109]);
            model.qUpdate(cd_s8[108]);
            model.pUpdate(cd_s8[107]);
            model.qUpdate(cd_s8[106]);
            model.pUpdate(cd_s8[105]);
            model.qUpdate(cd_s8[104]);
            model.pUpdate(cd_s8[103]);
            model.qUpdate(cd_s8[102]);
            model.pUpdate(cd_s8[101]);
            model.qUpdate(cd_s8[100]);
            model.pUpdate(cd_s8[99]);
            model.qUpdate(cd_s8[98]);
            model.pUpdate(cd_s8[97]);
            model.qUpdate(cd_s8[96]);
            model.pUpdate(cd_s8[95]);
            model.qUpdate(cd_s8[94]);
            model.pUpdate(cd_s8[93]);
            model.qUpdate(cd_s8[92]);
            model.pUpdate(cd_s8[91]);
            model.qUpdate(cd_s8[90]);
            model.pUpdate(cd_s8[89]);
            model.qUpdate(cd_s8[88]);
            model.pUpdate(cd_s8[87]);
            model.qUpdate(cd_s8[86]);
            model.pUpdate(cd_s8[85]);
            model.qUpdate(cd_s8[84]);
            model.pUpdate(cd_s8[83]);
            model.qUpdate(cd_s8[82]);
            model.pUpdate(cd_s8[81]);
            model.qUpdate(cd_s8[80]);
            model.pUpdate(cd_s8[79]);
            model.qUpdate(cd_s8[78]);
            model.pUpdate(cd_s8[77]);
            model.qUpdate(cd_s8[76]);
            model.pUpdate(cd_s8[75]);
            model.qUpdate(cd_s8[74]);
            model.pUpdate(cd_s8[73]);
            model.qUpdate(cd_s8[72]);
            model.pUpdate(cd_s8[71]);
            model.qUpdate(cd_s8[70]);
            model.pUpdate(cd_s8[69]);
            model.qUpdate(cd_s8[68]);
            model.pUpdate(cd_s8[67]);
            model.qUpdate(cd_s8[66]);
            model.pUpdate(cd_s8[65]);
            model.qUpdate(cd_s8[64]);
            model.pUpdate(cd_s8[63]);
            model.qUpdate(cd_s8[62]);
            model.pUpdate(cd_s8[61]);
            model.qUpdate(cd_s8[60]);
            model.pUpdate(cd_s8[59]);
            model.qUpdate(cd_s8[58]);
            model.pUpdate(cd_s8[57]);
            model.qUpdate(cd_s8[56]);
            model.pUpdate(cd_s8[55]);
            model.qUpdate(cd_s8[54]);
            model.pUpdate(cd_s8[53]);
            model.qUpdate(cd_s8[52]);
            model.pUpdate(cd_s8[51]);
            model.qUpdate(cd_s8[50]);
            model.pUpdate(cd_s8[49]);
            model.qUpdate(cd_s8[48]);
            model.pUpdate(cd_s8[47]);
            model.qUpdate(cd_s8[46]);
            model.pUpdate(cd_s8[45]);
            model.qUpdate(cd_s8[44]);
            model.pUpdate(cd_s8[43]);
            model.qUpdate(cd_s8[42]);
            model.pUpdate(cd_s8[41]);
            model.qUpdate(cd_s8[40]);
            model.pUpdate(cd_s8[39]);
            model.qUpdate(cd_s8[38]);
            model.pUpdate(cd_s8[37]);
            model.qUpdate(cd_s8[36]);
            model.pUpdate(cd_s8[35]);
            model.qUpdate(cd_s8[34]);
            model.pUpdate(cd_s8[33]);
            model.qUpdate(cd_s8[32]);
            model.pUpdate(cd_s8[31]);
            model.qUpdate(cd_s8[30]);
            model.pUpdate(cd_s8[29]);
            model.qUpdate(cd_s8[28]);
            model.pUpdate(cd_s8[27]);
            model.qUpdate(cd_s8[26]);
            model.pUpdate(cd_s8[25]);
            model.qUpdate(cd_s8[24]);
            model.pUpdate(cd_s8[23]);
            model.qUpdate(cd_s8[22]);
            model.pUpdate(cd_s8[21]);
            model.qUpdate(cd_s8[20]);
            model.pUpdate(cd_s8[19]);
            model.qUpdate(cd_s8[18]);
            model.pUpdate(cd_s8[17]);
            model.qUpdate(cd_s8[16]);
            model.pUpdate(cd_s8[15]);
            model.qUpdate(cd_s8[14]);
            model.pUpdate(cd_s8[13]);
            model.qUpdate(cd_s8[12]);
            model.pUpdate(cd_s8[11]);
            model.qUpdate(cd_s8[10]);
            model.pUpdate(cd_s8[9]);
            model.qUpdate(cd_s8[8]);
            model.pUpdate(cd_s8[7]);
            model.qUpdate(cd_s8[6]);
            model.pUpdate(cd_s8[5]);
            model.qUpdate(cd_s8[4]);
            model.pUpdate(cd_s8[3]);
            model.qUpdate(cd_s8[2]);
            model.pUpdate(cd_s8[1]);
            model.qUpdate(cd_s8[0]);
        }
    }
}
