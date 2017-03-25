/*
Copyright (c) 2014, 2015, 2016, 2017 Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
using Json;
using GLib.Math;

/**
 * Various simulations based on symplectic integration of separable Hamiltonian problems
 */
namespace Simulations {

    /**
     * Internal interface for the parameter generators
     */
    public interface IGenerator : GLib.Object {
        /**
         * Sets up and runs the solver
         *
         * @param input the JSON IC generation data
         */
        public abstract void generateInitialConditions (Json.Object input);

        /**
         * Writes the potential data to STDOUT for plotting
         *
         * @param input the generated JSON simulation data
         */
        public abstract void printPotentials (Json.Object input);
    }

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
         * @param start start time
         * @param end finish time
         * @param tr plot ratio i.e. only plot every tr-th point
         *
         * @return an array of iteration counters
         */
        public abstract int64[] solve (ISymplectic integrator, double start, double end, int64 tr);
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
         *
         * @param model class providing equations of motion for a separable hamiltonian problem
         */
        public abstract void step (IModel model);

        /**
         * Getter method for the time step.
         *
         * @return the baked-in time step
         */
        public abstract double getH ();
    }

    /**
     * Parse JSON initial conditions data from stdin
     *
     * @return an object containing the parsed data
     */
    private static Json.Object getJson () {
        var json = new StringBuilder();
        var line = stdin.read_line();
        while (line != null) {
            json.append_printf("%s\n", line);
            stderr.printf("%s\n", line);
            line = stdin.read_line();
        }
        unowned Json.Object obj;
        var p = new Parser();
        try {
            p.load_from_data(json.str);
            obj = p.get_root().get_object();
        } catch (Error e) {
            stderr.printf("Unable to parse the input data: %s\n", e.message);
            assert_not_reached();
        }
        return obj;
    }
}
