/*
Copyright (c) 2014, 2015, 2016, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
namespace Simulations {

    void add_test_solve_rk4_start_0 () {
        Test.add_func ("/KdS/test_solve_rk4_start_0", () => {
            var start = 0.0;
            var end = 100;
            var step = 0.001;
            var interval = 10;
            var counts = BhRk4.getInstanceFromArgs(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, start, end, step, interval, "rk4").solve();
            assert(counts[0] == 100000);
            assert(counts[1] == 10000);
        });
    }

    void add_test_solve_rk4_start_non_0 () {
        Test.add_func ("/KdS/test_solve_rk4_start_non_0", () => {
            var start = 50.0;
            var end = 100;
            var step = 0.001;
            var interval = 10;
            var counts = BhRk4.getInstanceFromArgs(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, start, end, step, interval, "rk4").solve();
            assert(counts[0] == 100000);
            assert(counts[1] == 5000);
        });
    }

    void add_test_solve_symplectic_start_0 () {
        Test.add_func ("/KdS/test_solve_symplectic_start_0", () => {
            var start = 0.0;
            var end = 100;
            var step = 0.0001;
            var interval = 1;
            var counts = BhSymp.getInstanceFromArgs(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, start, end, step, interval, "sb4").solve();
            assert(counts[0] == 42766);
            assert(counts[1] == 42766);
        });
    }

    void add_test_solve_symplectic_start_non_0 () {
        Test.add_func ("/KdS/test_solve_symplectic_start_non_0", () => {
            var start = 50.0;
            var end = 100;
            var step = 0.0001;
            var interval = 10;
            var counts = BhSymp.getInstanceFromArgs(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, start, end, step, interval, "sb4").solve();
            assert(counts[0] == 42766);
            assert(counts[1] == 16073);
        });
    }

    void main (string[] args) {
        Test.init(ref args);
        //Test.verbose();

        add_test_solve_rk4_start_0();
        add_test_solve_rk4_start_non_0();

        add_test_solve_symplectic_start_0();
        add_test_solve_symplectic_start_non_0();

        Test.run();
    }
}
