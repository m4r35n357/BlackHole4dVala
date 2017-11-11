/*
Copyright (c) 2014, 2015, 2016, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
void add_test_solve_symplectic_start_non_0_kl6() {
    Test.add_func ("/KdS/test_solve_symplectic_start_non_0_kl6", () => {
        var start = 50.0;
        var end = 100;
        var step = 0.00005;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "kl6", 0), step, start, end, interval);
        assert(counts[0] == 85532);
        assert(counts[1] == 3215);
    });
}

void add_test_solve_symplectic_start_non_0_kl8() {
    Test.add_func ("/KdS/test_solve_symplectic_start_non_0_kl8", () => {
        var start = 50.0;
        var end = 100;
        var step = 0.00005;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "kl8", 0), step, start, end, interval);
        assert(counts[0] == 85532);
        assert(counts[1] == 3215);
    });
}

void add_test_solve_symp_polar_b1() {
    Test.add_func ("/KdS/test_solve_symp_polar_b1", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.00001;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 1.0, 1.0, 0.96210432940242041, 5.6843449527674236e-13, 15.914691393798241, 12.0, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b1", 0), step, start, end, interval);
        assert(counts[0] == 69174);
        assert(counts[1] == 6918);
    });
}

void add_test_solve_symp_polar_b2() {
    Test.add_func ("/KdS/test_solve_symp_polar_b2", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.00001;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 1.0, 1.0, 0.96210432940242041, 5.6843449527674236e-13, 15.914691393798241, 12.0, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b2", 0), step, start, end, interval);
        assert(counts[0] == 69174);
        assert(counts[1] == 6918);
    });
}

void add_test_solve_symp_polar_b4() {
    Test.add_func ("/KdS/test_solve_symp_polar_b4", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.00001;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 1.0, 1.0, 0.96210432940242041, 5.6843449527674236e-13, 15.914691393798241, 12.0, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b4", 5), step, start, end, interval);
        assert(counts[0] == 69174);
        assert(counts[1] == 6918);
    });
}

void add_test_solve_symp_polar_b6() {
    Test.add_func ("/KdS/test_solve_symp_polar_b6", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.00001;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 1.0, 1.0, 0.96210432940242041, 5.6843449527674236e-13, 15.914691393798241, 12.0, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b6", 5), step, start, end, interval);
        assert(counts[0] == 69174);
        assert(counts[1] == 6918);
    });
}

void add_test_solve_symp_polar_b8() {
    Test.add_func ("/KdS/test_solve_symp_polar_b8", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.00001;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 1.0, 1.0, 0.96210432940242041, 5.6843449527674236e-13, 15.914691393798241, 12.0, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b8", 5), step, start, end, interval);
        assert(counts[0] == 69174);
        assert(counts[1] == 6918);
    });
}

void add_test_solve_symp_light_b1() {
    Test.add_func ("/KdS/test_solve_symp_light_b1", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.0001;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 1.0, 0.0, 1.0, -2.0, 27.0, 3.0, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b1", 0), step, start, end, interval);
        assert(counts[0] == 105951);
        assert(counts[1] == 10596);
    });
}

void add_test_solve_symp_light_b2() {
    Test.add_func ("/KdS/test_solve_symp_light_b2", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.0001;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 1.0, 0.0, 1.0, -2.0, 27.0, 3.0, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b2", 0), step, start, end, interval);
        assert(counts[0] == 105951);
        assert(counts[1] == 10596);
    });
}

void add_test_solve_symp_light_b4() {
    Test.add_func ("/KdS/test_solve_symp_light_b4", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.0001;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 1.0, 0.0, 1.0, -2.0, 27.0, 3.0, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b4", 5), step, start, end, interval);
        assert(counts[0] == 105951);
        assert(counts[1] == 10596);
    });
}

void add_test_solve_symp_light_b6() {
    Test.add_func ("/KdS/test_solve_symp_light_b6", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.0001;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 1.0, 0.0, 1.0, -2.0, 27.0, 3.0, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b6", 5), step, start, end, interval);
        assert(counts[0] == 105951);
        assert(counts[1] == 10596);
    });
}

void add_test_solve_symp_light_b8() {
    Test.add_func ("/KdS/test_solve_symp_light_b8", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.0001;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 1.0, 0.0, 1.0, -2.0, 27.0, 3.0, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b8", 5), step, start, end, interval);
        assert(counts[0] == 105951);
        assert(counts[1] == 10596);
    });
}

void add_test_solve_symplectic_start_0_b1() {
    Test.add_func ("/KdS/test_solve_symplectic_start_0_b1", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.00005;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b1", 0), step, start, end, interval);
        assert(counts[0] == 85862);
        assert(counts[1] == 8587);
    });
}

void add_test_solve_symplectic_start_0_b2() {
    Test.add_func ("/KdS/test_solve_symplectic_start_0_b2", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.00005;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b2", 0), step, start, end, interval);
        assert(counts[0] == 85532);
        assert(counts[1] == 8554);
    });
}

void add_test_solve_symplectic_start_0_b4() {
    Test.add_func ("/KdS/test_solve_symplectic_start_0_b4", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.00005;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b4", 5), step, start, end, interval);
        assert(counts[0] == 85532);
        assert(counts[1] == 8554);
    });
}

void add_test_solve_symplectic_start_0_b6() {
    Test.add_func ("/KdS/test_solve_symplectic_start_0_b6", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.00005;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b6", 5), step, start, end, interval);
        assert(counts[0] == 85532);
        assert(counts[1] == 8554);
    });
}

void add_test_solve_symplectic_start_0_b8() {
    Test.add_func ("/KdS/test_solve_symplectic_start_0_b8", () => {
        var start = 0.0;
        var end = 100;
        var step = 0.00005;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b8", 5), step, start, end, interval);
        assert(counts[0] == 85532);
        assert(counts[1] == 8554);
    });
}

void add_test_solve_symplectic_start_non_0_b1() {
    Test.add_func ("/KdS/test_solve_symplectic_start_non_0_b1", () => {
        var start = 50.0;
        var end = 100;
        var step = 0.00005;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b1", 0), step, start, end, interval);
        assert(counts[0] == 85862);
        assert(counts[1] == 3242);
    });
}

void add_test_solve_symplectic_start_non_0_b2() {
    Test.add_func ("/KdS/test_solve_symplectic_start_non_0_b2", () => {
        var start = 50.0;
        var end = 100;
        var step = 0.00005;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b2", 0), step, start, end, interval);
        assert(counts[0] == 85532);
        assert(counts[1] == 3215);
    });
}

void add_test_solve_symplectic_start_non_0_b4() {
    Test.add_func ("/KdS/test_solve_symplectic_start_non_0_b4", () => {
        var start = 50.0;
        var end = 100;
        var step = 0.00005;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b4", 5), step, start, end, interval);
        assert(counts[0] == 85532);
        assert(counts[1] == 3215);
    });
}

void add_test_solve_symplectic_start_non_0_b6() {
    Test.add_func ("/KdS/test_solve_symplectic_start_non_0_b6", () => {
        var start = 50.0;
        var end = 100;
        var step = 0.00005;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b6", 5), step, start, end, interval);
        assert(counts[0] == 85532);
        assert(counts[1] == 3215);
    });
}

void add_test_solve_symplectic_start_non_0_b8() {
    Test.add_func ("/KdS/test_solve_symplectic_start_non_0_b8", () => {
        var start = 50.0;
        var end = 100;
        var step = 0.00005;
        var interval = 10;
        var model = new Simulations.Bh3d(0.0, 0.8, 1.0, 0.94550509567490792, 1.4343745095317371, 7.9787599589278697, 7.5, 0.0, true);
        var counts = model.solve(Integrators.getIntegrator(model, step, "b8", 5), step, start, end, interval);
        assert(counts[0] == 85532);
        assert(counts[1] == 3215);
    });
}

void main (string[] args) {
    Test.init(ref args);

    add_test_solve_symplectic_start_non_0_kl6();
    add_test_solve_symplectic_start_non_0_kl8();

    add_test_solve_symp_polar_b1();
    add_test_solve_symp_polar_b2();
    add_test_solve_symp_polar_b4();
    add_test_solve_symp_polar_b6();
    add_test_solve_symp_polar_b8();

    add_test_solve_symp_light_b1();
    add_test_solve_symp_light_b2();
    add_test_solve_symp_light_b4();
    add_test_solve_symp_light_b6();
    add_test_solve_symp_light_b8();

    add_test_solve_symplectic_start_0_b1();
    add_test_solve_symplectic_start_0_b2();
    add_test_solve_symplectic_start_0_b4();
    add_test_solve_symplectic_start_0_b6();
    add_test_solve_symplectic_start_0_b8();

    add_test_solve_symplectic_start_non_0_b1();
    add_test_solve_symplectic_start_non_0_b2();
    add_test_solve_symplectic_start_non_0_b4();
    add_test_solve_symplectic_start_non_0_b6();
    add_test_solve_symplectic_start_non_0_b8();

    Test.run();
}
