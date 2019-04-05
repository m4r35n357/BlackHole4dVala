#!/usr/bin/env python3

import copy
from json import loads
from sys import argv, stderr
from gmpy2 import mpfr, get_context, sin, acos, cos
get_context().precision = 113  # Set this BEFORE importing any Taylor Series stuff!
from Symplectic import D1, D2
from dual import Dual, make_mpfr


def nelder_mead(f, x_0, x_δ, stuck_threshold=10e-12, stuck_break=10, max_iterations=0, α=1.0, γ=2.0, ρ=-0.5, σ=0.5):

    # initial simplex
    nt = nr = ne = nc = ns = 0
    dim = len(x_0)
    prev_best = f(x_0)
    no_improvement = 0
    results = [[x_0, prev_best]]
    for i in range(dim):
        x = copy.copy(x_0)
        x[i] += x_δ[i]
        score = f(x)
        results.append([x, score])

    iterations = 0
    while True:
        nt += 1

        # 1. order
        results.sort(key=lambda z: z[1])
        best = results[0][1]

        # break after max_iter
        if max_iterations and iterations >= max_iterations:
            return results[0], nt - 1, nr, ne, nc, ns
        iterations += 1

        # break after no_improvement iterations with no improvement
        # print('...best so far:', best)
        if best < prev_best - stuck_threshold:
            no_improvement = 0
            prev_best = best
        else:
            no_improvement += 1
        if no_improvement >= stuck_break:
            return results[0], nt - 1, nr, ne, nc, ns

        # 2. centroid
        x0 = [0.0] * dim
        for tup in results[:-1]:
            for i, c in enumerate(tup[0]):
                x0[i] += c / (len(results)-1)

        # 3. reflect
        worst = results[-1][0]
        xr = [0.0] * dim
        for i in range(dim):
            xr[i] = x0[i] + α * (x0[i] - worst[i])
        r_score = f(xr)
        if results[0][1] <= r_score < results[-2][1]:
            nr += 1
            del results[-1]
            results.append([xr, r_score])
            continue

        # 4. expand
        if r_score < results[0][1]:
            ne += 1
            worst = results[-1][0]
            xe = [0.0] * dim
            for i in range(dim):
                xe[i] = x0[i] + γ * (x0[i] - worst[i])
            e_score = f(xe)
            if e_score < r_score:
                del results[-1]
                results.append([xe, e_score])
                continue
            else:
                del results[-1]
                results.append([xr, r_score])
                continue

        # 5. contract
        worst = results[-1][0]
        xc = [0.0] * dim
        for i in range(dim):
            xc[i] = x0[i] + ρ * (x0[i] - worst[i])
        c_score = f(xc)
        if c_score < results[-1][1]:
            nc += 1
            del results[-1]
            results.append([xc, c_score])
            continue

        # 6. shrink
        ns += 1
        x1 = results[0][0]
        reduced = []
        xs = [0.0] * dim
        for tup in results:
            for i in range(dim):
                xs[i] = x1[i] + σ * (tup[0][i] - x1[i])
            score = f(xs)
            reduced.append([xs, score])
        results = reduced


class Potentials(object):
    def __init__(self, a, r_min, r_max, elevation):
        self.a = a
        self.a2 = a**2
        self.μ2 = make_mpfr(1)
        self.r_min = make_mpfr(r_min)
        self.r_max = make_mpfr(r_max)
        self.θ = make_mpfr((D1 - (make_mpfr(90) - elevation) / make_mpfr(180)) * acos(make_mpfr(-1)))

    def f_r(self, r, e, l, q):
        r2 =r.sqr
        ra2 = r2 + self.a2
        return (ra2 * e - self.a * l).sqr - (ra2 - D2 * r) * (self.μ2 * r2 + q + (l - self.a * e)**2)

    def f_θ(self, θ, e, l, q):
        return q - cos(θ)**2 * (self.a2 * (self.μ2 - e**2) + (l / sin(θ))**2)

    def f_spherical(self, x):
        E, L, Q = x
        r_potential = self.f_r(Dual.from_number(self.r_max, variable=True), E, L, Q)
        return r_potential.val**2 + r_potential.der**2 + self.f_θ(self.θ, E, L, Q)**2

    def f_nonspherical(self, x):
        E, L, Q = x
        r_min_potential = self.f_r(Dual.from_number(self.r_min), E, L, Q).val
        r_max_potential = self.f_r(Dual.from_number(self.r_max), E, L, Q).val
        return r_min_potential**2 + r_max_potential**2 + self.f_θ(self.θ, E, L, Q)**2


if __name__ == "__main__":
    # Example: ./NelderMead.py icgen-data.json 1.0 5.0 0.0 1.0 1.0 1.0
    print("Generator: {}".format(argv[0]), file=stderr)
    input_data = open(argv[1]).read()
    ic = loads(input_data, parse_float=mpfr)
    print(input_data, file=stderr)
    if ic.get('rMax'):
        potentials = Potentials(ic['spin'], ic['rMin'], ic['rMax'], ic['elevation'])
        func = potentials.f_nonspherical
    else:
        potentials = Potentials(ic['spin'], ic['r'], ic['r'], ic['elevation'])
        func = potentials.f_spherical
    print(nelder_mead(func, [make_mpfr(argv[2]), make_mpfr(argv[3]), make_mpfr(argv[4])], [make_mpfr(argv[5]), make_mpfr(argv[6]), make_mpfr(argv[7])]))

