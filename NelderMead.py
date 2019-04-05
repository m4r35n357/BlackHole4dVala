
import copy
from json import loads
from sys import argv, stdin, stderr
# import numpy as np
from gmpy2 import mpfr, get_context, sin, acos, cos
get_context().precision = 113  # Set this BEFORE importing any Taylor Series stuff!
from Symplectic import D1, D2
from dual import Dual, make_mpfr


def nelder_mead(f, x_start,
                step=0.1, no_improve_thr=10e-12,
                no_improv_break=10, max_iter=0,
                alpha=1., gamma=2., rho=-0.5, sigma=0.5):

    # init
    dim = len(x_start)
    prev_best = f(x_start)
    no_improvement = 0
    res = [[x_start, prev_best]]

    for i in range(dim):
        x = copy.copy(x_start)
        x[i] = x[i] + step
        score = f(x)
        res.append([x, score])

    # simplex iter
    iters = 0
    while 1:
        # order
        res.sort(key=lambda x: x[1])
        best = res[0][1]

        # break after max_iter
        if max_iter and iters >= max_iter:
            return res[0]
        iters += 1

        # break after no_improvement iterations with no improvement
        # print('...best so far:', best)

        if best < prev_best - no_improve_thr:
            no_improvement = 0
            prev_best = best
        else:
            no_improvement += 1

        if no_improvement >= no_improv_break:
            return res[0]

        # centroid
        x0 = [0.] * dim
        for tup in res[:-1]:
            for i, c in enumerate(tup[0]):
                x0[i] += c / (len(res)-1)

        # reflection
        worst = res[-1][0]
        xr = [0.0, 0.0, 0.0]
        for i in range(dim):
            xr[i] = x0[i] + alpha*(x0[i] - worst[i])
        rscore = f(xr)
        if res[0][1] <= rscore < res[-2][1]:
            del res[-1]
            res.append([xr, rscore])
            continue

        # expansion
        if rscore < res[0][1]:
            worst = res[-1][0]
            xe = [0.0, 0.0, 0.0]
            for i in range(dim):
                xe[i] = x0[i] + gamma * (x0[i] - worst[i])
            escore = f(xe)
            if escore < rscore:
                del res[-1]
                res.append([xe, escore])
                continue
            else:
                del res[-1]
                res.append([xr, rscore])
                continue

        # contraction
        worst = res[-1][0]
        xc = [0.0, 0.0, 0.0]
        for i in range(dim):
            xc[i] = x0[i] + rho*(x0[i] - worst[i])
        cscore = f(xc)
        if cscore < res[-1][1]:
            del res[-1]
            res.append([xc, cscore])
            continue

        # reduction
        x1 = res[0][0]
        nres = []
        for tup in res:
            redx = x1 + sigma*(tup[0] - x1)
            score = f(redx)
            nres.append([redx, score])
        res = nres


class Potentials(object):
    def __init__(self, a, r_min, r_max, elevation):
        self.a = a
        self.a2 = a**2
        self.μ2 = make_mpfr(1)
        self.r_min = make_mpfr(r_min)
        self.r_max = make_mpfr(r_max)
        self.θ = make_mpfr((D1 - (make_mpfr(90) - elevation) / make_mpfr(180)) * acos(make_mpfr(-1)))

    def f_r(self, r, E, L, Q):
        r2 =r.sqr
        ra2 = r2 + self.a2
        return (ra2 * E - self.a * L).sqr - (ra2 - D2 * r) * (self.μ2 * r2 + Q + (L - self.a * E)**2)

    def f_θ(self, θ, E, L, Q):
        return Q - cos(θ)**2 * (self.a2 * (self.μ2 - E**2) + (L / sin(θ))**2)

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
    print("Generator: {}".format(argv[0]), file=stderr)
    if argv[1]:
        input_data = open(argv[1]).read()
    else:
        input_data = stdin.read()
    ic = loads(input_data, parse_float=mpfr)
    print(input_data, file=stderr)
    if ic.get('rMax'):
        potentials = Potentials(ic['spin'], ic['rMin'], ic['rMax'], ic['elevation'])
        f = potentials.f_nonspherical
        print('NONSPHERICAL')
    else:
        potentials = Potentials(ic['spin'], ic['r'], ic['r'], ic['elevation'])
        f = potentials.f_spherical
        print('SPHERICAL')
    print(nelder_mead(f, [0.9, 1.5, 8.0]))

