#!/usr/bin/env python3

from json import loads
from sys import argv, stderr
from gmpy2 import mpfr, get_context, sin, acos, cos
get_context().precision = 53  # Set this BEFORE importing any mathematical stuff!
from Symplectic import D1, D2
from NelderMead import nelder_mead
from dual import Dual, make_mpfr


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
    # Example: ./Generator.py icgen-data.json 1.0 5.0 0.0 1.0 1.0 1.0
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
    print(nelder_mead(func, [make_mpfr(argv[2]), make_mpfr(argv[3]), make_mpfr(argv[4])],
                            [make_mpfr(argv[5]), make_mpfr(argv[6]), make_mpfr(argv[7])], 1.0e-9))

