
from sys import stderr
from copy import copy
from collections import namedtuple

Point = namedtuple('PointType', ['x', 'f'])

def replace_worst(data, new):
    del data[-1]
    data.append(new)

def nelder_mead(f, x_0, x_δ, ε=10e-12, stuck_break=100, max_iterations=1000, α=1.0, γ=2.0, ρ=-0.5, σ=0.5):

    dim = len(x_0)
    assert dim == len(x_δ)
    dimensions = range(dim)
    prev_best = f(x_0)
    simplex = [Point(x=x_0, f=prev_best)]
    for i in dimensions:
        x = copy(x_0)
        x[i] += x_δ[i]
        simplex.append(Point(x=x, f=f(x)))
    iterations = stuck_counter = nt = nr = ne = nc = ns = 0
    latest = ""

    while True:
        nt += 1

        simplex.sort(key=lambda z: z.f)
        best = simplex[0]
        worst = simplex[-1]
        second_worst = simplex[-2]
        best_value = best.f

        data = '{} {} {}'.format(iterations, simplex, latest)
        if best_value < prev_best:
            stuck_counter = 0
            prev_best = best_value
        else:
            stuck_counter += 1
        if stuck_counter >= stuck_break:
            raise RuntimeError("STUCK for {} steps! ".format(stuck_counter) + data)
        if max_iterations and iterations >= max_iterations:
            raise RuntimeError("UNFINISHED! " + data)
        print(data, file=stderr)
        if sum((best.x[i] - worst.x[i])**2 for i in dimensions) < ε and (best.f - worst.f)**2 < ε:
            return best, nt - 1, nr, ne, nc, ns
        iterations += 1

        centroid = [0.0] * dim
        for result in simplex[:-1]:
            for i, c in enumerate(result.x):
                centroid[i] += c / (len(simplex)-1)

        xr = [centroid[i] + α * (centroid[i] - worst.x[i]) for i in dimensions]
        r_score = f(xr)
        if best.f <= r_score < second_worst.f:
            nr += 1
            replace_worst(simplex, Point(x=xr, f=r_score))
            latest = "reflection"
            continue

        if r_score < best.f:
            xe = [centroid[i] + γ * (centroid[i] - worst.x[i]) for i in dimensions]
            e_score = f(xe)
            ne += 1
            replace_worst(simplex, Point(x=xe, f=e_score) if e_score < r_score else Point(x=xr, f=r_score))
            latest = "expansion" + ("(e)" if e_score < r_score else "(r)")
            continue

        xc = [centroid[i] + ρ * (centroid[i] - worst.x[i]) for i in dimensions]
        c_score = f(xc)
        if c_score < worst.f:
            nc += 1
            replace_worst(simplex, Point(x=xc, f=c_score))
            latest = "contraction"
            continue

        reduced = []
        for vertex in simplex[1:]:
            xs = [best.x[i] + σ * (vertex.x[i] - best.x[i]) for i in dimensions]
            ns += 1
            reduced.append(Point(x=xs, f=f(xs)))
        simplex = reduced
        latest = "reduction"


print(__name__ + " module loaded", file=stderr)
