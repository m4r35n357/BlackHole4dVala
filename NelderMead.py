
from sys import stderr
from copy import copy
from collections import namedtuple

Point = namedtuple('PointType', ['x', 'f'])

def replace_worst(data, new):
    del data[-1]
    data.append(new)

def nelder_mead(f, x_0, x_δ, ε=10e-12, stuck_break=100, max_iterations=1000, α=1.0, γ=2.0, ρ=-0.5, σ=0.5):

    # initial simplex
    dim = len(x_0)
    assert dim == len(x_δ)
    dimensions = range(dim)
    prev_best = f(x_0)
    results = [Point(x=x_0, f=prev_best)]
    for i in dimensions:
        x = copy(x_0)
        x[i] += x_δ[i]
        s_score = f(x)
        results.append(Point(x=x, f=s_score))
    iterations = stuck_counter = nt = nr = ne = nc = ns = 0
    latest = ""

    while True:
        nt += 1

        # 1. order
        results.sort(key=lambda z: z.f)
        best = results[0]
        worst = results[-1]
        second_worst = results[-2]
        best_value = best.f

        data = '{} {} {}'.format(iterations, results, latest)
        if max_iterations and iterations >= max_iterations:
            raise RuntimeWarning("UNFINISHED! " + data)
        print(data, file=stderr)
        if best_value < prev_best:
            stuck_counter = 0
            prev_best = best_value
        else:
            stuck_counter += 1
        if sum((best.x[i] - worst.x[i])**2 for i in dimensions) < ε and (best.f - worst.f)**2 < ε:
            return best, nt - 1, nr, ne, nc, ns
        if stuck_counter >= stuck_break:
            raise RuntimeWarning("STUCK for {} steps!".format(stuck_counter) + data)
        iterations += 1

        # 2. centroid
        centroid = [0.0] * dim
        for result in results[:-1]:
            for i, c in enumerate(result.x):
                centroid[i] += c / (len(results)-1)

        # 3. reflect
        xr = [centroid[i] + α * (centroid[i] - worst.x[i]) for i in dimensions]
        r_score = f(xr)
        if best[1] <= r_score < second_worst[1]:
            nr += 1
            replace_worst(results, Point(x=xr, f=r_score))
            latest = "reflection"
            continue

        # 4. expand
        if r_score < best[1]:
            xe = [centroid[i] + γ * (centroid[i] - worst.x[i]) for i in dimensions]
            e_score = f(xe)
            ne += 1
            replace_worst(results, Point(x=xe, f=e_score) if e_score < r_score else Point(x=xr, f=r_score))
            latest = "expansion" + ("(e)" if e_score < r_score else "(r)")
            continue

        # 5. contract
        xc = [centroid[i] + ρ * (centroid[i] - worst.x[i]) for i in dimensions]
        c_score = f(xc)
        if c_score < worst.f:
            nc += 1
            replace_worst(results, Point(x=xc, f=c_score))
            latest = "contraction"
            continue

        # 6. shrink
        reduced = []
        for vertex in results[1:]:
            xs = [best.vertices[i] + σ * (results[vertex][i] - best.vertices[i]) for i in dimensions]
            ns += 1
            reduced.append(Point(x=xs, f=f(xs)))
        results = reduced
        latest = "reduction"


print(__name__ + " module loaded", file=stderr)