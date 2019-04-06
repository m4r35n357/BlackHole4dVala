from copy import copy
from sys import stderr

def replace_worst(data, new):
    del data[-1]
    data.append(new)

def nelder_mead(f, x_0, x_δ, stuck_threshold=10e-12, stuck_break=10, max_iterations=0, α=1.0, γ=2.0, ρ=-0.5, σ=0.5):

    # initial simplex
    dim = len(x_0)
    assert dim == len(x_δ)
    dimensions = range(dim)
    prev_best = f(x_0)
    results = [[x_0, prev_best]]
    for i in dimensions:
        x = copy(x_0)
        x[i] += x_δ[i]
        score = f(x)
        results.append([x, score])
    iterations = no_improvement = nt = nr = ne = nc = ns = 0

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
        print('...best so far:', best, file=stderr)
        if best < prev_best - stuck_threshold:
            no_improvement = 0
            prev_best = best
        else:
            no_improvement += 1
        if no_improvement >= stuck_break:
            return results[0], nt - 1, nr, ne, nc, ns

        # 2. centroid
        x0 = [0.0] * dim
        for result in results[:-1]:
            for i, c in enumerate(result[0]):
                x0[i] += c / (len(results)-1)

        # 3. reflect
        xr = [x0[i] + α * (x0[i] - results[-1][0][i]) for i in dimensions]
        r_score = f(xr)
        if results[0][1] <= r_score < results[-2][1]:
            nr += 1
            replace_worst(results, [xr, r_score])
            continue

        # 4. expand
        if r_score < results[0][1]:
            xe = [x0[i] + γ * (x0[i] - results[-1][0][i]) for i in dimensions]
            e_score = f(xe)
            ne += 1
            replace_worst(results, [xe, e_score] if e_score < r_score else [xr, r_score])
            continue

        # 5. contract
        xc = [x0[i] + ρ * (x0[i] - results[-1][0][i]) for i in dimensions]
        c_score = f(xc)
        if c_score < results[-1][1]:
            nc += 1
            replace_worst(results, [xc, c_score])
            continue

        # 6. shrink
        x1 = results[0][0]
        reduced = []
        for result in results:
            xs = [x1[i] + σ * (result[0][i] - x1[i]) for i in dimensions]
            score = f(xs)
            ns += 1
            reduced.append([xs, score])
        results = reduced

print(__name__ + " module loaded", file=stderr)