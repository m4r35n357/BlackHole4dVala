from sys import stderr

def nelder_mead(f, x_0, x_δ, ε, stuck=100, limit=1000, α=1.0, γ=2.0, ρ=-0.5, σ=0.5):
    n = len(x_0)
    assert n == len(x_δ)
    dim = range(n)
    best = f(x_0)
    s = [[x_0, best]]
    for i in dim:
        v = [x for x in x_0]
        v[i] += x_δ[i]
        s.append([v, f(v)])
    count = stuck_count = nr = ne = nc = ns = 0
    latest = ""

    while True:
        s.sort(key=lambda z: z[1])
        c = [0.0] * n
        for v in s[:-1]:
            for i, x in enumerate(v[0]):
                c[i] += x / n

        data = '{} {} {}'.format(count, s, latest)
        if s[0][1] < best:
            stuck_count = 0
            best = s[0][1]
        else:
            stuck_count += 1
        if stuck and stuck_count > stuck:
            raise RuntimeError("STUCK for {} steps! ".format(stuck_count) + data)
        if limit and count > limit:
            raise RuntimeError("ABANDONED after {} steps! ".format(count) + data)
        if ε and max([abs(s[0][0][i] - c[i]) for i in dim]) < ε and abs(s[0][1] - s[-1][1]) < ε:
                return s[0], count, nr, ne, nc, ns
        else:
            print(s[0], count, nr, ne, nc, ns, file=stderr)
        count += 1

        xr = [c[i] + α * (c[i] - s[-1][0][i]) for i in dim]
        fr = f(xr)
        if s[0][1] <= fr < s[-2][1]:
            nr += 1
            del s[-1]
            s.append([xr, fr])
            latest = "reflection"
            continue

        if fr < s[0][1]:
            xe = [c[i] + γ * (c[i] - s[-1][0][i]) for i in dim]
            fe = f(xe)
            ne += 1
            del s[-1]
            s.append([xe, fe] if fe < fr else [xr, fr])
            latest = "expansion" + ("(e)" if fe < fr else "(r)")
            continue

        xc = [c[i] + ρ * (c[i] - s[-1][0][i]) for i in dim]
        fc = f(xc)
        if fc < s[-1][1]:
            nc += 1
            del s[-1]
            s.append([xc, fc])
            latest = "contraction"
            continue

        reduced = [s[0]]
        for v in s[1:]:
            xs = [s[0][0][i] + σ * (v[0][i] - s[0][0][i]) for i in dim]
            reduced.append([xs, f(xs)])
        ns += 1
        s = reduced
        latest = "reduction"

def secant(f, a, b, ε, limit=101):
    f_a = f(a)
    f_b = f(b)
    count = δx = c = f_c = 1
    while abs(f_c) > ε or abs(δx) > ε:
        if count == limit:
            raise RuntimeError("{}\n After {} iterations, current: {}, previous: {}".format(f, count - 1, b, a))
        c = (b * f_a - a * f_b) / (f_a - f_b)
        f_c = f(c)
        b = a
        f_b = f_a
        a = c
        f_a = f_c
        δx = b - a
        count += 1
    return c, f_c, δx, count - 1

def bisect(f, a, b, ε, limit=101):
    f_a = f(a)
    count = δx = c = f_c = 1
    while abs(f_c) > ε or abs(δx) > ε:
        if count == limit:
            raise RuntimeError("{}\n After {} iterations, a: {}, b: {}".format(f, count - 1, a, b))
        c = (a + b) / 2
        f_c = f(c)
        if f_a * f_c > 0.0:
            a = c
        else:
            b = c
        δx = b - a
        count += 1
    return c, f_c, δx, count - 1

print(__name__ + " module loaded", file=stderr)
