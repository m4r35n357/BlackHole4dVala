from sys import stderr

def nelder_mead(f, x_0, x_δ, ε, stuck_limit=100, limit=1000, α=1.0, γ=2.0, ρ=-0.5, σ=0.5):
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
        if stuck_count > stuck_limit:
            raise RuntimeError("STUCK for {} steps! ".format(stuck_count) + data)
        if count > limit:
            raise RuntimeError("ABANDONED after {} steps! ".format(count) + data)
        print(data, file=stderr)
        if max([abs(s[0][0][i] - c[i]) for i in dim]) < ε and abs(s[0][1] - s[-1][1]) < ε:
            return s[0], count, nr, ne, nc, ns
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

print(__name__ + " module loaded", file=stderr)
