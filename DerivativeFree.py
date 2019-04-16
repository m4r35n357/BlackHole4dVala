from sys import argv
from NelderMead import secant,  bisect
from timeit import default_timer as timer

def polynomial(x):
    return x * (x - 1) * (x - 2) * (x - 3)

if __name__ == "__main__":
    # Example: python DerivativeFree.py 1.0 0.5 2>/dev/null
    tol = 1.0e-12
    value = float(argv[1])
    variation = float(argv[2])
    n = 100000
    result = None
    try:
        print("SECANT")
        start = timer()
        for run in range(n):
            result = secant(polynomial, value - variation, value + variation, tol)
        end = timer()
        print("{:8.3f} {}".format(end - start, result))
    except RuntimeError as e:
        print(e)
    try:
        print("BISECTION")
        start = timer()
        for run in range(n):
            result = bisect(polynomial, value - variation, value + variation, tol)
        end = timer()
        print("{:8.3f} {}".format(end - start, result))
    except RuntimeError as e:
        print(e)
