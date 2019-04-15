from sys import argv
from NelderMead import secant, false_position, bisect


def polynomial(x):
    return x * (x - 1) * (x - 2) * (x - 3)

if __name__ == "__main__":
    # Example: python DerivativeFree.py 1.0 0.1
    tol = 1.0e-9
    value = float(argv[1])
    variation = float(argv[2])
    try:
        print(secant(polynomial, value - variation, value + variation, tol))
    except RuntimeError as e:
        print(e)
    try:
        print(false_position(polynomial, value - variation, value + variation, tol))
    except RuntimeError as e:
        print(e)
    try:
        print(bisect(polynomial, value - variation, value + variation, tol))
    except RuntimeError as e:
        print(e)
