from sys import argv, float_info
from NelderMead import nelder_mead

class Rosenbrock(object):

    def __init__(self, a, b):
        self.a = a
        self.b = b

    def func(self, z):
        x, y = z
        return (self.a - x)**2 + self.b * (y - x**2)**2

def himmelblau(z):
    x, y = z
    return (x**2 + y - 11)**2 + (x + y**2 - 7)**2

def basic(z):
    x, y = z
    f = x**2 + y**2 - 4.0
    c = 0.0 if x + y >= 2.0 else float_info.max
    return f + c

if __name__ == "__main__":
    # Example: python Rosenbrock.py 3.0 3.0 0.1 0.1
    init = [float(argv[1]), float(argv[2])]
    deltas = [float(argv[3]), float(argv[4])]
    print(nelder_mead(Rosenbrock(1.0, 100.0).func, init, deltas, 1.0e-9))
    print(nelder_mead(himmelblau, init, deltas, 1.0e-9))
    print(nelder_mead(basic, init, deltas, 1.0e-9))
