using GLib;
using GLib.Math;
using Gsl;

public class MultiRootSample : GLib.Object
{
    struct IcGenParams {
        public double mu2;
        public double rMin;
        public double rMax;
        public double thMin;
        public double a;
    }

    static double R (double r, double E, double L, double Q, void* params) {
        double a = ((IcGenParams*) params) -> a;
        double mu2 = ((IcGenParams*) params) -> mu2;
        return (E * (r*r + a*a) - a * L) * (E * (r*r + a*a) - a * L) - (r*r + a*a - 2.0*r) * (mu2 * r*r + Q + (L - a * E) * (L - a * E));
    }

    static double dR (double r, double E, double L, double Q, void* params) {
        double a = ((IcGenParams*) params) -> a;
        double mu2 = ((IcGenParams*) params) -> mu2;
        return - (2.0*r - 2.0) * (mu2 * r*r + Q + (L - a * E) * (L - a * E)) - 2.0*mu2 * r * (r*r + a*a - 2.0*r) + 4*r*E*(E * (r*r + a*a) - a * L);
    }

    static double THETA (double theta, double E, double L, double Q, void* params) {
        double a = ((IcGenParams*) params) -> a;
        double mu2 = ((IcGenParams*) params) -> mu2;
        return Q - cos(theta)*cos(theta) * (a*a * (mu2 - E*E) + L*L / (sin(theta)*sin(theta)));
    }

    static int spherical (Vector x, void* params, Vector f) {
        double E = x.get(0);
        double L = x.get(1);
        double Q = x.get(2);

        f.set(0, R(((IcGenParams*) params) -> rMin, E, L, Q, params));
        f.set(1, dR(((IcGenParams*) params) -> rMin, E, L, Q, params));
        f.set(2, THETA(((IcGenParams*) params) -> thMin, E, L, Q, params));

        return Status.SUCCESS;
    }

    static int nonSpherical (Vector x, void* params, Vector f) {
        double E = x.get(0);
        double L = x.get(1);
        double Q = x.get(2);

        f.set(0, R(((IcGenParams*) params) -> rMin, E, L, Q, params));
        f.set(1, R(((IcGenParams*) params) -> rMax, E, L, Q, params));
        f.set(2, THETA(((IcGenParams*) params) -> thMin, E, L, Q, params));

        return Status.SUCCESS;
    }

    static void print_state (size_t iter, MultirootFsolver s) {
        stdout.printf("iter = %3u x = % .3f % .3f % .3f f(x) = % .3e % .3e % .3e\n", (uint) iter, s.x.get(0), s.x.get(1), s.x.get(2), s.f.get(0), s.f.get(1), s.f.get(2));
    }

    public static void main (string[] args) {
        size_t nDim = 3;
        Vector x = new Vector(nDim);
        x.set(0, 1.0);
        x.set(1, 5.0);
        x.set(2, 0.0);

        MultirootFsolver solver;
        switch (args[1]) {
            case "dnewton":
                solver = new MultirootFsolver((MultirootFsolverType*) MultirootFsolverTypes.dnewton, nDim);
                break;
            case "broyden":
                solver = new MultirootFsolver((MultirootFsolverType*) MultirootFsolverTypes.broyden, nDim);
                break;
            case "hybrid":
                solver = new MultirootFsolver((MultirootFsolverType*) MultirootFsolverTypes.hybrid, nDim);
                break;
            case "hybrids":
                solver = new MultirootFsolver((MultirootFsolverType*) MultirootFsolverTypes.hybrids, nDim);
                break;
            default:
                return_if_reached();
        }
        IcGenParams params = { 1.0, 3.0, 12.0, 0.15 * PI, 1.0 };
        MultirootFunction f = MultirootFunction() { f = nonSpherical, n = nDim, params = &params };
        solver.set (&f, x);

        int status = 0;
        size_t iterations = 0;
        print_state(iterations, solver);
        do {
            iterations++;
            status = solver.iterate();
            print_state(iterations, solver);
            if ((bool) status)
                break;
            status = MultirootTest.residual(solver.f, 1.0e-9);
        } while (status == Status.CONTINUE && iterations < 1000);
    }
}

