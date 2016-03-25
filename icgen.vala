using GLib;
using GLib.Math;
using Gsl;

namespace Sim {

    public class IcGenerator : GLib.Object {

        private struct IcGenParams {
            public double mu2;
            public double rMin;
            public double rMax;
            public double thMin;
            public double a;
        }

        private static double R (double r, double E, double L, double Q, void* params) {
            double a = ((IcGenParams*) params) -> a;
            double mu2 = ((IcGenParams*) params) -> mu2;
            return (E * (r * r + a * a) - a * L) * (E * (r * r + a * a) - a * L)
                    - (r * r + a * a - 2.0 * r) * (mu2 * r * r + Q + (L - a * E) * (L - a * E));
        }

        private static double dRdr (double r, double E, double L, double Q, void* params) {
            double a = ((IcGenParams*) params) -> a;
            double mu2 = ((IcGenParams*) params) -> mu2;
            return 4.0 * r * E * (E * (r * r + a * a) - a * L)
                    - (2.0 * r - 2.0) * (mu2 * r * r + Q + (L - a * E) * (L - a * E)) - 2.0 * mu2 * r * (r * r + a * a - 2.0 * r);
        }

        private static double THETA (double theta, double E, double L, double Q, void* params) {
            double a = ((IcGenParams*) params) -> a;
            double mu2 = ((IcGenParams*) params) -> mu2;
            return Q - cos(theta) * cos(theta) * (a * a * (mu2 - E * E) + L * L / (sin(theta) * sin(theta)));
        }

        private static int spherical (Vector x, void* params, Vector f) {
            double E = x.get(0);
            double L = x.get(1);
            double Q = x.get(2);

            f.set(0, R(((IcGenParams*) params) -> rMax, E, L, Q, params));
            f.set(1, dRdr(((IcGenParams*) params) -> rMax, E, L, Q, params));
            f.set(2, THETA(((IcGenParams*) params) -> thMin, E, L, Q, params));

            return Status.SUCCESS;
        }

        private static int nonSpherical (Vector x, void* params, Vector f) {
            double E = x.get(0);
            double L = x.get(1);
            double Q = x.get(2);

            f.set(0, R(((IcGenParams*) params) -> rMin, E, L, Q, params));
            f.set(1, R(((IcGenParams*) params) -> rMax, E, L, Q, params));
            f.set(2, THETA(((IcGenParams*) params) -> thMin, E, L, Q, params));

            return Status.SUCCESS;
        }

        private void print_state (size_t iter, MultirootFsolver s) {
            stderr.printf("iter = %3u x = % .6f % .6f % .6f f(x) = % .3e % .3e % .3e\n",
                            (uint) iter, s.x.get(0), s.x.get(1), s.x.get(2), s.f.get(0), s.f.get(1), s.f.get(2));
        }

        private void print_inital_conditions (MultirootFsolver s, void* params) {
            stdout.printf("{ \"M\" : %.1f,\n", 1.0);
            stdout.printf("  \"a\" : %.1f,\n", ((IcGenParams*) params) -> a);
            stdout.printf("  \"mu\" : %.1f,\n", ((IcGenParams*) params) -> mu2);
            stdout.printf("  \"E\" : %.12f,\n", s.x.get(0));
            stdout.printf("  \"Lz\" : %.12f,\n", s.x.get(1));
            stdout.printf("  \"C\" : %.12f,\n", s.x.get(2));
            stdout.printf("  \"r\" : %.1f,\n", ((IcGenParams*) params) -> rMin);
            stdout.printf("  \"theta\" : %.9f,\n", 0.5 * PI);
            stdout.printf("  \"start\" : %.1f,\n", 0.0);
            stdout.printf("  \"duration\" : %.1f,\n", 5000.0);
            stdout.printf("  \"step\" : %.3f,\n", 0.001);
            stdout.printf("  \"plotratio\" : %.1d,\n", 500);
            stdout.printf("  \"integrator\" : \"%s\"\n", "sc4");
            stdout.printf("}\n");
        }

        public void generate (string[] args) {
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
            IcGenParams params;
            MultirootFunction f;
            switch (args.length) {
                case 6:
                    params = { 1.0, double.parse(args[2]), double.parse(args[3]), (1.0 - double.parse(args[4])) * PI, double.parse(args[5]) };
                    f = MultirootFunction() { f = nonSpherical, n = nDim, params = &params };
                    break;
                case 5:
                    params = { 1.0, double.parse(args[2]), double.parse(args[2]), (1.0 - double.parse(args[3])) * PI, double.parse(args[4]) };
                    f = MultirootFunction() { f = spherical, n = nDim, params = &params };
                    break;
                default:
                    return_if_reached();
            }
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

            print_inital_conditions(solver, &params);
        }
    }

    public static void main (string[] args) {
        new IcGenerator().generate(args);
    }
}
