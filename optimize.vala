using GLib;
using Gsl;

public class MultiRootSample : GLib.Object
{
    struct RParams {
        public double a;
        public double b;
    }

    static int rosenbrock_f (Vector x, void* params, Vector f) {
        double a = ((RParams*)params)->a;
        double b = ((RParams*)params)->b;

        double x0 = x.get (0);
        double x1 = x.get (1);

        double y0 = a * (1.0 - x0);
        double y1 = b * (x1 - x0 * x0);

        f.set (0, y0);
        f.set (1, y1);

        return Status.SUCCESS;
    }

    static void print_state (size_t iter, MultirootFsolver s) {
        stdout.printf ("iter = %3u x = % .3f % .3f f(x) = % .3e % .3e\n", (uint)iter, s.x.get (0), s.x.get (1), s.f.get (0), s.f.get (1));
    }

    public static void main (string[] args) {
        MultirootFsolverType* T;
        MultirootFsolver s;

        int status = 0;
        size_t iter=0;

        size_t n = 2;

        RParams p = { 1.0, 10.0 };
        MultirootFunction f = MultirootFunction () { f = rosenbrock_f, n = n, params = &p };

        double[] x_init = new double[] { -10.0, -5.0 };
        Vector x = new Vector (n);

        x.set (0, x_init[0]);
        x.set (1, x_init[1]);

        T = (MultirootFsolverType*)MultirootFsolverTypes.dnewton;
        s = new MultirootFsolver (T, n);
        s.set (&f, x);

        print_state (iter, s);

        do {
            iter++;
            status = s.iterate ();
            print_state (iter, s);
            if ((bool)status)
                break;
            status = MultirootTest.residual (s.f, 1.0e-7);
        } while (status==Status.CONTINUE && iter < 1000);
    }
}

