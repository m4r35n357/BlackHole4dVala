using GLib.Math;
using Json;

namespace Simulations {

    public class BhSymp : GLib.Object {
        private double l_3;
        private double a;
        private double a2;
        private double a2l_3;
        private double mu2;
        private double a2mu2;
        private double X2;
        private double E;
        private double L;
        private double K;
        private double aE;
        private double aL;
        private double start;
        private double end;
        private double h;
        private int64 tr;
        private double r2;
        private double ra2;
        private double sth2;
        private double cth2;
        private double D_r;
        private double D_th;
        private double S;
        private double P;
        private double R;
        private double T;
        private double TH;
        private double t = 0.0;
        private double r;
        private double th;
        private double ph = 0.0;
        private double Ut;
        private double Ur;
        private double Uth;
        private double Uph;
        private double[] y;

        public BhSymp (double lambda, double a, double mu2, double E, double L, double Q, double r0, double th0,
                            double tau0, double tauN, double tStep, int64 plotRatio) {
            this.l_3 = lambda / 3.0;
            this.a = a;
            this.mu2 = mu2;
            this.E = E;
            this.L = L;
            this.a2 = a * a;
            this.a2l_3 = a2 * l_3;
            this.a2mu2 = a2 * mu2;
            this.aE = a * E;
            this.aL = a * L;
            this.X2 = (1.0 + a2l_3) * (1.0 + a2l_3);
            this.K = Q + X2 * (L - aE) * (L - aE);
            this.start = tau0;
            this.end = tauN;
            this.h = tStep;
            this.tr = plotRatio;
            this.r = r0;
            this.th = (90.0 - th0) * PI / 180.0;
            var CBRT2 = pow(2.0, (1.0 / 3.0));
            this.y = { 0.5 / (2.0 - CBRT2), 1.0 / (2.0 - CBRT2), 0.5 * (1.0 - CBRT2) / (2.0 - CBRT2), - CBRT2 / (2.0 - CBRT2) };
            refresh(r, th);
            this.Ur = - sqrt(fabs(R));
            this.Uth = - sqrt(fabs(TH));
        }

        private void refresh (double radius, double theta) {
            r2 = radius * radius;
            ra2 = r2 + a2;
            P = ra2 * E - aL;
            D_r = (1.0 - l_3 * r2) * ra2 - 2.0 * radius;
            R = X2 * P * P - D_r * (mu2 * r2 + K);
            sth2 = sin(theta) * sin(theta);
            cth2 = 1.0 - sth2;
            T = aE * sth2 - L;
            D_th = 1.0 + a2l_3 * cth2;
            TH = D_th * (K - a2mu2 * cth2) - X2 * T * T / sth2;
            var P_Dr = P / D_r;
            var T_Dth = T / D_th;
            S = r2 + a2 * cth2;
            Ut = (P_Dr * ra2 - T_Dth * a) * X2;
            Uph = (P_Dr * a - T_Dth / sth2) * X2;
        }

        private void qUp (double c) {
            t += c * h * Ut;
            r += c * h * Ur;
            th += c * h * Uth;
            ph += c * h * Uph;
            refresh(r, th);
        }

        private void pUp (double d) {
            Ur += d * h * (2.0 * r * E * P * X2 - (r * (1.0 - l_3 * r2) - 1.0 - l_3 * r * ra2) * (K + mu2 * r2) - mu2 * r * D_r);
            var sth = sin(th);
            var cth = cos(th);
            Uth += d * h * (cth * sth * a2 * (mu2 * D_th - l_3 * (K - a2mu2 * cth2)) + cth * X2 * T / sth * (T / sth2 - 2.0 * aE));
        }

        public void solve () {
            var mino = 0.0;
            var tau = 0.0;
            var i = 0;
            var plotCount = 0;
            while (tau < end) {
                if ((tau >= start) && (i % tr == 0)) {
                    plot(mino, tau);
                    plotCount += 1;
                }
                qUp(y[0]); pUp(y[1]); qUp(y[2]); pUp(y[3]); qUp(y[2]); pUp(y[1]); qUp(y[0]);
                i += 1;
                mino = i * h;
                tau += h * S;
            }
            plot(mino, tau);
        }

        private void plot (double mino, double tau) {
            var U1 = a * Ut - ra2 * Uph;
            var U4 = Ut - a * sth2 * Uph;
            var S3X2 = S * S * S * X2;
            var e = fabs(mu2 + (sth2 * D_th / S3X2 * U1 * U1 + Ur * Ur / (S * D_r) + Uth * Uth / (S * D_th) - D_r / S3X2 * U4 * U4));
            stdout.printf("{\"mino\":%.9e, \"tau\":%.9e, \"v4e\":%.1f, \"D_r\":%.9e, \"R\":%.9e, \"TH\":%.9e, \"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, \"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e}\n",
                        mino, tau, 10.0 * log10(e > 1.0e-18 ? e : 1.0e-18), D_r, R, TH, t, r, th, ph, Ut / S, Ur / S, Uth / S, Uph / S);
        }
    }

    public static void main (string[] args) {
        stderr.printf("Executable: %s\n", args[0]);
        var json = new StringBuilder();
        var line = stdin.read_line();
        while (line != null) {
            json.append_printf("%s\n", line);
            stderr.printf("%s\n", line);
            line = stdin.read_line();
        }
        var parser = new Parser();
        try {
            parser.load_from_data(json.str);
        } catch (Error e) {
            stderr.printf("Unable to parse the input data: %s\n", e.message);
            return_if_reached();
        }
        var ic = parser.get_root().get_object().get_object_member("IC");
        new BhSymp(ic.get_double_member("lambda"), ic.get_double_member("a"), ic.get_double_member("mu"), ic.get_double_member("E"),
                ic.get_double_member("L"), ic.get_double_member("Q"), ic.get_double_member("r0"), ic.get_double_member("th0"),
                ic.get_double_member("start"), ic.get_double_member("end"), ic.get_double_member("step"), ic.get_int_member("plotratio")).solve();
    }
}
