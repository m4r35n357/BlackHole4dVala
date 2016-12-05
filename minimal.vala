using GLib.Math;
using Json;

namespace Simulations {

    public class BhRk4 : GLib.Object {
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
        private double timeStep;
        private int64 plotRatio;
        private double ra2;
        private double sth2;
        private double D_r;
        private double D_th;
        private double S;
        private double R;
        private double TH;
        private double t = 0.0;
        private double r;
        private double th;
        private double ph = 0.0;
        private double Ut;
        private double Ur;
        private double Uth;
        private double Uph;
        private double[] kt = { 0.0, 0.0, 0.0, 0.0 };
        private double[] kr = { 0.0, 0.0, 0.0, 0.0 };
        private double[] kth = { 0.0, 0.0, 0.0, 0.0 };
        private double[] kph = { 0.0, 0.0, 0.0, 0.0 };
        private int signR = -1;
        private int signTH = -1;

        public BhRk4 (double l, double a, double mu2, double E, double L, double Q, double r0, double th0, double t0, double tN, double h, int64 pr) {
            this.l_3 = l / 3.0;
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
            this.start = t0;
            this.end = tN;
            this.timeStep = h;
            this.plotRatio = pr;
            this.r = r0;
            this.th = (90.0 - th0) * PI / 180.0;
            evaluate(r, th, 0);
        }

        private void evaluate (double radius, double theta, int stage) {
            var r2 = radius * radius;
            ra2 = r2 + a2;
            var P = ra2 * E - aL;
            D_r = (1.0 - l_3 * r2) * ra2 - 2.0 * radius;
            R = X2 * P * P - D_r * (mu2 * r2 + K);
            sth2 = sin(theta) * sin(theta);
            var cth2 = 1.0 - sth2;
            var T = aE * sth2 - L;
            D_th = 1.0 + a2l_3 * cth2;
            TH = D_th * (K - a2mu2 * cth2) - X2 * T * T / sth2;
            var P_Dr = P / D_r;
            var T_Dth = T / D_th;
            S = r2 + a2 * cth2;
            Ut = (P_Dr * ra2 - T_Dth * a) * X2 / S;
            Ur = sqrt(R > 0.0 ? R : -R) / S;
            Uth = sqrt(TH > 0.0 ? TH : -TH) / S;
            Uph = (P_Dr * a - T_Dth / sth2) * X2 / S;
            kt[stage] = timeStep * Ut;
            kr[stage] = timeStep * Ur;
            kth[stage] = timeStep * Uth;
            kph[stage] = timeStep * Uph;
        }

        private void plot (double tau) {
            var U1 = a * Ut - ra2 * Uph;
            var U4 = Ut - a * sth2 * Uph;
            var SX2 = S * X2;
            var e = fabs(mu2 + sth2 * D_th / SX2 * U1 * U1 + S / D_r * Ur * Ur + S / D_th * Uth * Uth - D_r / SX2 * U4 * U4);
            stdout.printf("{\"tau\":%.9e, \"v4e\":%.1f, \"D_r\":%.9e, \"R\":%.9e, \"TH\":%.9e, \"t\":%.9e, \"r\":%.9e, \"th\":%.9e, \"ph\":%.9e, \"tP\":%.9e, \"rP\":%.9e, \"thP\":%.9e, \"phP\":%.9e}\n",
                            tau, 10.0 * log10(e > 1.0e-18 ? e : 1.0e-18), D_r, R, TH, t, r, th, ph, Ut, Ur, Uth, Uph);
        }

        public void solve () {
            var tau = 0.0;
            var iterationCount = 0;
            while (tau < end) {
                if ((tau >= start) && (iterationCount % plotRatio == 0)) {
                    plot(tau);
                }
                signR = R > 0.0 ? signR : -signR;
                signTH = TH > 0.0 ? signTH : -signTH;
                evaluate(r + 0.5 * kr[0], th + 0.5 * kth[0], 1);
                evaluate(r + 0.5 * kr[1], th + 0.5 * kth[1], 2);
                evaluate(r + kr[2], th + kth[2], 3);
                t += (kt[0] + 2.0 * (kt[1] + kt[2]) + kt[3]) / 6.0;
                r += (kr[0] + 2.0 * (kr[1] + kr[2]) + kr[3]) / 6.0 * signR;
                th += (kth[0] + 2.0 * (kth[1] + kth[2]) + kth[3]) / 6.0 * signTH;
                ph += (kph[0] + 2.0 * (kph[1] + kph[2]) + kph[3]) / 6.0;
                evaluate(r, th, 0);
                iterationCount += 1;
                tau = iterationCount * timeStep;
            }
            plot(tau);
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
        new BhRk4(ic.get_double_member("lambda"), ic.get_double_member("a"), ic.get_double_member("mu"), ic.get_double_member("E"),
                ic.get_double_member("L"), ic.get_double_member("Q"), ic.get_double_member("r0"), ic.get_double_member("th0"),
                ic.get_double_member("start"), ic.get_double_member("end"), ic.get_double_member("step"), ic.get_int_member("plotratio")).solve();
    }
}
