#!/usr/bin/env pypy
from json import loads
from math import sqrt, sin, pi, fabs, log10, cos
from sys import stdin, stderr, stdout, argv


class BhSymp(object):
    def __init__(self, Lambda, a, mu2, E, L, C, r0, thetaMin, start, end, timestep):
        self.l_3 = Lambda / 3.0
        self.a = a
        self.mu2 = mu2
        self.E = E
        self.L = L
        self.a2 = a**2
        self.a2l_3 = self.a2 * self.l_3
        self.a2mu2 = self.a2 * mu2
        self.aE = a * E
        self.aL = a * L
        self.X2 = (1.0 + self.a2l_3)**2
        self.K = C + self.X2 * (L - self.aE)**2
        self.start_time = start
        self.end_time = end
        self.h = timestep
        self.sgnR = self.sgnTH = -1
        self.t = self.ph = 0.0
        self.r = r0
        self.th = (90.0 - thetaMin) * pi / 180.0
        self.cbrt2 = 2.0**(1.0 / 3.0)
        self.y = [0.5 * 1.0 / (2.0 - self.cbrt2), 1.0 / (2.0 - self.cbrt2), 0.5 * (1.0 - self.cbrt2) * 1.0 / (2.0 - self.cbrt2), - self.cbrt2 * 1.0 / (2.0 - self.cbrt2)]
        self.refresh()
        self.Ur = - sqrt(fabs(self.X2 * self.P**2 - self.D_r * (self.mu2 * self.r2 + self.K)))
        self.Uth = - sqrt(fabs(self.D_th * (self.K - self.a2mu2 * self.cth2) - self.X2 * self.T**2 / self.sth2))

    def refresh(self):
        self.r2 = self.r**2
        self.sth2 = sin(self.th)**2
        self.cth2 = 1.0 - self.sth2
        self.ra2 = self.r2 + self.a2
        self.P = self.ra2 * self.E - self.aL
        self.D_r = (1.0 - self.l_3 * self.r2) * self.ra2 - 2.0 * self.r
        self.T = self.aE * self.sth2 - self.L
        self.D_th = 1.0 + self.a2l_3 * self.cth2
        P_Dr = self.P / self.D_r
        T_Dth = self.T / self.D_th
        self.S = self.r2 + self.a2 * self.cth2
        self.Ut = (P_Dr * self.ra2 - T_Dth * self.a) * self.X2
        self.Uph = (P_Dr * self.a - T_Dth / self.sth2) * self.X2

    def qUp(self, d):
        self.t += d * self.h * self.Ut
        self.r += d * self.h * self.Ur
        self.th += d * self.h * self.Uth
        self.ph += d * self.h * self.Uph
        self.refresh()

    def pUp(self, c):
        self.Ur += c * self.h * (self.r * (2.0 * self.E * self.P * self.X2- self.mu2 * self.D_r)
                                 - (self.r * (1.0 - self.l_3 * self.r2) - 1.0 - self.l_3 * self.r * self.ra2) * (self.K + self.mu2 * self.r2))
        sth = sin(self.th)
        self.Uth += c * self.h * cos(self.th) * (sth * self.a2 * (self.mu2 * self.D_th - self.l_3 * (self.K - self.a2mu2 * self.cth2))
                                  + self.X2 * self.T / sth * (self.T / self.sth2 - 2.0 * self.aE))

    def plot(self, mino, tau):
        S3X2 = self.S**2 * self.S * self.X2
        e = fabs(self.mu2 + self.sth2 * self.D_th / S3X2 * (self.a * self.Ut - self.ra2 * self.Uph)**2
                    + self.Ur**2 / (self.S * self.D_r)  + self.Uth**2 / (self.S * self.D_th)
                    - self.D_r / S3X2 * (self.Ut - self.a * self.sth2 * self.Uph)**2)
        print >> stdout, '{"mino":%.9e, "tau":%.9e, "v4e":%.1f, "t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "tP":%.9e, "rP":%.9e, "thP":%.9e, "phP":%.9e}' \
                % (mino, tau, 10.0 * log10(e if e > 1.0e-18 else 1.0e-18), self.t, self.r, self.th, self.ph, self.Ut / self.S, self.Ur / self.S, self.Uth / self.S, self.Uph / self.S)

    def solve(self):
        mino = tau = 0.0
        i = 0
        while tau <= self.end_time:
            if tau >= self.start_time:
                self.plot(mino, tau)
            self.qUp(self.y[0]); self.pUp(self.y[1]); self.qUp(self.y[2]); self.pUp(self.y[3]); self.qUp(self.y[2]); self.pUp(self.y[1]); self.qUp(self.y[0])
            i += 1
            mino = i * self.h
            tau += self.h * self.S
        self.plot(mino, tau)


def main():
    print >> stderr, "Simulator: {}".format(argv[0])
    input_data = stdin.read()
    ic = loads(input_data)['IC']
    print >> stderr, input_data
    BhSymp(ic['lambda'], ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'], ic['start'], ic['end'], ic['step']).solve()


if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"
