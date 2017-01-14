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
        self.f2 = 1.0 / (2.0 - self.cbrt2)
        self.coefficients = [0.5 * self.f2, self.f2, 0.5 * (1.0 - self.cbrt2) * self.f2, - self.cbrt2 * self.f2]
        self.refresh(self.r, self.th)
        self.Ur = - sqrt(fabs(self.R))
        self.Uth = - sqrt(fabs(self.TH))

    def refresh(self, radius, theta):
        self.r2 = radius**2
        self.sth = sin(theta)
        self.cth = cos(theta)
        self.sth2 = self.sth**2
        self.cth2 = 1.0 - self.sth2
        self.ra2 = self.r2 + self.a2
        self.P = self.ra2 * self.E - self.aL
        self.D_r = (1.0 - self.l_3 * self.r2) * self.ra2 - 2.0 * radius
        self.R = self.X2 * self.P**2 - self.D_r * (self.mu2 * self.r2 + self.K)
        self.T = self.aE * self.sth2 - self.L
        self.D_th = 1.0 + self.a2l_3 * self.cth2
        self.TH = self.D_th * (self.K - self.a2mu2 * self.cth2) - self.X2 * self.T**2 / self.sth2
        P_Dr = self.P / self.D_r
        T_Dth = self.T / self.D_th
        self.S = self.r2 + self.a2 * self.cth2
        self.Ut = (P_Dr * self.ra2 - T_Dth * self.a) * self.X2
        self.Uph = (P_Dr * self.a - T_Dth / self.sth2) * self.X2

    def v4_error(self, Ut, Ur, Uth, Uph):  # norm squared, xDot means dx/dTau !!!
        SX2 = self.S * self.X2
        return fabs(self.mu2 + self.sth2 * self.D_th / SX2 * (self.a * Ut - self.ra2 * Uph)**2
                    + self.S / self.D_r * Ur**2 + self.S / self.D_th * Uth**2
                    - self.D_r / SX2 * (Ut - self.a * self.sth2 * Uph)**2)

    @staticmethod
    def log_error(e):
        return 10.0 * log10(e if e > 1.0e-18 else 1.0e-18)


    @staticmethod
    def modH(xdot, x):
        return 0.5 * fabs(xdot**2 - x)

    def qUp(self, d):
        self.t += d * self.h * self.Ut
        self.r += d * self.h * self.Ur
        self.th += d * self.h * self.Uth
        self.ph += d * self.h * self.Uph
        self.refresh(self.r, self.th)

    def pUp(self, c):
        self.Ur += c * self.h * (2.0 * self.r * self.E * self.P * self.X2
                                 - (self.r * (1.0 - self.l_3 * self.r2) - 1.0 - self.l_3 * self.r * self.ra2) * (self.K + self.mu2 * self.r2)
                                 - self.mu2 * self.r * self.D_r)
        self.Uth += c * self.h * (self.cth * self.sth * self.a2 * (self.mu2 * self.D_th - self.l_3 * (self.K - self.a2mu2 * self.cth2))
                                  + self.cth * self.X2 * self.T / self.sth * (self.T / self.sth2 - 2.0 * self.aE))

    def plot(self, mino, tau):
        eR = self.log_error(self.modH(self.Ur, self.R))
        eTh = self.log_error(self.modH(self.Uth, self.TH))
        v4e = self.log_error(self.v4_error(self.Ut / self.S, self.Ur / self.S, self.Uth / self.S, self.Uph / self.S))  # d/dTau = 1/sigma * d/dLambda !!!
        print >> stdout, '{"mino":%.9e, "tau":%.9e, "v4e":%.1f, "v4c":%.1f, "ER":%.1f, "ETh":%.1f, "t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "tP":%.9e, "rP":%.9e, "thP":%.9e, "phP":%.9e}' \
                % (mino, tau, v4e, -180.0, eR, eTh, self.t, self.r, self.th, self.ph, self.Ut / self.S, self.Ur / self.S, self.Uth / self.S, self.Uph / self.S)  # Log data,  d/dTau = 1/sigma * d/dLambda !!!

    def solve(self):
        mino = tau = 0.0
        iterationCount = plotCount = 0
        while tau <= self.end_time:
            if tau >= self.start_time:
                self.plot(mino, tau)
                plotCount += 1
            self.qUp(self.coefficients[0])  # w * c1
            self.pUp(self.coefficients[1])  # w * d1
            self.qUp(self.coefficients[2])  # w * c2
            self.pUp(self.coefficients[3])  # w * d2
            self.qUp(self.coefficients[2])  # w * c3
            self.pUp(self.coefficients[1])  # w * d3
            self.qUp(self.coefficients[0])  # w * c4
            iterationCount += 1
            mino = iterationCount * self.h
            tau += self.h * self.S  # dTau = sigma * dlambda  - NB lambda is affine parameter here, not the cc !!!
        self.plot(mino, tau)
        return iterationCount, plotCount


def main():
    print >> stderr, "Simulator: {}".format(argv[0])
    input_data = stdin.read()
    ic = loads(input_data)['IC']
    print >> stderr, input_data
    BhSymp(ic['lambda'], ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'],
           ic['start'], ic['end'], ic['step']).solve()


if __name__ == "__main__":
    main()
else:
    print >> stderr, __name__ + " module loaded"
