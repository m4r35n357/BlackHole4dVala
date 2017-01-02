#!/usr/bin/env pypy
from json import loads
from math import sqrt, sin, pi, fabs, log10
from sys import stdin, stderr, stdout, argv

class BhRk4(object):
    def __init__(self, Lambda, a, mu2, E, L, C, r0, thetaMin, start, end, time_step, plot_ratio):
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
        self.start = start
        self.end = end
        self.time_step = time_step
        self.plot_ratio = plot_ratio
        self.t = self.ph = 0.0
        self.r = r0
        self.th = (90.0 - thetaMin) * pi / 180.0
        self.kt = [0.0, 0.0, 0.0, 0.0]
        self.kr = [0.0, 0.0, 0.0, 0.0]
        self.kth = [0.0, 0.0, 0.0, 0.0]
        self.kph = [0.0, 0.0, 0.0, 0.0]
        self.sign_r = self.sign_th = -1
        self.evaluate(self.r, self.th, 0)

    def evaluate(self, radius, theta, stage):
        r2 = radius**2
        self.sth2 = sin(theta)**2
        cth2 = 1.0 - self.sth2
        self.ra2 = r2 + self.a2
        P = self.ra2 * self.E - self.aL
        self.D_r = (1.0 - self.l_3 * r2) * self.ra2 - 2.0 * radius
        self.R = self.X2 * P**2 - self.D_r * (self.mu2 * r2 + self.K)
        T = self.aE * self.sth2 - self.L
        self.D_th = 1.0 + self.a2l_3 * cth2
        self.TH = self.D_th * (self.K - self.a2mu2 * cth2) - self.X2 * T**2 / self.sth2
        P_Dr = P / self.D_r
        T_Dth = T / self.D_th
        self.S = r2 + self.a2 * cth2
        self.kt[stage] = self.time_step * (P_Dr * self.ra2 - T_Dth * self.a) * self.X2 / self.S
        self.kr[stage] = self.time_step * sqrt(self.R if self.R >= 0.0 else -self.R) / self.S
        self.kth[stage] = self.time_step * sqrt(self.TH if self.TH >= 0.0 else -self.TH) / self.S
        self.kph[stage] = self.time_step * (P_Dr * self.a - T_Dth / self.sth2) * self.X2 / self.S

    def solve(self):
        tau = 0.0
        iteration_count = 0
        while tau < self.end:
            if tau >= self.start and iteration_count % self.plot_ratio == 0:
                self.plot(tau, self.kt[0] / self.time_step, self.kr[0] / self.time_step, self.kth[0] / self.time_step, self.kph[0] / self.time_step)
            self.sign_r = self.sign_r if self.R > 0.0 else -self.sign_r
            self.sign_th = self.sign_th if self.TH > 0.0 else -self.sign_th
            self.evaluate(self.r + 0.5 * self.kr[0], self.th + 0.5 * self.kth[0], 1)
            self.evaluate(self.r + 0.5 * self.kr[1], self.th + 0.5 * self.kth[1], 2)
            self.evaluate(self.r + self.kr[2], self.th + self.kth[2], 3)
            self.t += (self.kt[0] + 2.0 * (self.kt[1] + self.kt[2]) + self.kt[3]) / 6.0
            self.r += (self.kr[0] + 2.0 * (self.kr[1] + self.kr[2]) + self.kr[3]) / 6.0 * self.sign_r
            self.th += (self.kth[0] + 2.0 * (self.kth[1] + self.kth[2]) + self.kth[3]) / 6.0 * self.sign_th
            self.ph += (self.kph[0] + 2.0 * (self.kph[1] + self.kph[2]) + self.kph[3]) / 6.0
            self.evaluate(self.r, self.th, 0)
            iteration_count += 1
            tau = iteration_count * self.time_step
        self.plot(tau, self.kt[0] / self.time_step, self.kr[0] / self.time_step, self.kth[0] / self.time_step, self.kph[0] / self.time_step)

    def plot(self, tau, Ut, Ur, Uth, Uph):
        SX2 = self.S * self.X2
        e = fabs(self.mu2 + self.sth2 * self.D_th / SX2 * (self.a * Ut - self.ra2 * Uph)**2
                    + self.S / self.D_r * Ur**2 + self.S / self.D_th * Uth**2
                    - self.D_r / SX2 * (Ut - self.a * self.sth2 * Uph)**2)
        print >> stdout, '{"tau":%.9e, "v4e":%.1f, "D_r":%.9e, "D_th":%.9e, "S":%.9e,' \
                         % (tau, 10.0 * log10(e if e > 1.0e-18 else 1.0e-18), self.D_r, self.D_th, self.S),
        print >> stdout, '"t":%.9e, "r":%.9e, "th":%.9e, "ph":%.9e, "tP":%.9e, "rP":%.9e, "thP":%.9e, "phP":%.9e}' \
                         % (self.t, self.r, self.th, self.ph, Ut, Ur, Uth, Uph)

if __name__ == "__main__":
    print >> stderr, "Executable: {}".format(argv[0])
    input_data = stdin.read()
    ic = loads(input_data)['IC']
    print >> stderr, input_data
    BhRk4(ic['lambda'], ic['a'], ic['mu'], ic['E'], ic['L'], ic['Q'], ic['r0'], ic['th0'],
          ic['start'], ic['end'], ic['step'], ic['plotratio']).solve()
else:
    print >> stderr, __name__ + " module loaded"
