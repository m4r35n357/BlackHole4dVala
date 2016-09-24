/*
Copyright (c) 2014, 2015, 2016, Ian Smith (m4r35n357)
All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
using GLib.Math;

namespace Simulations {

    protected class KerrBase : GLib.Object {
        /**
         * All fields are protected
         */
        // Constants from IC file
        protected double l_3;
        protected double a;
        protected double a2;
        protected double a2l_3;
        protected double mu2;
        protected double a2mu2;
        protected double X2;
        protected double E;
        protected double L;
        protected double K;
        protected double aE;
        protected double aL;
        protected double start;
        protected double end;
        protected double h;
        protected int64 tr;
        protected double horizon;
        // variables
        protected double sth2;
        protected double ra2;
        protected double D_r;
        protected double D_th;
        protected double S;
        protected double R;
        protected double TH;
        // state
        protected double tau = 0.0;
        protected double t = 0.0;
        protected double r;
        protected double th;
        protected double ph = 0.0;
        protected double Ut;
        protected double Ur;
        protected double Uth;
        protected double Uph;

        /**
         * Protected constructor, use the static factory in subclass
         */
        protected KerrBase (double lambda, double a, double mu2, double E, double L, double Q, double r0, double th0,
                    double tau0, double deltaTau, double tStep, int64 tRatio) {
            this.l_3 = lambda / 3.0;
            this.a = a;
            this.mu2 = mu2;
            this.E = E;
            this.L = L;
            this.a2 = a * a;
            this.a2l_3 = a2 * l_3;
            this.a2mu2 = a2 * mu2;
            this.horizon = 1.0 + sqrt(1.0 - a2);
            this.aE = a * E;
            this.aL = a * L;
            this.X2 = (1.0 + a2l_3) * (1.0 + a2l_3);
            this.K = Q + X2 * (L - aE) * (L - aE);
            this.start = tau0;
            this.end = tau0 + deltaTau;
            this.h = tStep;
            this.tr = tRatio;
            this.r = r0;
            this.th = (90.0 - th0) * PI / 180.0;
        }
    }
}
