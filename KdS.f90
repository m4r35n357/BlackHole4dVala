!Copyright (c) 2014, 2015, 2016, 2017 Ian Smith (m4r35n357)
!All rights reserved.
!
!Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
!
!1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
!
!2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
!
!3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
!
!THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

program KdS
    implicit none
    real(kind=16) :: l_3, a, a2, a2l_3, mu2, a2mu2, X2, E, L, CC_K, aE, two_EX2, two_aE, aL, step, start, finish, pi = acos(-1.0)
    integer :: plotratio
    logical :: cross
    real(kind=16) :: r2, ra2, sth, cth, sth2, cth2, D_r, D_th, SIGMA, R_POT, R_tmp, TH_tmp, TH_POT  ! INTERMEDIATE VARIABLES
    real(kind=16) :: t = 0.0, r, th, ph = 0.0, Ut, Ur, Uth, Uph  ! VARIABLES
    call init_vars()
contains
    subroutine init_vars()
        real(kind=16) :: lambda, spin, pMass2, energy, angMom, CC_Q, r0, th0
        read(*,*) lambda, spin, pMass2, energy, angMom, CC_Q, r0, th0, step, start, finish, plotratio
        l_3 = lambda / 3.0
        a = spin
        mu2 = pMass2
        E = energy
        L = angMom
        a2 = a**2
        a2l_3 = a2 * l_3
        a2mu2 = a2 * mu2
        aE = a * E
        aL = a * L
        X2 = (1.0 + a2l_3)**2
        two_EX2 = 2.0 * E * X2
        two_aE = 2.0 * aE
        CC_K = CC_Q + X2 * (L - aE)**2
        r = r0
        th = (90.0 - th0) * pi / 180.0
        cross = .false.
        call solve()
    end subroutine init_vars

    subroutine refresh()
        real(kind=16) :: P_Dr, T_Dth
        r2 = r**2
        ra2 = r2 + a2
        R_tmp = ra2 * E - aL
        D_r = (1.0 - l_3 * r2) * ra2 - 2.0 * r
        R_POT = X2 * R_tmp**2 - D_r * (mu2 * r2 + CC_K)
        sth = sin(th)
        cth = cos(th)
        sth2 = sth**2
        cth2 = 1.0 - sth2
        TH_tmp = aE * sth2 - L
        D_th = 1.0 + a2l_3 * cth2
        TH_POT = D_th * (CC_K - a2mu2 * cth2) - X2 * TH_tmp**2 / sth2
        P_Dr = R_tmp / D_r
        T_Dth = TH_tmp / D_th
        Ut = (P_Dr * ra2 - T_Dth * a) * X2
        Uph = (P_Dr * a - T_Dth / sth2) * X2
        SIGMA = r2 + a2 * cth2
    end subroutine refresh

    subroutine qUpdate(c)
        real(kind=16) :: c
        t = t + c * Ut
        r = r + c * Ur
        th = th + c * Uth
        ph = ph + c * Uph
        call refresh()
    end subroutine qUpdate

    subroutine pUpdate(d)
        real(kind=16) :: d
        Ur = Ur + d * (r * (two_EX2 * R_tmp - mu2 * D_r) - (r * (1.0 - l_3 * (r2 + ra2)) - 1.0) * (CC_K + mu2 * r2))
        Uth = Uth + d * cth * (sth * a2 * (mu2 * D_th - l_3 * (CC_K - a2mu2 * cth2)) + X2 * TH_tmp / sth * (TH_tmp / sth2 - two_aE))
    end subroutine pUpdate

    subroutine solve()
        real(kind=16) :: mino = 0.0, tau = 0.0, theta
        real(kind=16), dimension(4) :: c_d
        integer :: counter = 0
        theta = 1.0 / (2.0 - 2.0**(1.0 / 3.0));
        c_d = (/ 0.5 * step * theta, step * theta, 0.5 * step * (1.0 - theta), step * (1.0 - 2.0 * theta) /)
        call refresh()
        Ur = - sqrt(merge(R_POT, -R_POT, R_POT >= 0.0))
        Uth = - sqrt(merge(TH_POT, -TH_POT, TH_POT >= 0.0))
        do while ((tau < finish) .and. (cross .or. D_r > 0.0))
            if ((tau >= start) .and. (mod(counter, plotratio) == 0)) then
                call plot(mino, tau, Ut / SIGMA, Ur / SIGMA, Uth / SIGMA, Uph / SIGMA)
            end if
            call qUpdate(c_d(1))
            call pUpdate(c_d(2))
            call qUpdate(c_d(3))
            call pUpdate(c_d(4))
            call qUpdate(c_d(3))
            call pUpdate(c_d(2))
            call qUpdate(c_d(1))
            counter = counter + 1
            mino = step * counter
            tau = tau + step * SIGMA
        end do
        call plot(mino, tau, Ut / SIGMA, Ur / SIGMA, Uth / SIGMA, Uph / SIGMA)
    end subroutine solve

    subroutine plot(mino, tau, Vt, Vr, Vth, Vph)
        real(kind=16) :: mino, tau, Vt, Vr, Vth, Vph
        write (*, '(A, 13(ES16.9, A))') '{"mino":', mino, ',"tau":', tau,&
                    ',"v4e":',mu2 + sth2 * D_th / (SIGMA * X2) * (a * Vt - ra2 * Vph)**2 + SIGMA / D_r * Vr**2&
                                  + SIGMA / D_th * Vth**2 - D_r / (SIGMA * X2) * (Vt - a * sth2 * Vph)**2,&
                    ',"ER":', 0.5 * (Vr**2 - R_POT / SIGMA**2), ',"ETh":', 0.5 * (Vth**2 - TH_POT / SIGMA**2),&
                    ',"t":', t, ',"r":', r, ',"th":', th, ',"ph":', ph,&
                    ',"tP":', Vt, ',"rP":', Vr, ',"thP":', Vth, ',"phP":', Vph, '}'
    end subroutine plot
end program KdS

