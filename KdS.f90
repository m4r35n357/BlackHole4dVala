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
    integer :: plotratio ! CONSTANTS
    logical :: cross = .false.
    real(16), parameter :: D0 = 0.0, D05 = 0.5, D1 = 1.0, D2 = 2.0, D3 = 3.0
    real(16) :: l_3, a, a2, a2l_3, mu2, a2mu2, X2, E, L, ccK, aE, two_EX2, two_aE, aL, step, start, finish, pi = acos(-D1)
    real(16) :: r2, ra2, sth, cth, sth2, cth2, Dr, Dth, Sigma, Rpot, Rint, THint, THpot  ! INTERMEDIATE VARIABLES
    real(16) :: t = D0, r, th, ph = D0, Ut, Ur, Uth, Uph  ! VARIABLES
    call init_vars()
contains
    subroutine init_vars()
        real(16) :: lambda, spin, pMass2, energy, angMom, ccQ, r0, th0
        read(*,*) lambda, spin, pMass2, energy, angMom, ccQ, r0, th0, step, start, finish, plotratio
        l_3 = lambda / D3
        a = spin
        mu2 = pMass2
        E = energy
        L = angMom
        a2 = a**2
        a2l_3 = a2 * l_3
        a2mu2 = a2 * mu2
        aE = a * E
        aL = a * L
        X2 = (D1 + a2l_3)**2
        two_EX2 = D2 * E * X2
        two_aE = D2 * aE
        ccK = ccQ + X2 * (L - aE)**2
        r = r0
        th = (90.0 - th0) * pi / 180.0
        call solve()
    end subroutine init_vars

    subroutine refresh()
        real(16) :: P_Dr, T_Dth
        r2 = r**2
        ra2 = r2 + a2
        Rint = ra2 * E - aL
        Dr = (D1 - l_3 * r2) * ra2 - D2 * r
        Rpot = X2 * Rint**2 - Dr * (mu2 * r2 + ccK)
        sth = sin(th)
        cth = cos(th)
        sth2 = sth**2
        cth2 = D1 - sth2
        THint = aE * sth2 - L
        Dth = D1 + a2l_3 * cth2
        THpot = Dth * (ccK - a2mu2 * cth2) - X2 * THint**2 / sth2
        P_Dr = Rint / Dr
        T_Dth = THint / Dth
        Ut = (P_Dr * ra2 - T_Dth * a) * X2
        Uph = (P_Dr * a - T_Dth / sth2) * X2
        Sigma = r2 + a2 * cth2
    end subroutine refresh

    subroutine qUpdate(c)
        real(16) :: c
        t = t + c * Ut
        r = r + c * Ur
        th = th + c * Uth
        ph = ph + c * Uph
        call refresh()
    end subroutine qUpdate

    subroutine pUpdate(d)
        real(16) :: d
        Ur = Ur + d * (r * (two_EX2 * Rint - mu2 * Dr) - (r * (D1 - l_3 * (r2 + ra2)) - D1) * (ccK + mu2 * r2))
        Uth = Uth + d * cth * (sth * a2 * (mu2 * Dth - l_3 * (ccK - a2mu2 * cth2)) + X2 * THint / sth * (THint / sth2 - two_aE))
    end subroutine pUpdate

    subroutine solve()
        real(16) :: mino = D0, tau = D0, theta = D1 / (D2 - D2**(D1 / D3))
        real(16), dimension(4) :: cd
        integer :: counter = 0
        cd = (/ D05 * step * theta, step * theta, D05 * step * (D1 - theta), step * (D1 - D2 * theta) /)
        call refresh()
        Ur = - sqrt(merge(Rpot, -Rpot, Rpot >= D0))
        Uth = - sqrt(merge(THpot, -THpot, THpot >= D0))
        do while ((tau < finish) .and. (cross .or. Dr > D0))
            if ((tau >= start) .and. (mod(counter, plotratio) == 0)) then
                call plot(mino, tau, Ut / Sigma, Ur / Sigma, Uth / Sigma, Uph / Sigma)
            end if
            call qUpdate(cd(1))
            call pUpdate(cd(2))
            call qUpdate(cd(3))
            call pUpdate(cd(4))
            call qUpdate(cd(3))
            call pUpdate(cd(2))
            call qUpdate(cd(1))
            counter = counter + 1
            mino = step * counter
            tau = tau + step * Sigma
        end do
        call plot(mino, tau, Ut / Sigma, Ur / Sigma, Uth / Sigma, Uph / Sigma)
    end subroutine solve

    subroutine plot(mino, tau, Vt, Vr, Vth, Vph)
        real(16) :: mino, tau, Vt, Vr, Vth, Vph
        write (*, '(A, 13(ES16.9, A))') '{"mino":', mino, ',"tau":', tau,&
                    ',"v4e":',mu2 + sth2 * Dth / (Sigma * X2) * (a * Vt - ra2 * Vph)**2 + Sigma / Dr * Vr**2&
                                  + Sigma / Dth * Vth**2 - Dr / (Sigma * X2) * (Vt - a * sth2 * Vph)**2,&
                    ',"ER":', D05 * (Vr**2 - Rpot / Sigma**2), ',"ETh":', D05 * (Vth**2 - THpot / Sigma**2),&
                    ',"t":', t, ',"r":', r, ',"th":', th, ',"ph":', ph,&
                    ',"tP":', Vt, ',"rP":', Vr, ',"thP":', Vth, ',"phP":', Vph, '}'
    end subroutine plot
end program KdS

