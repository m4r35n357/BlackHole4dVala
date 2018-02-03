!Copyright (c) 2014-2018 Ian Smith (m4r35n357)
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
module Model
    use iso_fortran_env
    implicit none
    double precision, parameter :: MD0=0.0, MD1=1.0, MD2=2.0, MD3=3.0  ! CONSTANTS
    double precision :: l_3, a, a2, a2l_3, mu2, a2mu2, X2, E, L, ccK, aE, twoEX2, twoEa, aL  ! IMMUTABLES
    double precision :: r2, ra2, sth, cth, sth2, cth2, Dr, Dth, Sigma, Rpot, Rint, THint, THpot  ! INTERMEDIATE VARIABLES
    double precision :: mino, t=MD0, r, th, ph=MD0, Ut, Ur, Uth, Uph  ! PARTICLE VARIABLES (proper time, coordinates, and velocities))
    logical :: cross, carry_on = .true.
contains
    subroutine init_model()
        double precision :: lambda, ccQ, th0
        write (error_unit, *) "Kerr-deSitter Geodesic"
        read (input_unit, *) cross, lambda, a, mu2, E, L, ccQ, r, th0
        l_3 = lambda / MD3
        a2 = a**2
        a2l_3 = a2 * l_3
        a2mu2 = a2 * mu2
        aE = a * E
        aL = a * L
        X2 = (MD1 + a2l_3)**2
        twoEX2 = MD2 * E * X2
        twoEa = MD2 * aE
        ccK = ccQ + X2 * (L - aE)**2
        th = (90.0 - th0) * acos(-MD1) / 180.0
        call refresh()
        Ur = - sqrt(merge(Rpot, -Rpot, Rpot >= MD0))
        Uth = - sqrt(merge(THpot, -THpot, THpot >= MD0))
    end subroutine init_model

    subroutine refresh()
        double precision :: P_Dr, T_Dth
        r2 = r**2
        ra2 = r2 + a2
        Rint = ra2 * E - aL
        Dr = (MD1 - l_3 * r2) * ra2 - MD2 * r
        Rpot = X2 * Rint**2 - Dr * (mu2 * r2 + ccK)
        sth = sin(th)
        cth = cos(th)
        sth2 = sth**2
        cth2 = MD1 - sth2
        THint = aE * sth2 - L
        Dth = MD1 + a2l_3 * cth2
        THpot = Dth * (ccK - a2mu2 * cth2) - X2 * THint**2 / sth2
        P_Dr = Rint / Dr
        T_Dth = THint / Dth
        Ut = (P_Dr * ra2 - T_Dth * a) * X2
        Uph = (P_Dr * a - T_Dth / sth2) * X2
        Sigma = r2 + a2 * cth2
    end subroutine refresh

    subroutine q_update(c)
        double precision, intent(in) :: c
        t = t + c * Ut
        r = r + c * Ur
        th = th + c * Uth
        ph = ph + c * Uph
        call refresh()
    end subroutine q_update

    subroutine p_update(d)
        double precision, intent(in) :: d
        Ur = Ur + d * (r * (twoEX2 * Rint - mu2 * Dr) - (r * (MD1 - l_3 * (r2 + ra2)) - MD1) * (ccK + mu2 * r2))
        Uth = Uth + d * cth * (sth * a2 * (mu2 * Dth - l_3 * (ccK - a2mu2 * cth2)) + X2 * THint / sth * (THint / sth2 - twoEa))
    end subroutine p_update

    double precision function t_update(tau, step, counter)
        double precision, intent(in) :: tau, step
        integer, intent(in) :: counter
        carry_on = cross .or. Dr > MD0
        mino = step * counter
        t_update = tau + step * Sigma
    end function t_update

    subroutine plot(tau)
        double precision, intent(in) :: tau
        double precision :: Vt, Vr, Vth, Vph
        Vt = Ut / Sigma
        Vr = Ur / Sigma
        Vth = Uth / Sigma
        Vph = Uph / Sigma
        write (output_unit, '(A, 13(ES16.9, A))') '{"mino":',mino,',"tau":',tau,&
                    ',"v4e":',mu2 + sth2 * Dth / (Sigma * X2) * (a * Vt - ra2 * Vph)**2 + Sigma / Dr * Vr**2&
                                  + Sigma / Dth * Vth**2 - Dr / (Sigma * X2) * (Vt - a * sth2 * Vph)**2,&
                    ',"ER":',Vr**2 - Rpot / Sigma**2,',"ETh":',Vth**2 - THpot / Sigma**2,&
                    ',"t":', t,',"r":',r,',"th":',th,',"ph":',ph,',"tP":',Vt,',"rP":',Vr,',"thP":',Vth,',"phP":',Vph,'}'
    end subroutine plot
end module Model

