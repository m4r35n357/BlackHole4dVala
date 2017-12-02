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
    real(16), parameter :: D0=0.0_16, D05=0.5_16, D1=1.0_16, D2=2.0_16, D3=3.0_16, D5=5.0_16, D7=7.0_16, D9=9.0_16  ! CONSTANTS
    real(16) :: l_3, a, a2, a2l_3, mu2, a2mu2, X2, E, L, ccK, aE, twoEX2, twoEa, aL, step, start, finish  ! IMMUTABLES
    real(16) :: root, w1, w3, x1, x3, y1, y3, z1, z3
    integer :: plotratio, stages, outer
    logical :: cross
    character (len=3) :: integrator
    real(16) :: r2, ra2, sth, cth, sth2, cth2, Dr, Dth, Sigma, Rpot, Rint, THint, THpot  ! INTERMEDIATE VARIABLES
    real(16) :: t = D0, r, th, ph = D0, Ut, Ur, Uth, Uph, mino = D0, tau = D0  ! PARTICLE VARIABLES
    character(len=32) :: arg
    call get_command_argument(0, arg)
    write (0, *) "Executable: ", trim(arg)
    call init_vars()
    if ((stages < 3) .or. (mod(stages, 2) == 0)) then
        error stop "'stages' should be odd and at least 3"
    end if
    select case (integrator)
        case ("b1")
            write (0, *) "1st Order Symplectic Integrator"
            call solve(first_order)
        case ("b2")
            write (0, *) "2nd Order Symplectic Integrator"
            call solve(second_order)
        case ("b4")
            write (0, *) "4th Order Symplectic Integrator (using explicit composition)"
            call solve(fourth_order)
        case ("b6")
            write (0, *) "6th Order Symplectic Integrator (using explicit composition)"
            call solve(sixth_order)
        case ("b8")
            write (0, *) "8th Order Symplectic Integrator (using explicit composition)"
            call solve(eightth_order)
        case ("b10")
            write (0, *) "10th Order Symplectic Integrator (using explicit composition)"
            call solve(tenth_order)
        case default
            error stop "Invalid integrator method"
    end select
contains
    subroutine init_vars()
        real(16) :: lambda, spin, pMass2, energy, angMom, ccQ, r0, th0
        write (0, *) "Kerr-deSitter Geodesic"
        read(*,*) lambda, spin, pMass2, energy, angMom, ccQ, r0, th0, step, start, finish, plotratio, cross, integrator, stages
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
        twoEX2 = D2 * E * X2
        twoEa = D2 * aE
        ccK = ccQ + X2 * (L - aE)**2
        r = r0
        th = (90.0_16 - th0) * acos(-D1) / 180.0_16
        root = stages - D1;
        outer = (stages - 1) / 2;
        w1 = D1 / (root - root**(D1 / D9))
        w3 = D1 - root * w1
        x1 = D1 / (root - root**(D1 / D7))
        x3 = D1 - root * x1
        y1 = D1 / (root - root**(D1 / D5))
        y3 = D1 - root * y1
        z1 = D1 / (root - root**(D1 / D3))
        z3 = D1 - root * z1
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
        real(16), intent(in) :: c
        t = t + c * Ut
        r = r + c * Ur
        th = th + c * Uth
        ph = ph + c * Uph
        call refresh()
    end subroutine qUpdate

    subroutine pUpdate(d)
        real(16), intent(in) :: d
        Ur = Ur + d * (r * (twoEX2 * Rint - mu2 * Dr) - (r * (D1 - l_3 * (r2 + ra2)) - D1) * (ccK + mu2 * r2))
        Uth = Uth + d * cth * (sth * a2 * (mu2 * Dth - l_3 * (ccK - a2mu2 * cth2)) + X2 * THint / sth * (THint / sth2 - twoEa))
    end subroutine pUpdate

    subroutine solve(method)
        integer :: counter = 0
        call refresh()
        Ur = - sqrt(merge(Rpot, -Rpot, Rpot >= D0))
        Uth = - sqrt(merge(THpot, -THpot, THpot >= D0))
        do while ((tau < finish) .and. (cross .or. Dr > D0))
            if ((tau >= start) .and. (mod(counter, plotratio) == 0)) then
                call plot(Ut / Sigma, Ur / Sigma, Uth / Sigma, Uph / Sigma)
            end if
            call method()
            counter = counter + 1
            mino = step * counter
            tau = tau + step * Sigma
        end do
        call plot(Ut / Sigma, Ur / Sigma, Uth / Sigma, Uph / Sigma)
    end subroutine solve

    subroutine plot(Vt, Vr, Vth, Vph)
        real(16), intent(in) :: Vt, Vr, Vth, Vph
        write (*, '(A, 13(ES16.9, A))') '{"mino":',mino,',"tau":',tau,&
                    ',"v4e":',mu2 + sth2 * Dth / (Sigma * X2) * (a * Vt - ra2 * Vph)**2 + Sigma / Dr * Vr**2&
                                  + Sigma / Dth * Vth**2 - Dr / (Sigma * X2) * (Vt - a * sth2 * Vph)**2,&
                    ',"ER":',Vr**2 - Rpot / Sigma**2,',"ETh":',Vth**2 - THpot / Sigma**2,&
                    ',"t":', t,',"r":',r,',"th":',th,',"ph":',ph,',"tP":',Vt,',"rP":',Vr,',"thP":',Vth,',"phP":',Vph,'}'
    end subroutine plot

    subroutine first_order()
        call qUpdate(step)
        call pUpdate(step)
    end subroutine first_order

    subroutine base2(s)
        real(16), intent(in) :: s
        call qUpdate(step * s * D05)
        call pUpdate(step * s)
        call qUpdate(step * s * D05)
    end subroutine base2

    subroutine second_order()
        call base2(D1)
    end subroutine second_order

    subroutine base4(s)
        real(16), intent(in) :: s
        integer i
        do i = 1, outer
            call base2(s * z1)
        end do
        call base2(s * z3)
        do i = 1, outer
            call base2(s * z1)
        end do
    end subroutine base4

    subroutine fourth_order()
        call base4(D1)
    end subroutine fourth_order

    subroutine base6(s)
        real(16), intent(in) :: s
        integer j
        do j = 1, outer
            call base4(s * y1)
        end do
        call base4(s * y3)
        do j = 1, outer
            call base4(s * y1)
        end do
    end subroutine base6

    subroutine sixth_order()
        call base6(D1)
    end subroutine sixth_order

    subroutine base8(s)
        real(16), intent(in) :: s
        integer k
        do k = 1, outer
            call base6(s * x1)
        end do
        call base6(s * x3)
        do k = 1, outer
            call base6(s * x1)
        end do
    end subroutine base8

    subroutine eightth_order()
        call base8(D1)
    end subroutine eightth_order

    subroutine base10(s)
        real(16), intent(in) :: s
        integer l
        do l = 1, outer
            call base8(s * w1)
        end do
        call base8(s * w3)
        do l = 1, outer
            call base8(s * w1)
        end do
    end subroutine base10

    subroutine tenth_order()
        call base10(D1)
    end subroutine tenth_order
end program KdS

