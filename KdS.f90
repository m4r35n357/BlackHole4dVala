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
program KdS
    use Model
    implicit none
    real(16), parameter :: D0=0.0_16, D05=0.5_16, D1=1.0_16, D2=2.0_16, D3=3.0_16, D5=5.0_16, D7=7.0_16, D9=9.0_16, D11=11.0_16  ! CONSTANTS
    real(16) :: root, v1, v3, w1, w3, x1, x3, y1, y3, z1, z3, time = D0, h, start, finish
    integer :: plotratio, stages, outer
    character (len=3) :: integrator
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
        case ("b12")
            write (0, *) "12th Order Symplectic Integrator (using explicit composition)"
            call solve(twelfth_order)
        case default
            error stop "Invalid integrator method"
    end select
contains
    subroutine init_vars()
        read(*,*) h, start, finish, plotratio, integrator, stages
        root = stages - D1;
        outer = (stages - 1) / 2;
        v1 = D1 / (root - root**(D1 / D11))
        v3 = D1 - root * v1
        w1 = D1 / (root - root**(D1 / D9))
        w3 = D1 - root * w1
        x1 = D1 / (root - root**(D1 / D7))
        x3 = D1 - root * x1
        y1 = D1 / (root - root**(D1 / D5))
        y3 = D1 - root * y1
        z1 = D1 / (root - root**(D1 / D3))
        z3 = D1 - root * z1
        call init_model_vars()
    end subroutine init_vars

    subroutine solve(method)
        integer :: counter = 0
        do while ((tau < finish) .and. (cross .or. Dr > D0))
            if ((tau >= start) .and. (mod(counter, plotratio) == 0)) then
                call plotModel(time)
            end if
            call method()
            counter = counter + 1
            time = h * counter
            call postLoop(h)
        end do
        call plotModel(time)
    end subroutine solve

    subroutine first_order()
        call qUpdate(h)
        call pUpdate(h)
    end subroutine first_order

    subroutine base2(s)
        real(16), intent(in) :: s
        call qUpdate(h * s * D05)
        call pUpdate(h * s)
        call qUpdate(h * s * D05)
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

    subroutine base12(s)
        real(16), intent(in) :: s
        integer m
        do m = 1, outer
            call base10(s * v1)
        end do
        call base10(s * v3)
        do m = 1, outer
            call base10(s * v1)
        end do
    end subroutine base12

    subroutine twelfth_order()
        call base12(D1)
    end subroutine twelfth_order
end program KdS

