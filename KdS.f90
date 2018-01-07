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
    real(16) :: root, v_outer, v_mid, w_outer, w_mid, x_outer, x_mid, y_outer, y_mid, z_outer, z_mid, time = D0, step, start, finish
    integer :: plot_ratio, stages, outer
    character (len=3) :: integrator
    character(len=32) :: arg
    call get_command_argument(0, arg)
    write (0, *) "Executable: ", trim(arg)
    call init_integrator_vars()
    if ((stages < 3) .or. (mod(stages, 2) == 0)) then
        error stop "'stages' should be odd and at least 3"
    end if
    select case (integrator)
        case ("b1")
            write (0, *) "1st Order Symplectic Integrator"
            call evolve(first_order)
        case ("b2")
            write (0, *) "2nd Order Symplectic Integrator"
            call evolve(second_order)
        case ("b4")
            write (0, *) "4th Order Symplectic Integrator (using explicit composition)"
            call evolve(fourth_order)
        case ("b6")
            write (0, *) "6th Order Symplectic Integrator (using explicit composition)"
            call evolve(sixth_order)
        case ("b8")
            write (0, *) "8th Order Symplectic Integrator (using explicit composition)"
            call evolve(eightth_order)
        case ("b10")
            write (0, *) "10th Order Symplectic Integrator (using explicit composition)"
            call evolve(tenth_order)
        case ("b12")
            write (0, *) "12th Order Symplectic Integrator (using explicit composition)"
            call evolve(twelfth_order)
        case default
            error stop "Invalid integrator method"
    end select
contains
    subroutine init_integrator_vars()
        read(*,*) step, start, finish, plot_ratio, integrator, stages
        root = stages - D1;
        outer = (stages - 1) / 2;
        v_outer = D1 / (root - root**(D1 / D11))
        v_mid = D1 - root * v_outer
        w_outer = D1 / (root - root**(D1 / D9))
        w_mid = D1 - root * w_outer
        x_outer = D1 / (root - root**(D1 / D7))
        x_mid = D1 - root * x_outer
        y_outer = D1 / (root - root**(D1 / D5))
        y_mid = D1 - root * y_outer
        z_outer = D1 / (root - root**(D1 / D3))
        z_mid = D1 - root * z_outer
        call init_model_vars()
    end subroutine init_integrator_vars

    subroutine evolve(nth_order)
        integer :: counter = 0
        do while ((time < finish) .and. carry_on)
            if ((time >= start) .and. (mod(counter, plot_ratio) == 0)) then
                call plot(time)
            end if
            call nth_order()
            counter = counter + 1
            time = t_update(time, step, counter)
        end do
        call plot(time)
    end subroutine evolve

    subroutine first_order()
        call q_update(step)
        call p_update(step)
    end subroutine first_order

    subroutine base_2(s)
        real(16), intent(in) :: s
        call q_update(step * s * D05)
        call p_update(step * s)
        call q_update(step * s * D05)
    end subroutine base_2

    subroutine second_order()
        call base_2(D1)
    end subroutine second_order

    subroutine base_4(s)
        real(16), intent(in) :: s
        integer i
        do i = 1, outer
            call base_2(s * z_outer)
        end do
        call base_2(s * z_mid)
        do i = 1, outer
            call base_2(s * z_outer)
        end do
    end subroutine base_4

    subroutine fourth_order()
        call base_4(D1)
    end subroutine fourth_order

    subroutine base_6(s)
        real(16), intent(in) :: s
        integer j
        do j = 1, outer
            call base_4(s * y_outer)
        end do
        call base_4(s * y_mid)
        do j = 1, outer
            call base_4(s * y_outer)
        end do
    end subroutine base_6

    subroutine sixth_order()
        call base_6(D1)
    end subroutine sixth_order

    subroutine base_8(s)
        real(16), intent(in) :: s
        integer k
        do k = 1, outer
            call base_6(s * x_outer)
        end do
        call base_6(s * x_mid)
        do k = 1, outer
            call base_6(s * x_outer)
        end do
    end subroutine base_8

    subroutine eightth_order()
        call base_8(D1)
    end subroutine eightth_order

    subroutine base_10(s)
        real(16), intent(in) :: s
        integer l
        do l = 1, outer
            call base_8(s * w_outer)
        end do
        call base_8(s * w_mid)
        do l = 1, outer
            call base_8(s * w_outer)
        end do
    end subroutine base_10

    subroutine tenth_order()
        call base_10(D1)
    end subroutine tenth_order

    subroutine base_12(s)
        real(16), intent(in) :: s
        integer m
        do m = 1, outer
            call base_10(s * v_outer)
        end do
        call base_10(s * v_mid)
        do m = 1, outer
            call base_10(s * v_outer)
        end do
    end subroutine base_12

    subroutine twelfth_order()
        call base_12(D1)
    end subroutine twelfth_order
end program KdS

