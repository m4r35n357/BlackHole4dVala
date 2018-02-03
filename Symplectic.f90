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
program Symplectic
    use iso_fortran_env
    use Model
    implicit none
    double precision, parameter :: D0=0.0, D05=0.5, D1=1.0, D4=4.0, D3=3.0, D5=5.0, D7=7.0, D9=9.0, D11=11.0
    double precision :: v_fwd, v_back, w_fwd, w_back, x_fwd, x_back, y_fwd, y_back, z_fwd, z_back, step, start, finish
    integer :: plot_ratio
    character (len=3) :: integrator
    character(len=32) :: arg
    call get_command_argument(0, arg)
    write (error_unit, *) "Executable: ", trim(arg)
    write (error_unit, *) 'double precision is:', precision(D0), ' decimal places'
    v_fwd = D1 / (D4 - D4**(D1 / D11))
    w_fwd = D1 / (D4 - D4**(D1 / D9))
    x_fwd = D1 / (D4 - D4**(D1 / D7))
    y_fwd = D1 / (D4 - D4**(D1 / D5))
    z_fwd = D1 / (D4 - D4**(D1 / D3))
    v_back = D1 - D4 * v_fwd
    w_back = D1 - D4 * w_fwd
    x_back = D1 - D4 * x_fwd
    y_back = D1 - D4 * y_fwd
    z_back = D1 - D4 * z_fwd
    read (input_unit, *) step, start, finish, plot_ratio, integrator
    call init_model()
    select case (integrator)
        case ("b2")
            write (error_unit, *) "2nd Order Integrator Base (Stormer-Verlet)"
            call evolve(second_order)
        case ("b4")
            write (error_unit, *) "4th Order Integrator (using Suzuki composition)"
            call evolve(fourth_order)
        case ("b6")
            write (error_unit, *) "6th Order Integrator (using Suzuki composition)"
            call evolve(sixth_order)
        case ("b8")
            write (error_unit, *) "8th Order Integrator (using Suzuki composition)"
            call evolve(eightth_order)
        case ("b10")
            write (error_unit, *) "10th Order Integrator (using Suzuki composition)"
            call evolve(tenth_order)
        case ("b12")
            write (error_unit, *) "12th Order Integrator (using Suzuki composition)"
            call evolve(twelfth_order)
        case default
            error stop "Invalid integrator method"
    end select
contains
    subroutine evolve (nth_order)
        double precision :: time = D0
        integer :: counter = 0
        do while ((time < finish) .and. carry_on)
            if ((time >= start) .and. (mod(counter, plot_ratio) == 0)) then
                call plot(time)
            end if
            call nth_order()
            counter = counter + 1
            time = t_update(time, step, counter)
        end do
        !call plot(time)
    end subroutine evolve

    subroutine suzuki (base, s, forward, back)
        double precision, intent(in) :: s, forward, back
        call base(s * forward)
        call base(s * forward)
        call base(s * back)
        call base(s * forward)
        call base(s * forward)
    end subroutine suzuki

    subroutine base_2 (s)
        double precision, intent(in) :: s
        call q_update(step * s * D05)
        call p_update(step * s)
        call q_update(step * s * D05)
    end subroutine base_2

    subroutine second_order ()
        call base_2(D1)
    end subroutine second_order

    subroutine base_4 (s)
        double precision, intent(in) :: s
        call suzuki(base_2, s, z_fwd, z_back)
    end subroutine base_4

    subroutine fourth_order ()
        call base_4(D1)
    end subroutine fourth_order

    subroutine base_6 (s)
        double precision, intent(in) :: s
        call suzuki(base_4, s, y_fwd, y_back)
    end subroutine base_6

    subroutine sixth_order ()
        call base_6(D1)
    end subroutine sixth_order

    subroutine base_8 (s)
        double precision, intent(in) :: s
        call suzuki(base_6, s, x_fwd, x_back)
    end subroutine base_8

    subroutine eightth_order ()
        call base_8(D1)
    end subroutine eightth_order

    subroutine base_10 (s)
        double precision, intent(in) :: s
        call suzuki(base_8, s, w_fwd, w_back)
    end subroutine base_10

    subroutine tenth_order ()
        call base_10(D1)
    end subroutine tenth_order

    subroutine base_12 (s)
        double precision, intent(in) :: s
        call suzuki(base_10, s, v_fwd, v_back)
    end subroutine base_12

    subroutine twelfth_order ()
        call base_12(D1)
    end subroutine twelfth_order
end program Symplectic
