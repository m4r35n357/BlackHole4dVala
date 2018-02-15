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
    use Model
    implicit none
    real(16), parameter :: D0=0.0, D05=0.5, D1=1.0, D4=4.0, D3=3.0, D5=5.0, D7=7.0, D9=9.0
    real(16) :: w_fwd, w_back, x_fwd, x_back, y_fwd, y_back, z_fwd, z_back, step, start, finish
    integer :: plot_ratio
    character (len=3) :: integrator
    character(len=32) :: arg
    call get_command_argument(0, arg)
    write (error_unit, *) "Executable: ", trim(arg)
    write (error_unit, *) 'double precision is:', precision(D0), ' decimal places'
    read (input_unit, *) step, start, finish, plot_ratio, integrator
    call init_model()
    w_fwd = D1 / (D4 - D4 ** (D1 / D9))
    x_fwd = D1 / (D4 - D4 ** (D1 / D7))
    y_fwd = D1 / (D4 - D4 ** (D1 / D5))
    z_fwd = D1 / (D4 - D4 ** (D1 / D3))
    w_back = D1 - D4 * w_fwd
    x_back = D1 - D4 * x_fwd
    y_back = D1 - D4 * y_fwd
    z_back = D1 - D4 * z_fwd
    select case (integrator)
        case ("b2")
            write (error_unit, *) "2nd Order Base (Stormer-Verlet)"
            call evolve(second_order_integrator)
        case ("b4")
            write (error_unit, *) "4th Order (Suzuki composition)"
            call evolve(fourth_order_integrator)
        case ("b6")
            write (error_unit, *) "6th Order (Suzuki composition)"
            call evolve(sixth_order_integrator)
        case ("b8")
            write (error_unit, *) "8th Order (Suzuki composition)"
            call evolve(eightth_order_integrator)
        case ("b10")
            write (error_unit, *) "10th Order (Suzuki composition)"
            call evolve(tenth_order_integrator)
        case default
            error stop "Invalid integrator method"
    end select
contains
    subroutine evolve (nth_order_integrator)
        real(16) :: time = D0
        integer(16) :: counter = 0
        do
            if ((time >= start) .and. (mod(counter, plot_ratio) == 0)) then
                call plot(time)
            end if
            call nth_order_integrator()
            counter = counter + 1
            time = t_update(time, step, counter)
            if (.not. (time < finish) .or. .not. carry_on) exit
        end do
    end subroutine evolve

    subroutine compose_suzuki (base_method, s, forward, back)
        real(16), intent(in) :: s, forward, back
        call base_method(s * forward)
        call base_method(s * forward)
        call base_method(s * back)
        call base_method(s * forward)
        call base_method(s * forward)
    end subroutine compose_suzuki

    subroutine base_2 (s)
        real(16), intent(in) :: s
        call q_update(s * step * D05)
        call p_update(s * step)
        call q_update(s * step * D05)
    end subroutine base_2

    subroutine second_order_integrator ()
        call base_2(D1)
    end subroutine second_order_integrator

    subroutine base_4 (s)
        real(16), intent(in) :: s
        call compose_suzuki(base_2, s, z_fwd, z_back)
    end subroutine base_4

    subroutine fourth_order_integrator ()
        call base_4(D1)
    end subroutine fourth_order_integrator

    subroutine base_6 (s)
        real(16), intent(in) :: s
        call compose_suzuki(base_4, s, y_fwd, y_back)
    end subroutine base_6

    subroutine sixth_order_integrator ()
        call base_6(D1)
    end subroutine sixth_order_integrator

    subroutine base_8 (s)
        real(16), intent(in) :: s
        call compose_suzuki(base_6, s, x_fwd, x_back)
    end subroutine base_8

    subroutine eightth_order_integrator ()
        call base_8(D1)
    end subroutine eightth_order_integrator

    subroutine tenth_order_integrator ()
        call compose_suzuki(base_8, D1, w_fwd, w_back)
    end subroutine tenth_order_integrator
end program Symplectic
