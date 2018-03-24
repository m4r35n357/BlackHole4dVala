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
    real(kind=16), parameter :: D0 = 0.0, D05 = 0.5, D1 = 1.0, D4 = 4.0, D3 = 3.0, D5 = 5.0, D7 = 7.0, D9 =9.0
    real(kind=16) :: w1, w0, x1, x0, y1, y0, z1, z0, h, start, finish
    real(kind=16), dimension(2) :: cd_2
    real(kind=16), dimension(6) :: cd_4
    real(kind=16), dimension(26) :: cd_6
    integer :: plot_ratio
    character (len=3) :: integrator
    character(len=32) :: arg
    call get_command_argument(0, arg)
    write (error_unit, *) "Executable: ", trim(arg)
    write (error_unit, *) 'float precision is:', precision(D0), ' decimal places'
    read (input_unit, *) h, start, finish, plot_ratio, integrator
    call init_model()
    w1 = D1 / (D4 - D4 ** (D1 / D9))
    x1 = D1 / (D4 - D4 ** (D1 / D7))
    y1 = D1 / (D4 - D4 ** (D1 / D5))
    z1 = D1 / (D4 - D4 ** (D1 / D3))
    w0 = D1 - D4 * w1
    x0 = D1 - D4 * x1
    y0 = D1 - D4 * y1
    z0 = D1 - D4 * z1
    cd_2(1) = h * D05
    cd_2(2) = h
    cd_4(1) = h * z1 * D05  !!
    cd_4(2) = h * z1
    cd_4(3) = h * z1
    cd_4(4) = h * z1
    cd_4(5) = h * (z0 + z1) * D05
    cd_4(6) = h * z0  !
    cd_6(1) = h * z1 * y1 * D05  !!
    cd_6(2) = h * z1 * y1
    cd_6(3) = h * z1 * y1
    cd_6(4) = h * z1 * y1
    cd_6(5) = h * (z0 + z1) * y1 * D05
    cd_6(6) = h * z0 * y1  !
    cd_6(7) = h * (z0 + z1) * y1 * D05
    cd_6(8) = h * z1 * y1
    cd_6(9) = h * z1 * y1
    cd_6(10) = h * z1 * y1
    cd_6(11) = h * z1 * y1  !!
    cd_6(12) = h * z1 * y1
    cd_6(13) = h * z1 * y1
    cd_6(14) = h * z1 * y1
    cd_6(15) = h * (z0 + z1) * y1 * D05
    cd_6(16) = h * z0 * y1  !
    cd_6(17) = h * (z0 + z1) * y1 * D05
    cd_6(18) = h * z1 * y1
    cd_6(19) = h * z1 * y1
    cd_6(20) = h * z1 * y1
    cd_6(21) = h * z1 * (y1 + y0) * D05  !!
    cd_6(22) = h * z1 * y0
    cd_6(23) = h * z1 * y0
    cd_6(24) = h * z1 * y0
    cd_6(25) = h * (z0 + z1) * y0 * D05
    cd_6(26) = h * z0 * y0  !
    select case (integrator)
        case ("b1")
            write (error_unit, *) "1st Order (Euler-Cromer)"
            call evolve(first_order)
        case ("b2")
            write (error_unit, *) "2nd Order (Stormer-Verlet)"
            call evolve(second_order)
        case ("b4")
            write (error_unit, *) "4th Order (Smith4)"
            call evolve(fourth_order)
        case ("b6")
            write (error_unit, *) "6th Order (Smith6)"
            call evolve(sixth_order)
        case ("b8")
            write (error_unit, *) "8th Order (Suzuki composition)"
            call evolve(eightth_order)
        case ("b10")
            write (error_unit, *) "10th Order (Suzuki composition)"
            call evolve(tenth_order)
        case default
            error stop "Invalid integrator method"
    end select
contains
    subroutine evolve (nth_order)
        real(kind=16) :: time = D0
        integer(kind=16) :: counter = 0
        do
            if ((time >= start) .and. (mod(counter, plot_ratio) == 0)) then
                call plot(time)
            end if
            call nth_order()
            counter = counter + 1
            time = t_update(time, h, counter)
            if (.not. (time < finish) .or. .not. carry_on) exit
        end do
    end subroutine evolve

    subroutine first_order()
        call q_update(h)
        call p_update(h)
    end subroutine first_order

    subroutine second_order ()
        call q_update(cd_2(1))
        call p_update(cd_2(2))
        call q_update(cd_2(1))
    end subroutine second_order

    subroutine fourth_order ()
        call q_update(cd_4(1))
        call p_update(cd_4(2))
        call q_update(cd_4(3))
        call p_update(cd_4(4))
        call q_update(cd_4(5))
        call p_update(cd_4(6))
        call q_update(cd_4(5))
        call p_update(cd_4(4))
        call q_update(cd_4(3))
        call p_update(cd_4(2))
        call q_update(cd_4(1))
    end subroutine fourth_order

    subroutine base_6 (s)
        real(kind=16), intent(in) :: s
        call q_update(s * cd_6(1))
        call p_update(s * cd_6(2))
        call q_update(s * cd_6(3))
        call p_update(s * cd_6(4))
        call q_update(s * cd_6(5))
        call p_update(s * cd_6(6))
        call q_update(s * cd_6(7))
        call p_update(s * cd_6(8))
        call q_update(s * cd_6(9))
        call p_update(s * cd_6(10))
        call q_update(s * cd_6(11))
        call p_update(s * cd_6(12))
        call q_update(s * cd_6(13))
        call p_update(s * cd_6(14))
        call q_update(s * cd_6(15))
        call p_update(s * cd_6(16))
        call q_update(s * cd_6(17))
        call p_update(s * cd_6(18))
        call q_update(s * cd_6(19))
        call p_update(s * cd_6(20))
        call q_update(s * cd_6(21))
        call p_update(s * cd_6(22))
        call q_update(s * cd_6(23))
        call p_update(s * cd_6(24))
        call q_update(s * cd_6(25))
        call p_update(s * cd_6(26))
        call q_update(s * cd_6(25))
        call p_update(s * cd_6(24))
        call q_update(s * cd_6(23))
        call p_update(s * cd_6(22))
        call q_update(s * cd_6(21))
        call p_update(s * cd_6(20))
        call q_update(s * cd_6(19))
        call p_update(s * cd_6(18))
        call q_update(s * cd_6(17))
        call p_update(s * cd_6(16))
        call q_update(s * cd_6(15))
        call p_update(s * cd_6(14))
        call q_update(s * cd_6(13))
        call p_update(s * cd_6(12))
        call q_update(s * cd_6(11))
        call p_update(s * cd_6(10))
        call q_update(s * cd_6(9))
        call p_update(s * cd_6(8))
        call q_update(s * cd_6(7))
        call p_update(s * cd_6(6))
        call q_update(s * cd_6(5))
        call p_update(s * cd_6(4))
        call q_update(s * cd_6(3))
        call p_update(s * cd_6(2))
        call q_update(s * cd_6(1))
    end subroutine base_6

    subroutine sixth_order ()
        call base_6(D1)
    end subroutine sixth_order

    subroutine compose (base, s, forward, back)
        real(kind=16), intent(in) :: s, forward, back
        call base(s * forward)
        call base(s * forward)
        call base(s * back)
        call base(s * forward)
        call base(s * forward)
    end subroutine compose

    subroutine base_8 (s)
        real(kind=16), intent(in) :: s
        call compose(base_6, s, x1, x0)
    end subroutine base_8

    subroutine eightth_order ()
        call base_8(D1)
    end subroutine eightth_order

    subroutine tenth_order ()
        call compose(base_8, D1, w1, w0)
    end subroutine tenth_order
end program Symplectic
