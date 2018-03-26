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
    real(kind=16), parameter :: D05 = 0.5, D1 = 1.0, D4 = 4.0, D3 = 3.0, D5 = 5.0, D7 = 7.0, D9 = 9.0
    real(kind=16) :: w1, w0, x1, x0, y1, y0, z1, z0, h, start, finish
    real(kind=16), dimension(2) :: cd_2
    real(kind=16), dimension(6) :: cd_4
    real(kind=16), dimension(26) :: cd_6
    real(kind=16), dimension(126) :: cd_8
    integer :: plot_ratio
    character (len=3) :: integrator
    character(len=32) :: arg
    call get_command_argument(0, arg)
    write (error_unit, *) "Executable: ", trim(arg)
    write (error_unit, *) 'float precision is:', precision(D1), ' decimal places'
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
    cd_8(1) = h * z1 * y1 * x1 * D05  !!
    cd_8(2) = h * z1 * y1 * x1
    cd_8(3) = h * z1 * y1 * x1
    cd_8(4) = h * z1 * y1 * x1
    cd_8(5) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(6) = h * z0 * y1 * x1  !
    cd_8(7) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(8) = h * z1 * y1 * x1
    cd_8(9) = h * z1 * y1 * x1
    cd_8(10) = h * z1 * y1 * x1
    cd_8(11) = h * z1 * y1 * x1  !!
    cd_8(12) = h * z1 * y1 * x1
    cd_8(13) = h * z1 * y1 * x1
    cd_8(14) = h * z1 * y1 * x1
    cd_8(15) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(16) = h * z0 * y1 * x1  !
    cd_8(17) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(18) = h * z1 * y1 * x1
    cd_8(19) = h * z1 * y1 * x1
    cd_8(20) = h * z1 * y1 * x1
    cd_8(21) = h * z1 * (y1 + y0) * x1 * D05  !!
    cd_8(22) = h * z1 * y0 * x1
    cd_8(23) = h * z1 * y0 * x1
    cd_8(24) = h * z1 * y0 * x1
    cd_8(25) = h * (z0 + z1) * y0 * x1 * D05
    cd_8(26) = h * z0 * y0 * x1  !
    cd_8(27) = h * (z0 + z1) * y0 * x1 * D05
    cd_8(28) = h * z1 * y0 * x1
    cd_8(29) = h * z1 * y0 * x1
    cd_8(30) = h * z1 * y0 * x1
    cd_8(31) = h * z1 * (y1 + y0) * x1 * D05  !!
    cd_8(32) = h * z1 * y1 * x1
    cd_8(33) = h * z1 * y1 * x1
    cd_8(34) = h * z1 * y1 * x1
    cd_8(35) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(36) = h * z0 * y1 * x1  !
    cd_8(37) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(38) = h * z1 * y1 * x1
    cd_8(39) = h * z1 * y1 * x1
    cd_8(40) = h * z1 * y1 * x1
    cd_8(41) = h * z1 * y1 * x1  !!
    cd_8(42) = h * z1 * y1 * x1
    cd_8(43) = h * z1 * y1 * x1
    cd_8(44) = h * z1 * y1 * x1
    cd_8(45) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(46) = h * z0 * y1 * x1  !
    cd_8(47) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(48) = h * z1 * y1 * x1
    cd_8(49) = h * z1 * y1 * x1
    cd_8(50) = h * z1 * y1 * x1
    cd_8(51) = h * z1 * y1 * x1  !!
    cd_8(52) = h * z1 * y1 * x1
    cd_8(53) = h * z1 * y1 * x1
    cd_8(54) = h * z1 * y1 * x1
    cd_8(55) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(56) = h * z0 * y1 * x1  !
    cd_8(57) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(58) = h * z1 * y1 * x1
    cd_8(59) = h * z1 * y1 * x1
    cd_8(60) = h * z1 * y1 * x1
    cd_8(61) = h * z1 * y1 * x1  !!
    cd_8(62) = h * z1 * y1 * x1
    cd_8(63) = h * z1 * y1 * x1
    cd_8(64) = h * z1 * y1 * x1
    cd_8(65) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(66) = h * z0 * y1 * x1  !
    cd_8(67) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(68) = h * z1 * y1 * x1
    cd_8(69) = h * z1 * y1 * x1
    cd_8(70) = h * z1 * y1 * x1
    cd_8(71) = h * z1 * (y1 + y0) * x1 * D05  !!
    cd_8(72) = h * z1 * y0 * x1
    cd_8(73) = h * z1 * y0 * x1
    cd_8(74) = h * z1 * y0 * x1
    cd_8(75) = h * (z0 + z1) * y0 * x1 * D05
    cd_8(76) = h * z0 * y0 * x1  !
    cd_8(77) = h * (z0 + z1) * y0 * x1 * D05
    cd_8(78) = h * z1 * y0 * x1
    cd_8(79) = h * z1 * y0 * x1
    cd_8(80) = h * z1 * y0 * x1
    cd_8(81) = h * z1 * (y1 + y0) * x1 * D05  !!
    cd_8(82) = h * z1 * y1 * x1
    cd_8(83) = h * z1 * y1 * x1
    cd_8(84) = h * z1 * y1 * x1
    cd_8(85) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(86) = h * z0 * y1 * x1  !
    cd_8(87) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(88) = h * z1 * y1 * x1
    cd_8(89) = h * z1 * y1 * x1
    cd_8(90) = h * z1 * y1 * x1
    cd_8(91) = h * z1 * y1 * x1  !!
    cd_8(92) = h * z1 * y1 * x1
    cd_8(93) = h * z1 * y1 * x1
    cd_8(94) = h * z1 * y1 * x1
    cd_8(95) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(96) = h * z0 * y1 * x1  !
    cd_8(97) = h * (z0 + z1) * y1 * x1 * D05
    cd_8(98) = h * z1 * y1 * x1
    cd_8(99) = h * z1 * y1 * x1
    cd_8(100) = h * z1 * y1 * x1
    cd_8(101) = h * z1 * y1 * (x1 + x0) * D05  !!
    cd_8(102) = h * z1 * y1 * x0
    cd_8(103) = h * z1 * y1 * x0
    cd_8(104) = h * z1 * y1 * x0
    cd_8(105) = h * (z0 + z1) * y1 * x0 * D05
    cd_8(106) = h * z0 * y1 * x0  !
    cd_8(107) = h * (z0 + z1) * y1 * x0 * D05
    cd_8(108) = h * z1 * y1 * x0
    cd_8(109) = h * z1 * y1 * x0
    cd_8(110) = h * z1 * y1 * x0
    cd_8(111) = h * z1 * y1 * x0  !!
    cd_8(112) = h * z1 * y1 * x0
    cd_8(113) = h * z1 * y1 * x0
    cd_8(114) = h * z1 * y1 * x0
    cd_8(115) = h * (z0 + z1) * y1 * x0 * D05
    cd_8(116) = h * z0 * y1 * x0  !
    cd_8(117) = h * (z0 + z1) * y1 * x0 * D05
    cd_8(118) = h * z1 * y1 * x0
    cd_8(119) = h * z1 * y1 * x0
    cd_8(120) = h * z1 * y1 * x0
    cd_8(121) = h * z1 * (y1 + y0) * x0 * D05  !!
    cd_8(122) = h * z1 * y0 * x0
    cd_8(123) = h * z1 * y0 * x0
    cd_8(124) = h * z1 * y0 * x0
    cd_8(125) = h * (z0 + z1) * y0 * x0 * D05
    cd_8(126) = h * z0 * y0 * x0  !
    select case (integrator)
        case ("b2")
            write (error_unit, *) "2nd Order (Stormer-Verlet)"
            call evolve(second_order)
        case ("b4")
            write (error_unit, *) "4th Order (Smith)"
            call evolve(fourth_order)
        case ("b6")
            write (error_unit, *) "6th Order (Smith)"
            call evolve(sixth_order)
        case ("b8")
            write (error_unit, *) "8th Order (Smith)"
            call evolve(eightth_order)
        case ("b10")
            write (error_unit, *) "10th Order (Suzuki composition)"
            call evolve(tenth_order)
        case default
            error stop "Invalid integrator method"
    end select
contains
    subroutine evolve (nth_order)
        real(kind=16) :: time = 0.0
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

    subroutine second_order ()
        call q_update(cd_2(1))
        call p_update(cd_2(2))
        call q_update(cd_2(1))
    end subroutine second_order

    subroutine smith (cd, s)
        real(kind=16), dimension(:) :: cd
        real(kind=16), intent(in) :: s
        integer :: array_size
        integer :: i
        array_size = size(cd)
        do i = 1, array_size
            if (mod(i, 2) == 0) then
                call p_update(cd(i) * s)
            else
                call q_update(cd(i) * s)
            end if
        end do
        do i = array_size - 1, 1, -1
            if (mod(i, 2) == 0) then
                call p_update(cd(i) * s)
            else
                call q_update(cd(i) * s)
            end if
        end do
    end subroutine smith

    subroutine fourth_order ()
        call smith(cd_4, D1)
    end subroutine fourth_order

    subroutine sixth_order ()
        call smith(cd_6, D1)
    end subroutine sixth_order

    subroutine eightth_order ()
        call smith(cd_8, D1)
    end subroutine eightth_order

    subroutine tenth_order ()
        call smith(cd_8, w1)
        call smith(cd_8, w1)
        call smith(cd_8, w0)
        call smith(cd_8, w1)
        call smith(cd_8, w1)
    end subroutine tenth_order
end program Symplectic
