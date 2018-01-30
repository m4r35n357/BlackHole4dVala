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
    real(16), parameter :: MD0=0.0, MD05=0.5
    real(16) :: x, Ux = MD0, k, m, h0
    logical :: carry_on = .true.
contains
    subroutine init_model()
        write (error_unit, *) "Harmonic Oscillator"
        read (input_unit, *) k, m, x
        h0 = hamiltonian()
    end subroutine init_model

    real(16) function hamiltonian ()
        hamiltonian = MD05 * (m * Ux**2 + k * x**2)
    end function hamiltonian

    subroutine q_update(c)
        real(16), intent(in) :: c
        x = x + c * m * Ux
    end subroutine q_update

    subroutine p_update(d)
        real(16), intent(in) :: d
        Ux = Ux - d * k * x
    end subroutine p_update

    real(16) function t_update(time, step, counter)
        real(16), intent(in) :: time, step
        integer(16), intent(in) :: counter
        t_update = step * counter
    end function t_update

    subroutine plot(time)
        real(16), intent(in) :: time
        write (output_unit, '(A, 13(ES16.9, A))') '{"tau":',time,',"v4e":',hamiltonian()-h0,',"t":',time,',"x":',x,',"xP":',Ux,'}'
    end subroutine plot
end module Model

