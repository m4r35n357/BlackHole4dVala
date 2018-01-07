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
    implicit none
    real(16), parameter :: MD0=0.0_16, MD05=0.5_16, MD1=1.0_16, MD2=2.0_16, PI_2=MD05*acos(-MD1)  ! CONSTANTS
    real(16) :: L, L2, E0, V, H0  ! INTERMEDIATE VARIABLES
    real(16) :: r, ph = MD0, rDot, phDot  ! PARTICLE VARIABLES (coordinates and velocities))
    logical :: carry_on = .true.
contains
    subroutine init_model_vars()
        real(16) :: lFac, r0, LC
        write (0, *) "Newtonian Central-Body Problem"
        read(*,*) lFac, r0
        r = r0
        LC = sqrt(r0)
        E0 = MD05 * LC**2 / (r * r) - MD1 / r
        L = lFac * LC
        L2 = L**2
        V = MD05 * L2 / (r * r) - MD1 / r
        rDot = - sqrt(MD2 * merge(E0 - V, V - E0, E0 > V))
        H0 = hamiltonian()
    end subroutine init_model_vars

    real(16) function hamiltonian ()
        hamiltonian = MD05 * (rDot * rDot + L2 / (r * r)) - MD1 / r
    end function hamiltonian

    subroutine q_update(c)
        real(16), intent(in) :: c
        r = r + c * rDot
        phDot = L / (r * r)
        ph = ph + c * phDot
    end subroutine q_update

    subroutine p_update(d)
        real(16), intent(in) :: d
        rDot = rDot - d * (MD1 / (r * r) - L2 / (r * r * r))
    end subroutine p_update

    real(16) function t_update(time, step, counter)
        real(16), intent(in) :: time, step
        integer, intent(in) :: counter
        carry_on = r > MD2
        t_update = step * counter
    end function t_update

    subroutine plot(time)
        real(16), intent(in) :: time
        write (*, '(A, 13(ES16.9, A))') '{"tau":',time, ',"v4e":',hamiltonian() - H0,&
            ',"t":', time,',"r":',r,',"th":',PI_2,',"ph":',ph, ',"tP":',1.0,',"rDot":',rDot,',"thP":',0.0,',"phP":',phDot,'}'
    end subroutine plot
end module Model

