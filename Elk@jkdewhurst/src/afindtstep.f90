
! Copyright (C) 2020 Peter Elliott, J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: afindtstep
! !INTERFACE:
subroutine afindtstep
! !USES:
use modmain
use modtddft
use modmpi
! !DESCRIPTION:
!   Performs a time step of the macroscopic Maxwell equation and updates the
!   induced vector potential ${\bf A}(t)$. In practice, a more general damped
!   Proca equation is solved:
!   $$ p_0{\bf A}+p_1\dot{\bf A}+p_2\ddot{\bf A}=
!    \frac{4\pi c}{\Omega}{\bf J}, $$
!   where $\Omega$ is the unit cell volume, ${\bf J}$ is the total current
!   across the unit cell, and the parameters $p_i$, $i=0,1,2$ are stored in the
!   array {\tt afindpm}. This generalisation allows for both a mass and damping
!   term, however the default values of $p_0=p_1=0$ and $p_2=1$ recover the
!   physical Maxwell equation.
!
! !REVISION HISTORY:
!   Created January 2020 (P. Elliott)
!   Added mass and damping terms, December 2022 (JKD)
!EOP
!BOC
implicit none
! local variables
integer i
real(8) dt,t1,t2,t3
! time step length
dt=times(itimes+1)-times(itimes)
! add to the time derivative of the induced A-field
t1=fourpi*solsc/omega
t2=dt/afindpm(2)
do i=1,3
  t3=t1*jtot(i)-afindpm(1)*afindt(i,1)-afindpm(0)*afindt(i,0)
  afindt(i,1)=afindt(i,1)+t2*t3
end do
! add to the induced A-field
afindt(:,0)=afindt(:,0)+afindt(:,1)*dt
! add to the total A-field
afieldt(:,itimes+1)=afieldt(:,itimes+1)+afindt(:,0)
! write the induced A-field and its time derivative to file
if (mp_mpi) then
  open(50,file='AFINDT.OUT',form='FORMATTED')
  write(50,'(6G18.10)') afindt(:,:)
  close(50)
end if
end subroutine
!EOC

