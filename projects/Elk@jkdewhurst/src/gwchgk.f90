
! Copyright (C) 2018 P. Elliott, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gwchgk(ik,chgk)
use modmain
use modgw
implicit none
! arguments
integer, intent(in) :: ik
real(8), intent(out) :: chgk
! local variables
integer ist,iw
real(8) e
complex(8) z1
! allocatable arrays
complex(8), allocatable :: gs(:),g(:,:),gf(:),gl(:)
complex(8), allocatable :: se(:,:,:)
! external functions
complex(8), external :: gtwsum
! read the self-energy from file
allocate(se(nstsv,nstsv,0:nwfm))
call getgwsefm(ik,se)
! allocate local arrays
allocate(gs(nstsv),g(nstsv,nstsv))
allocate(gf(nstsv),gl(nstsv))
chgk=0.d0
do iw=0,nwfm
! compute the diagonal matrix G_s
  do ist=1,nstsv
    e=evalsv(ist,ik)-efermi
    gs(ist)=1.d0/(wfm(iw)-e)
  end do
! compute 1 - G_s Sigma
  do ist=1,nstsv
    z1=-gs(ist)
    g(ist,:)=z1*se(ist,:,iw)
    g(ist,ist)=g(ist,ist)+1.d0
  end do
! invert this matrix
  call zminv(nstsv,g)
! take the trace of G = (1 - G_s Sigma)^(-1) G_s
  do ist=1,nstsv
    g(ist,ist)=g(ist,ist)*gs(ist)
    chgk=chgk+dble(g(ist,ist))
  end do
! store the Green's function at the first and last frequencies
  if (iw.eq.0) then
    do ist=1,nstsv
      gf(ist)=g(ist,ist)
    end do
  end if
  if (iw.eq.nwfm) then
    do ist=1,nstsv
      gl(ist)=g(ist,ist)
    end do
  end if
end do
! add the Matsubara tails analytically
do ist=1,nstsv
  chgk=chgk+dble(gtwsum(gf(ist),gl(ist)))
end do
chgk=chgk*wkpt(ik)*occmax*kboltz*tempk
deallocate(se,gs,g,gf,gl)
end subroutine

