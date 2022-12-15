
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gwdmatk(ik)
use modmain
use modgw
implicit none
! arguments
integer, intent(in) :: ik
! local variables
integer ist,jst,iw
real(8) e,t1
complex(8) z1
! allocatable arrays
complex(8), allocatable :: evecsv(:,:),se(:,:,:)
complex(8), allocatable :: gs(:),g(:,:),gf(:,:),gl(:,:)
complex(8), allocatable :: d(:,:),a(:,:)
! external functions
complex(8), external :: gtwsum
! read the self-energy from file
allocate(se(nstsv,nstsv,0:nwfm))
call getgwsefm(ik,se)
! allocate local arrays
allocate(gs(nstsv),g(nstsv,nstsv))
allocate(gf(nstsv,nstsv),gl(nstsv,nstsv))
allocate(d(nstsv,nstsv))
! zero the density matrix
d(:,:)=0.d0
! loop over fermionic Matsubara frequencies
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
! compute G = (1 - G_s Sigma)^(-1) G_s
  do jst=1,nstsv
    z1=gs(jst)
    g(:,jst)=g(:,jst)*z1
  end do
! add to the density matrix
  d(:,:)=d(:,:)+g(:,:)
! store the Green's function at the first and last frequencies
  if (iw.eq.0) gf(:,:)=g(:,:)
  if (iw.eq.nwfm) gl(:,:)=g(:,:)
end do
! add the Matsubara tails analytically
do jst=1,nstsv
  do ist=1,nstsv
    d(ist,jst)=d(ist,jst)+gtwsum(gf(ist,jst),gl(ist,jst))
  end do
end do
! multiply by 1/beta
t1=kboltz*tempk
d(:,:)=t1*d(:,:)
deallocate(se,gs,g,gf,gl)
! diagonalise the density matrix for the natural orbitals and occupation numbers
call eveqnzh(nstsv,nstsv,d,occsv(:,ik))
occsv(:,ik)=occsv(:,ik)*occmax
! get the second-variational eigenvectors from file
allocate(evecsv(nstsv,nstsv))
call getevecsv(filext,ik,vkl(:,ik),evecsv)
! apply unitary transformation to the third-variational states so that they
! refer to the first-variational basis
allocate(a(nstsv,nstsv))
call zgemm('N','N',nstsv,nstsv,nstsv,zone,evecsv,nstsv,d,nstsv,zzero,a,nstsv)
! write the density matrix to file as second-variational eigenvectors
call putevecsv(filext,ik,a)
deallocate(evecsv,d,a)
end subroutine

