
! Copyright (C) 2020 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine initeph
use modmain
use modphonon
use modbog
implicit none
! local variables
integer iq,ik,n,i
! automatic arrays
complex(8) u(nstsv,nstsv),v(nstsv,nstsv)
complex(8) w(nbph,nbph),x(nbph,nbph),y(nbph)
! allocatable arrays
complex(8), allocatable :: dynq(:,:,:),ev(:,:)
complex(8), allocatable :: ephmat(:,:,:)

! combined target array for fermionic and bosonic density matrices
if (allocated(duvwx)) deallocate(duvwx)
n=2*nstsv*nstsv*nkpt+2*nbph*nbph*nqpt
allocate(duvwx(n))

!------------------------------!
!     electronic variables     !
!------------------------------!
if (allocated(evaluv)) deallocate(evaluv)
allocate(evaluv(nstsv,nkpt))
if (allocated(vnorm)) deallocate(vnorm)
allocate(vnorm(nstsv,nkpt))
! associate the electronic density matrices with target
dvv(1:nstsv,1:nstsv,1:nkpt)=>duvwx(1:)
i=nstsv*nstsv*nkpt+1
duv(1:nstsv,1:nstsv,1:nkpt)=>duvwx(i:)
i=i+nstsv*nstsv*nkpt
if (task.eq.270) then
! initialise the density matrices to random numbers
  dvv(:,:,:)=0.d0
  duv(:,:,:)=0.d0
  do ik=1,nkpt
    call rndevsv(1.d0,dvv(:,:,ik))
    call rndevsv(1.d0,duv(:,:,ik))
  end do
else
  do ik=1,nkpt
! get the eigenvalues from file
    call getevaluv(ik,evaluv(:,ik))
! get the eigenvectors from file
    call getevecuv(ik,vkl(:,ik),u,v)
! calculate the density matrices
    call dmatuv(nstsv,efermi,evalsv(:,ik),u,v,dvv(:,:,ik),duv(:,:,ik), &
     vnorm(:,ik))
  end do
end if

!----------------------------!
!     phononic variables     !
!----------------------------!
allocate(dynq(nbph,nbph,nqpt),ev(nbph,nbph))
! read in the dynamical matrices
call readdyn(dynq)
! apply the acoustic sum rule
call sumrule(dynq)
! find the eigenvalues of the dynamical matrices and store in global array
do iq=1,nqpt
  call dynev(dynq(:,:,iq),wphq(:,iq),ev)
end do
deallocate(dynq,ev)
if (allocated(evalwx)) deallocate(evalwx)
allocate(evalwx(nbph,nqpt))
if (allocated(xnorm)) deallocate(xnorm)
allocate(xnorm(nbph,nqpt))
! associate the phononic density matrices with target
dxx(1:nbph,1:nbph,1:nqpt)=>duvwx(i:)
i=i+nbph*nbph*nqpt
dwx(1:nbph,1:nbph,1:nqpt)=>duvwx(i:)
if (task.eq.270) then
! zero the density matrices
  dxx(:,:,:)=0.d0
  dwx(:,:,:)=0.d0
else
  do iq=1,nqpt
! get the eigenvalues from file
    call getevalwx(iq,evalwx(:,iq))
! get the eigenvectors from file
    call getevecwxy(iq,w,x,y)
! calculate the density matrices
    call dmatwx(nbph,w,x,dxx(:,:,iq),dwx(:,:,iq),xnorm(:,iq))
  end do
end if

!-----------------------------------!
!     electron-phonon variables     !
!-----------------------------------!
if (any(task.eq.[270,271])) then
! allocate the electron-phonon matrix elements array
  if (allocated(ephmkq)) deallocate(ephmkq)
  allocate(ephmkq(nstsv,nstsv,nbph,nkptnr,nqpt))
! read the matrix elements from file and store in global array
  allocate(ephmat(nstsv,nstsv,nbph))
  do iq=1,nqpt
    do ik=1,nkptnr
      call getephmat(iq,ik,ephmat)
      ephmkq(:,:,:,ik,iq)=ephmat(:,:,:)
    end do
! zero the electron-phonon coupling for phonon small phonon frequencies
    do i=1,nbph
      if (wphq(i,iq).lt.wphcut) ephmkq(:,:,i,:,iq)=0.d0
    end do
  end do
  deallocate(ephmat)
end if

end subroutine

