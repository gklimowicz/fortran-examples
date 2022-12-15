
! Copyright (C) 2018 P. Elliott, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine ksinvert
use modmain
use modmpi
implicit none
! local variables
integer idm,it
real(8) tau
! allocatable arrays
real(8), allocatable :: rhomt_(:,:),rhoir_(:)
real(8), allocatable :: magmt_(:,:,:),magir_(:,:)
real(8), allocatable :: rfmt(:,:),rfir(:)
! copy the existing density and magnetisation
allocate(rhomt_(npmtmax,natmtot),rhoir_(ngtot))
call rfcopy(rhomt,rhoir,rhomt_,rhoir_)
allocate(rfmt(npmtmax,natmtot),rfir(ngtot))
if (spinpol) then
  allocate(magmt_(npmtmax,natmtot,ndmag),magir_(ngtot,ndmag))
  do idm=1,ndmag
    call rfcopy(magmt(:,:,idm),magir(:,idm),magmt_(:,:,idm),magir_(:,idm))
  end do
end if
! initialise the exchange-correlation potential and magnetic field
tau=tauoep(2)
vxcmt(:,:)=tau*vxcmt(:,:)
vxcir(:)=tau*vxcir(:)
if (spinpol) then
  bxcmt(:,:,:)=tau*bxcmt(:,:,:)
  bxcir(:,:)=tau*bxcir(:,:)
end if
! iterative method step size
tau=tauoep(1)
if (mp_mpi) then
  write(*,*)
  write(*,'("Info(ksinvert): inverting the Kohn-Sham equations")')
end if
do it=1,maxitoep
  if (mp_mpi.and.(mod(it,10).eq.0)) then
    write(*,'("Info(ksinvert): done ",I4," iterations of ",I4)') it,maxitoep
  end if
  call hmlrad
  call genevfsv
  call occupy
  call rhomag
! determine the residual and add it to the exchange-correlation potential
  call residual
! set the constant part of V_xc equal to zero
  call rfint0(0.d0,vxcmt,vxcir)
! add the external and Hartree potentials
  call potks(.false.)
  call genvsig
end do
deallocate(rhomt_,rhoir_,rfmt,rfir)
if (spinpol) deallocate(magmt_,magir_)
! synchronise MPI processes
call mpi_barrier(mpicom,ierror)
return

contains

subroutine residual
implicit none
! local variables
integer is,ias,np
! external functions
real(8), external :: rfinp
do ias=1,natmtot
  is=idxis(ias)
  np=npmt(is)
  rfmt(1:np,ias)=rhomt(1:np,ias)-rhomt_(1:np,ias)
  call daxpy(np,tau,rfmt(:,ias),1,vxcmt(:,ias),1)
end do
rfir(:)=rhoir(:)-rhoir_(:)
call daxpy(ngtot,tau,rfir,1,vxcir,1)
resoep=rfinp(rfmt,rfir,rfmt,rfir)
if (spinpol) then
  do idm=1,ndmag
    do ias=1,natmtot
      is=idxis(ias)
      np=npmt(is)
      rfmt(1:np,ias)=magmt(1:np,ias,idm)-magmt_(1:np,ias,idm)
      call daxpy(np,tau,rfmt(:,ias),1,bxcmt(:,ias,idm),1)
    end do
    rfir(:)=magir(:,idm)-magir_(:,idm)
    call daxpy(ngtot,tau,rfir,1,bxcir(:,idm),1)
    resoep=resoep+rfinp(rfmt,rfir,rfmt,rfir)
  end do
end if
resoep=sqrt(resoep)/omega
end subroutine

end subroutine

