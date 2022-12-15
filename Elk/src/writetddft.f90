
! Copyright (C) 2014 K. Krieger, J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetddft
use modmain
use modtddft
use moddftu
implicit none
! local variables
integer is,ia,ias
real(8) vc(3),vl(3)
character(256) fext
! allocatable arrays
real(8), allocatable :: rvfmt(:,:,:),rvfir(:,:)
! file extension
write(fext,'("_TS",I8.8,".OUT")') itimes
! delete all files at first time step
if (itimes.le.1) then
  open(50,file='CHARGEMT_TD.OUT')
  close(50,status='DELETE')
  open(50,file='CHARGEIR_TD.OUT')
  close(50,status='DELETE')
  open(50,file='MOMENT_TD.OUT')
  close(50,status='DELETE')
  open(50,file='MOMENTM_TD.OUT')
  close(50,status='DELETE')
  open(50,file='MOMENTMT_TD.OUT')
  close(50,status='DELETE')
  open(50,file='MOMENTIR_TD.OUT')
  close(50,status='DELETE')
  open(50,file='JTOT_TD.OUT')
  close(50,status='DELETE')
  open(50,file='JTOTM_TD.OUT')
  close(50,status='DELETE')
  open(50,file='TOTENERGY_TD.OUT')
  close(50,status='DELETE')
  if (tddos) then
    open(50,file='TDTEMP.OUT')
    close(50,status='DELETE')
  end if
  if (tdlsj) then
    open(50,file='TDLSJ.OUT')
    close(50,status='DELETE')
  end if
  if (tforce) then
    open(50,file='FORCETOT_TD.OUT')
    close(50,status='DELETE')
  end if
  if (tdatpos) then
    open(50,file='ATPOSL_TD.OUT')
    close(50,status='DELETE')
  end if
  if (tafindt) then
    open(50,file='AFIND_TD.OUT')
    close(50,status='DELETE')
  end if
end if
! muffin-tin charges
open(50,file='CHARGEMT_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(G18.10)') times(itimes)
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,'(2I4,G18.10)') is,ia,chgmt(ias)
  end do
end do
write(50,*)
close(50)
! interstitial charge
open(50,file='CHARGEIR_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(2G18.10)') times(itimes),chgir
close(50)
! spin moment
open(50,file='MOMENT_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(4G18.10)') times(itimes),momtot(1:ndmag)
close(50)
open(50,file='MOMENTM_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(2G18.10)') times(itimes),momtotm
close(50)
! muffin-tin moments
open(50,file='MOMENTMT_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(G18.10)') times(itimes)
do is=1,nspecies
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,'(2I4,3G18.10)') is,ia,mommt(1:ndmag,ias)
  end do
end do
write(50,*)
close(50)
! interstitial moment
open(50,file='MOMENTIR_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(4G18.10)') times(itimes),momir(1:ndmag)
close(50)
! total current
open(50,file='JTOT_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(4G18.10)') times(itimes),jtot(:)
close(50)
! total current magnitude
open(50,file='JTOTM_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(4G18.10)') times(itimes),jtotm
close(50)
! total energy
open(50,file='TOTENERGY_TD.OUT',form='FORMATTED',position='APPEND')
write(50,'(2G18.10)') times(itimes),engytot
close(50)
! total force
if (tforce) then
  open(50,file='FORCETOT_TD.OUT',form='FORMATTED',position='APPEND')
  write(50,'(I8,G18.10)') itimes,times(itimes)
  do is=1,nspecies
    do ia=1,natoms(is)
      ias=idxas(ia,is)
      write(50,'(2I4,3G18.10)') is,ia,forcetot(:,ias)
    end do
  end do
  close(50)
end if
! write the time-dependent atomic positions in lattice coordinates
if (tdatpos) then
  open(50,file='ATPOSL_TD.OUT',form='FORMATTED',position='APPEND')
  write(50,'(I8,G18.10)') itimes,times(itimes)
  do is=1,nspecies
    do ia=1,natoms(is)
      vc(:)=atposc(:,ia,is)+datposc(:,0,ia,is)
      call r3mv(ainv,vc,vl)
      write(50,'(2I4,3G18.10)') is,ia,vl(:)
    end do
  end do
  close(50)
end if
if (tafindt) then
  open(50,file='AFIND_TD.OUT',form='FORMATTED',position='APPEND')
  write(50,'(4G18.10)') times(itimes),afindt(:,0)
  close(50)
end if
! write optional quantities
if (ntswrite(1).le.0) return
if (mod(itimes-1,ntswrite(1)).ne.0) return
if ((itimes.gt.1).and.(itimes.lt.ntswrite(2))) return
! charge density in 1D
if (tdrho1d) then
  open(50,file='RHO1D'//trim(fext),form='FORMATTED',action='WRITE')
  open(51,file='RHOLINES.OUT',form='FORMATTED',action='WRITE')
  call plot1d(50,51,1,rhomt,rhoir)
  close(50)
  close(51)
end if
! charge density in 2D
if (tdrho2d) then
  open(50,file='RHO2D'//trim(fext),form='FORMATTED',action='WRITE')
  call plot2d(.false.,50,1,rhomt,rhoir)
  close(50)
end if
! charge density in 3D
if (tdrho3d) then
  open(50,file='RHO3D'//trim(fext),form='FORMATTED',action='WRITE')
  call plot3d(50,1,rhomt,rhoir)
  close(50)
end if
! magnetisation in 1D, 2D or 3D
if ((tdmag1d.or.tdmag2d.or.tdmag3d).and.(spinpol)) then
  allocate(rvfmt(npmtmax,natmtot,3),rvfir(ngtot,3))
  if (ncmag) then
! non-collinear
    rvfmt(:,:,:)=magmt(:,:,:)
    rvfir(:,:)=magir(:,:)
  else
! collinear
    rvfmt(:,:,1:2)=0.d0
    rvfir(:,1:2)=0.d0
    rvfmt(:,:,3)=magmt(:,:,1)
    rvfir(:,3)=magir(:,1)
  end if
  if (tdmag1d) then
    open(50,file='MAG1D'//trim(fext),form='FORMATTED',action='WRITE')
    open(51,file='MAGLINES.OUT',form='FORMATTED',action='WRITE')
    call plot1d(50,51,3,rvfmt,rvfir)
    close(50)
    close(51)
  end if
  if (tdmag2d) then
    open(50,file='MAG2D'//trim(fext),form='FORMATTED',action='WRITE')
    call plot2d(.true.,50,3,rvfmt,rvfir)
    close(50)
  end if
  if (tdmag3d) then
    open(50,file='MAG3D'//trim(fext),form='FORMATTED',action='WRITE')
    call plot3d(50,3,rvfmt,rvfir)
    close(50)
  end if
  deallocate(rvfmt,rvfir)
end if
! gauge-invariant current density in 1D
if (tdjr1d) then
  open(50,file='JR1D'//trim(fext),form='FORMATTED',action='WRITE')
  open(51,file='JRLINES.OUT',form='FORMATTED',action='WRITE')
  call plot1d(50,51,3,jrmt,jrir)
  close(50)
  close(51)
end if
! gauge-invariant current density in 2D
if (tdjr2d) then
  open(50,file='JR2D'//trim(fext),form='FORMATTED',action='WRITE')
  call plot2d(.true.,50,3,jrmt,jrir)
  close(50)
end if
! gauge-invariant current density in 3D
if (tdjr3d) then
  open(50,file='JR3D'//trim(fext),form='FORMATTED',action='WRITE')
  call plot3d(50,3,jrmt,jrir)
  close(50)
end if
! calculate and write tensor moments
if (dftu.ne.0) then
  if (tmwrite) call writetm3td
end if
! write time-dependent DOS
if (tddos) then
  call writetddos(fext)
end if
! write the atomic forces
if (tforce) then
  open(50,file='FORCES'//trim(fext),form='FORMATTED',action='WRITE')
  call writeforces(50)
  close(50)
end if
end subroutine

