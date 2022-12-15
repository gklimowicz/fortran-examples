
! Copyright (C) 2011 S. Sharma, J. K. Dewhurst and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine genwfsvp(tsh,tgp,nst,idx,ngdg,igf,vpl,ngp,igpig,wfmt,ld,wfir)
use modmain
implicit none
! arguments
logical, intent(in) :: tsh,tgp
integer, intent(in) :: nst,idx(*),ngdg(3),igf(*)
real(8), intent(in) :: vpl(3)
integer, intent(out) :: ngp(nspnfv),igpig(ngkmax,nspnfv)
complex(8), intent(out) :: wfmt(npcmtmax,natmtot,nspinor,nst)
integer, intent(in) :: ld
complex(8), intent(out) :: wfir(ld,nspinor,nst)
! local variables
integer ispn
real(8) vl(3),vc(3)
! automatic arrays
real(8) vgpl(3,ngkmax,nspnfv),vgpc(3,ngkmax),gpc(ngkmax)
! allocatable arrays
complex(8), allocatable :: sfacgp(:,:),apwalm(:,:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
allocate(sfacgp(ngkmax,natmtot))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot,nspnfv))
! loop over first-variational spins
do ispn=1,nspnfv
  vl(:)=vpl(:)
  vc(:)=bvec(:,1)*vpl(1)+bvec(:,2)*vpl(2)+bvec(:,3)*vpl(3)
! spin-spiral case
  if (spinsprl) then
    if (ispn.eq.1) then
      vl(:)=vl(:)+0.5d0*vqlss(:)
      vc(:)=vc(:)+0.5d0*vqcss(:)
    else
      vl(:)=vl(:)-0.5d0*vqlss(:)
      vc(:)=vc(:)-0.5d0*vqcss(:)
    end if
  end if
! generate the G+p-vectors
  call gengkvec(ngvc,ivg,vgc,vl,vc,gkmax,ngkmax,ngp(ispn),igpig(:,ispn), &
   vgpl(:,:,ispn),vgpc,gpc)
! generate structure factors for G+p-vectors
  call gensfacgp(ngp(ispn),vgpc,ngkmax,sfacgp)
! find the matching coefficients
  call match(ngp(ispn),vgpc,gpc,sfacgp,apwalm(:,:,:,:,ispn))
end do
deallocate(sfacgp)
! get the first- and second-variational eigenvectors from file
allocate(evecfv(nmatmax,nstfv,nspnfv),evecsv(nstsv,nstsv))
call getevecfv(filext,0,vpl,vgpl,evecfv)
call getevecsv(filext,0,vpl,evecsv)
! calculate the second-variational wavefunctions
call genwfsv(tsh,tgp,nst,idx,ngdg,igf,ngp,igpig,apwalm,evecfv,evecsv,wfmt,ld, &
 wfir)
deallocate(apwalm,evecfv,evecsv)
end subroutine

