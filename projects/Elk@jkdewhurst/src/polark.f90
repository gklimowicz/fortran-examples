
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine polark(ik,l,expqmt,pl)
use modmain
implicit none
! arguments
integer, intent(in) :: ik,l
complex(8), intent(in) :: expqmt(npcmtmax,natmtot)
real(8), intent(inout) :: pl
! local variables
integer jk,nst,ist
real(8) vkql(3)
complex(8) z1
! automatic arrays
integer idx(nstsv),ngp(nspnfv)
! allocatable arrays
integer, allocatable :: igpig(:,:)
complex(8), allocatable :: wfmt(:,:,:,:),wfir(:,:,:)
complex(8), allocatable :: wfmtq(:,:,:,:),wfgkq(:,:,:)
complex(8), allocatable :: oq(:,:)
! external functions
complex(8), external :: zmdet
! equivalent reduced k-point
jk=ivkik(ivk(1,ik),ivk(2,ik),ivk(3,ik))
! find the adjacent k-point in lattice coordinates
vkql(:)=vkl(:,ik)
vkql(l)=vkql(l)+1.d0/dble(ngridk(l))
! count and index the occupied states
nst=0
do ist=1,nstsv
  if (evalsv(ist,jk).gt.efermi) cycle
  nst=nst+1
  idx(nst)=ist
end do
allocate(igpig(ngkmax,nspnfv))
! generate the wavefunctions for occupied states at k
allocate(wfmt(npcmtmax,natmtot,nspinor,nst),wfir(ngtot,nspinor,nst))
call genwfsvp(.false.,.false.,nst,idx,ngridg,igfft,vkl(:,ik),ngp,igpig,wfmt, &
 ngtot,wfir)
! generate the wavefunctions for occupied states at k+q
allocate(wfmtq(npcmtmax,natmtot,nspinor,nst),wfgkq(ngkmax,nspinor,nst))
call genwfsvp(.false.,.true.,nst,idx,ngridg,igfft,vkql,ngp,igpig,wfmtq,ngkmax, &
 wfgkq)
! determine the overlap matrix for all occupied states
allocate(oq(nst,nst))
call genolpq(nst,expqmt,ngp,igpig,wfmt,wfir,wfmtq,wfgkq,oq)
! compute the determinant of the matrix
z1=zmdet(nst,oq)
! determine the phase of the determinant and add to total polarisation
pl=pl+atan2(aimag(z1),dble(z1))
deallocate(igpig,wfmt,wfir,wfmtq,wfgkq,oq)
end subroutine

