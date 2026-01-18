
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine vclcore(wfmt,vmat)
use modmain
use modomp
implicit none
! arguments
complex(4), intent(in) :: wfmt(npcmtmax,natmtot,nspinor,nstsv)
complex(8), intent(inout) :: vmat(nstsv,nstsv)
! local variables
integer ist1,ist2,ist3
integer is,ia,ias,m,nthd
integer nrc,nrci,npc
! automatic arrays
complex(4) wfcr(npcmtmax,2)
complex(8) zfmt(npcmtmax),v(nstsv,nstsv)
! allocatable arrays
complex(8), allocatable :: zrhomt(:,:)
! external functions
complex(8), external :: zdotc
allocate(zrhomt(npcmtmax,nstsv))
call holdthd(nstsv,nthd)
v(:,:)=0.d0
do is=1,nspecies
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  npc=npcmt(is)
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    do ist3=1,nstsp(is)
      if (spcore(ist3,is)) then
        do m=-ksp(ist3,is),ksp(ist3,is)-1
! generate the core wavefunction in spherical coordinates (pass in m-1/2)
          call wavefcr(.false.,lradstp,is,ia,ist3,m,npcmtmax,wfcr)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(zfmt,ist1,ist2) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
          do ist1=1,nstsv
! calculate the complex overlap density in spherical harmonics
            if (spinpol) then
              call zrho2(npc,wfcr,wfcr(:,2),wfmt(:,ias,1,ist1), &
               wfmt(:,ias,2,ist1),zfmt)
            else
              call zrho1(npc,wfcr,wfmt(:,ias,1,ist1),zfmt)
            end if
            call zfsht(nrc,nrci,zfmt,zrhomt(:,ist1))
          end do
!$OMP END DO
!$OMP DO
          do ist2=1,nstsv
            call zpotclmt(nrc,nrci,nrcmtmax,rlcmt(:,:,is),wprcmt(:,:,is), &
             zrhomt(:,ist2),zfmt)
            call zfmtwr(nrc,nrci,wrcmt(:,is),zfmt)
            do ist1=1,ist2
              v(ist1,ist2)=v(ist1,ist2)-zdotc(npc,zrhomt(:,ist1),1,zfmt,1)
            end do
          end do
!$OMP END DO
!$OMP END PARALLEL
        end do
      end if
    end do
  end do
end do
call freethd(nthd)
! set the lower triangular part of the matrix
do ist1=1,nstsv
  do ist2=1,ist1-1
    v(ist1,ist2)=conjg(v(ist2,ist1))
  end do
end do
! scale the Coulomb matrix elements in the case of a hybrid functional
if (hybrid) v(:,:)=hybridc*v(:,:)
! add to input matrix
vmat(:,:)=vmat(:,:)+v(:,:)
deallocate(zrhomt)
return

contains

pure subroutine zrho1(n,wf1,wf2,zrho)
implicit none
integer, intent(in) :: n
complex(4), intent(in) :: wf1(n),wf2(n)
complex(8), intent(out) :: zrho(n)
zrho(:)=conjg(wf1(:))*wf2(:)
end subroutine

pure subroutine zrho2(n,wf11,wf12,wf21,wf22,zrho)
implicit none
integer, intent(in) :: n
complex(4), intent(in) :: wf11(n),wf12(n),wf21(n),wf22(n)
complex(8), intent(out) :: zrho(n)
zrho(:)=conjg(wf11(:))*wf21(:)+conjg(wf12(:))*wf22(:)
end subroutine

end subroutine

