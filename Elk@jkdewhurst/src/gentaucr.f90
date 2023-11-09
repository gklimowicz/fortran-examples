
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine gentaucr
use modmain
use modomp
implicit none
! local variables
integer ist,ispn,jspn
integer is,ia,ias,nthd
integer nr,nri,np,i,m
! automatic arrays
real(8) rfmt(npmtmax)
! allocatable arrays
complex(4), allocatable :: wfcr(:,:)
complex(8), allocatable :: zfmt(:),gzfmt(:,:)
taucr(:,:,:)=0.d0
call holdthd(natmtot,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP PRIVATE(wfcr,zfmt,gzfmt) &
!$OMP PRIVATE(is,ia,nr,nri,np) &
!$OMP PRIVATE(ist,m,ispn,jspn,i) &
!$OMP NUM_THREADS(nthd)
allocate(wfcr(npmtmax,2),zfmt(npmtmax),gzfmt(npmtmax,3))
!$OMP DO
do ias=1,natmtot
  is=idxis(ias)
  ia=idxia(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  np=npmt(is)
  do ist=1,nstsp(is)
    if (spcore(ist,is)) then
      do m=-ksp(ist,is),ksp(ist,is)-1
! generate the core wavefunction in spherical harmonics (pass in m-1/2)
        call wavefcr(.true.,1,is,ia,ist,m,npmtmax,wfcr)
        do ispn=1,2
          if (spinpol) then
            jspn=ispn
          else
            jspn=1
          end if
! compute the gradient of the wavefunction
          zfmt(1:np)=wfcr(1:np,ispn)
          call gradzfmt(nr,nri,rlmt(:,-1,is),wcrmt(:,:,is),zfmt,npmtmax,gzfmt)
          do i=1,3
! convert gradient to spherical coordinates
            call zbsht(nr,nri,gzfmt(:,i),zfmt)
! add to total in muffin-tin
            taucr(1:np,ias,jspn)=taucr(1:np,ias,jspn) &
             +0.5d0*(dble(zfmt(1:np))**2+aimag(zfmt(1:np))**2)
          end do
        end do
      end do
    end if
  end do
end do
!$OMP END DO
deallocate(wfcr,zfmt,gzfmt)
!$OMP END PARALLEL
call freethd(nthd)
! convert core tau to spherical harmonics
do ispn=1,nspinor
  do ias=1,natmtot
    is=idxis(ias)
    rfmt(1:npmt(is))=taucr(1:npmt(is),ias,ispn)
    call rfsht(nrmt(is),nrmti(is),rfmt,taucr(:,ias,ispn))
  end do
end do
end subroutine

