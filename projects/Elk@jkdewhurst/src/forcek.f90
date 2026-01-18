
! Copyright (C) 2002-2006 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: forcek
! !INTERFACE:
subroutine forcek(ik)
! !USES:
use modmain
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   ik : reduced k-point number (in,integer)
! !DESCRIPTION:
!   Computes the {\bf k}-dependent contribution to the incomplete basis set
!   (IBS) force. See the calling routine {\tt force} for a full description.
!
! !REVISION HISTORY:
!   Created June 2006 (JKD)
!   Updated for spin-spiral case, May 2007 (Francesco Cricchio and JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: ik
! local variables
integer ispn0,ispn1,ispn,jspn
integer n,nm,nm2,is,ias,ist,jst
integer j1,j2,j3,ig,i,j,k,l
integer nthd
real(8) v1,v2,v3,sm,t1
complex(8) z1,z2
! automatic arrays
real(8) evalfv(nstfv,nspnfv)
complex(8) vh(nmatmax),vo(nmatmax)
complex(8) ffv(nstfv,nstfv),y(nstfv)
! allocatable arrays
integer, allocatable :: ijg(:)
real(8), allocatable :: dp(:)
complex(8), allocatable :: apwalm(:,:,:,:)
complex(8), allocatable :: evecfv(:,:,:),evecsv(:,:)
complex(8), allocatable :: h(:),o(:),dlh(:),dlo(:)
! external functions
complex(8), external :: zdotc
nm2=nmatmax**2
! allocate local arrays
allocate(ijg(nm2),dp(nm2))
allocate(apwalm(ngkmax,apwordmax,lmmaxapw,natmtot))
allocate(evecfv(nmatmax,nstfv,nspnfv))
allocate(h(nm2),o(nm2),dlh(nm2),dlo(nm2))
! get the eigenvalues/vectors from file
call getevalfv(filext,ik,vkl(:,ik),evalfv)
call getevecfv(filext,ik,vkl(:,ik),vgkl(:,:,:,ik),evecfv)
if (tevecsv) then
  allocate(evecsv(nstsv,nstsv))
  call getevecsv(filext,ik,vkl(:,ik),evecsv)
end if
! loop over first-variational spin components
do jspn=1,nspnfv
  if (spinsprl) then
    ispn0=jspn; ispn1=jspn
  else
    ispn0=1; ispn1=nspinor
  end if
  n=ngk(jspn,ik)
  nm=nmat(jspn,ik)
  do j=1,n
    k=(j-1)*nm
    ig=igkig(j,jspn,ik)
    j1=ivg(1,ig); j2=ivg(2,ig); j3=ivg(3,ig)
    v1=0.5d0*vgkc(1,j,jspn,ik)
    v2=0.5d0*vgkc(2,j,jspn,ik)
    v3=0.5d0*vgkc(3,j,jspn,ik)
    do i=1,j
      k=k+1
      ig=igkig(i,jspn,ik)
      ijg(k)=ivgig(ivg(1,ig)-j1,ivg(2,ig)-j2,ivg(3,ig)-j3)
      dp(k)=vgkc(1,i,jspn,ik)*v1 &
           +vgkc(2,i,jspn,ik)*v2 &
           +vgkc(3,i,jspn,ik)*v3
    end do
  end do
! find the matching coefficients
  call match(n,vgkc(:,:,jspn,ik),gkc(:,jspn,ik),sfacgk(:,:,jspn,ik),apwalm)
! zero the local-orbital-local-orbital contribution
  do j=n+1,nm
    k=(j-1)*nm+n
    do i=n+1,j
      k=k+1
      dlh(k)=0.d0
      dlo(k)=0.d0
    end do
  end do
! loop over species and atoms
  do ias=1,natmtot
    is=idxis(ias)
! Hamiltonian and overlap matrices
    h(:)=0.d0
    call hmlaa(.false.,is,ias,n,apwalm(:,:,:,ias),nm,h)
    call hmlalo(is,ias,n,apwalm(:,:,:,ias),nm,h)
    o(:)=0.d0
    call olpaa(.false.,is,n,apwalm(:,:,:,ias),nm,o)
    call olpalo(is,ias,n,apwalm(:,:,:,ias),nm,o)
! loop over Cartesian directions
    do l=1,3
! APW-APW contribution
      do j=1,n
        k=(j-1)*nm
        do i=1,j
          k=k+1
          ig=ijg(k)
          t1=vgc(l,ig)
          z1=-ffacg(ig,is)*conjg(sfacg(ig,ias))
          z2=t1*(dp(k)*z1+h(k))
          dlh(k)=cmplx(-aimag(z2),dble(z2),8)
          z2=t1*(z1+o(k))
          dlo(k)=cmplx(-aimag(z2),dble(z2),8)
        end do
      end do
! APW-local-orbital contribution
      do j=n+1,nm
        k=(j-1)*nm
        do i=1,n
          k=k+1
          t1=vgkc(l,i,jspn,ik)
          z1=t1*h(k)
          dlh(k)=cmplx(-aimag(z1),dble(z1),8)
          z1=t1*o(k)
          dlo(k)=cmplx(-aimag(z1),dble(z1),8)
        end do
      end do
! compute the force matrix elements in the first-variational basis
      call holdthd(nstfv,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(vh,vo,t1,ist,z1,z2) &
!$OMP NUM_THREADS(nthd)
      do jst=1,nstfv
        call zhemv('U',nm,zone,dlh,nm,evecfv(:,jst,jspn),1,zzero,vh,1)
        call zhemv('U',nm,zone,dlo,nm,evecfv(:,jst,jspn),1,zzero,vo,1)
        t1=evalfv(jst,jspn)
        do ist=1,nstfv
          z1=zdotc(nm,evecfv(:,ist,jspn),1,vh,1)
          z2=zdotc(nm,evecfv(:,ist,jspn),1,vo,1)
          ffv(ist,jst)=z1-t1*z2
        end do
      end do
!$OMP END PARALLEL DO
      call freethd(nthd)
! compute the force using the second-variational coefficients if required
      sm=0.d0
      if (tevecsv) then
! spin-polarised case
        do j=1,nstsv
          do ispn=ispn0,ispn1
            i=(ispn-1)*nstfv+1
            call zgemv('N',nstfv,nstfv,zone,ffv,nstfv,evecsv(i,j),1,zzero,y,1)
            z1=zdotc(nstfv,evecsv(i,j),1,y,1)
            sm=sm+occsv(j,ik)*dble(z1)
          end do
        end do
      else
! spin-unpolarised case
        do j=1,nstsv
          sm=sm+occsv(j,ik)*dble(ffv(j,j))
        end do
      end if
!$OMP ATOMIC
      forceibs(l,ias)=forceibs(l,ias)+wkpt(ik)*sm
! end loop over Cartesian components
    end do
! end loop over atoms and species
  end do
! end loop over first-variational spins
end do
deallocate(ijg,dp,apwalm,evecfv)
deallocate(h,o,dlh,dlo)
if (tevecsv) deallocate(evecsv)
end subroutine
!EOC

