
! Copyright (C) 2008 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: writetm3
! !INTERFACE:
subroutine writetm3
! !USES:
use modmain
use moddftu
use modtest
use modvars
! !DESCRIPTION:
!   Decompose the density matrix into 3-index tensor moments and write to
!   {\tt TENSMOM.OUT}. See {\it Phys. Rev. B} {\bf 80}, 035121 (2009) and
!   {\it J. Phys.: Condens. Matter} {\bf 7} 9947 (1995). See also the routines
!   {\tt tm2todm} and {\tt tm3todm}.
!
! !REVISION HISTORY:
!   Created April 2008 (F. Cricchio and L. Nordstrom)
!   Updated, December 2021 (JKD)
!EOP
!BOC
implicit none
! local variables
integer is,ia,ias,idu
integer l,p,k,r,t
real(8) ehx,sm,t0,t1
! automatic arrays
real(8) wkpr(-lmmaxdm:lmmaxdm)
complex(8) zkpr(-lmmaxdm:lmmaxdm)
complex(8) gamma(lmmaxdm,2,lmmaxdm,2)
! external functions
real(8), external :: trzhmm
wkpr(:)=0.d0
! open TENSMOM.OUT file
open(50,file='TENSMOM'//trim(filext),form='FORMATTED',action='WRITE')
write(50,'("Density matrix decomposition in coupled tensor moments")')
write(50,'("Components are in the spherical basis")')
write(50,'(" (see Phys. Rev. B. 80, 035121 (2009))")')
! loop over DFT+U entries
do idu=1,ndftu
  is=isldu(1,idu)
  l=isldu(2,idu)
! scale factor for conventional normalisation
  t0=sqrt(dble((2*l+1)*2))
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,*)
    write(50,'("Species : ",I4," (",A,"), atom : ",I4)') is,trim(spsymb(is)),ia
    write(50,'(" l = ",I1)') l
    do k=0,2*l
      do p=0,1
        do r=abs(k-p),k+p
! decompose density matrix in 3-index tensor moment components
          call dmtotm3(l,k,p,r,lmmaxdm,dmatmt(:,:,:,:,ias),wkpr)
! determine the contribution to muffin-tin Hartree + exchange energy
          call tm3todm(l,k,p,r,lmmaxdm,wkpr,gamma)
          ehx=0.5d0*trzhmm(lmmaxdm*2,gamma,vmatmt(:,:,:,:,ias))
! write out tensor moment components and vector magnitude
          write(50,*)
          write(50,'("  k = ",I1,", p = ",I1,", r = ",I1)') k,p,r
          if (tm3old) then
! old complex convention
            call tm3rtoz(l,k,p,r,lmmaxdm,wkpr,zkpr)
            do t=-r,r
              write(50,'("   t = ",I2," : ",2F14.8)') t,zkpr(t)
            end do
          else
! new real convention
            sm=0.d0
            do t=-r,r
              t1=t0*wkpr(t)
              write(50,'("   t = ",I2," : ",F14.8)') t,t1
              sm=sm+t1**2
            end do
            sm=sqrt(sm)
            write(50,'("  magnitude : ",F14.8)') sm
          end if
          write(50,'("  Hartree + exchange energy : ",F14.8)') ehx
! write to VARIABLES.OUT if required, but only on the last iteration
          if (wrtvars.and.tlast) then
            call writevars('wkpr',n1=is,n2=ia,n3=l,n4=k,n5=p,n6=r,nv=2*r+1, &
             rva=wkpr(-r:r))
          end if
        end do
      end do
    end do
! end loop over atoms and species
  end do
end do
close(50)
! write last entry of tensor moment components to test file if required
if (test) then
  call writetest(820,'Coupled tensor moments',nv=size(wkpr),tol=5.d-4,rva=wkpr)
end if
end subroutine
!EOC

