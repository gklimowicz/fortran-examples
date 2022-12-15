
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine dhmlrad
use modmain
use modphonon
implicit none
! local variables
integer is,ias
integer nr,nri,iro,i0,i1
integer l1,l2,l3,lm2
integer io,jo,ilo,jlo
complex(8) zsm
! begin loops over atoms and species
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  iro=nri+1
!---------------------------!
!     APW-APW integrals     !
!---------------------------!
  do l1=0,lmaxapw
    do io=1,apword(l1,is)
      do l3=0,lmaxapw
        do jo=1,apword(l3,is)
          do l2=0,lmaxi
            do lm2=l2**2+1,(l2+1)**2
              i1=lmmaxi*(nri-1)+lm2
              zsm=sum(apwfr(1:nri,1,io,l1,ias)*apwfr(1:nri,1,jo,l3,ias) &
               *wrmt(1:nri,is)*dvsmt(lm2:i1:lmmaxi,ias))
              i0=i1+lmmaxi
              i1=lmmaxo*(nr-iro)+i0
              zsm=zsm+sum(apwfr(iro:nr,1,io,l1,ias)*apwfr(iro:nr,1,jo,l3,ias) &
               *wrmt(iro:nr,is)*dvsmt(i0:i1:lmmaxo,ias))
              dhaa(lm2,jo,l3,io,l1,ias)=zsm
            end do
          end do
          do l2=lmaxi+1,lmaxo
            do lm2=l2**2+1,(l2+1)**2
              i0=lmmaxi*nri+lm2
              i1=lmmaxo*(nr-iro)+i0
              zsm=sum(apwfr(iro:nr,1,io,l1,ias)*apwfr(iro:nr,1,jo,l3,ias) &
               *wrmt(iro:nr,is)*dvsmt(i0:i1:lmmaxo,ias))
              dhaa(lm2,jo,l3,io,l1,ias)=zsm
            end do
          end do
        end do
      end do
    end do
  end do
!-------------------------------------!
!     local-orbital-APW integrals     !
!-------------------------------------!
  do ilo=1,nlorb(is)
    l1=lorbl(ilo,is)
    do l3=0,lmaxapw
      do io=1,apword(l3,is)
        do l2=0,lmaxi
          do lm2=l2**2+1,(l2+1)**2
            i1=lmmaxi*(nri-1)+lm2
            zsm=sum(lofr(1:nri,1,ilo,ias)*apwfr(1:nri,1,io,l3,ias) &
             *wrmt(1:nri,is)*dvsmt(lm2:i1:lmmaxi,ias))
            i0=i1+lmmaxi
            i1=lmmaxo*(nr-iro)+i0
            zsm=zsm+sum(lofr(iro:nr,1,ilo,ias)*apwfr(iro:nr,1,io,l3,ias) &
             *wrmt(iro:nr,is)*dvsmt(i0:i1:lmmaxo,ias))
            dhloa(lm2,io,l3,ilo,ias)=zsm
          end do
        end do
        do l2=lmaxi+1,lmaxo
          do lm2=l2**2+1,(l2+1)**2
            i0=lmmaxi*nri+lm2
            i1=lmmaxo*(nr-iro)+i0
            zsm=sum(lofr(iro:nr,1,ilo,ias)*apwfr(iro:nr,1,io,l3,ias) &
             *wrmt(iro:nr,is)*dvsmt(i0:i1:lmmaxo,ias))
            dhloa(lm2,io,l3,ilo,ias)=zsm
          end do
        end do
      end do
    end do
  end do
!-----------------------------------------------!
!     local-orbital-local-orbital integrals     !
!-----------------------------------------------!
  do ilo=1,nlorb(is)
    l1=lorbl(ilo,is)
    do jlo=1,nlorb(is)
      l3=lorbl(jlo,is)
      do l2=0,lmaxi
        do lm2=l2**2+1,(l2+1)**2
          i1=lmmaxi*(nri-1)+lm2
          zsm=sum(lofr(1:nri,1,ilo,ias)*lofr(1:nri,1,jlo,ias)*wrmt(1:nri,is) &
           *dvsmt(lm2:i1:lmmaxi,ias))
          i0=i1+lmmaxi
          i1=lmmaxo*(nr-iro)+i0
          zsm=zsm+sum(lofr(iro:nr,1,ilo,ias)*lofr(iro:nr,1,jlo,ias) &
           *wrmt(iro:nr,is)*dvsmt(i0:i1:lmmaxo,ias))
          dhlolo(lm2,jlo,ilo,ias)=zsm
        end do
      end do
      do l2=lmaxi+1,lmaxo
        do lm2=l2**2+1,(l2+1)**2
          i0=lmmaxi*nri+lm2
          i1=lmmaxo*(nr-iro)+i0
          zsm=sum(lofr(iro:nr,1,ilo,ias)*lofr(iro:nr,1,jlo,ias)*wrmt(iro:nr,is)&
           *dvsmt(i0:i1:lmmaxo,ias))
          dhlolo(lm2,jlo,ilo,ias)=zsm
        end do
      end do
    end do
  end do
! end loops over atoms and species
end do
end subroutine

