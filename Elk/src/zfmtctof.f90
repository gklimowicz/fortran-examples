
! Copyright (C) 2013 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine zfmtctof(zfmt)
use modmain
use modomp
implicit none
! arguments
complex(8), intent(inout) :: zfmt(npmtmax,natmtot)
! local variables
integer is,ias,lm
integer nr,nri,nro,iro
integer nrc,nrci,nrco,irco
integer i0,i1,nthd
! automatic arrays
real(8) fi1(nrcmtmax),fi2(nrcmtmax)
real(8) fo1(nrmtmax),fo2(nrmtmax)
complex(8) zfmt1(npcmtmax)
if (lradstp.eq.1) return
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(zfmt1,fi1,fi2,fo1,fo2,is) &
!$OMP PRIVATE(nr,nri,nro,iro) &
!$OMP PRIVATE(nrc,nrci,nrco,irco) &
!$OMP PRIVATE(lm,i0,i1) &
!$OMP NUM_THREADS(nthd)
do ias=1,natmtot
  is=idxis(ias)
  nr=nrmt(is)
  nri=nrmti(is)
  nro=nr-nri
  iro=nri+1
  nrc=nrcmt(is)
  nrci=nrcmti(is)
  nrco=nrc-nrci
  irco=nrci+1
! copy the input function
  zfmt1(1:npcmt(is))=zfmt(1:npcmt(is),ias)
! interpolate up to lmaxi over entire muffin-tin
  do lm=1,lmmaxi
    i1=lmmaxi*(nrci-1)+lm
    fi1(1:nrci)=dble(zfmt1(lm:i1:lmmaxi))
    fi2(1:nrci)=aimag(zfmt1(lm:i1:lmmaxi))
    i0=i1+lmmaxi
    i1=lmmaxo*(nrc-irco)+i0
    fi1(irco:nrc)=dble(zfmt1(i0:i1:lmmaxo))
    fi2(irco:nrc)=aimag(zfmt1(i0:i1:lmmaxo))
    call rfinterp(nrc,rcmt(:,is),wcrcmt(:,:,is),fi1,nr,rlmt(:,1,is),fo1)
    call rfinterp(nrc,rcmt(:,is),wcrcmt(:,:,is),fi2,nr,rlmt(:,1,is),fo2)
    i1=lmmaxi*(nri-1)+lm
    zfmt(lm:i1:lmmaxi,ias)=cmplx(fo1(1:nri),fo2(1:nri),8)
    i0=i1+lmmaxi
    i1=lmmaxo*(nr-iro)+i0
    zfmt(i0:i1:lmmaxo,ias)=cmplx(fo1(iro:nr),fo2(iro:nr),8)
  end do
! interpolate up to lmaxo on outer part of muffin-tin
  do lm=lmmaxi+1,lmmaxo
    i0=lmmaxi*nrci+lm
    i1=lmmaxo*(nrc-irco)+i0
    fi1(irco:nrc)=dble(zfmt1(i0:i1:lmmaxo))
    fi2(irco:nrc)=aimag(zfmt1(i0:i1:lmmaxo))
    call rfinterp(nrco,rcmt(irco,is),wcrcmt(:,irco,is),fi1(irco),nro, &
     rsp(iro,is),fo1(iro))
    call rfinterp(nrco,rcmt(irco,is),wcrcmt(:,irco,is),fi2(irco),nro, &
     rsp(iro,is),fo2(iro))
    i0=lmmaxi*nri+lm
    i1=lmmaxo*(nr-iro)+i0
    zfmt(i0:i1:lmmaxo,ias)=cmplx(fo1(iro:nr),fo2(iro:nr),8)
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine

