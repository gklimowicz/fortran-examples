
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: rfmtctof
! !INTERFACE:
subroutine rfmtctof(rfmt)
! !USES:
use modmain
use modomp
! !INPUT/OUTPUT PARAMETERS:
!   rfmt : real muffin-tin function (in,real(npmtmax,natmtot))
! !DESCRIPTION:
!   Converts a real muffin-tin function from a coarse to a fine radial mesh by
!   using cubic spline interpolation. See {\tt rfinterp} and {\tt spline}.
!
! !REVISION HISTORY:
!   Created October 2003 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(inout) :: rfmt(npmtmax,natmtot)
! local variables
integer is,ias,lm
integer nr,nri,nro,iro
integer nrc,nrci,nrco,irco
integer i0,i1,nthd
! automatic arrays
real(8) rfmt1(npcmtmax),fi(nrcmtmax),fo(nrmtmax)
if (lradstp.eq.1) return
call holdthd(natmtot,nthd)
!$OMP PARALLEL DO DEFAULT(SHARED) &
!$OMP PRIVATE(rfmt1,fi,fo,is) &
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
  rfmt1(1:npcmt(is))=rfmt(1:npcmt(is),ias)
! interpolate up to lmaxi over entire muffin-tin
  do lm=1,lmmaxi
    i1=lmmaxi*(nrci-1)+lm
    fi(1:nrci)=rfmt1(lm:i1:lmmaxi)
    i0=i1+lmmaxi
    i1=lmmaxo*(nrc-irco)+i0
    fi(irco:nrc)=rfmt1(i0:i1:lmmaxo)
    call rfinterp(nrc,rcmt(:,is),wcrcmt(:,:,is),fi,nr,rlmt(:,1,is),fo)
    i1=lmmaxi*(nri-1)+lm
    rfmt(lm:i1:lmmaxi,ias)=fo(1:nri)
    i0=i1+lmmaxi
    i1=lmmaxo*(nr-iro)+i0
    rfmt(i0:i1:lmmaxo,ias)=fo(iro:nr)
  end do
! interpolate up to lmaxo on outer part of muffin-tin
  do lm=lmmaxi+1,lmmaxo
    i0=lmmaxi*nrci+lm
    i1=lmmaxo*(nrc-irco)+i0
    fi(irco:nrc)=rfmt1(i0:i1:lmmaxo)
    call rfinterp(nrco,rcmt(irco,is),wcrcmt(:,irco,is),fi(irco),nro, &
     rsp(iro,is),fo(iro))
    i0=lmmaxi*nri+lm
    i1=lmmaxo*(nr-iro)+i0
    rfmt(i0:i1:lmmaxo,ias)=fo(iro:nr)
  end do
end do
!$OMP END PARALLEL DO
call freethd(nthd)
end subroutine
!EOC

