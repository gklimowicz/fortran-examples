
! Copyright (C) 2009 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: grad2rfmt
! !INTERFACE:
subroutine grad2rfmt(nr,nri,ri,ri2,wcr,rfmt,g2rfmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr     : number of radial mesh points (in,integer)
!   nri    : number of points on the inner part of the muffin-tin (in,integer)
!   ri     : 1/r on the radial mesh (in,real(nr))
!   ri2    : 1/r^2 on the radial mesh (in,real(nr))
!   wcr    : weights for spline coefficients on radial mesh (in,real(12,nr))
!   rfmt   : real muffin-tin function (in,real(*))
!   g2rfmt : laplacian of the input function (out,real(*))
! !DESCRIPTION:
!   Calculates the Laplacian of a real muffin-tin function. In other words,
!   given the real spherical harmonic expansion coefficients $f_{lm}(r)$ of a
!   function $f({\bf r})$, the routine returns
!   $$ F_{lm}(r)=\frac{1}{r}\frac{\partial^2}{\partial r^2}\big(rf_{lm}(r)\big)
!    -\frac{l(l+1)}{r^2}f_{lm}(r) $$
!   which yields
!   $$ \nabla^2 f({\bf r})=\sum_{lm} F_{lm}(r)R_{lm}(\hat{\bf r}), $$
!   where $R_{lm}$ is a real spherical harmonic function.
!
! !REVISION HISTORY:
!   Created July 2009 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: ri(nr),ri2(nr),wcr(12,nr),rfmt(*)
real(8), intent(out) :: g2rfmt(*)
! local variables
integer nro,iro,npi
integer l,lm,i1,j0,j1
real(8) t1
! automatic arrays
real(8) fr(nr),cf(3,nr)
nro=nr-nri
iro=nri+1
npi=lmmaxi*nri
do l=0,lmaxi
  t1=-dble(l*(l+1))
  do lm=l**2+1,(l+1)**2
! use a cubic spline to compute radial derivatives
    i1=lmmaxi*(nri-1)+lm
    j0=i1+lmmaxi
    j1=lmmaxo*(nr-iro)+j0
    fr(1:nri)=rfmt(lm:i1:lmmaxi)
    fr(iro:nr)=rfmt(j0:j1:lmmaxo)
    call splinew(nr,wcr,fr,cf)
! apply Laplacian
    g2rfmt(lm:i1:lmmaxi)=2.d0*(ri(1:nri)*cf(1,1:nri)+cf(2,1:nri)) &
     +t1*ri2(1:nri)*rfmt(lm:i1:lmmaxi)
    g2rfmt(j0:j1:lmmaxo)=2.d0*(ri(iro:nr)*cf(1,iro:nr)+cf(2,iro:nr)) &
     +t1*ri2(iro:nr)*rfmt(j0:j1:lmmaxo)
  end do
end do
do l=lmaxi+1,lmaxo
  t1=-dble(l*(l+1))
  do lm=l**2+1,(l+1)**2
    j0=lmmaxi*nri+lm
    j1=lmmaxo*(nr-iro)+j0
    fr(iro:nr)=rfmt(j0:j1:lmmaxo)
    call splinew(nro,wcr(:,iro),fr(iro),cf(1,iro))
    g2rfmt(j0:j1:lmmaxo)=2.d0*(ri(iro:nr)*cf(1,iro:nr)+cf(2,iro:nr)) &
     +t1*ri2(iro:nr)*rfmt(j0:j1:lmmaxo)
  end do
end do
! apply smoothing if required
call rfmtsm(msmooth,nr,nri,g2rfmt)
end subroutine
!EOC

