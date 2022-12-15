
! Copyright (C) 2002-2009 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: gradzfmt
! !INTERFACE:
subroutine gradzfmt(nr,nri,ri,wcr,zfmt,ld,gzfmt)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   nr    : number of radial mesh points (in,integer)
!   nri   : number of points on inner part of muffin-tin (in,integer)
!   ri    : 1/r on the radial mesh (in,real(nr))
!   wcr   : weights for spline coefficients on radial mesh (in,real(12,nr))
!   zfmt  : complex muffin-tin function (in,complex(*))
!   ld    : leading dimension (in,integer)
!   gzfmt : gradient of zfmt (out,complex(ld,3))
! !DESCRIPTION:
!   Calculates the gradient of a complex muffin-tin function. In other words,
!   given the spherical harmonic expansion coefficients, $f_{lm}(r)$, of a
!   function $f({\bf r})$, the routine returns ${\bf F}_{lm}$ where
!   $$ \sum_{lm}{\bf F}_{lm}(r)Y_{lm}(\hat{\bf r})=\nabla f({\bf r}). $$
!   This is done using the gradient formula (see, for example, V. Devanathan,
!   {\em Angular Momentum Techniques In Quantum Mechanics})
!   \begin{align*}
!    \nabla f_{lm}(r)Y_{lm}(\hat{\bf r})&=-\sqrt{\frac{l+1}{2l+1}}
!    \left(\frac{d}{dr}-\frac{l}{r}\right)f_{lm}(r)
!    {\bf Y}_{lm}^{l+1}(\hat{\bf r})\\
!    &+\sqrt{\frac{l}{2l+1}}\left(\frac{d}{dr}+\frac{l+1}{r}\right)f_{lm}(r)
!    {\bf Y}_{lm}^{l-1}(\hat{\bf r}),
!   \end{align*}
!   where the vector spherical harmonics are determined from Clebsch-Gordan
!   coefficients as follows:
!   $$ {\bf Y}_{lm}^{l'}(\hat{\bf r})=\sum_{m'\mu}
!    \begin{bmatrix} l' & 1 & l \\ m' & \mu & m \end{bmatrix}
!    Y_{lm}(\hat{\bf r})\hat{\bf e}^{\mu} $$
!   and the (contravariant) spherical unit vectors are given by
!   $$ \hat{\bf e}_{+1}=-\frac{\hat{\bf x}+i\hat{\bf y}}{\sqrt{2}},
!    \qquad\hat{\bf e}_0=\hat{\bf z},\qquad
!    \hat{\bf e}_{-1}=\frac{\hat{\bf x}-i\hat{\bf y}}{\sqrt{2}}. $$
!
! !REVISION HISTORY:
!   Rewritten May 2009 (JKD)
!   Modified, February 2020 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: nr,nri
real(8), intent(in) :: ri(nr),wcr(12,nr)
complex(8), intent(in) :: zfmt(*)
integer, intent(in) :: ld
complex(8), intent(out) :: gzfmt(ld,3)
! local variables
integer nro,iro,ir,mu
integer np,npi,i,i1,j
integer l,m,lm,lm1
! real constant 1/sqrt(2)
real(8), parameter :: c1=0.7071067811865475244d0
real(8) t1,t2,t3
complex(8) z1
! automatic arrays
real(8) f1(nr),f2(nr),g1(nr),g2(nr)
complex(8) drmt(ld)
! external functions
real(8), external :: clebgor
nro=nr-nri
iro=nri+1
npi=lmmaxi*nri
np=npi+lmmaxo*nro
!----------------------------------------!
!     compute the radial derivatives     !
!----------------------------------------!
do lm=1,lmmaxi
  i=lm
  do ir=1,nri
    f1(ir)=dble(zfmt(i))
    f2(ir)=aimag(zfmt(i))
    i=i+lmmaxi
  end do
  do ir=iro,nr
    f1(ir)=dble(zfmt(i))
    f2(ir)=aimag(zfmt(i))
    i=i+lmmaxo
  end do
  call splined(nr,wcr,f1,g1)
  call splined(nr,wcr,f2,g2)
  i=lm
  do ir=1,nri
    drmt(i)=cmplx(g1(ir),g2(ir),8)
    i=i+lmmaxi
  end do
  do ir=iro,nr
    drmt(i)=cmplx(g1(ir),g2(ir),8)
    i=i+lmmaxo
  end do
end do
do lm=lmmaxi+1,lmmaxo
  i=npi+lm
  do ir=iro,nr
    f1(ir)=dble(zfmt(i))
    f2(ir)=aimag(zfmt(i))
    i=i+lmmaxo
  end do
  call splined(nro,wcr(1,iro),f1(iro),g1(iro))
  call splined(nro,wcr(1,iro),f2(iro),g2(iro))
  i=npi+lm
  do ir=iro,nr
    drmt(i)=cmplx(g1(ir),g2(ir),8)
    i=i+lmmaxo
  end do
end do
!-----------------------------------------------------!
!     compute the gradient in the spherical basis     !
!-----------------------------------------------------!
! zero the gradient array
gzfmt(1:np,:)=0.d0
! inner part of muffin-tin
lm=0
do l=0,lmaxi
  t1=-sqrt(dble(l+1)/dble(2*l+1))
  if (l.gt.0) then
    t2=sqrt(dble(l)/dble(2*l+1))
  else
    t2=0.d0
  end if
  do m=-l,l
    lm=lm+1
    j=1
    do mu=-1,1
      if (mu.eq.0) j=3
      if (mu.eq.1) j=2
      if (l+1.le.lmaxi) then
! index to (l,m) is l*(l+1)+m+1, therefore index to (l+1,m-mu) is
        lm1=(l+1)*(l+2)+(m-mu)+1
        t3=t1*clebgor(l+1,1,l,m-mu,mu,m)
        i=lm; i1=lm1
        do ir=1,nri
          gzfmt(i1,j)=gzfmt(i1,j)+t3*(drmt(i)-dble(l)*ri(ir)*zfmt(i))
          i=i+lmmaxi; i1=i1+lmmaxi
        end do
      end if
      if (abs(m-mu).le.l-1) then
! index to (l-1,m-mu)
        lm1=(l-1)*l+(m-mu)+1
        t3=t2*clebgor(l-1,1,l,m-mu,mu,m)
        i=lm; i1=lm1
        do ir=1,nri
          gzfmt(i1,j)=gzfmt(i1,j)+t3*(drmt(i)+dble(l+1)*ri(ir)*zfmt(i))
          i=i+lmmaxi; i1=i1+lmmaxi
        end do
      end if
    end do
  end do
end do
! outer part of muffin-tin
lm=0
do l=0,lmaxo
  t1=-sqrt(dble(l+1)/dble(2*l+1))
  if (l.gt.0) then
    t2=sqrt(dble(l)/dble(2*l+1))
  else
    t2=0.d0
  end if
  do m=-l,l
    lm=lm+1
    j=1
    do mu=-1,1
      if (mu.eq.0) j=3
      if (mu.eq.1) j=2
      if (l+1.le.lmaxo) then
        lm1=(l+1)*(l+2)+(m-mu)+1
        t3=t1*clebgor(l+1,1,l,m-mu,mu,m)
        i=npi+lm; i1=npi+lm1
        do ir=iro,nr
          gzfmt(i1,j)=gzfmt(i1,j)+t3*(drmt(i)-dble(l)*ri(ir)*zfmt(i))
          i=i+lmmaxo; i1=i1+lmmaxo
        end do
      end if
      if (abs(m-mu).le.l-1) then
        lm1=(l-1)*l+(m-mu)+1
        t3=t2*clebgor(l-1,1,l,m-mu,mu,m)
        i=npi+lm; i1=npi+lm1
        do ir=iro,nr
          gzfmt(i1,j)=gzfmt(i1,j)+t3*(drmt(i)+dble(l+1)*ri(ir)*zfmt(i))
          i=i+lmmaxo; i1=i1+lmmaxo
        end do
      end if
    end do
  end do
end do
!---------------------------------------------------!
!     convert from spherical to Cartesian basis     !
!---------------------------------------------------!
! note that the gradient transforms as a covariant vector, i.e. y -> -y
i=0
do ir=1,nri
  do lm=1,lmmaxi
    i=i+1
    z1=gzfmt(i,1)
    gzfmt(i,1)=c1*(z1-gzfmt(i,2))
    z1=c1*(z1+gzfmt(i,2))
    gzfmt(i,2)=cmplx(aimag(z1),-dble(z1),8)
  end do
end do
do ir=iro,nr
  do lm=1,lmmaxo
    i=i+1
    z1=gzfmt(i,1)
    gzfmt(i,1)=c1*(z1-gzfmt(i,2))
    z1=c1*(z1+gzfmt(i,2))
    gzfmt(i,2)=cmplx(aimag(z1),-dble(z1),8)
  end do
end do
end subroutine
!EOC

