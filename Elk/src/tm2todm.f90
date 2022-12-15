
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: tm2todm
! !INTERFACE:
subroutine tm2todm(l,k,p,ld,wkp,dm)
! !INPUT/OUTPUT PARAMETERS:
!   l   : angular momentum quantum number (in,integer)
!   k   : angular momentum tensor moment label (in,integer)
!   p   : spin tensor moment label (in,integer)
!   ld  : leading dimension (in,integer)
!   wkp : 2-index tensor moment components (in,real(-ld:ld,-1:1))
!   dm  : real density matrix (out,real(ld,2,ld,2))
! !DESCRIPTION:
!   Calculates the real density matrix
!   $$ D=\sum_{y=-p}^p\sum_{x=-k}^k w_{xy}^{kp}\,\Gamma_{xy}^{kp} $$
!   from the real 2-index coefficients $w_{xy}^{kp}$ and the uncoupled tensor
!   moment matrices given by
!   $$ \Gamma_{xy}^{kp}(m_1\sigma_1,m_2\sigma_2)=
!    (-1)^{l-m_2+s-\sigma_2}\sqrt{(2k+1)(2p+1)}
!    \begin{pmatrix} l & k & l \\ -m_2 & x & m_1 \end{pmatrix}
!    \begin{pmatrix} s & p & s \\ -\sigma_2 & y & \sigma_1 \end{pmatrix}, $$
!   where $l$ is the angular momentum quantum number, $s=\frac{1}{2}$ and the
!   irreducible representations are labeled by $k\in\{0,\ldots,2l\}$ and
!   $p\in\{0,1\}$. The variables $x\in\{-k,\ldots, k\}$ and $y\in\{-1,0,1\}$
!   index the components in the array {\tt wkp}. These matrices are real and
!   orthonormal in the sense
!   $$ \tr\big(\Gamma_{xy}^{kp}\Gamma_{x'y'}^{k'p'}\big)=
!    \delta_{kk'}\delta_{pp'}\delta_{xx'}\delta_{yy'}. $$
!   For a detailed derivation see {\it Phys. Rev. B} {\bf 80}, 035121 (2009) and
!   {\it J. Phys.: Condens. Matter} {\bf 7}, 9947 (1995). See also the routine
!   {\tt tm3todm}.
!
! !REVISION HISTORY:
!   Created 2007 (Francesco Cricchio and Lars Nordstrom)
!   Changed normalisation and decoupled loops, January 2022 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: l,k,p,ld
real(8), intent(in) :: wkp(-ld:ld,-1:1)
real(8), intent(out) :: dm(ld,2,ld,2)
! local variables
integer ispn,jspn
integer m1,m2,n,x,y
integer lm0,lm1,lm2
real(8) t0,t1
! automatic arrays
real(8) dlm(2*l+1,2*l+1,-k:k),dsp(2,2,-p:p)
! external functions
real(8), external :: wigner3j,wigner3jf
if (l.lt.0) then
  write(*,*)
  write(*,'("Error(tm2todm): l < 0 : ",I8)') l
  write(*,*)
  stop
end if
if (k.lt.0) then
  write(*,*)
  write(*,'("Error(tm2todm): k < 0 : ",I8)') k
  write(*,*)
  stop
end if
if (k.gt.2*l) then
  write(*,*)
  write(*,'("Error(tm2todm): k > 2*l : ",2I8)') k,2*l
  write(*,*)
  stop
end if
if ((p.lt.0).or.(p.gt.1)) then
  write(*,*)
  write(*,'("Error(tm2todm): p should be 0 or 1 : ",I8)') p
  write(*,*)
  stop
end if
! calculate the angular momentum matrices
t0=sqrt(dble(2*k+1))
do x=-k,k
  dlm(:,:,x)=0.d0
  lm2=0
  do m2=-l,l
    lm2=lm2+1
    if (mod(l-m2,2).eq.0) then
      t1=t0
    else
      t1=-t0
    end if
    lm1=0
    do m1=-l,l
      lm1=lm1+1
      dlm(lm1,lm2,x)=t1*wigner3j(l,k,l,-m2,x,m1)
    end do
  end do
end do
! calculate the spin matrices
t0=sqrt(dble(2*p+1))
do y=-p,p
  dsp(:,:,y)=0.d0
  do jspn=1,2
    if (jspn.eq.1) then
      t1=t0
    else
      t1=-t0
    end if
    do ispn=1,2
      dsp(ispn,jspn,y)=t1*wigner3jf(1,2*p,1,2*jspn-3,2*y,3-2*ispn)
    end do
  end do
end do
! determine the full matrix from the Kronecker product of dlm and dsp
dm(:,:,:,:)=0.d0
lm0=l**2
n=2*l+1
do y=-p,p
  do x=-k,k
    t1=wkp(x,y)
    if (abs(t1).lt.1.d-8) cycle
    do jspn=1,2
      do lm2=1,n
        do ispn=1,2
          do lm1=1,n
            dm(lm0+lm1,ispn,lm0+lm2,jspn)=dm(lm0+lm1,ispn,lm0+lm2,jspn) &
             +t1*dlm(lm1,lm2,x)*dsp(ispn,jspn,y)
          end do
        end do
      end do
    end do
  end do
end do
end subroutine
!EOC

