
! Copyright (C) 2008 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

!BOP
! !ROUTINE: ylmroty
! !INTERFACE:
subroutine ylmroty(beta,lmax,ld,dy)
! !INPUT/OUTPUT PARAMETERS:
!   beta : rotation angle about y-axis (in,real)
!   lmax : maximum angular momentum (in,integer)
!   ld   : leading dimension (in,integer)
!   dy   : rotation matrix for complex spherical harmonics (out,real(ld,*))
! !DESCRIPTION:
!   Returns the rotation matrix in the basis of complex spherical harmonics for
!   a rotation of angle $\beta$ about the $y$-axis. This matrix is real and is
!   given by the formula
!   \begin{align*}
!    d^l_{m_1m_2}(\beta)=&[(l+m_1)!(l-m_1)!(l+m_2)!(l-m_2)!]^{1/2}\\
!    &\times\sum_k(-1)^k\frac{\left(\cos\frac{\beta}{2}\right)^{2(l-k)-m_2+m_1}
!    \left(\sin\frac{\beta}{2}\right)^{2k+m_2-m_1}}
!    {k!(l+m_1-k)!(l-m_2-k)!(m_2-m_1+k)!},
!   \end{align*}
!   where $k$ runs through all integer values for which the factorials exist.
!
! !REVISION HISTORY:
!   Created December 2008 (JKD)
!EOP
!BOC
implicit none
! arguments
real(8), intent(in) :: beta
integer, intent(in) :: lmax,ld
real(8), intent(out) :: dy(ld,*)
! local variables
integer j,k,l,m1,m2,lm1,lm2
real(8) cb,sb,sm,t1,t2
! external functions
real(8), external :: factn
t1=0.5d0*beta
cb=cos(t1)
sb=sin(t1)
lm1=0
do l=0,lmax
! generate rotation operator for m-components of current l
  do m1=-l,l
    lm1=lm1+1
    t1=factn(l+m1)*factn(l-m1)
    lm2=l**2
    do m2=-l,l
      lm2=lm2+1
      sm=0.d0
      do k=0,min(l+m1,l-m2)
        if (((l+m1-k).ge.0).and.((l-m2-k).ge.0).and.((m2-m1+k).ge.0)) then
          j=2*(l-k)+m1-m2
          if (j.eq.0) then
            t2=1.d0
          else
            t2=cb**j
          end if
          j=2*k+m2-m1
          if (j.ne.0) t2=t2*sb**j
          t2=t2/(factn(k)*factn(l+m1-k)*factn(l-m2-k)*factn(m2-m1+k))
          if (mod(k,2).ne.0) t2=-t2
          sm=sm+t2
        end if
      end do
      dy(lm1,lm2)=sqrt(t1*factn(l+m2)*factn(l-m2))*sm
    end do
  end do
end do
end subroutine
!EOC

