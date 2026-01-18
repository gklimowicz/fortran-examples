
! Copyright (C) 2007 F. Bultmark, F. Cricchio, L. Nordstrom and J. K. Dewhurst.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: tm3todm
! !INTERFACE:
subroutine tm3todm(l,k,p,r,ld,wkpr,dm)
! !INPUT/OUTPUT PARAMETERS:
!   l    : angular momentum quantum number (in,integer)
!   k    : k-index of tensor moment (in,integer)
!   p    : p-index of tensor moment (in,integer)
!   r    : r-index of tensor moment (in,integer)
!   ld   : leading dimension (in,integer)
!   wkpr : 3-index tensor moment components (in,real(-ld:ld))
!   dm   : complex Hermitian density matrix (out,complex(ld,2,ld,2))
! !DESCRIPTION:
!   The 3-index coupled tensor moment matrices are given by
!   $$ \Gamma_t^{kpr}=
!    \sqrt{2r+1}\sum_{x=-k}^k\sum_{y=-p}^p
!    \begin{pmatrix} k & r & p \\ -x & t & -y \end{pmatrix}
!    \Gamma_{xy}^{kp}, $$
!   where the irreducible representations are labeled by $k\in\{0,\ldots,2l\}$,
!   $p\in\{0,1\}$, $r\in\{|k-p|,\ldots,k+p\}$ and $\Gamma_{xy}^{kp}$ are the
!   uncoupled tensor moments (note that the phase $(-1)^{x+y}$ in the original
!   formula has been removed because of the Wigner $3j$ condition $x+y=t$). The
!   coupled tensor moment matrices are real and orthonormal in the sense
!   $$ \tr\big(\Gamma_t^{kpr}\Gamma_{t'}^{k'p'r'}\big)=
!    \delta_{kk'}\delta_{pp'}\delta_{rr'}\delta_{tt'}. $$
!   It can also be shown that the matrices are complete, thus any general
!   complex matrix $D$ of dimension $2(2l+1)$ can be expanded as
!   $$ D=\sum_{k=0}^{2l}\sum_{p=0}^1\sum_{r=|k-p|}^{k+p}\sum_{t=-r}^r
!    z_t^{kpr}\Gamma_t^{kpr} $$
!   where $z_t^{kpr}$ are complex numbers. Likewise, any real matrix can be
!   expanded in real tensor moments $w_t^{kpr}$. Using the the symmetry
!   properties of the Wigner $3j$-symbols, one can show that the transpose
!   $$ \big(\Gamma_t^{kpr}\big)^t=(-1)^{k+p+r+t}\,\Gamma_{-t}^{kpr} $$
!   and thus both the symmetric and antisymmetric parts of $\Gamma_t^{kpr}$
!   transform under rotation within the same irreducible representation.
!   Consequently, any complex Hermitian matrix $D$ can be written as
!   $$ D=\sum_{k,p,r,t} w_t^{kpr}\big[(\Gamma_t^{kpr})_{\rm S}
!    +i(\Gamma_t^{kpr})_{\rm A}\big], $$
!   where the subscripts S and A refer to the symmetric and antisymmetric parts
!   of the matrix, respectively. This routine generates the Hermitian density
!   matrix $D$ as described above from the real tensor moments $w_t^{kpr}$. For
!   a detailed derivation see {\it Phys. Rev. B} {\bf 80}, 035121 (2009) and
!   {\it J. Phys.: Condens. Matter} {\bf 7}, 9947 (1995). See also the routines
!   {\tt tm2todm} and {\tt tm3rtoz}.
!
! !REVISION HISTORY:
!   Created 2007 (Francesco Cricchio and Lars Nordstrom)
!   Changed normalisation, made the moments real and the matrix Hermitian,
!    January 2022 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: l,k,p,r,ld
real(8), intent(in) :: wkpr(-ld:ld)
complex(8), intent(out) :: dm(ld,2,ld,2)
! local variables
integer x,y,t
real(8) t0,t1
! automatic arrays
real(8) wkp(-ld:ld,-1:1),dmr(ld,2,ld,2)
! external functions
real(8), external :: wigner3j
if (l.lt.0) then
  write(*,*)
  write(*,'("Error(tm3todm): l < 0 : ",I8)') l
  write(*,*)
  stop
end if
if (k.lt.0) then
  write(*,*)
  write(*,'("Error(tm3todm): k < 0 : ",I8)') k
  write(*,*)
  stop
end if
if (k.gt.2*l) then
  write(*,*)
  write(*,'("Error(tm3todm): k > 2*l : ",2I8)') k,2*l
  write(*,*)
  stop
end if
if ((p.lt.0).or.(p.gt.1)) then
  write(*,*)
  write(*,'("Error(tm3todm): p should be 0 or 1 : ",I8)') p
  write(*,*)
  stop
end if
if (r.lt.abs(k-p)) then
  write(*,*)
  write(*,'("Error(tm3todm): r < |k-p| : ",2I8)') r,abs(k-p)
  write(*,*)
  stop
end if
if (r.gt.(k+p)) then
  write(*,*)
  write(*,'("Error(tm3todm): r > k+p : ",2I8)') r,k+p
  write(*,*)
  stop
end if
! compute 2-index tensor moment from 3-index tensor moment
wkp(:,:)=0.d0
t0=sqrt(dble(2*r+1))
do t=-r,r
  t1=wkpr(t)
  if (abs(t1).lt.1.d-8) cycle
  t1=t0*t1
  do x=-k,k
    do y=-p,p
      wkp(x,y)=wkp(x,y)+t1*wigner3j(k,r,p,-x,t,-y)
    end do
  end do
end do
! compute the real matrix from the 2-index tensor moment
call tm2todm(l,k,p,ld,wkp,dmr)
! convert to complex Hermitian matrix
call dmrtoz(ld*2,dmr,dm)
return

contains

pure subroutine dmrtoz(n,dmr,dm)
implicit none
! arguments
integer, intent(in) :: n
real(8), intent(in) :: dmr(n,n)
complex(8), intent(out) :: dm(n,n)
! local variables
integer i,j
real(8) a,b
do j=1,n
  do i=1,j-1
! symmetric part
    a=0.5d0*(dmr(i,j)+dmr(j,i))
! antisymmetric part
    b=0.5d0*(dmr(i,j)-dmr(j,i))
    dm(i,j)=cmplx(a,b,8)
    dm(j,i)=cmplx(a,-b,8)
  end do
  dm(j,j)=dmr(j,j)
end do
end subroutine

end subroutine
!EOC

