
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine mixbroyden(iscl,n,msd,alpha,w0,nu,mu,f,df,u,a,d)
use modomp
implicit none
! arguments
integer, intent(in) :: iscl,n,msd
real(8), intent(in) :: alpha,w0
real(8), intent(inout) :: nu(n),mu(n,2)
real(8), intent(inout) :: f(n,2),df(n,msd)
real(8), intent(inout) :: u(n,msd),a(msd,msd)
real(8), intent(out) :: d
! local variables
integer jc,kp,kc
integer k,l,m
integer info,nthd
real(8) t1
! automatic arrays
integer ipiv(msd)
real(8) c(msd),beta(msd,msd),gamma(msd),work(msd)
if (n.lt.1) then
  write(*,*)
  write(*,'("Error(mixbroyden): n < 1 : ",I8)') n
  write(*,*)
  stop
end if
if (msd.lt.2) then
  write(*,*)
  write(*,'("Error(mixbroyden): msd < 2 : ",I8)') msd
  write(*,*)
  stop
end if
! initialise mixer
if (iscl.le.0) then
  mu(:,1)=nu(:)
  mu(:,2)=nu(:)
  f(:,1)=0.d0
  df(:,1)=0.d0
  u(:,1)=0.d0
  a(:,:)=0.d0
  d=1.d0
  return
end if
! current subspace dimension
m=min(iscl+1,msd)
! current index modulo m
jc=mod(iscl,m)+1
! previous index modulo 2
kp=mod(iscl-1,2)+1
! current index modulo 2
kc=mod(iscl,2)+1
f(:,kc)=nu(:)-mu(:,kp)
d=sum(f(:,kc)**2)
d=sqrt(d/dble(n))
df(:,jc)=f(:,kc)-f(:,kp)
t1=norm2(df(:,jc))
if (t1.gt.1.d-8) t1=1.d0/t1
df(:,jc)=t1*df(:,jc)
u(:,jc)=alpha*df(:,jc)+t1*(mu(:,kp)-mu(:,kc))
call holdthd(m,nthd)
!$OMP PARALLEL DEFAULT(SHARED) &
!$OMP NUM_THREADS(nthd)
!$OMP DO
do k=1,m
  c(k)=dot_product(df(:,k),f(:,kc))
end do
!$OMP END DO NOWAIT
!$OMP DO
do k=1,m
  a(k,jc)=dot_product(df(:,jc),df(:,k))
  a(jc,k)=a(k,jc)
end do
!$OMP END DO
!$OMP END PARALLEL
call freethd(nthd)
beta(:,:)=a(:,:)
do k=1,m
  beta(k,k)=beta(k,k)+w0**2
end do
! invert beta
call dgetrf(m,m,beta,msd,ipiv,info)
if (info.eq.0) call dgetri(m,beta,msd,ipiv,work,m,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(mixbroyden): could not invert matrix")')
  write(*,*)
  stop
end if
do l=1,m
  gamma(l)=0.d0
  do k=1,m
    gamma(l)=gamma(l)+c(k)*beta(k,l)
  end do
end do
nu(:)=mu(:,kp)+alpha*f(:,kc)
do l=1,m
  t1=-gamma(l)
  nu(:)=nu(:)+t1*u(:,l)
end do
mu(:,kc)=nu(:)
end subroutine

