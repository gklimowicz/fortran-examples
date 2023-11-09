
! Copyright (C) 2019 Chung-Yu Wang, J. K. Dewhurst, S. Sharma and
! E. K. U. Gross. This file is distributed under the terms of the GNU General
! Public License. See the file COPYING for license details.

subroutine eveqnwxy(n,p,dw,ex,fy,w)
implicit none
! arguments
integer, intent(in) :: n,p
complex(8), intent(inout) :: dw(n,n),ex(n,n),fy(n)
real(8), intent(out) :: w(n)
! local variables
integer n2,i,j
real(8) t1,t2
! allocatable arrays
integer, allocatable :: idx(:)
real(8), allocatable :: r(:)
complex(8), allocatable :: w2(:),h(:,:)
complex(8), allocatable :: x(:),a(:,:)
! external functions
real(8), external :: dznrm2
n2=2*n
! setup the bosonic Bogoliubov Hamiltonian
allocate(w2(n2),h(n2,n2))
do j=1,n
  do i=1,n
    h(i,j)=dw(i,j)
    h(n+i,n+j)=-dw(i,j)
    h(i,n+j)=-ex(i,j)
    h(n+i,j)=ex(i,j)
  end do
end do
! find the eigenvalues and right eigenvectors
call eveqnzg(n2,n2,h,w2)
! select the eigenpairs corresponding to W†W - X†X = I
allocate(idx(n2),r(n2))
do j=1,n2
  t1=dznrm2(n,h(1,j),1)**2
  t2=dznrm2(n,h(n+1,j),1)**2
  r(j)=t1-t2
end do
call sortidx(n2,r,idx)
! pseudo-normalise the eigenvectors and store in output arrays
do i=1,n
  j=idx(n+i)
  t1=abs(r(j))+1.d-8
  t1=(1.d0-(1.d0-t1)**p)/sqrt(t1)
  w(i)=dble(w2(j))
  dw(1:n,i)=t1*h(1:n,j)
  ex(1:n,i)=t1*h(n+1:n2,j)
end do
deallocate(idx,r,w2,h)
! solve for the vector y
allocate(x(n),a(n,n))
a(:,:)=dw(:,:)-ex(:,:)
x(:)=fy(:)
call zgemv('T',n,n,(1.d0,0.d0),a,n,x,1,(0.d0,0.d0),fy,1)
do i=1,n
  if (w(i).gt.1.d-6) then
    fy(i)=fy(i)/w(i)
  else
    fy(i)=0.d0
  end if
end do
deallocate(x,a)
end subroutine

