
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: unitary
! !INTERFACE:
subroutine unitary(n,a)
! !INPUT/OUTPUT PARAMETERS:
!   n : order of matrix (in,integer)
!   a : complex square matrix (inout,complex(n,n))
! !DESCRIPTION:
!   Finds the closest unitary matrix (in terms of the Frobenius norm) to a
!   given matrix $A$. Let $U\Sigma V^{\dag}$ be the singular value
!   decomposition of $A$. Then it can be shown that $UV^{\dag}$ is the closest
!   unitary matrix to $A$. The input matrix is overwritten by this matrix.
!
! !REVISION HISTORY:
!   Created January 2020 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(in) :: n
complex(8), intent(inout) :: a(n,n)
! local variables
integer lwork,info
complex(8), parameter :: zzero=(0.d0,0.d0),zone=(1.d0,0.d0)
! automatic arrays
real(8) s(n),rwork(5*n)
! allocatable arrays
complex(8), allocatable :: u(:,:),vt(:,:),work(:)
! perform singular value decomposition on matrix
lwork=3*n
allocate(u(n,n),vt(n,n),work(lwork))
call zgesvd('A','A',n,n,a,n,s,u,n,vt,n,work,lwork,rwork,info)
if (info.ne.0) then
  write(*,*)
  write(*,'("Error(unitary): singular value decomposition failed")')
  write(*,'(" ZGESVD returned INFO = ",I8)') info
  write(*,*)
  stop
end if
! multiply the two unitary matrices together and store in the input matrix
call zgemm('N','N',n,n,n,zone,u,n,vt,n,zzero,a,n)
deallocate(u,vt,work)
end subroutine
!EOC

