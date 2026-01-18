
! Copyright (C) 2002-2005 J. K. Dewhurst, S. Sharma and C. Ambrosch-Draxl.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine zfftifc(nd,n,sgn,z)
use mkl_dfti
implicit none
! arguments
integer, intent(in) :: nd,n(nd),sgn
complex(8), intent(inout) :: z(*)
! local variables
integer status,p
real(8) t1
type(DFTI_DESCRIPTOR), pointer :: handle
! interface to the Intel MKL advanced Discreet Fourier Transform (DFT) routines
! (with thanks to Torbjorn Bjorkman)
p=product(n(:))
t1=1.d0/dble(p)
status=DftiCreateDescriptor(handle,DFTI_DOUBLE,DFTI_COMPLEX,nd,n)
status=DftiSetValue(handle,DFTI_FORWARD_SCALE,t1)
status=DftiCommitDescriptor(handle)
if (sgn.eq.-1) then
  status=DftiComputeForward(handle,z)
else
  status=DftiComputeBackward(handle,z)
end if
status=DftiFreeDescriptor(handle)
end subroutine

subroutine rzfftifc(nd,n,sgn,r,z)
use mkl_dfti
implicit none
! arguments
integer, intent(in) :: nd,n(nd),sgn
real(8), intent(inout) :: r(*)
complex(8), intent(inout) :: z(*)
! local variables
integer status,p,i
real(8) t1
type(DFTI_DESCRIPTOR), pointer :: handle
! automatic arrays
integer strides(0:nd)
p=product(n(:))
t1=1.d0/dble(p)
status=DftiCreateDescriptor(handle,DFTI_DOUBLE,DFTI_REAL,nd,n)
status=DftiSetValue(handle,DFTI_CONJUGATE_EVEN_STORAGE,DFTI_COMPLEX_COMPLEX)
status=DftiSetValue(handle,DFTI_PLACEMENT,DFTI_NOT_INPLACE)
status=DftiSetValue(handle,DFTI_FORWARD_SCALE,t1)
strides(0)=0
strides(1)=1
if (nd.gt.1) then
  strides(2)=n(1)/2+1
  do i=2,nd-1
    strides(i+1)=strides(i)*n(i)
  end do
end if
if (sgn.eq.-1) then
  status=DftiSetValue(handle,DFTI_OUTPUT_STRIDES,strides)
  status=DftiCommitDescriptor(handle)
  status=DftiComputeForward(handle,r,z)
else
  status=DftiSetValue(handle,DFTI_INPUT_STRIDES,strides)
  status=DftiCommitDescriptor(handle)
  status=DftiComputeBackward(handle,z,r)
end if
status=DftiFreeDescriptor(handle)
end subroutine

