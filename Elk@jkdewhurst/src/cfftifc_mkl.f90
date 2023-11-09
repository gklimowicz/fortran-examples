
! Copyright (C) 2022 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine cfftifc(nd,n,sgn,c)
use mkl_dfti
implicit none
! arguments
integer, intent(in) :: nd,n(nd),sgn
complex(4), intent(inout) :: c(*)
! local variables
! local variables
integer status,p
real(4) t1
type(DFTI_DESCRIPTOR), pointer :: handle
! interface to the Intel MKL advanced Discreet Fourier Transform (DFT) routines
! (with thanks to Torbjorn Bjorkman)
p=product(n(:))
t1=1.e0/real(p)
status=DftiCreateDescriptor(handle,DFTI_SINGLE,DFTI_COMPLEX,nd,n)
status=DftiSetValue(handle,DFTI_FORWARD_SCALE,t1)
status=DftiCommitDescriptor(handle)
if (sgn.eq.-1) then
  status=DftiComputeForward(handle,c)
else
  status=DftiComputeBackward(handle,c)
end if
status=DftiFreeDescriptor(handle)
end subroutine

subroutine rcfftifc(nd,n,sgn,r,c)
use mkl_dfti
implicit none
! arguments
integer, intent(in) :: nd,n(nd),sgn
real(4), intent(inout) :: r(*)
complex(4), intent(inout) :: c(*)
! local variables
integer status,p,i
real(4) t1
type(DFTI_DESCRIPTOR), pointer :: handle
! automatic arrays
integer strides(0:nd)
p=product(n(:))
t1=1.e0/real(p)
status=DftiCreateDescriptor(handle,DFTI_SINGLE,DFTI_REAL,nd,n)
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
  status=DftiComputeForward(handle,r,c)
else
  status=DftiSetValue(handle,DFTI_INPUT_STRIDES,strides)
  status=DftiCommitDescriptor(handle)
  status=DftiComputeBackward(handle,c,r)
end if
status=DftiFreeDescriptor(handle)
end subroutine
