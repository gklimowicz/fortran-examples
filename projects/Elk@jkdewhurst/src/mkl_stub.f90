
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! Stub routines for Intel MKL

subroutine mkl_set_num_threads(num_threads)
implicit none
integer, intent(in) :: num_threads
end subroutine

integer function mkl_set_num_threads_local(num_threads)
implicit none
integer, intent(in) :: num_threads
mkl_set_num_threads_local=0
end function

subroutine mkl_set_dynamic(dynamic)
implicit none
logical, intent(in) :: dynamic
end subroutine

