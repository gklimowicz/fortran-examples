
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

! Stub routines for OpenMP.

integer function omp_get_num_procs()
implicit none
omp_get_num_procs=1
end function

integer function omp_get_max_threads()
implicit none
omp_get_max_threads=1
end function

integer function omp_get_level()
implicit none
omp_get_level=0
end function

subroutine omp_set_num_threads(num_threads)
implicit none
integer, intent(in) :: num_threads
end subroutine

subroutine omp_set_nested(nested)
implicit none
logical, intent(in) :: nested
end subroutine

subroutine omp_set_max_active_levels(max_levels)
implicit none
integer, intent(in) :: max_levels
end subroutine

subroutine omp_set_dynamic(dynamic)
implicit none
logical, intent(in) :: dynamic
end subroutine

subroutine omp_init_lock(lock)
implicit none
integer(8), intent(in) :: lock
end subroutine

subroutine omp_destroy_lock(lock)
implicit none
integer(8), intent(in) :: lock
end subroutine

subroutine omp_set_lock(lock)
implicit none
integer(8), intent(in) :: lock
end subroutine

subroutine omp_unset_lock(lock)
implicit none
integer(8), intent(in) :: lock
end subroutine

