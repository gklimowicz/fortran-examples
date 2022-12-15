
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine blis_init
implicit none
! set the initial number of BLIS threads equal to one
call bli_thread_set_num_threads(1)
end subroutine

