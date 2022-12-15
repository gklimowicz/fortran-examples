
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine oblas_init
implicit none
! set the initial number of OpenBLAS threads equal to one
call openblas_set_num_threads(1)
end subroutine

