
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine mkl_init
use modomp
implicit none
! set the initial global number of MKL threads equal to one
call mkl_set_num_threads(1)
! set the maximum number of threads available to MKL
select case(maxthdmkl)
case(:-1)
  maxthdmkl=maxthd/abs(maxthdmkl)
  maxthdmkl=max(maxthdmkl,1)
case(0)
  maxthdmkl=maxthd
case default
  maxthdmkl=min(maxthdmkl,maxthd)
end select
! enable dynamic thread allocation
call mkl_set_dynamic(.true.)
end subroutine

