
! Copyright (C) 2015 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modomp

! maximum number of OpenMP threads available
integer maxthd
! maximum number of OpenMP threads for the first nesting level
integer maxthd1
! maximum number of threads available to MKL
integer maxthdmkl
! maximum OpenMP nesting level
integer maxlvl
! number of active OpenMP threads for each nesting level
integer, allocatable, protected :: nathd(:)

interface

integer function omp_get_num_procs()
end function

integer function omp_get_max_threads()
end function

integer function omp_get_level()
end function

subroutine omp_set_num_threads(num_threads)
integer, intent(in) :: num_threads
end subroutine

integer function mkl_set_num_threads_local(num_threads)
integer, intent(in) :: num_threads
end function

subroutine omp_set_nested(nested)
logical, intent(in) :: nested
end subroutine

subroutine omp_set_max_active_levels(max_levels)
integer, intent(in) :: max_levels
end subroutine

subroutine omp_set_dynamic(dynamic_threads)
logical, intent(in) :: dynamic_threads
end subroutine

end interface

contains

subroutine omp_init
implicit none
! determine the maximum number of available threads
select case(maxthd)
case(:-1)
! set the number of threads equal to a fraction of the number of processors
  maxthd=omp_get_num_procs()/abs(maxthd)
  maxthd=max(maxthd,1)
  call omp_set_num_threads(maxthd)
case(0)
! use the system default number of threads
  maxthd=omp_get_max_threads()
case default
! use the number of threads specified in the input file
  call omp_set_num_threads(maxthd)
end select
! determine the maximum number of threads available at first nesting level
select case(maxthd1)
case(:-1)
  maxthd1=maxthd/abs(maxthd1)
  maxthd1=max(maxthd1,1)
case(0)
  maxthd1=maxthd
case default
  maxthd1=min(maxthd1,maxthd)
end select
! switch off dynamic allocation of threads
call omp_set_dynamic(.false.)
! allow nested parallelism (deprecated in OpenMP version 5)
call omp_set_nested(.true.)
! set the maximum nesting level
call omp_set_max_active_levels(maxlvl)
! allocate the number of active threads array
if (allocated(nathd)) deallocate(nathd)
allocate(nathd(0:maxlvl))
! initialise the number of active threads
call omp_reset
end subroutine

subroutine omp_reset
implicit none
! number of active threads at each nesting level
nathd(0)=1
nathd(1:)=0
end subroutine

subroutine holdthd(nloop,nthd)
implicit none
! arguments
integer, intent(in) :: nloop
integer, intent(out) :: nthd
! local variables
integer lvl,na,n
! current nesting level
lvl=omp_get_level()
if ((lvl.lt.0).or.(lvl.ge.maxlvl)) then
  nthd=1
  return
end if
! determine number of active threads at the current nesting level
na=nathd(lvl)
na=max(min(na,maxthd),1)
! number of threads allowed for this loop
nthd=maxthd/na
if (mod(maxthd,na).gt.0) nthd=nthd+1
if (lvl.eq.0) nthd=min(nthd,maxthd1)
nthd=max(min(nthd,maxthd,nloop),1)
! add to number of active threads in next nesting level
n=nathd(lvl+1)+nthd
n=max(min(n,maxthd),0)
!$OMP ATOMIC WRITE
nathd(lvl+1)=n
end subroutine

subroutine freethd(nthd)
implicit none
! arguments
integer, intent(in) :: nthd
! local variables
integer lvl,n
! current nesting level
lvl=omp_get_level()
if ((lvl.lt.0).or.(lvl.ge.maxlvl)) return
! subtract from the number of active threads in next nesting level
n=nathd(lvl+1)-nthd
n=max(min(n,maxthd),0)
!$OMP ATOMIC WRITE
nathd(lvl+1)=n
end subroutine

end module

