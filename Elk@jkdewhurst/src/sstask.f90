
! Copyright (C) 2012 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine sstask(fnum,fext)
use modmain
use modmpi
implicit none
! arguments
integer, intent(in) :: fnum
character(*), intent(out) :: fext
! local variables
logical exist
! only master process should search for file
if (.not.mp_mpi) goto 10
do iqss=1,nqpt
! construct the spin-spiral file extension
  call ssfext(iqss,fext)
! determine if the SS file exists
  inquire(file='SS'//trim(fext),exist=exist)
  if (.not.exist) then
    open(fnum,file='SS'//trim(fext),form='FORMATTED')
    goto 10
  end if
end do
iqss=0
write(*,*)
write(*,'("Info(sstask): nothing more to do")')
10 continue
! broadcast to all other processes
call mpi_bcast(iqss,1,mpi_integer,0,mpicom,ierror)
if (iqss.eq.0) then
  fext='.OUT'
else
  call ssfext(iqss,fext)
end if
end subroutine

