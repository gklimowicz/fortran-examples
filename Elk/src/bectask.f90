
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine bectask(fnum,fext)
use modmain
use modphonon
use modmpi
implicit none
! arguments
integer, intent(in) :: fnum
character(*), intent(out) :: fext
! local variables
logical exist
! only master process should search for file
if (.not.mp_mpi) goto 10
do ipph=1,3
  do isph=1,nspecies
    do iaph=1,natoms(isph)
! Born effective charge file extension
      call becfext(isph,iaph,ipph,fext)
! determine if the BEC file with this extension exists
      inquire(file='BEC'//trim(fext),exist=exist)
      if (.not.exist) then
        open(fnum,file='BEC'//trim(fext),form='FORMATTED')
        iasph=idxas(iaph,isph)
        goto 10
      end if
    end do
  end do
end do
isph=0; iaph=0; iasph=0; ipph=0
write(*,*)
write(*,'("Info(bectask): nothing more to do")')
10 continue
! broadcast to all other MPI processes
call mpi_bcast(isph,1,mpi_integer,0,mpicom,ierror)
call mpi_bcast(iaph,1,mpi_integer,0,mpicom,ierror)
call mpi_bcast(iasph,1,mpi_integer,0,mpicom,ierror)
call mpi_bcast(ipph,1,mpi_integer,0,mpicom,ierror)
if (isph.eq.0) then
  fext='.OUT'
else
  call becfext(isph,iaph,ipph,fext)
end if
end subroutine

