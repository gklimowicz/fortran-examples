
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modvars

use modmain
use moddftu

! if wrtvars is .true. then variables are written to VARIABLES.OUT
logical wrtvars
! if batch is .true. then Elk will run in batch mode
logical batch

contains

subroutine delvars
implicit none
! delete existing variables file
open(95,file='VARIABLES.OUT')
close(95,status='DELETE')
end subroutine

subroutine writevars(vname,n1,n2,n3,n4,n5,n6,nv,iv,iva,rv,rva,zv,zva,sv,sva)
implicit none
! arguments
character(*), intent(in) :: vname
integer, optional, intent(in) :: n1,n2,n3,n4,n5,n6
integer, optional, intent(in) :: nv,iv,iva(*)
real(8), optional, intent(in) :: rv,rva(*)
complex(8), optional, intent(in) :: zv,zva(*)
character(*), optional, intent(in) :: sv,sva(*)
! local variables
integer i
if ((present(iva)).or.(present(rva)).or.(present(zva)).or.(present(sva))) then
  if (.not.present(nv)) then
    write(*,*)
    write(*,'("Error(writevars): missing argument nv")')
    write(*,*)
    stop
  else
    if (nv.lt.0) then
      write(*,*)
      write(*,'("Error(writevars): nv < 0 : ",I8)') nv
      write(*,*)
      stop
    end if
  end if
end if
open(95,file='VARIABLES.OUT',position='APPEND',form='FORMATTED',action='WRITE')
write(95,*)
write(95,'(A," ")',advance='NO') trim(vname)
if (present(n1)) write(95,'(I8)',advance='NO') n1
if (present(n2)) write(95,'(I8)',advance='NO') n2
if (present(n3)) write(95,'(I8)',advance='NO') n3
if (present(n4)) write(95,'(I8)',advance='NO') n4
if (present(n5)) write(95,'(I8)',advance='NO') n5
if (present(n6)) write(95,'(I8)',advance='NO') n6
write(95,*)
if (present(iv)) then
  write(95,'(2I8)') 1,1
  write(95,'(I8)') iv
else if (present(rv)) then
  write(95,'(2I8)') 2,1
  write(95,'(G22.12)') rv
else if (present(zv)) then
  write(95,'(2I8)') 3,1
  write(95,'(2G22.12)') dble(zv),aimag(zv)
else if (present(sv)) then
  write(95,'(2I8)') 4,1
  write(95,'(A)') trim(sv)
else if (present(iva)) then
  write(95,'(2I8)') 1,nv
  do i=1,nv
    write(95,'(I8)') iva(i)
  end do
else if (present(rva)) then
  write(95,'(2I8)') 2,nv
  do i=1,nv
    write(95,'(G22.12)') rva(i)
  end do
else if (present(zva)) then
  write(95,'(2I8)') 3,nv
  do i=1,nv
    write(95,'(2G22.12)') dble(zva(i)),aimag(zva(i))
  end do
else if (present(sva)) then
  write(95,'(2I8)') 4,nv
  do i=1,nv
    write(95,'(A)') trim(sva(i))
  end do
end if
close(95)
end subroutine

end module

