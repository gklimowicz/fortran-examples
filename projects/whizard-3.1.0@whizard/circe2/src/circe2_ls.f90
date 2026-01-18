! circe2_ls.f90 -- beam spectra for linear colliders and photon colliders
! Copyright (C) 2001-2022 by Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!
! Circe2 is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! Circe2 is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!-----------------------------------------------------------------------
program circe2_ls
   use circe2
   use kinds
   implicit none
   integer :: i, lun
   character(len=132) :: buffer
   character(len=60) :: design, polspt
   integer :: pid1, hel1, pid2, hel2, nc
   real(kind=default) :: roots, lumi
   integer :: status
   logical :: exists, isopen
   character(len=1024) :: filename
   scan: do lun = 10, 99
      inquire (unit = lun, exist = exists, opened = isopen, iostat = status)
      if (status == 0 .and. exists .and. .not.isopen) exit scan
   end do scan
   if (lun > 99) lun = -1
   if (lun < 0) then
      write (*, '(A)') 'circe2_ls: no free unit'
      stop
   end if
   files: do i = 1, command_argument_count ()
      call get_command_argument (i, value = filename, status = status)
      if (status /= 0) then
         exit files
      else
         open (unit = lun, file = filename, status = 'old', iostat = status)
         if (status /= 0) then
            write (*, "(A,1X,A)") "circe2: can't open", trim(filename)
         else
            write (*, "(A,1X,A)") "file:", trim(filename)
            lines: do
               read (lun, '(A)', iostat = status) buffer
               if (status /= 0) exit lines
               if (buffer(1:7) == 'design,') then
                  read (lun, *) design, roots
                  read (lun, *)
                  read (lun, *) nc, polspt
                  write (*, '(A,1X,A)')    '        design:', trim(design)
                  write (*, '(A,1X,F7.1)') '       sqrt(s):', roots
                  write (*, '(A,1X,I3)')   '     #channels:', nc
                  write (*, '(A,1X,A)')      '  polarization:', trim(polspt)
                  write (*, '(4X,4(A5,2X),A)') &
                       'pid#1', 'hel#1', 'pid#2', 'hel#2', 'luminosity / (10^32cm^-2sec^-1)'
               else if (buffer(1:5) == 'pid1,') then
                  read (lun, *) pid1, hel1, pid2, hel2, lumi
                  write (*, '(4X,4(I5,2X),F10.2)') pid1, hel1, pid2, hel2, lumi
               end if
            end do lines
         end if
         close (unit = lun)
      end if
   end do files
end program circe2_ls
!-----------------------------------------------------------------------
