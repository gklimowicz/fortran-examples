! WHIZARD 3.1.0 Dec 14 2022
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!
!     with contributions from
!     cf. main AUTHORS file
!
! WHIZARD is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! WHIZARD is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file has been stripped of most comments.  For documentation, refer
! to the source 'whizard.nw'

submodule (prc_core_def) prc_core_def_s

  use io_units
  use diagnostics

  implicit none

contains

  module subroutine prc_core_def_set_md5sum (core_def, md5sum)
    class(prc_core_def_t), intent(inout) :: core_def
    character(32) :: md5sum
    if (allocated (core_def%writer))  core_def%writer%md5sum = md5sum
  end subroutine prc_core_def_set_md5sum

  module function prc_core_def_needs_code () result (flag)
    logical :: flag
    flag = .false.
  end function prc_core_def_needs_code

  module subroutine allocate_core_def (template, name, core_def)
    type(prc_template_t), dimension(:), intent(in) :: template
    type(string_t), intent(in) :: name
    class(prc_core_def_t), allocatable :: core_def
    integer :: i
    do i = 1, size (template)
       if (template(i)%core_def%type_string () == name) then
          allocate (core_def, source = template(i)%core_def)
          return
       end if
    end do
  end subroutine allocate_core_def


end submodule prc_core_def_s

