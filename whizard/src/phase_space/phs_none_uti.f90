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

module phs_none_uti

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use flavors
  use lorentz
  use model_data
  use process_constants
  use phs_base

  use phs_none

  use phs_base_ut, only: init_test_process_data, init_test_decay_data

  implicit none
  private

  public :: phs_none_1

contains

  subroutine phs_none_1 (u)
    integer, intent(in) :: u
    type(model_data_t), target :: model
    type(process_constants_t) :: process_data
    class(phs_config_t), allocatable :: phs_data
    real(default) :: sqrts

    write (u, "(A)")  "* Test output: phs_none_1"
    write (u, "(A)")  "*   Purpose: display &
         &phase-space configuration data"
    write (u, "(A)")

    allocate (phs_none_config_t :: phs_data)
    call phs_data%init (process_data, model)

    sqrts = 1000._default
    call phs_data%configure (sqrts, azimuthal_dependence=.false.)

    call phs_data%write (u)

    call phs_data%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_none_1"

  end subroutine phs_none_1


end module phs_none_uti
