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

submodule (dispatch_me_methods) dispatch_me_methods_s

  implicit none

contains

  module subroutine dispatch_core_update &
       (core, model, helicity_selection, qcd, saved_core)

    class(prc_core_t), allocatable, intent(inout) :: core
    class(model_data_t), intent(in), optional, target :: model
    type(helicity_selection_t), intent(in), optional :: helicity_selection
    type(qcd_t), intent(in), optional :: qcd
    class(prc_core_t), allocatable, intent(inout), optional :: saved_core

    if (present (saved_core)) then
       allocate (saved_core, source = core)
    end if
    select type (core)
    type is (test_t)
    type is (prc_omega_t)
       call core%set_parameters (model, helicity_selection, qcd)
       call core%activate_parameters ()
    class is (prc_external_t)
      call msg_message ("Updating user defined cores is not implemented yet.")
    class default
       call msg_bug ("Process core update: unexpected process definition type")
    end select
  end subroutine dispatch_core_update

  module subroutine dispatch_core_restore (core, saved_core)

    class(prc_core_t), allocatable, intent(inout) :: core
    class(prc_core_t), allocatable, intent(inout) :: saved_core

    call move_alloc (from = saved_core, to = core)
    select type (core)
    type is (test_t)
    type is (prc_omega_t)
       call core%activate_parameters ()
    class default
       call msg_bug ("Process core restore: unexpected process definition type")
    end select
  end subroutine dispatch_core_restore


end submodule dispatch_me_methods_s

