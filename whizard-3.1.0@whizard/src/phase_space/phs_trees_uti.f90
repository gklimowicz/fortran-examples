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

module phs_trees_uti

!!!  use kinds, only: default
  use kinds, only: TC
  use iso_varying_string, string_t => varying_string
  use flavors, only: flavor_t
  use model_data, only: model_data_t

  use resonances, only: resonance_history_t
  use mappings, only: mapping_defaults_t

  use phs_trees

  implicit none
  private

  public :: phs_tree_1
  public :: phs_tree_2

contains

  subroutine phs_tree_1 (u)
    integer, intent(in) :: u
    type(phs_tree_t) :: tree
    type(model_data_t), target :: model
    type(flavor_t), dimension(5) :: flv
    integer :: i

    write (u, "(A)")  "* Test output: phs_tree_1"
    write (u, "(A)")  "*   Purpose: test PHS tree routines"
    write (u, "(A)")

    write (u, "(A)")  "* Read model file"

    call model%init_sm_test ()

    write (u, "(A)")
    write (u, "(A)")  "* Set up flavors"
    write (u, "(A)")

    call flv%init ([1, -2, 24, 5, -5], model)
    do i = 1, 5
       write (u, "(1x)", advance="no")
       call flv(i)%write (u)
    end do
    write (u, *)

    write (u, "(A)")
    write (u, "(A)")  "* Create tree"
    write (u, "(A)")

    call tree%init (2, 3, 0, 0)
    call tree%from_array ([integer(TC) :: 1, 2, 3, 4, 7, 8, 16])
    call tree%set_mass_sum (flv)
    call tree%set_effective_masses ()

    call tree%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call tree%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_tree_1"

  end subroutine phs_tree_1

  subroutine phs_tree_2 (u)
    integer, intent(in) :: u
    type(phs_tree_t) :: tree
    type(model_data_t), target :: model
    type(mapping_defaults_t) :: mapping_defaults
    type(flavor_t), dimension(5) :: flv
    type(resonance_history_t) :: res_history
    integer :: i

    write (u, "(A)")  "* Test output: phs_tree_2"
    write (u, "(A)")  "*   Purpose: test PHS tree with resonances"
    write (u, "(A)")

    write (u, "(A)")  "* Read model file"

    call model%init_sm_test ()

    write (u, "(A)")
    write (u, "(A)")  "* Set up flavors"
    write (u, "(A)")

    call flv%init ([1, -2, 24, 5, -5], model)
    do i = 1, 5
       write (u, "(1x)", advance="no")
       call flv(i)%write (u)
    end do
    write (u, *)

    write (u, "(A)")
    write (u, "(A)")  "* Create tree with mappings"
    write (u, "(A)")

    call tree%init (2, 3, 0, 0)
    call tree%from_array ([integer(TC) :: 1, 2, 3, 4, 7, 8, 16])
    call tree%set_mass_sum (flv)

    call tree%init_mapping (3_TC, var_str ("s_channel"), -24, model)
    call tree%init_mapping (7_TC, var_str ("s_channel"), 23, model)

    call tree%set_mapping_parameters (mapping_defaults, variable_limits=.false.)
    call tree%set_effective_masses ()

    call tree%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Extract resonances from mappings"
    write (u, "(A)")

    call tree%extract_resonance_history (res_history)
    call res_history%write (u)

    write (u, "(A)")
    write (u, "(A)")  "* Cleanup"

    call tree%final ()
    call model%final ()

    write (u, "(A)")
    write (u, "(A)")  "* Test output end: phs_tree_2"

  end subroutine phs_tree_2


end module phs_trees_uti
