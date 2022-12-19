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

module prc_test_core

  use kinds, only: default
  use lorentz
  use interactions
  use prc_test
  use prc_core

  implicit none
  private

  public :: test_t

  type, extends (prc_core_t) :: test_t
   contains
     procedure :: write => test_write
     procedure :: write_name => test_write_name
     procedure :: needs_mcset => test_needs_mcset
     procedure :: get_n_terms => test_get_n_terms
     procedure :: is_allowed => test_is_allowed
     procedure :: compute_hard_kinematics => test_compute_hard_kinematics
     procedure :: compute_eff_kinematics => test_compute_eff_kinematics
     procedure :: recover_kinematics => test_recover_kinematics
     procedure :: compute_amplitude => test_compute_amplitude
  end type test_t


  interface
    module subroutine test_write (object, unit)
      class(test_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine test_write
    module subroutine test_write_name (object, unit)
      class(test_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine test_write_name
    module function test_needs_mcset (object) result (flag)
      class(test_t), intent(in) :: object
      logical :: flag
    end function test_needs_mcset
    module function test_get_n_terms (object) result (n)
      class(test_t), intent(in) :: object
      integer :: n
    end function test_get_n_terms
    module function test_is_allowed (object, i_term, f, h, c) result (flag)
      class(test_t), intent(in) :: object
      integer, intent(in) :: i_term, f, h, c
      logical :: flag
    end function test_is_allowed
    module subroutine test_compute_hard_kinematics &
         (object, p_seed, i_term, int_hard, core_state)
      class(test_t), intent(in) :: object
      type(vector4_t), dimension(:), intent(in) :: p_seed
      integer, intent(in) :: i_term
      type(interaction_t), intent(inout) :: int_hard
      class(prc_core_state_t), intent(inout), allocatable :: core_state
    end subroutine test_compute_hard_kinematics
    module subroutine test_compute_eff_kinematics &
         (object, i_term, int_hard, int_eff, core_state)
      class(test_t), intent(in) :: object
      integer, intent(in) :: i_term
      type(interaction_t), intent(in) :: int_hard
      type(interaction_t), intent(inout) :: int_eff
      class(prc_core_state_t), intent(inout), allocatable :: core_state
    end subroutine test_compute_eff_kinematics
    module subroutine test_recover_kinematics &
         (object, p_seed, int_hard, int_eff, core_state)
      class(test_t), intent(in) :: object
      type(vector4_t), dimension(:), intent(inout) :: p_seed
      type(interaction_t), intent(inout) :: int_hard
      type(interaction_t), intent(inout) :: int_eff
      class(prc_core_state_t), intent(inout), allocatable :: core_state
    end subroutine test_recover_kinematics
    module function test_compute_amplitude &
         (object, j, p, f, h, c, fac_scale, ren_scale, alpha_qcd_forced, &
         core_state) result (amp)
      class(test_t), intent(in) :: object
      integer, intent(in) :: j
      type(vector4_t), dimension(:), intent(in) :: p
      integer, intent(in) :: f, h, c
      real(default), intent(in) :: fac_scale, ren_scale
      real(default), intent(in), allocatable :: alpha_qcd_forced
      class(prc_core_state_t), intent(inout), allocatable, optional :: &
           core_state
      complex(default) :: amp
    end function test_compute_amplitude
  end interface

end module prc_test_core
