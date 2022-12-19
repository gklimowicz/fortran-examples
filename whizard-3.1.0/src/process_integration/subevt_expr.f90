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
module subevt_expr

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use lorentz
  use subevents
  use variables
  use flavors
  use quantum_numbers
  use interactions
  use particles
  use expr_base

  implicit none
  private

  public :: parton_expr_t
  public :: event_expr_t

  type, extends (subevt_t), abstract :: subevt_expr_t
     logical :: subevt_filled = .false.
     type(var_list_t) :: var_list
     real(default) :: sqrts_hat = 0
     integer :: n_in = 0
     integer :: n_out = 0
     integer :: n_tot = 0
     logical :: has_selection = .false.
     class(expr_t), allocatable :: selection
     logical :: colorize_subevt = .false.
   contains
     procedure :: base_write => subevt_expr_write
     procedure (subevt_expr_final), deferred :: final
     procedure :: base_final => subevt_expr_final
     procedure (subevt_expr_setup_vars), deferred :: setup_vars
     procedure :: base_setup_vars => subevt_expr_setup_vars
     procedure :: setup_var_self => subevt_expr_setup_var_self
     procedure :: link_var_list => subevt_expr_link_var_list
     procedure :: setup_selection => subevt_expr_setup_selection
     procedure :: colorize => subevt_expr_colorize
     procedure :: reset_contents => subevt_expr_reset_contents
     procedure :: base_reset_contents => subevt_expr_reset_contents
     procedure :: base_evaluate => subevt_expr_evaluate
  end type subevt_expr_t

  type, extends (subevt_expr_t) :: parton_expr_t
     integer, dimension(:), allocatable :: i_beam
     integer, dimension(:), allocatable :: i_in
     integer, dimension(:), allocatable :: i_out
     logical :: has_scale = .false.
     logical :: has_fac_scale = .false.
     logical :: has_ren_scale = .false.
     logical :: has_weight = .false.
     class(expr_t), allocatable :: scale
     class(expr_t), allocatable :: fac_scale
     class(expr_t), allocatable :: ren_scale
     class(expr_t), allocatable :: weight
   contains
     procedure :: final => parton_expr_final
     procedure :: write => parton_expr_write
     procedure :: setup_vars => parton_expr_setup_vars
     procedure :: setup_scale => parton_expr_setup_scale
     procedure :: setup_fac_scale => parton_expr_setup_fac_scale
     procedure :: setup_ren_scale => parton_expr_setup_ren_scale
     procedure :: setup_weight => parton_expr_setup_weight
     procedure :: setup_subevt => parton_expr_setup_subevt
     procedure :: renew_flv_content_subevt => parton_expr_renew_flv_content_subevt
     procedure :: fill_subevt => parton_expr_fill_subevt
     procedure :: evaluate => parton_expr_evaluate
     procedure :: get_beam_index => parton_expr_get_beam_index
     procedure :: get_in_index => parton_expr_get_in_index
  end type parton_expr_t

  type, extends (subevt_expr_t) :: event_expr_t
     logical :: has_reweight = .false.
     logical :: has_analysis = .false.
     class(expr_t), allocatable :: reweight
     class(expr_t), allocatable :: analysis
     logical :: has_id = .false.
     type(string_t) :: id
     logical :: has_num_id = .false.
     integer :: num_id = 0
     logical :: has_index = .false.
     integer :: index = 0
     logical :: has_sqme_ref = .false.
     real(default) :: sqme_ref = 0
     logical :: has_sqme_prc = .false.
     real(default) :: sqme_prc = 0
     logical :: has_weight_ref = .false.
     real(default) :: weight_ref = 0
     logical :: has_weight_prc = .false.
     real(default) :: weight_prc = 0
     logical :: has_excess_prc = .false.
     real(default) :: excess_prc = 0
     integer :: n_alt = 0
     logical :: has_sqme_alt = .false.
     real(default), dimension(:), allocatable :: sqme_alt
     logical :: has_weight_alt = .false.
     real(default), dimension(:), allocatable :: weight_alt
   contains
     procedure :: final => event_expr_final
     procedure :: write => event_expr_write
     procedure :: init => event_expr_init
     procedure :: setup_vars => event_expr_setup_vars
     procedure :: setup_analysis => event_expr_setup_analysis
     procedure :: setup_reweight => event_expr_setup_reweight
     procedure :: set_process_id => event_expr_set_process_id
     procedure :: set_process_num_id => event_expr_set_process_num_id
     procedure :: reset_contents => event_expr_reset_contents
     procedure :: set => event_expr_set
     procedure :: has_event_index => event_expr_has_event_index
     procedure :: get_event_index => event_expr_get_event_index
     procedure :: set_event_index => event_expr_set_event_index
     procedure :: reset_event_index => event_expr_reset_event_index
     procedure :: increment_event_index => event_expr_increment_event_index
     procedure :: fill_subevt => event_expr_fill_subevt
     procedure :: evaluate => event_expr_evaluate
  end type event_expr_t


  interface interaction_momenta_to_subevt
     module procedure interaction_momenta_to_subevt_id
     module procedure interaction_momenta_to_subevt_tr
  end interface


  interface
    module subroutine subevt_expr_write (object, unit, pacified)
      class(subevt_expr_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: pacified
    end subroutine subevt_expr_write
    module subroutine subevt_expr_final (object)
      class(subevt_expr_t), intent(inout) :: object
    end subroutine subevt_expr_final
    module subroutine subevt_expr_setup_vars (expr, sqrts)
      class(subevt_expr_t), intent(inout), target :: expr
      real(default), intent(in) :: sqrts
    end subroutine subevt_expr_setup_vars
    module subroutine subevt_expr_setup_var_self (expr)
      class(subevt_expr_t), intent(inout), target :: expr
    end subroutine subevt_expr_setup_var_self
    module subroutine subevt_expr_link_var_list (expr, var_list)
      class(subevt_expr_t), intent(inout) :: expr
      type(var_list_t), intent(in), target :: var_list
    end subroutine subevt_expr_link_var_list
    module subroutine subevt_expr_setup_selection (expr, ef_cuts)
      class(subevt_expr_t), intent(inout), target :: expr
      class(expr_factory_t), intent(in) :: ef_cuts
    end subroutine subevt_expr_setup_selection
    module subroutine subevt_expr_colorize (expr, colorize_subevt)
      class(subevt_expr_t), intent(inout), target :: expr
      logical, intent(in) :: colorize_subevt
    end subroutine subevt_expr_colorize
    module subroutine subevt_expr_reset_contents (expr)
      class(subevt_expr_t), intent(inout) :: expr
    end subroutine subevt_expr_reset_contents
    module subroutine subevt_expr_evaluate (expr, passed)
      class(subevt_expr_t), intent(inout) :: expr
      logical, intent(out) :: passed
    end subroutine subevt_expr_evaluate
    module subroutine parton_expr_final (object)
      class(parton_expr_t), intent(inout) :: object
    end subroutine parton_expr_final
    module subroutine parton_expr_write (object, unit, prefix, pacified)
      class(parton_expr_t), intent(in) :: object
      integer, intent(in), optional :: unit
      character(*), intent(in), optional :: prefix
      logical, intent(in), optional :: pacified
    end subroutine parton_expr_write
    module subroutine parton_expr_setup_vars (expr, sqrts)
      class(parton_expr_t), intent(inout), target :: expr
      real(default), intent(in) :: sqrts
    end subroutine parton_expr_setup_vars
    module subroutine parton_expr_setup_scale (expr, ef_scale)
      class(parton_expr_t), intent(inout), target :: expr
      class(expr_factory_t), intent(in) :: ef_scale
    end subroutine parton_expr_setup_scale
    module subroutine parton_expr_setup_fac_scale (expr, ef_fac_scale)
      class(parton_expr_t), intent(inout), target :: expr
      class(expr_factory_t), intent(in) :: ef_fac_scale
    end subroutine parton_expr_setup_fac_scale
    module subroutine parton_expr_setup_ren_scale (expr, ef_ren_scale)
      class(parton_expr_t), intent(inout), target :: expr
      class(expr_factory_t), intent(in) :: ef_ren_scale
    end subroutine parton_expr_setup_ren_scale
    module subroutine parton_expr_setup_weight (expr, ef_weight)
      class(parton_expr_t), intent(inout), target :: expr
      class(expr_factory_t), intent(in) :: ef_weight
    end subroutine parton_expr_setup_weight
    module subroutine parton_expr_setup_subevt (expr, int, &
         i_beam, i_in, i_out, f_beam, f_in, f_out)
      class(parton_expr_t), intent(inout) :: expr
      type(interaction_t), intent(in), target :: int
      integer, dimension(:), intent(in) :: i_beam, i_in, i_out
      type(flavor_t), dimension(:), intent(in) :: f_beam, f_in, f_out
    end subroutine parton_expr_setup_subevt
    module subroutine parton_expr_renew_flv_content_subevt (expr, int, &
         i_beam, i_in, i_out, f_beam, f_in, f_out)
      class(parton_expr_t), intent(inout) :: expr
      type(interaction_t), intent(in), target :: int
      integer, dimension(:), intent(in) :: i_beam, i_in, i_out
      type(flavor_t), dimension(:), intent(in) :: f_beam, f_in, f_out
    end subroutine parton_expr_renew_flv_content_subevt
    module subroutine interaction_momenta_to_subevt_id &
         (int, j_beam, j_in, j_out, subevt)
      type(interaction_t), intent(in) :: int
      integer, dimension(:), intent(in) :: j_beam, j_in, j_out
      type(subevt_t), intent(inout) :: subevt
    end subroutine interaction_momenta_to_subevt_id
    module subroutine interaction_momenta_to_subevt_tr &
         (int, j_beam, j_in, j_out, lt, subevt)
      type(interaction_t), intent(in) :: int
      integer, dimension(:), intent(in) :: j_beam, j_in, j_out
      type(subevt_t), intent(inout) :: subevt
      type(lorentz_transformation_t), intent(in) :: lt
    end subroutine interaction_momenta_to_subevt_tr
    module subroutine parton_expr_fill_subevt (expr, int)
      class(parton_expr_t), intent(inout) :: expr
      type(interaction_t), intent(in), target :: int
    end subroutine parton_expr_fill_subevt
    module subroutine parton_expr_evaluate (expr, passed, scale, fac_scale, &
         ren_scale, weight, scale_forced, force_evaluation)
      class(parton_expr_t), intent(inout) :: expr
      logical, intent(out) :: passed
      real(default), intent(out) :: scale
      real(default), allocatable, intent(out) :: fac_scale
      real(default), allocatable, intent(out) :: ren_scale
      real(default), intent(out) :: weight
      real(default), intent(in), allocatable, optional :: scale_forced
      logical, intent(in), optional :: force_evaluation
    end subroutine parton_expr_evaluate
    module subroutine parton_expr_get_beam_index (expr, i_beam)
      class(parton_expr_t), intent(in) :: expr
      integer, dimension(:), intent(out) :: i_beam
    end subroutine parton_expr_get_beam_index
    module subroutine parton_expr_get_in_index (expr, i_in)
      class(parton_expr_t), intent(in) :: expr
      integer, dimension(:), intent(out) :: i_in
    end subroutine parton_expr_get_in_index
    module subroutine event_expr_final (object)
      class(event_expr_t), intent(inout) :: object
    end subroutine event_expr_final
    module subroutine event_expr_write (object, unit, prefix, pacified)
      class(event_expr_t), intent(in) :: object
      integer, intent(in), optional :: unit
      character(*), intent(in), optional :: prefix
      logical, intent(in), optional :: pacified
    end subroutine event_expr_write
    module subroutine event_expr_init (expr, n_alt)
      class(event_expr_t), intent(out) :: expr
      integer, intent(in), optional :: n_alt
    end subroutine event_expr_init
    module subroutine event_expr_setup_vars (expr, sqrts)
      class(event_expr_t), intent(inout), target :: expr
      real(default), intent(in) :: sqrts
    end subroutine event_expr_setup_vars
    module subroutine event_expr_setup_analysis (expr, ef_analysis)
      class(event_expr_t), intent(inout), target :: expr
      class(expr_factory_t), intent(in) :: ef_analysis
    end subroutine event_expr_setup_analysis
    module subroutine event_expr_setup_reweight (expr, ef_reweight)
      class(event_expr_t), intent(inout), target :: expr
      class(expr_factory_t), intent(in) :: ef_reweight
    end subroutine event_expr_setup_reweight
    module subroutine event_expr_set_process_id (expr, id)
      class(event_expr_t), intent(inout) :: expr
      type(string_t), intent(in) :: id
    end subroutine event_expr_set_process_id
    module subroutine event_expr_set_process_num_id (expr, num_id)
      class(event_expr_t), intent(inout) :: expr
      integer, intent(in) :: num_id
    end subroutine event_expr_set_process_num_id
    module subroutine event_expr_reset_contents (expr)
      class(event_expr_t), intent(inout) :: expr
    end subroutine event_expr_reset_contents
    module subroutine event_expr_set (expr, &
         weight_ref, weight_prc, weight_alt, &
         excess_prc, &
         sqme_ref, sqme_prc, sqme_alt)
      class(event_expr_t), intent(inout) :: expr
      real(default), intent(in), optional :: weight_ref, weight_prc
      real(default), intent(in), optional :: excess_prc
      real(default), intent(in), optional :: sqme_ref, sqme_prc
      real(default), dimension(:), intent(in), optional :: sqme_alt, weight_alt
    end subroutine event_expr_set
    module function event_expr_has_event_index (expr) result (flag)
      class(event_expr_t), intent(in) :: expr
      logical :: flag
    end function event_expr_has_event_index
    module function event_expr_get_event_index (expr) result (index)
      class(event_expr_t), intent(in) :: expr
      integer :: index
    end function event_expr_get_event_index
    module subroutine event_expr_set_event_index (expr, index)
      class(event_expr_t), intent(inout) :: expr
      integer, intent(in) :: index
    end subroutine event_expr_set_event_index
    module subroutine event_expr_reset_event_index (expr)
      class(event_expr_t), intent(inout) :: expr
    end subroutine event_expr_reset_event_index
    module subroutine event_expr_increment_event_index (expr, offset)
      class(event_expr_t), intent(inout) :: expr
      integer, intent(in), optional :: offset
    end subroutine event_expr_increment_event_index
    module subroutine event_expr_fill_subevt (expr, particle_set)
      class(event_expr_t), intent(inout) :: expr
      type(particle_set_t), intent(in) :: particle_set
    end subroutine event_expr_fill_subevt
    module subroutine event_expr_evaluate &
         (expr, passed, reweight, analysis_flag)
      class(event_expr_t), intent(inout) :: expr
      logical, intent(out) :: passed
      real(default), intent(out) :: reweight
      logical, intent(out) :: analysis_flag
    end subroutine event_expr_evaluate
  end interface

end module subevt_expr
