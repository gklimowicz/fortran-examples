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

module event_base

  use kinds, only: default
  use kinds, only: i64
  use iso_varying_string, string_t => varying_string
  use model_data
  use particles

  implicit none
  private

  public :: generic_event_t
  public :: event_normalization_mode
  public :: event_normalization_string
  public :: event_normalization_update
  public :: event_callback_t
  public :: event_callback_nop_t

  integer, parameter, public :: NORM_UNDEFINED = 0
  integer, parameter, public :: NORM_UNIT = 1
  integer, parameter, public :: NORM_N_EVT = 2
  integer, parameter, public :: NORM_SIGMA = 3
  integer, parameter, public :: NORM_S_N = 4


  type, abstract :: generic_event_t
     !private
     logical :: particle_set_is_valid = .false.
     type(particle_set_t), pointer :: particle_set => null ()
     logical :: sqme_ref_known = .false.
     real(default) :: sqme_ref = 0
     logical :: sqme_prc_known = .false.
     real(default) :: sqme_prc = 0
     logical :: weight_ref_known = .false.
     real(default) :: weight_ref = 0
     logical :: weight_prc_known = .false.
     real(default) :: weight_prc = 0
     logical :: excess_prc_known = .false.
     real(default) :: excess_prc = 0
     logical :: n_dropped_known = .false.
     integer :: n_dropped = 0
     integer :: n_alt = 0
     logical :: sqme_alt_known = .false.
     real(default), dimension(:), allocatable :: sqme_alt
     logical :: weight_alt_known = .false.
     real(default), dimension(:), allocatable :: weight_alt
   contains
     procedure :: base_init => generic_event_init
     procedure :: has_valid_particle_set => generic_event_has_valid_particle_set
     procedure :: accept_particle_set => generic_event_accept_particle_set
     procedure :: discard_particle_set => generic_event_discard_particle_set
     procedure :: get_particle_set_ptr => generic_event_get_particle_set_ptr
     procedure :: link_particle_set => generic_event_link_particle_set
     procedure :: sqme_prc_is_known => generic_event_sqme_prc_is_known
     procedure :: sqme_ref_is_known => generic_event_sqme_ref_is_known
     procedure :: sqme_alt_is_known => generic_event_sqme_alt_is_known
     procedure :: weight_prc_is_known => generic_event_weight_prc_is_known
     procedure :: weight_ref_is_known => generic_event_weight_ref_is_known
     procedure :: weight_alt_is_known => generic_event_weight_alt_is_known
     procedure :: excess_prc_is_known => generic_event_excess_prc_is_known
     procedure :: get_n_alt => generic_event_get_n_alt
     procedure :: get_sqme_prc => generic_event_get_sqme_prc
     procedure :: get_sqme_ref => generic_event_get_sqme_ref
     generic :: get_sqme_alt => &
          generic_event_get_sqme_alt_0, generic_event_get_sqme_alt_1
     procedure :: generic_event_get_sqme_alt_0
     procedure :: generic_event_get_sqme_alt_1
     procedure :: get_weight_prc => generic_event_get_weight_prc
     procedure :: get_weight_ref => generic_event_get_weight_ref
     generic :: get_weight_alt => &
          generic_event_get_weight_alt_0, generic_event_get_weight_alt_1
     procedure :: generic_event_get_weight_alt_0
     procedure :: generic_event_get_weight_alt_1
     procedure :: get_n_dropped => generic_event_get_n_dropped
     procedure :: get_excess_prc => generic_event_get_excess_prc
     procedure :: set_sqme_prc => generic_event_set_sqme_prc
     procedure :: set_sqme_ref => generic_event_set_sqme_ref
     procedure :: set_sqme_alt => generic_event_set_sqme_alt
     procedure :: set_weight_prc => generic_event_set_weight_prc
     procedure :: set_weight_ref => generic_event_set_weight_ref
     procedure :: set_weight_alt => generic_event_set_weight_alt
     procedure :: set_excess_prc => generic_event_set_excess_prc
     procedure :: set_n_dropped => generic_event_set_n_dropped
     procedure :: set => generic_event_set
     procedure (generic_event_write), deferred :: write
     procedure (generic_event_generate), deferred :: generate
     procedure (generic_event_set_hard_particle_set), deferred :: &
          set_hard_particle_set
     procedure (generic_event_set_index), deferred :: set_index
     procedure (generic_event_handler), deferred :: reset_index
     procedure (generic_event_increment_index), deferred :: increment_index
     procedure (generic_event_handler), deferred :: evaluate_expressions
     procedure (generic_event_select), deferred :: select
     procedure (generic_event_get_model_ptr), deferred :: get_model_ptr
     procedure (generic_event_has_index), deferred :: has_index
     procedure (generic_event_get_index), deferred :: get_index
     procedure (generic_event_get_fac_scale), deferred :: get_fac_scale
     procedure (generic_event_get_alpha_s), deferred :: get_alpha_s
     procedure (generic_event_get_sqrts), deferred :: get_sqrts
     procedure (generic_event_get_polarization), deferred :: get_polarization
     procedure (generic_event_get_beam_file), deferred :: get_beam_file
     procedure (generic_event_get_process_name), deferred :: &
          get_process_name
     procedure (generic_event_set_alpha_qcd_forced), deferred :: &
          set_alpha_qcd_forced
     procedure (generic_event_set_scale_forced), deferred :: &
          set_scale_forced
     procedure :: reset_contents => generic_event_reset_contents
     procedure :: base_reset_contents => generic_event_reset_contents
     procedure :: pacify_particle_set => generic_event_pacify_particle_set
  end type generic_event_t

  type, abstract :: event_callback_t
     private
   contains
     procedure(event_callback_write), deferred :: write
     procedure(event_callback_proc), deferred :: proc
  end type event_callback_t

  type, extends (event_callback_t) :: event_callback_nop_t
     private
   contains
     procedure :: write => event_callback_nop_write
     procedure :: proc => event_callback_nop
  end type event_callback_nop_t


  abstract interface
     subroutine generic_event_write (object, unit, &
          show_process, show_transforms, &
          show_decay, verbose, testflag)
       import
       class(generic_event_t), intent(in) :: object
       integer, intent(in), optional :: unit
       logical, intent(in), optional :: show_process
       logical, intent(in), optional :: show_transforms
       logical, intent(in), optional :: show_decay
       logical, intent(in), optional :: verbose
       logical, intent(in), optional :: testflag
     end subroutine generic_event_write
  end interface

  abstract interface
     subroutine generic_event_generate (event, i_mci, r, i_nlo)
       import
       class(generic_event_t), intent(inout) :: event
       integer, intent(in) :: i_mci
       real(default), dimension(:), intent(in), optional :: r
       integer, intent(in), optional :: i_nlo
     end subroutine generic_event_generate
  end interface

  abstract interface
     subroutine generic_event_set_hard_particle_set (event, particle_set)
       import
       class(generic_event_t), intent(inout) :: event
       type(particle_set_t), intent(in) :: particle_set
     end subroutine generic_event_set_hard_particle_set
  end interface

  abstract interface
     subroutine generic_event_set_index (event, index)
       import
       class(generic_event_t), intent(inout) :: event
       integer, intent(in) :: index
     end subroutine generic_event_set_index
  end interface

  abstract interface
     subroutine generic_event_handler (event)
       import
       class(generic_event_t), intent(inout) :: event
     end subroutine generic_event_handler
  end interface

  abstract interface
     subroutine generic_event_increment_index (event, offset)
       import
       class(generic_event_t), intent(inout) :: event
       integer, intent(in), optional :: offset
     end subroutine generic_event_increment_index
  end interface

  abstract interface
     subroutine generic_event_select (event,  i_mci, i_term, channel)
       import
       class(generic_event_t), intent(inout) :: event
       integer, intent(in) :: i_mci, i_term, channel
     end subroutine generic_event_select
  end interface

  abstract interface
     function generic_event_get_model_ptr (event) result (model)
       import
       class(generic_event_t), intent(in) :: event
       class(model_data_t), pointer :: model
     end function generic_event_get_model_ptr
  end interface

  abstract interface
     function generic_event_has_index (event) result (flag)
       import
       class(generic_event_t), intent(in) :: event
       logical :: flag
     end function generic_event_has_index
  end interface

  abstract interface
     function generic_event_get_index (event) result (index)
       import
       class(generic_event_t), intent(in) :: event
       integer :: index
     end function generic_event_get_index
  end interface

  abstract interface
     function generic_event_get_fac_scale (event) result (fac_scale)
       import
       class(generic_event_t), intent(in) :: event
       real(default) :: fac_scale
     end function generic_event_get_fac_scale
  end interface

  abstract interface
     function generic_event_get_alpha_s (event) result (alpha_s)
       import
       class(generic_event_t), intent(in) :: event
       real(default) :: alpha_s
     end function generic_event_get_alpha_s
  end interface

  abstract interface
     function generic_event_get_sqrts (event) result (sqrts)
       import
       class(generic_event_t), intent(in) :: event
       real(default) :: sqrts
     end function generic_event_get_sqrts
  end interface

  abstract interface
     function generic_event_get_polarization (event) result (pol)
       import
       class(generic_event_t), intent(in) :: event
       real(default), dimension(:), allocatable :: pol
     end function generic_event_get_polarization
  end interface

  abstract interface
     function generic_event_get_beam_file (event) result (file)
       import
       class(generic_event_t), intent(in) :: event
       type(string_t) :: file
     end function generic_event_get_beam_file
  end interface

  abstract interface
     function generic_event_get_process_name (event) result (name)
       import
       class(generic_event_t), intent(in) :: event
       type(string_t) :: name
     end function generic_event_get_process_name
  end interface

  abstract interface
     subroutine generic_event_set_alpha_qcd_forced (event, alpha_qcd)
       import
       class(generic_event_t), intent(inout) :: event
       real(default), intent(in) :: alpha_qcd
     end subroutine generic_event_set_alpha_qcd_forced
  end interface

  abstract interface
     subroutine generic_event_set_scale_forced (event, scale)
       import
       class(generic_event_t), intent(inout) :: event
       real(default), intent(in) :: scale
     end subroutine generic_event_set_scale_forced
  end interface

  abstract interface
     subroutine event_callback_write (event_callback, unit)
       import
       class(event_callback_t), intent(in) :: event_callback
       integer, intent(in), optional :: unit
     end subroutine event_callback_write
  end interface

  abstract interface
     subroutine event_callback_proc (event_callback, i, event)
       import
       class(event_callback_t), intent(in) :: event_callback
       integer(i64), intent(in) :: i
       class(generic_event_t), intent(in) :: event
     end subroutine event_callback_proc
  end interface


  interface
    module subroutine generic_event_init (event, n_alt)
      class(generic_event_t), intent(out) :: event
      integer, intent(in) :: n_alt
    end subroutine generic_event_init
    module function generic_event_has_valid_particle_set (event) result (flag)
      class(generic_event_t), intent(in) :: event
      logical :: flag
    end function generic_event_has_valid_particle_set

    module subroutine generic_event_accept_particle_set (event)
      class(generic_event_t), intent(inout) :: event
    end subroutine generic_event_accept_particle_set

    module subroutine generic_event_discard_particle_set (event)
      class(generic_event_t), intent(inout) :: event
    end subroutine generic_event_discard_particle_set
    module function generic_event_get_particle_set_ptr (event) result (ptr)
      class(generic_event_t), intent(in) :: event
      type(particle_set_t), pointer :: ptr
    end function generic_event_get_particle_set_ptr
    module subroutine generic_event_link_particle_set (event, particle_set)
      class(generic_event_t), intent(inout) :: event
      type(particle_set_t), intent(in), target :: particle_set
    end subroutine generic_event_link_particle_set
    module function generic_event_sqme_prc_is_known (event) result (flag)
      class(generic_event_t), intent(in) :: event
      logical :: flag
    end function generic_event_sqme_prc_is_known
    module function generic_event_sqme_ref_is_known (event) result (flag)
      class(generic_event_t), intent(in) :: event
      logical :: flag
    end function generic_event_sqme_ref_is_known
    module function generic_event_sqme_alt_is_known (event) result (flag)
      class(generic_event_t), intent(in) :: event
      logical :: flag
    end function generic_event_sqme_alt_is_known
    module function generic_event_weight_prc_is_known (event) result (flag)
      class(generic_event_t), intent(in) :: event
      logical :: flag
    end function generic_event_weight_prc_is_known
    module function generic_event_weight_ref_is_known (event) result (flag)
      class(generic_event_t), intent(in) :: event
      logical :: flag
    end function generic_event_weight_ref_is_known
    module function generic_event_weight_alt_is_known (event) result (flag)
      class(generic_event_t), intent(in) :: event
      logical :: flag
    end function generic_event_weight_alt_is_known
    module function generic_event_excess_prc_is_known (event) result (flag)
      class(generic_event_t), intent(in) :: event
      logical :: flag
    end function generic_event_excess_prc_is_known
    module function generic_event_get_n_alt (event) result (n)
      class(generic_event_t), intent(in) :: event
      integer :: n
    end function generic_event_get_n_alt
    module function generic_event_get_sqme_prc (event) result (sqme)
      class(generic_event_t), intent(in) :: event
      real(default) :: sqme
    end function generic_event_get_sqme_prc
    module function generic_event_get_sqme_ref (event) result (sqme)
      class(generic_event_t), intent(in) :: event
      real(default) :: sqme
    end function generic_event_get_sqme_ref
    module function generic_event_get_sqme_alt_0 (event, i) result (sqme)
      class(generic_event_t), intent(in) :: event
      integer, intent(in) :: i
      real(default) :: sqme
    end function generic_event_get_sqme_alt_0
    module function generic_event_get_sqme_alt_1 (event) result (sqme)
      class(generic_event_t), intent(in) :: event
      real(default), dimension(event%n_alt) :: sqme
    end function generic_event_get_sqme_alt_1
    module function generic_event_get_weight_prc (event) result (weight)
      class(generic_event_t), intent(in) :: event
      real(default) :: weight
    end function generic_event_get_weight_prc
    module function generic_event_get_weight_ref (event) result (weight)
      class(generic_event_t), intent(in) :: event
      real(default) :: weight
    end function generic_event_get_weight_ref
    module function generic_event_get_weight_alt_0 (event, i) result (weight)
      class(generic_event_t), intent(in) :: event
      integer, intent(in) :: i
      real(default) :: weight
    end function generic_event_get_weight_alt_0
    module function generic_event_get_weight_alt_1 (event) result (weight)
      class(generic_event_t), intent(in) :: event
      real(default), dimension(event%n_alt) :: weight
    end function generic_event_get_weight_alt_1
    module function generic_event_get_excess_prc (event) result (excess)
      class(generic_event_t), intent(in) :: event
      real(default) :: excess
    end function generic_event_get_excess_prc
    module function generic_event_get_n_dropped (event) result (n_dropped)
      class(generic_event_t), intent(in) :: event
      integer :: n_dropped
    end function generic_event_get_n_dropped
    module subroutine generic_event_set_sqme_prc (event, sqme)
      class(generic_event_t), intent(inout) :: event
      real(default), intent(in) :: sqme
    end subroutine generic_event_set_sqme_prc
    module subroutine generic_event_set_sqme_ref (event, sqme)
      class(generic_event_t), intent(inout) :: event
      real(default), intent(in) :: sqme
    end subroutine generic_event_set_sqme_ref
    module subroutine generic_event_set_sqme_alt (event, sqme)
      class(generic_event_t), intent(inout) :: event
      real(default), dimension(:), intent(in) :: sqme
    end subroutine generic_event_set_sqme_alt
    module subroutine generic_event_set_weight_prc (event, weight)
      class(generic_event_t), intent(inout) :: event
      real(default), intent(in) :: weight
    end subroutine generic_event_set_weight_prc
    module subroutine generic_event_set_weight_ref (event, weight)
      class(generic_event_t), intent(inout) :: event
      real(default), intent(in) :: weight
    end subroutine generic_event_set_weight_ref
    module subroutine generic_event_set_weight_alt (event, weight)
      class(generic_event_t), intent(inout) :: event
      real(default), dimension(:), intent(in) :: weight
    end subroutine generic_event_set_weight_alt
    module subroutine generic_event_set_excess_prc (event, excess)
      class(generic_event_t), intent(inout) :: event
      real(default), intent(in) :: excess
    end subroutine generic_event_set_excess_prc
    module subroutine generic_event_set_n_dropped (event, n_dropped)
      class(generic_event_t), intent(inout) :: event
      integer, intent(in) :: n_dropped
    end subroutine generic_event_set_n_dropped
    module subroutine generic_event_set (event, &
         weight_ref, weight_prc, weight_alt, &
         excess_prc, n_dropped, &
         sqme_ref, sqme_prc, sqme_alt)
      class(generic_event_t), intent(inout) :: event
      real(default), intent(in), optional :: weight_ref, weight_prc
      real(default), intent(in), optional :: sqme_ref, sqme_prc
      real(default), dimension(:), intent(in), optional :: sqme_alt, weight_alt
      real(default), intent(in), optional :: excess_prc
      integer, intent(in), optional :: n_dropped
    end subroutine generic_event_set
    module subroutine generic_event_reset_contents (event)
      class(generic_event_t), intent(inout) :: event
    end subroutine generic_event_reset_contents
    module subroutine generic_event_pacify_particle_set (event)
      class(generic_event_t), intent(inout) :: event
    end subroutine generic_event_pacify_particle_set
    module function event_normalization_mode (string, unweighted) result (mode)
      integer :: mode
      type(string_t), intent(in) :: string
      logical, intent(in) :: unweighted
    end function event_normalization_mode
    module function event_normalization_string (norm_mode) result (string)
      integer, intent(in) :: norm_mode
      type(string_t) :: string
    end function event_normalization_string
    module subroutine event_normalization_update &
         (weight, sigma, n, mode_new, mode_old)
      real(default), intent(inout) :: weight
      real(default), intent(in) :: sigma
      integer, intent(in) :: n
      integer, intent(in) :: mode_new, mode_old
    end subroutine event_normalization_update
    module subroutine event_callback_nop_write (event_callback, unit)
      class(event_callback_nop_t), intent(in) :: event_callback
      integer, intent(in), optional :: unit
    end subroutine event_callback_nop_write
    module subroutine event_callback_nop (event_callback, i, event)
      class(event_callback_nop_t), intent(in) :: event_callback
      integer(i64), intent(in) :: i
      class(generic_event_t), intent(in) :: event
    end subroutine event_callback_nop
  end interface

end module event_base
