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

module variables

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use numeric_utils, only: pacify
  use os_interface, only: paths_t
  use pdg_arrays
  use subevents
  use var_base

  implicit none
  private

  public :: obs_unary_int
  public :: obs_unary_real
  public :: obs_binary_int
  public :: obs_binary_real
  public :: obs_sev_int
  public :: obs_sev_real
  public :: var_list_t

  integer, parameter, public :: V_NONE = 0, V_LOG = 1, V_INT = 2, V_REAL = 3
  integer, parameter, public :: V_CMPLX = 4, V_SEV = 5, V_PDG = 6, V_STR = 7
  integer, parameter, public :: V_OBS1_INT = 11, V_OBS2_INT = 12
  integer, parameter, public :: V_OBS1_REAL = 21, V_OBS2_REAL = 22
  integer, parameter, public :: V_OBSEV_INT = 13, V_OBSEV_REAL = 23
  integer, parameter, public :: V_UOBS1_INT = 31, V_UOBS2_INT = 32
  integer, parameter, public :: V_UOBS1_REAL = 41, V_UOBS2_REAL = 42


  type :: var_entry_t
     private
     integer :: type = V_NONE
     type(string_t) :: name
     logical :: is_allocated = .false.
     logical :: is_defined = .false.
     logical :: is_locked = .false.
     logical :: is_intrinsic = .false.
     logical :: is_user_var = .false.
     logical, pointer :: is_known => null ()
     logical,           pointer :: lval => null ()
     integer,           pointer :: ival => null ()
     real(default),     pointer :: rval => null ()
     complex(default), pointer :: cval => null ()
     type(subevt_t),  pointer :: pval => null ()
     type(pdg_array_t), pointer :: aval => null ()
     type(string_t),    pointer :: sval => null ()
     procedure(obs_unary_int),   nopass, pointer :: obs1_int  => null ()
     procedure(obs_unary_real),  nopass, pointer :: obs1_real => null ()
     procedure(obs_binary_int),  nopass, pointer :: obs2_int  => null ()
     procedure(obs_binary_real), nopass, pointer :: obs2_real => null ()
     procedure(obs_sev_int),  nopass, pointer :: obsev_int  => null ()
     procedure(obs_sev_real), nopass, pointer :: obsev_real => null ()
     type(prt_t), pointer :: prt1 => null ()
     type(prt_t), pointer :: prt2 => null ()
     type(var_entry_t), pointer :: next => null ()
     type(var_entry_t), pointer :: previous => null ()
     type(string_t) :: description
  end type var_entry_t

  type, extends (vars_t) :: var_list_t
     private
     type(var_entry_t), pointer :: first => null ()
     type(var_entry_t), pointer :: last => null ()
     type(var_list_t), pointer :: next => null ()
   contains
     procedure :: link => var_list_link
     procedure :: sort => var_list_sort
     procedure :: get_previous => var_list_get_previous
     procedure :: swap_with_next => var_list_swap_with_next
     generic :: append_log => var_list_append_log_s, var_list_append_log_c
     procedure, private :: var_list_append_log_s
     procedure, private :: var_list_append_log_c
     generic :: append_int => var_list_append_int_s, var_list_append_int_c
     procedure, private :: var_list_append_int_s
     procedure, private :: var_list_append_int_c
     generic :: append_real => var_list_append_real_s, var_list_append_real_c
     procedure, private :: var_list_append_real_s
     procedure, private :: var_list_append_real_c
     generic :: append_cmplx => var_list_append_cmplx_s, var_list_append_cmplx_c
     procedure, private :: var_list_append_cmplx_s
     procedure, private :: var_list_append_cmplx_c
     generic :: append_subevt => var_list_append_subevt_s, var_list_append_subevt_c
     procedure, private :: var_list_append_subevt_s
     procedure, private :: var_list_append_subevt_c
     generic :: append_pdg_array => var_list_append_pdg_array_s, var_list_append_pdg_array_c
     procedure, private :: var_list_append_pdg_array_s
     procedure, private :: var_list_append_pdg_array_c
     generic :: append_string => var_list_append_string_s, var_list_append_string_c
     procedure, private :: var_list_append_string_s
     procedure, private :: var_list_append_string_c
     procedure :: append_log_ptr => var_list_append_log_ptr
     procedure :: append_int_ptr => var_list_append_int_ptr
     procedure :: append_real_ptr => var_list_append_real_ptr
     procedure :: append_cmplx_ptr => var_list_append_cmplx_ptr
     procedure :: append_pdg_array_ptr => var_list_append_pdg_array_ptr
     procedure :: append_subevt_ptr => var_list_append_subevt_ptr
     procedure :: append_string_ptr => var_list_append_string_ptr
     procedure :: final => var_list_final
     procedure :: write => var_list_write
     procedure :: write_var => var_list_write_var
     procedure :: get_type => var_list_get_type
     procedure :: contains => var_list_exists
     procedure :: is_intrinsic => var_list_is_intrinsic
     procedure :: is_known => var_list_is_known
     procedure :: is_locked => var_list_is_locked
     procedure :: get_var_properties => var_list_get_var_properties
     procedure :: get_lval => var_list_get_lval
     procedure :: get_ival => var_list_get_ival
     procedure :: get_rval => var_list_get_rval
     procedure :: get_cval => var_list_get_cval
     procedure :: get_pval => var_list_get_pval
     procedure :: get_aval => var_list_get_aval
     procedure :: get_sval => var_list_get_sval
     procedure :: get_lptr => var_list_get_lptr
     procedure :: get_iptr => var_list_get_iptr
     procedure :: get_rptr => var_list_get_rptr
     procedure :: get_cptr => var_list_get_cptr
     procedure :: get_aptr => var_list_get_aptr
     procedure :: get_pptr => var_list_get_pptr
     procedure :: get_sptr => var_list_get_sptr
     procedure :: get_obs1_iptr => var_list_get_obs1_iptr
     procedure :: get_obs2_iptr => var_list_get_obs2_iptr
     procedure :: get_obsev_iptr => var_list_get_obsev_iptr
     procedure :: get_obs1_rptr => var_list_get_obs1_rptr
     procedure :: get_obs2_rptr => var_list_get_obs2_rptr
     procedure :: get_obsev_rptr => var_list_get_obsev_rptr
     procedure :: set_procvar_int => var_list_set_procvar_int
     procedure :: set_procvar_real => var_list_set_procvar_real
     procedure :: append_obs1_iptr => var_list_append_obs1_iptr
     procedure :: append_obs2_iptr => var_list_append_obs2_iptr
     procedure :: append_obs1_rptr => var_list_append_obs1_rptr
     procedure :: append_obs2_rptr => var_list_append_obs2_rptr
     procedure :: append_obsev_iptr => var_list_append_obsev_iptr
     procedure :: append_obsev_rptr => var_list_append_obsev_rptr
     procedure :: append_uobs_int => var_list_append_uobs_int
     procedure :: append_uobs_real => var_list_append_uobs_real
     procedure :: unset => var_list_clear
     procedure :: set_ival => var_list_set_ival
     procedure :: set_rval => var_list_set_rval
     procedure :: set_cval => var_list_set_cval
     procedure :: set_lval => var_list_set_lval
     procedure :: set_sval => var_list_set_sval
     procedure :: set_log => var_list_set_log
     procedure :: set_int => var_list_set_int
     procedure :: set_real => var_list_set_real
     procedure :: set_cmplx => var_list_set_cmplx
     procedure :: set_subevt => var_list_set_subevt
     procedure :: set_pdg_array => var_list_set_pdg_array
     procedure :: set_string => var_list_set_string
     procedure :: import => var_list_import
     procedure :: undefine => var_list_undefine
     procedure :: init_snapshot => var_list_init_snapshot
     procedure :: check_user_var => var_list_check_user_var
     procedure :: init_defaults => var_list_init_defaults
     procedure :: set_beams_defaults => var_list_set_beams_defaults
     procedure :: set_core_defaults => var_list_set_core_defaults
     procedure :: set_integration_defaults => var_list_set_integration_defaults
     procedure :: set_phase_space_defaults => var_list_set_phase_space_defaults
     procedure :: set_gamelan_defaults => var_list_set_gamelan_defaults
     procedure :: set_clustering_defaults => var_list_set_clustering_defaults
     procedure :: set_isolation_recomb_defaults => &
          var_list_set_isolation_recomb_defaults
     procedure :: set_eio_defaults => var_list_set_eio_defaults
     procedure :: set_shower_defaults => var_list_set_shower_defaults
     procedure :: set_hadronization_defaults => var_list_set_hadronization_defaults
     procedure :: set_tauola_defaults => var_list_set_tauola_defaults
     procedure :: set_mlm_matching_defaults => var_list_set_mlm_matching_defaults
     procedure :: set_powheg_matching_defaults => &
          var_list_set_powheg_matching_defaults
     procedure :: set_openmp_defaults => var_list_set_openmp_defaults
     procedure :: set_mpi_defaults => var_list_set_mpi_defaults
     procedure :: set_nlo_defaults => var_list_set_nlo_defaults
  end type var_list_t


  abstract interface
     function obs_unary_int (prt1) result (ival)
       import
       integer :: ival
       type(prt_t), intent(in) :: prt1
     end function obs_unary_int
  end interface
  abstract interface
     function obs_unary_real (prt1) result (rval)
       import
       real(default) :: rval
       type(prt_t), intent(in) :: prt1
     end function obs_unary_real
  end interface
  abstract interface
     function obs_binary_int (prt1, prt2) result (ival)
       import
       integer :: ival
       type(prt_t), intent(in) :: prt1, prt2
     end function obs_binary_int
  end interface
  abstract interface
     function obs_binary_real (prt1, prt2) result (rval)
       import
       real(default) :: rval
       type(prt_t), intent(in) :: prt1, prt2
     end function obs_binary_real
  end interface
  abstract interface
     function obs_sev_int (sev) result (ival)
       import
       integer :: ival
       type(subevt_t), intent(in) :: sev
     end function obs_sev_int
  end interface
  abstract interface
     function obs_sev_real (sev) result (rval)
       import
       real(default) :: rval
       type(subevt_t), intent(in) :: sev
     end function obs_sev_real
  end interface

  interface var_list_append_log
     module procedure var_list_append_log_s
     module procedure var_list_append_log_c
  end interface
  interface var_list_append_int
     module procedure var_list_append_int_s
     module procedure var_list_append_int_c
  end interface
  interface var_list_append_real
     module procedure var_list_append_real_s
     module procedure var_list_append_real_c
  end interface
  interface var_list_append_cmplx
     module procedure var_list_append_cmplx_s
     module procedure var_list_append_cmplx_c
  end interface
  interface var_list_append_subevt
     module procedure var_list_append_subevt_s
     module procedure var_list_append_subevt_c
  end interface
  interface var_list_append_pdg_array
     module procedure var_list_append_pdg_array_s
     module procedure var_list_append_pdg_array_c
  end interface
  interface var_list_append_string
     module procedure var_list_append_string_s
     module procedure var_list_append_string_c
  end interface

  interface
    module subroutine var_list_link (vars, target_vars)
      class(var_list_t), intent(inout) :: vars
      class(vars_t), intent(in), target :: target_vars
    end subroutine var_list_link
    module subroutine var_list_sort (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_sort
    module function var_list_get_previous &
         (var_list, var_entry) result (previous)
      type(var_entry_t), pointer :: previous
      class(var_list_t), intent(in) :: var_list
      type(var_entry_t), intent(in) :: var_entry
    end function var_list_get_previous
    module subroutine var_list_swap_with_next (var_list, var_entry)
      class(var_list_t), intent(inout) :: var_list
      type(var_entry_t), intent(in) :: var_entry
    end subroutine var_list_swap_with_next
    module subroutine var_list_append_log_s &
         (var_list, name, lval, locked, verbose, intrinsic, user, description)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      logical, intent(in), optional :: lval
      logical, intent(in), optional :: locked, verbose, intrinsic, user
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_log_s
    module subroutine var_list_append_int_s &
         (var_list, name, ival, locked, verbose, intrinsic, user, description)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      integer, intent(in), optional :: ival
      logical, intent(in), optional :: locked, verbose, intrinsic, user
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_int_s
    module subroutine var_list_append_real_s &
         (var_list, name, rval, locked, verbose, intrinsic, user, description)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      real(default), intent(in), optional :: rval
      logical, intent(in), optional :: locked, verbose, intrinsic, user
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_real_s
    module subroutine var_list_append_cmplx_s &
         (var_list, name, cval, locked, verbose, intrinsic, user, description)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      complex(default), intent(in), optional :: cval
      logical, intent(in), optional :: locked, verbose, intrinsic, user
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_cmplx_s
    module subroutine var_list_append_subevt_s &
         (var_list, name, pval, locked, verbose, intrinsic, user, description)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      type(subevt_t), intent(in), optional :: pval
      logical, intent(in), optional :: locked, verbose, intrinsic, user
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_subevt_s
    module subroutine var_list_append_pdg_array_s &
         (var_list, name, aval, locked, verbose, intrinsic, user, description)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      type(pdg_array_t), intent(in), optional :: aval
      logical, intent(in), optional :: locked, verbose, intrinsic, user
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_pdg_array_s
    module subroutine var_list_append_string_s &
         (var_list, name, sval, locked, verbose, intrinsic, user, description)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      type(string_t), intent(in), optional :: sval
      logical, intent(in), optional :: locked, verbose, intrinsic, user
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_string_s
    module subroutine var_list_append_log_c &
         (var_list, name, lval, locked, verbose, intrinsic, user, description)
      class(var_list_t), intent(inout) :: var_list
      character(*), intent(in) :: name
      logical, intent(in), optional :: lval
      logical, intent(in), optional :: locked, verbose, intrinsic, user
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_log_c
    module subroutine var_list_append_int_c &
         (var_list, name, ival, locked, verbose, intrinsic, user, description)
      class(var_list_t), intent(inout) :: var_list
      character(*), intent(in) :: name
      integer, intent(in), optional :: ival
      logical, intent(in), optional :: locked, verbose, intrinsic, user
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_int_c
    module subroutine var_list_append_real_c &
         (var_list, name, rval, locked, verbose, intrinsic, user, description)
      class(var_list_t), intent(inout) :: var_list
      character(*), intent(in) :: name
      real(default), intent(in), optional :: rval
      logical, intent(in), optional :: locked, verbose, intrinsic, user
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_real_c
    module subroutine var_list_append_cmplx_c &
         (var_list, name, cval, locked, verbose, intrinsic, user, description)
      class(var_list_t), intent(inout) :: var_list
      character(*), intent(in) :: name
      complex(default), intent(in), optional :: cval
      logical, intent(in), optional :: locked, verbose, intrinsic, user
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_cmplx_c
    module subroutine var_list_append_subevt_c &
         (var_list, name, pval, locked, verbose, intrinsic, user, description)
      class(var_list_t), intent(inout) :: var_list
      character(*), intent(in) :: name
      type(subevt_t), intent(in), optional :: pval
      logical, intent(in), optional :: locked, verbose, intrinsic, user
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_subevt_c
    module subroutine var_list_append_pdg_array_c &
         (var_list, name, aval, locked, verbose, intrinsic, user, description)
      class(var_list_t), intent(inout) :: var_list
      character(*), intent(in) :: name
      type(pdg_array_t), intent(in), optional :: aval
      logical, intent(in), optional :: locked, verbose, intrinsic, user
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_pdg_array_c
    module subroutine var_list_append_string_c &
         (var_list, name, sval, locked, verbose, intrinsic, user, description)
      class(var_list_t), intent(inout) :: var_list
      character(*), intent(in) :: name
      character(*), intent(in), optional :: sval
      logical, intent(in), optional :: locked, verbose, intrinsic, user
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_string_c
    module subroutine var_list_append_log_ptr &
         (var_list, name, lval, is_known, locked, verbose, &
          intrinsic, description)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      logical, intent(in), target :: lval
      logical, intent(in), target :: is_known
      logical, intent(in), optional :: locked, verbose, intrinsic
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_log_ptr
    module subroutine var_list_append_int_ptr &
         (var_list, name, ival, is_known, locked, verbose, &
          intrinsic, description)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      integer, intent(in), target :: ival
      logical, intent(in), target :: is_known
      logical, intent(in), optional :: locked, verbose, intrinsic
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_int_ptr
    module subroutine var_list_append_real_ptr &
         (var_list, name, rval, is_known, locked, verbose, &
          intrinsic, description)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      real(default), intent(in), target :: rval
      logical, intent(in), target :: is_known
      logical, intent(in), optional :: locked, verbose, intrinsic
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_real_ptr
    module subroutine var_list_append_cmplx_ptr &
         (var_list, name, cval, is_known, locked, verbose, &
          intrinsic, description)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      complex(default), intent(in), target :: cval
      logical, intent(in), target :: is_known
      logical, intent(in), optional :: locked, verbose, intrinsic
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_cmplx_ptr
    module subroutine var_list_append_pdg_array_ptr &
         (var_list, name, aval, is_known, locked, verbose, &
          intrinsic, description)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      type(pdg_array_t), intent(in), target :: aval
      logical, intent(in), target :: is_known
      logical, intent(in), optional :: locked, verbose, intrinsic
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_pdg_array_ptr
    module subroutine var_list_append_subevt_ptr &
         (var_list, name, pval, is_known, locked, verbose, &
          intrinsic, description)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      type(subevt_t), intent(in), target :: pval
      logical, intent(in), target :: is_known
      logical, intent(in), optional :: locked, verbose, intrinsic
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_subevt_ptr
    module subroutine var_list_append_string_ptr &
         (var_list, name, sval, is_known, locked, verbose, &
          intrinsic, description)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      type(string_t), intent(in), target :: sval
      logical, intent(in), target :: is_known
      logical, intent(in), optional :: locked, verbose, intrinsic
      type(string_t), intent(in), optional :: description
    end subroutine var_list_append_string_ptr
    recursive module subroutine var_list_final (vars, follow_link)
      class(var_list_t), intent(inout) :: vars
      logical, intent(in), optional :: follow_link
    end subroutine var_list_final
    recursive module subroutine var_list_write &
         (var_list, unit, follow_link, only_type, prefix, model_name, &
          intrinsic, pacified, descriptions, ascii_output)
      class(var_list_t), intent(in), target :: var_list
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: follow_link
      integer, intent(in), optional :: only_type
      character(*), intent(in), optional :: prefix
      type(string_t), intent(in), optional :: model_name
      logical, intent(in), optional :: intrinsic
      logical, intent(in), optional :: pacified
      logical, intent(in), optional :: descriptions
      logical, intent(in), optional :: ascii_output
    end subroutine var_list_write
    recursive module subroutine var_list_write_var &
         (var_list, name, unit, type, follow_link, &
         model_name, pacified, defined, descriptions, ascii_output)
      class(var_list_t), intent(in), target :: var_list
      type(string_t), intent(in) :: name
      integer, intent(in), optional :: unit
      integer, intent(in), optional :: type
      logical, intent(in), optional :: follow_link
      type(string_t), intent(in), optional :: model_name
      logical, intent(in), optional :: pacified
      logical, intent(in), optional :: defined
      logical, intent(in), optional :: descriptions
      logical, intent(in), optional :: ascii_output
    end subroutine var_list_write_var
    module function var_list_get_type &
         (var_list, name, follow_link) result (type)
      class(var_list_t), intent(in), target :: var_list
      type(string_t), intent(in) :: name
      logical, intent(in), optional :: follow_link
      integer :: type
    end function var_list_get_type
    module function var_list_exists (vars, name, follow_link) result (lval)
      logical :: lval
      type(string_t), intent(in) :: name
      class(var_list_t), intent(in) :: vars
      logical, intent(in), optional :: follow_link
    end function var_list_exists
    module function var_list_is_intrinsic &
         (vars, name, follow_link) result (lval)
      logical :: lval
      type(string_t), intent(in) :: name
      class(var_list_t), intent(in) :: vars
      logical, intent(in), optional :: follow_link
    end function var_list_is_intrinsic
    module function var_list_is_known (vars, name, follow_link) result (lval)
      logical :: lval
      type(string_t), intent(in) :: name
      class(var_list_t), intent(in) :: vars
      logical, intent(in), optional :: follow_link
    end function var_list_is_known
    module function var_list_is_locked (vars, name, follow_link) result (lval)
      logical :: lval
      type(string_t), intent(in) :: name
      class(var_list_t), intent(in) :: vars
      logical, intent(in), optional :: follow_link
    end function var_list_is_locked
    module subroutine var_list_get_var_properties &
         (vars, name, req_type, follow_link, &
          type, is_defined, is_known, is_locked)
      class(var_list_t), intent(in) :: vars
      type(string_t), intent(in) :: name
      integer, intent(in), optional :: req_type
      logical, intent(in), optional :: follow_link
      integer, intent(out), optional :: type
      logical, intent(out), optional :: is_defined, is_known, is_locked
    end subroutine var_list_get_var_properties
    module function var_list_get_lval (vars, name, follow_link) result (lval)
      logical :: lval
      type(string_t), intent(in) :: name
      class(var_list_t), intent(in) :: vars
      logical, intent(in), optional :: follow_link
    end function var_list_get_lval
    module function var_list_get_ival (vars, name, follow_link) result (ival)
      integer :: ival
      type(string_t), intent(in) :: name
      class(var_list_t), intent(in) :: vars
      logical, intent(in), optional :: follow_link
    end function var_list_get_ival
    module function var_list_get_rval (vars, name, follow_link) result (rval)
      real(default) :: rval
      type(string_t), intent(in) :: name
      class(var_list_t), intent(in) :: vars
      logical, intent(in), optional :: follow_link
    end function var_list_get_rval
    module function var_list_get_cval (vars, name, follow_link) result (cval)
      complex(default) :: cval
      type(string_t), intent(in) :: name
      class(var_list_t), intent(in) :: vars
      logical, intent(in), optional :: follow_link
    end function var_list_get_cval
    module function var_list_get_aval (vars, name, follow_link) result (aval)
      type(pdg_array_t) :: aval
      type(string_t), intent(in) :: name
      class(var_list_t), intent(in) :: vars
      logical, intent(in), optional :: follow_link
    end function var_list_get_aval
    module function var_list_get_pval (vars, name, follow_link) result (pval)
      type(subevt_t) :: pval
      type(string_t), intent(in) :: name
      class(var_list_t), intent(in) :: vars
      logical, intent(in), optional :: follow_link
    end function var_list_get_pval
    module function var_list_get_sval (vars, name, follow_link) result (sval)
      type(string_t) :: sval
      type(string_t), intent(in) :: name
      class(var_list_t), intent(in) :: vars
      logical, intent(in), optional :: follow_link
    end function var_list_get_sval
    module subroutine var_list_get_lptr (var_list, name, lptr, known)
      class(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: name
      logical, pointer, intent(out) :: lptr
      logical, pointer, intent(out), optional :: known
    end subroutine var_list_get_lptr
    module subroutine var_list_get_iptr (var_list, name, iptr, known)
      class(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: name
      integer, pointer, intent(out) :: iptr
      logical, pointer, intent(out), optional :: known
    end subroutine var_list_get_iptr
    module subroutine var_list_get_rptr (var_list, name, rptr, known)
      class(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: name
      real(default), pointer, intent(out) :: rptr
      logical, pointer, intent(out), optional :: known
    end subroutine var_list_get_rptr
    module subroutine var_list_get_cptr (var_list, name, cptr, known)
      class(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: name
      complex(default), pointer, intent(out) :: cptr
      logical, pointer, intent(out), optional :: known
    end subroutine var_list_get_cptr
    module subroutine var_list_get_aptr (var_list, name, aptr, known)
      class(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: name
      type(pdg_array_t), pointer, intent(out) :: aptr
      logical, pointer, intent(out), optional :: known
      type(var_entry_t), pointer :: var
    end subroutine var_list_get_aptr
    module subroutine var_list_get_pptr (var_list, name, pptr, known)
      class(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: name
      type(subevt_t), pointer, intent(out) :: pptr
      logical, pointer, intent(out), optional :: known
    end subroutine var_list_get_pptr
    module subroutine var_list_get_sptr (var_list, name, sptr, known)
      class(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: name
      type(string_t), pointer, intent(out) :: sptr
      logical, pointer, intent(out), optional :: known
    end subroutine var_list_get_sptr
    module subroutine var_list_get_obs1_iptr (var_list, name, obs1_iptr, p1)
      class(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: name
      procedure(obs_unary_int), pointer, intent(out) :: obs1_iptr
      type(prt_t), pointer, intent(out) :: p1
    end subroutine var_list_get_obs1_iptr
    module subroutine var_list_get_obs2_iptr (var_list, name, obs2_iptr, p1, p2)
      class(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: name
      procedure(obs_binary_int), pointer, intent(out) :: obs2_iptr
      type(prt_t), pointer, intent(out) :: p1, p2
    end subroutine var_list_get_obs2_iptr
    module subroutine var_list_get_obsev_iptr (var_list, name, obsev_iptr, pval)
      class(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: name
      procedure(obs_sev_int), pointer, intent(out) :: obsev_iptr
      type(subevt_t), pointer, intent(out) :: pval
    end subroutine var_list_get_obsev_iptr
    module subroutine var_list_get_obs1_rptr (var_list, name, obs1_rptr, p1)
      class(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: name
      procedure(obs_unary_real), pointer, intent(out) :: obs1_rptr
      type(prt_t), pointer, intent(out) :: p1
    end subroutine var_list_get_obs1_rptr
    module subroutine var_list_get_obs2_rptr (var_list, name, obs2_rptr, p1, p2)
      class(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: name
      procedure(obs_binary_real), pointer, intent(out) :: obs2_rptr
      type(prt_t), pointer, intent(out) :: p1, p2
    end subroutine var_list_get_obs2_rptr
    module subroutine var_list_get_obsev_rptr (var_list, name, obsev_rptr, pval)
      class(var_list_t), intent(in) :: var_list
      type(string_t), intent(in) :: name
      procedure(obs_sev_real), pointer, intent(out) :: obsev_rptr
      type(subevt_t), pointer, intent(out) :: pval
    end subroutine var_list_get_obsev_rptr
    module subroutine var_list_set_procvar_int (var_list, proc_id, name, ival)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: proc_id
      type(string_t), intent(in) :: name
      integer, intent(in), optional :: ival
    end subroutine var_list_set_procvar_int
    module subroutine var_list_set_procvar_real (var_list, proc_id, name, rval)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: proc_id
      type(string_t), intent(in) :: name
      real(default), intent(in), optional :: rval
    end subroutine var_list_set_procvar_real
    module subroutine var_list_append_obs1_iptr &
         (var_list, name, obs1_iptr, p1)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      procedure(obs_unary_int) :: obs1_iptr
      type(prt_t), intent(in), target :: p1
    end subroutine var_list_append_obs1_iptr
    module subroutine var_list_append_obs2_iptr &
         (var_list, name, obs2_iptr, p1, p2)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      procedure(obs_binary_int) :: obs2_iptr
      type(prt_t), intent(in), target :: p1, p2
    end subroutine var_list_append_obs2_iptr
    module subroutine var_list_append_obsev_iptr &
         (var_list, name, obsev_iptr, sev)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      procedure(obs_sev_int) :: obsev_iptr
      type(subevt_t), intent(in), target :: sev
    end subroutine var_list_append_obsev_iptr
    module subroutine var_list_append_obs1_rptr &
         (var_list, name, obs1_rptr, p1)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      procedure(obs_unary_real) :: obs1_rptr
      type(prt_t), intent(in), target :: p1
    end subroutine var_list_append_obs1_rptr
    module subroutine var_list_append_obs2_rptr &
         (var_list, name, obs2_rptr, p1, p2)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      procedure(obs_binary_real) :: obs2_rptr
      type(prt_t), intent(in), target :: p1, p2
    end subroutine var_list_append_obs2_rptr
    module subroutine var_list_append_obsev_rptr &
         (var_list, name, obsev_rptr, sev)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      procedure(obs_sev_real) :: obsev_rptr
      type(subevt_t), intent(in), target :: sev
    end subroutine var_list_append_obsev_rptr
    module subroutine var_list_append_uobs_int (var_list, name, p1, p2)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      type(prt_t), intent(in), target :: p1
      type(prt_t), intent(in), target, optional :: p2
    end subroutine var_list_append_uobs_int
    module subroutine var_list_append_uobs_real (var_list, name, p1, p2)
      class(var_list_t), intent(inout) :: var_list
      type(string_t), intent(in) :: name
      type(prt_t), intent(in), target :: p1
      type(prt_t), intent(in), target, optional :: p2
    end subroutine var_list_append_uobs_real
    module subroutine var_list_clear (vars, name, follow_link)
      class(var_list_t), intent(inout) :: vars
      type(string_t), intent(in) :: name
      logical, intent(in), optional :: follow_link
    end subroutine var_list_clear
    module subroutine var_list_set_ival (vars, name, ival, follow_link)
      class(var_list_t), intent(inout) :: vars
      type(string_t), intent(in) :: name
      integer, intent(in) :: ival
      logical, intent(in), optional :: follow_link
    end subroutine var_list_set_ival
    module subroutine var_list_set_rval (vars, name, rval, follow_link)
      class(var_list_t), intent(inout) :: vars
      type(string_t), intent(in) :: name
      real(default), intent(in) :: rval
      logical, intent(in), optional :: follow_link
    end subroutine var_list_set_rval
    module subroutine var_list_set_cval (vars, name, cval, follow_link)
      class(var_list_t), intent(inout) :: vars
      type(string_t), intent(in) :: name
      complex(default), intent(in) :: cval
      logical, intent(in), optional :: follow_link
    end subroutine var_list_set_cval
    module subroutine var_list_set_lval (vars, name, lval, follow_link)
      class(var_list_t), intent(inout) :: vars
      type(string_t), intent(in) :: name
      logical, intent(in) :: lval
      logical, intent(in), optional :: follow_link
    end subroutine var_list_set_lval
    module subroutine var_list_set_sval (vars, name, sval, follow_link)
      class(var_list_t), intent(inout) :: vars
      type(string_t), intent(in) :: name
      type(string_t), intent(in) :: sval
      logical, intent(in), optional :: follow_link
    end subroutine var_list_set_sval
    module subroutine var_list_set_log &
         (var_list, name, lval, is_known, ignore, force, verbose, model_name)
      class(var_list_t), intent(inout), target :: var_list
      type(string_t), intent(in) :: name
      logical, intent(in) :: lval
      logical, intent(in) :: is_known
      logical, intent(in), optional :: ignore, force, verbose
      type(string_t), intent(in), optional :: model_name
    end subroutine var_list_set_log
    module subroutine var_list_set_int &
         (var_list, name, ival, is_known, ignore, force, verbose, model_name)
      class(var_list_t), intent(inout), target :: var_list
      type(string_t), intent(in) :: name
      integer, intent(in) :: ival
      logical, intent(in) :: is_known
      logical, intent(in), optional :: ignore, force, verbose
      type(string_t), intent(in), optional :: model_name
    end subroutine var_list_set_int
    module subroutine var_list_set_real &
         (var_list, name, rval, is_known, ignore, force, &
          verbose, model_name, pacified)
      class(var_list_t), intent(inout), target :: var_list
      type(string_t), intent(in) :: name
      real(default), intent(in) :: rval
      logical, intent(in) :: is_known
      logical, intent(in), optional :: ignore, force, verbose, pacified
      type(string_t), intent(in), optional :: model_name
    end subroutine var_list_set_real
    module subroutine var_list_set_cmplx &
         (var_list, name, cval, is_known, ignore, force, &
          verbose, model_name, pacified)
      class(var_list_t), intent(inout), target :: var_list
      type(string_t), intent(in) :: name
      complex(default), intent(in) :: cval
      logical, intent(in) :: is_known
      logical, intent(in), optional :: ignore, force, verbose, pacified
      type(string_t), intent(in), optional :: model_name
    end subroutine var_list_set_cmplx
    module subroutine var_list_set_pdg_array &
         (var_list, name, aval, is_known, ignore, force, verbose, model_name)
      class(var_list_t), intent(inout), target :: var_list
      type(string_t), intent(in) :: name
      type(pdg_array_t), intent(in) :: aval
      logical, intent(in) :: is_known
      logical, intent(in), optional :: ignore, force, verbose
      type(string_t), intent(in), optional :: model_name
    end subroutine var_list_set_pdg_array
    module subroutine var_list_set_subevt &
         (var_list, name, pval, is_known, ignore, force, verbose, model_name)
      class(var_list_t), intent(inout), target :: var_list
      type(string_t), intent(in) :: name
      type(subevt_t), intent(in) :: pval
      logical, intent(in) :: is_known
      logical, intent(in), optional :: ignore, force, verbose
      type(string_t), intent(in), optional :: model_name
    end subroutine var_list_set_subevt
    module subroutine var_list_set_string &
         (var_list, name, sval, is_known, ignore, force, verbose, model_name)
      class(var_list_t), intent(inout), target :: var_list
      type(string_t), intent(in) :: name
      type(string_t), intent(in) :: sval
      logical, intent(in) :: is_known
      logical, intent(in), optional :: ignore, force, verbose
      type(string_t), intent(in), optional :: model_name
    end subroutine var_list_set_string
    module subroutine var_list_import (var_list, src_list)
      class(var_list_t), intent(inout) :: var_list
      type(var_list_t), intent(in) :: src_list
    end subroutine var_list_import
    recursive module subroutine var_list_undefine (var_list, follow_link)
      class(var_list_t), intent(inout) :: var_list
      logical, intent(in), optional :: follow_link
    end subroutine var_list_undefine
    recursive module subroutine var_list_init_snapshot &
         (var_list, vars_in, follow_link)
      class(var_list_t), intent(out) :: var_list
      type(var_list_t), intent(in) :: vars_in
      logical, intent(in), optional :: follow_link
    end subroutine var_list_init_snapshot
    module subroutine var_list_check_user_var (var_list, name, type, new)
      class(var_list_t), intent(in), target :: var_list
      type(string_t), intent(in) :: name
      integer, intent(inout) :: type
      logical, intent(in) :: new
    end subroutine var_list_check_user_var
    module subroutine var_list_init_defaults (var_list, seed, paths)
      class(var_list_t), intent(out) :: var_list
      integer, intent(in) :: seed
      type(paths_t), intent(in), optional :: paths
    end subroutine var_list_init_defaults
    module subroutine var_list_set_beams_defaults (var_list, paths)
      type(paths_t), intent(in), optional :: paths
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_beams_defaults
    module subroutine var_list_set_core_defaults (var_list, seed)
      class(var_list_t), intent(inout) :: var_list
      integer, intent(in) :: seed
    end subroutine var_list_set_core_defaults
    module subroutine var_list_set_integration_defaults (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_integration_defaults
    module subroutine var_list_set_phase_space_defaults (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_phase_space_defaults
    module subroutine var_list_set_gamelan_defaults (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_gamelan_defaults
    module subroutine var_list_set_clustering_defaults (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_clustering_defaults
    module subroutine var_list_set_isolation_recomb_defaults (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_isolation_recomb_defaults
    module subroutine var_list_set_eio_defaults (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_eio_defaults
    module subroutine var_list_set_shower_defaults (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_shower_defaults
    module subroutine var_list_set_hadronization_defaults (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_hadronization_defaults
    module subroutine var_list_set_tauola_defaults (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_tauola_defaults
    module subroutine var_list_set_mlm_matching_defaults (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_mlm_matching_defaults
    module subroutine var_list_set_powheg_matching_defaults (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_powheg_matching_defaults
    module subroutine var_list_set_openmp_defaults (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_openmp_defaults
    module subroutine var_list_set_mpi_defaults (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_mpi_defaults
    module subroutine var_list_set_nlo_defaults (var_list)
      class(var_list_t), intent(inout) :: var_list
    end subroutine var_list_set_nlo_defaults
  end interface

end module variables
