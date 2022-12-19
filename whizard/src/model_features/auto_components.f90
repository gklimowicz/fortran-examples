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

module auto_components

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use model_data
  use pdg_arrays

  implicit none
  private

  public :: split_constraints_t
  public :: constrain_n_tot
  public :: constrain_n_loop
  public :: constrain_splittings
  public :: constrain_insert
  public :: constrain_require
  public :: constrain_radiation
  public :: constrain_mass_sum
  public :: constrain_in_state
  public :: constrain_photon_induced_processes
  public :: constrain_couplings
  public :: ps_table_t
  public :: ds_table_t
  public :: fs_table_t
  public :: if_table_t

  type, abstract :: split_constraint_t
  contains
    procedure :: check_before_split  => split_constraint_check_before_split
    procedure :: check_before_insert => split_constraint_check_before_insert
    procedure :: check_before_record => split_constraint_check_before_record
  end type split_constraint_t

  type :: split_constraint_wrap_t
     class(split_constraint_t), allocatable :: c
  end type split_constraint_wrap_t

  type :: split_constraints_t
     class(split_constraint_wrap_t), dimension(:), allocatable :: cc
   contains
     procedure :: init => split_constraints_init
     procedure :: set => split_constraints_set
     procedure :: check_before_split  => split_constraints_check_before_split
     procedure :: check_before_insert => split_constraints_check_before_insert
     procedure :: check_before_record => split_constraints_check_before_record
  end type split_constraints_t

  type, extends (split_constraint_t) :: constraint_n_tot
     private
     integer :: n_max = 0
   contains
     procedure :: check_before_split => constraint_n_tot_check_before_split
     procedure :: check_before_record => constraint_n_tot_check_before_record
  end type constraint_n_tot

  type, extends (split_constraint_t) :: constraint_n_loop
     private
     integer :: n_loop_max = 0
   contains
     procedure :: check_before_record => constraint_n_loop_check_before_record
  end type constraint_n_loop

  type, extends (split_constraint_t) :: constraint_splittings
     private
     type(pdg_list_t) :: pl_match, pl_excluded_gauge_splittings
   contains
     procedure :: check_before_insert => constraint_splittings_check_before_insert
  end type constraint_splittings

  type, extends (split_constraint_t) :: constraint_insert
     private
     type(pdg_list_t) :: pl_match
   contains
     procedure :: check_before_insert => constraint_insert_check_before_insert
  end type constraint_insert

  type, extends (split_constraint_t) :: constraint_require
     private
     type(pdg_list_t) :: pl
   contains
     procedure :: check_before_record => constraint_require_check_before_record
  end type constraint_require

  type, extends (split_constraint_t) :: constraint_radiation
     private
   contains
     procedure :: check_before_insert => &
          constraint_radiation_check_before_insert
  end type constraint_radiation

  type, extends (split_constraint_t) :: constraint_mass_sum
     private
     real(default) :: mass_limit = 0
     logical :: strictly_less = .false.
     real(default) :: margin = 0
   contains
     procedure :: check_before_record => constraint_mass_sum_check_before_record
  end type constraint_mass_sum

  type, extends (split_constraint_t) :: constraint_in_state
     private
     type(pdg_list_t) :: pl
   contains
     procedure :: check_before_record => constraint_in_state_check_before_record
  end type constraint_in_state

  type, extends (split_constraint_t) :: constraint_photon_induced_processes
     private
     integer :: n_in
   contains
     procedure :: check_before_record => &
          constraint_photon_induced_processes_check_before_record
  end type constraint_photon_induced_processes

  type, extends (split_constraint_t) :: constraint_coupling_t
    private
    logical :: qed = .false.
    logical :: qcd = .true.
    logical :: ew = .false.
    integer :: n_nlo_correction_types
  contains
    procedure :: check_before_insert => constraint_coupling_check_before_insert
  end type constraint_coupling_t

  type, extends (pdg_list_t) :: ps_entry_t
     integer :: n_loop = 0
     integer :: n_rad = 0
     type(ps_entry_t), pointer :: previous => null ()
     type(ps_entry_t), pointer :: next => null ()
  end type ps_entry_t

  type, abstract :: ps_table_t
     private
     class(model_data_t), pointer :: model => null ()
     logical :: loops = .false.
     type(ps_entry_t), pointer :: first => null ()
     type(ps_entry_t), pointer :: last => null ()
     integer :: proc_type
   contains
     procedure :: final => ps_table_final
     procedure :: base_write => ps_table_base_write
     procedure (ps_table_write), deferred :: write
     procedure :: get_particle_string => ps_table_get_particle_string
     generic :: init => ps_table_init
     procedure, private :: ps_table_init
     procedure :: enable_loops => ps_table_enable_loops
     procedure :: split => ps_table_split
     procedure :: insert => ps_table_insert
     procedure :: record_sorted => ps_table_record_sorted
     procedure :: record => ps_table_record
     procedure :: get_length => ps_table_get_length
     procedure :: get_emitters => ps_table_get_emitters
     procedure :: get_pdg_out => ps_table_get_pdg_out
  end type ps_table_t

  type, extends (ps_table_t) :: ds_table_t
     private
     integer :: pdg_in = 0
   contains
     procedure :: write => ds_table_write
     procedure :: make => ds_table_make
  end type ds_table_t

  type, extends (ps_table_t) :: fs_table_t
   contains
     procedure :: write => fs_table_write
     procedure :: radiate => fs_table_radiate
  end type fs_table_t

  type, extends (fs_table_t) :: if_table_t
   contains
     procedure :: write => if_table_write
     generic :: init => if_table_init
     procedure, private :: if_table_init
     procedure :: insert => if_table_insert
     procedure :: record_sorted => if_table_record_sorted
  end type if_table_t


  interface
     subroutine ps_table_write (object, unit)
       import
       class(ps_table_t), intent(in) :: object
       integer, intent(in), optional :: unit
     end subroutine ps_table_write
  end interface

  interface
    module subroutine split_constraint_check_before_split (c, table, pl, k, passed)
      class(split_constraint_t), intent(in) :: c
      class(ps_table_t), intent(in) :: table
      type(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: k
      logical, intent(out) :: passed
    end subroutine split_constraint_check_before_split
    module subroutine split_constraint_check_before_insert (c, table, pa, pl, passed)
      class(split_constraint_t), intent(in) :: c
      class(ps_table_t), intent(in) :: table
      type(pdg_array_t), intent(in) :: pa
      type(pdg_list_t), intent(inout) :: pl
      logical, intent(out) :: passed
    end subroutine split_constraint_check_before_insert
    module subroutine split_constraint_check_before_record (c, table, pl, n_loop, passed)
      class(split_constraint_t), intent(in) :: c
      class(ps_table_t), intent(in) :: table
      type(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: n_loop
      logical, intent(out) :: passed
    end subroutine split_constraint_check_before_record
    module subroutine split_constraints_set (constraints, i, c)
      class(split_constraints_t), intent(inout) :: constraints
      integer, intent(in) :: i
      class(split_constraint_t), intent(in) :: c
    end subroutine split_constraints_set
    module subroutine split_constraints_check_before_split &
         (constraints, table, pl, k, passed)
      class(split_constraints_t), intent(in) :: constraints
      class(ps_table_t), intent(in) :: table
      type(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: k
      logical, intent(out) :: passed
    end subroutine split_constraints_check_before_split
    module subroutine split_constraints_check_before_insert &
         (constraints, table, pa, pl, passed)
      class(split_constraints_t), intent(in) :: constraints
      class(ps_table_t), intent(in) :: table
      type(pdg_array_t), intent(in) :: pa
      type(pdg_list_t), intent(inout) :: pl
      logical, intent(out) :: passed
    end subroutine split_constraints_check_before_insert
    module subroutine split_constraints_check_before_record &
         (constraints, table, pl, n_loop, passed)
      class(split_constraints_t), intent(in) :: constraints
      class(ps_table_t), intent(in) :: table
      type(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: n_loop
      logical, intent(out) :: passed
    end subroutine split_constraints_check_before_record
    module function constrain_n_tot (n_max) result (c)
      integer, intent(in) :: n_max
      type(constraint_n_tot) :: c
    end function constrain_n_tot
    module subroutine constraint_n_tot_check_before_split (c, table, pl, k, passed)
      class(constraint_n_tot), intent(in) :: c
      class(ps_table_t), intent(in) :: table
      type(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: k
      logical, intent(out) :: passed
    end subroutine constraint_n_tot_check_before_split
    module subroutine constraint_n_tot_check_before_record (c, table, pl, n_loop, passed)
      class(constraint_n_tot), intent(in) :: c
      class(ps_table_t), intent(in) :: table
      type(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: n_loop
      logical, intent(out) :: passed
    end subroutine constraint_n_tot_check_before_record
    module function constrain_n_loop (n_loop_max) result (c)
      integer, intent(in) :: n_loop_max
      type(constraint_n_loop) :: c
    end function constrain_n_loop
    module subroutine constraint_n_loop_check_before_record &
         (c, table, pl, n_loop, passed)
      class(constraint_n_loop), intent(in) :: c
      class(ps_table_t), intent(in) :: table
      type(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: n_loop
      logical, intent(out) :: passed
    end subroutine constraint_n_loop_check_before_record
    module function constrain_splittings &
         (pl_match, pl_excluded_gauge_splittings) result (c)
      type(pdg_list_t), intent(in) :: pl_match
      type(pdg_list_t), intent(in) :: pl_excluded_gauge_splittings
      type(constraint_splittings) :: c
    end function constrain_splittings
    module subroutine constraint_splittings_check_before_insert &
         (c, table, pa, pl, passed)
      class(constraint_splittings), intent(in) :: c
      class(ps_table_t), intent(in) :: table
      type(pdg_array_t), intent(in) :: pa
      type(pdg_list_t), intent(inout) :: pl
      logical, intent(out) :: passed
    end subroutine constraint_splittings_check_before_insert
    module function constrain_insert (pl_match) result (c)
      type(pdg_list_t), intent(in) :: pl_match
      type(constraint_insert) :: c
    end function constrain_insert
    module subroutine constraint_insert_check_before_insert (c, table, pa, pl, passed)
      class(constraint_insert), intent(in) :: c
      class(ps_table_t), intent(in) :: table
      type(pdg_array_t), intent(in) :: pa
      type(pdg_list_t), intent(inout) :: pl
      logical, intent(out) :: passed
    end subroutine constraint_insert_check_before_insert
    module function constrain_require (pl) result (c)
      type(pdg_list_t), intent(in) :: pl
      type(constraint_require) :: c
    end function constrain_require
    module subroutine constraint_require_check_before_record &
         (c, table, pl, n_loop, passed)
      class(constraint_require), intent(in) :: c
      class(ps_table_t), intent(in) :: table
      type(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: n_loop
      logical, intent(out) :: passed
    end subroutine constraint_require_check_before_record
    module function constrain_radiation () result (c)
      type(constraint_radiation) :: c
    end function constrain_radiation
    module subroutine constraint_radiation_check_before_insert &
         (c, table, pa, pl, passed)
      class(constraint_radiation), intent(in) :: c
      class(ps_table_t), intent(in) :: table
      type(pdg_array_t), intent(in) :: pa
      type(pdg_list_t), intent(inout) :: pl
      logical, intent(out) :: passed
    end subroutine constraint_radiation_check_before_insert
    module function constrain_mass_sum (mass_limit, margin) result (c)
      real(default), intent(in) :: mass_limit
      real(default), intent(in), optional :: margin
      type(constraint_mass_sum) :: c
    end function constrain_mass_sum
    module subroutine constraint_mass_sum_check_before_record &
         (c, table, pl, n_loop, passed)
      class(constraint_mass_sum), intent(in) :: c
      class(ps_table_t), intent(in) :: table
      type(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: n_loop
      logical, intent(out) :: passed
    end subroutine constraint_mass_sum_check_before_record
    module function constrain_in_state (pl) result (c)
      type(pdg_list_t), intent(in) :: pl
      type(constraint_in_state) :: c
    end function constrain_in_state
    module subroutine constraint_in_state_check_before_record &
         (c, table, pl, n_loop, passed)
      class(constraint_in_state), intent(in) :: c
      class(ps_table_t), intent(in) :: table
      type(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: n_loop
      logical, intent(out) :: passed
    end subroutine constraint_in_state_check_before_record
    module function constrain_photon_induced_processes (n_in) result (c)
      integer, intent(in) :: n_in
      type(constraint_photon_induced_processes) :: c
    end function constrain_photon_induced_processes
    module subroutine constraint_photon_induced_processes_check_before_record &
         (c, table, pl, n_loop, passed)
      class(constraint_photon_induced_processes), intent(in) :: c
      class(ps_table_t), intent(in) :: table
      type(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: n_loop
      logical, intent(out) :: passed
    end subroutine constraint_photon_induced_processes_check_before_record
    module function constrain_couplings (qcd, qed, n_nlo_correction_types) result (c)
      type(constraint_coupling_t) :: c
      logical, intent(in) :: qcd, qed
      integer, intent(in) :: n_nlo_correction_types
    end function constrain_couplings
    module subroutine constraint_coupling_check_before_insert (c, table, pa, pl, passed)
      class(constraint_coupling_t), intent(in) :: c
      class(ps_table_t), intent(in) :: table
      type(pdg_array_t), intent(in) :: pa
      type(pdg_list_t), intent(inout) :: pl
      logical, intent(out) :: passed
    end subroutine constraint_coupling_check_before_insert
    module subroutine ps_table_final (object)
      class(ps_table_t), intent(inout) :: object
    end subroutine ps_table_final
    module subroutine ps_table_base_write (object, unit, n_in)
      class(ps_table_t), intent(in) :: object
      integer, intent(in), optional :: unit
      integer, intent(in), optional :: n_in
    end subroutine ps_table_base_write
    module subroutine ds_table_write (object, unit)
      class(ds_table_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine ds_table_write
    module subroutine fs_table_write (object, unit)
      class(fs_table_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine fs_table_write
    module subroutine if_table_write (object, unit)
      class(if_table_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine if_table_write
    module subroutine ps_table_get_particle_string (object, index, prt_in, prt_out)
      class(ps_table_t), intent(in) :: object
      integer, intent(in) :: index
      type(string_t), intent(out), dimension(:), allocatable :: prt_in, prt_out
    end subroutine ps_table_get_particle_string
    module subroutine ps_table_init &
         (table, model, pl, constraints, n_in, do_not_check_regular)
      class(ps_table_t), intent(out) :: table
      class(model_data_t), intent(in), target :: model
      type(pdg_list_t), dimension(:), intent(in) :: pl
      type(split_constraints_t), intent(in) :: constraints
      integer, intent(in), optional :: n_in
      logical, intent(in), optional :: do_not_check_regular
    end subroutine ps_table_init
    module subroutine if_table_init (table, model, pl_in, pl_out, constraints)
      class(if_table_t), intent(out) :: table
      class(model_data_t), intent(in), target :: model
      type(pdg_list_t), dimension(:), intent(in) :: pl_in, pl_out
      type(split_constraints_t), intent(in) :: constraints
    end subroutine if_table_init
    module subroutine ps_table_enable_loops (table)
      class(ps_table_t), intent(inout) :: table
    end subroutine ps_table_enable_loops
    module subroutine ds_table_make (table, model, pdg_in, constraints)
      class(ds_table_t), intent(out) :: table
      class(model_data_t), intent(in), target :: model
      integer, intent(in) :: pdg_in
      type(split_constraints_t), intent(in) :: constraints
    end subroutine ds_table_make
    module subroutine fs_table_radiate (table, constraints, do_not_check_regular)
      class(fs_table_t), intent(inout) :: table
      type(split_constraints_t) :: constraints
      logical, intent(in), optional :: do_not_check_regular
    end subroutine fs_table_radiate
    recursive module subroutine ps_table_split (table, pl, n_rad, constraints, &
          record, do_not_check_regular)
      class(ps_table_t), intent(inout) :: table
      class(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: n_rad
      type(split_constraints_t), intent(in) :: constraints
      logical, intent(in), optional :: record, do_not_check_regular
    end subroutine ps_table_split
    recursive module subroutine ps_table_insert &
         (table, pl, n_rad, i, pdg, constraints, n_in, do_not_check_regular)
      class(ps_table_t), intent(inout) :: table
      class(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: n_rad, i
      integer, dimension(:), intent(in) :: pdg
      type(split_constraints_t), intent(in) :: constraints
      integer, intent(in), optional :: n_in
      logical, intent(in), optional :: do_not_check_regular
    end subroutine ps_table_insert
    recursive module subroutine if_table_insert  &
         (table, pl, n_rad, i, pdg, constraints, n_in, do_not_check_regular)
      class(if_table_t), intent(inout) :: table
      class(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: n_rad, i
      integer, dimension(:), intent(in) :: pdg
      type(split_constraints_t), intent(in) :: constraints
      integer, intent(in), optional :: n_in
      logical, intent(in), optional :: do_not_check_regular
    end subroutine if_table_insert
    module subroutine ps_table_record_sorted &
         (table, pl, n_loop, n_rad, constraints, do_not_check_regular, passed)
      class(ps_table_t), intent(inout) :: table
      type(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: n_loop, n_rad
      type(split_constraints_t), intent(in) :: constraints
      logical, intent(in), optional :: do_not_check_regular
      logical, intent(out) :: passed
    end subroutine ps_table_record_sorted
    module subroutine if_table_record_sorted &
         (table, pl, n_loop, n_rad, constraints, do_not_check_regular, passed)
      class(if_table_t), intent(inout) :: table
      type(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: n_loop, n_rad
      type(split_constraints_t), intent(in) :: constraints
      logical, intent(in), optional :: do_not_check_regular
      logical, intent(out) :: passed
    end subroutine if_table_record_sorted
    module subroutine ps_table_record (table, pl, n_loop, n_rad, constraints, &
         do_not_check_regular, passed)
      class(ps_table_t), intent(inout) :: table
      type(pdg_list_t), intent(in) :: pl
      integer, intent(in) :: n_loop, n_rad
      type(split_constraints_t), intent(in) :: constraints
      logical, intent(in), optional :: do_not_check_regular
      logical, intent(out) :: passed
    end subroutine ps_table_record
    module function ps_table_get_length (ps_table) result (n)
      class(ps_table_t), intent(in) :: ps_table
      integer :: n
    end function ps_table_get_length
    module subroutine ps_table_get_emitters (table, constraints, emitters)
      class(ps_table_t), intent(in) :: table
      type(split_constraints_t), intent(in) :: constraints
      integer, dimension(:), allocatable, intent(out) :: emitters
    end subroutine ps_table_get_emitters
    module subroutine ps_table_get_pdg_out (ps_table, i, pa_out, n_loop, n_rad)
      class(ps_table_t), intent(in) :: ps_table
      integer, intent(in) :: i
      type(pdg_array_t), dimension(:), allocatable, intent(out) :: pa_out
      integer, intent(out), optional :: n_loop, n_rad
    end subroutine ps_table_get_pdg_out
  end interface

contains

  subroutine split_constraints_init (constraints, n)
    class(split_constraints_t), intent(out) :: constraints
    integer, intent(in) :: n
    allocate (constraints%cc (n))
  end subroutine split_constraints_init


end module auto_components
