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

module pythia8

  use kinds, only: default
  use debug_master, only: debug_on

  use, intrinsic :: iso_c_binding
  use constants, only: tiny_10, tiny_07
  use diagnostics
  use iso_varying_string, string_t => varying_string
  use lorentz, only: assignment(=), operator(/=), &
       vector3_moving, vector4_t, vector4_moving, vector4_null, &
       vector4_write
  use numeric_utils, only: vanishes, nearly_equal
  use model_data, only: model_data_t, find_model
  use pdg_arrays, only: is_elementary, is_colored, is_gluon
  use colors, only: color_t
  use flavors, only: flavor_t
  use helicities, only: helicity_t
  use particles, only: particle_set_t, particle_t, &
       PRT_DEFINITE_HELICITY, PRT_GENERIC_POLARIZATION, PRT_UNPOLARIZED
  use event_base, only: generic_event_t
  use subevents, only: PRT_BEAM, PRT_INCOMING, PRT_OUTGOING, &
       PRT_UNDEFINED, PRT_VIRTUAL, PRT_RESONANT, PRT_BEAM_REMNANT
  use rng_base, only: rng_t
  use whizard_lha

  implicit none
  private

  public :: pythia8_t

  integer, parameter :: C_MAX_STR_LEN = 100



  type :: whizard_rndm_t
     class(rng_t), pointer :: rng
  end type whizard_rndm_t

  type :: pythia8_t
     private
     type(c_ptr) :: cptr
     type(whizard_rndm_t) :: rndm
   contains
     procedure :: init => whizard_pythia8_init
     procedure :: final => whizard_pythia8_final
     procedure :: set_lhaup_ptr => whizard_pythia8_set_lhaup_ptr
     procedure :: import_rng => whizard_pythia8_import_rng
     procedure :: set_rng_seed => pythia8_set_rng_seed
     procedure :: read_file => whizard_pythia8_read_file
     procedure :: read_string => whizard_pythia8_read_string
     procedure :: parse_and_set_config => whizard_pythia8_parse_and_set_config
     procedure :: init_pythia => whizard_pythia8_init_pythia
     procedure :: next => whizard_pythia8_next
     procedure :: list_lha_event => whizard_pythia8_list_lha_event
     procedure :: list_event => whizard_pythia8_list_event
     procedure :: get_event_size => whizard_pythia8_get_event_size
     procedure, private :: get_single_event => whizard_pythia8_get_single_event
     procedure, private :: get_particle_status => &
          whizard_pythia8_get_particle_status
     procedure, private :: get_particle_id => whizard_pythia8_get_particle_id
     procedure, private :: get_particle_momentum => &
          whizard_pythia8_get_particle_momentum
     procedure :: get_final_colored_ME_momenta => &
          whizard_pythia8_get_final_colored_ME_momenta
     procedure, private :: get_shower_mask => pythia8_get_shower_mask
     procedure, private :: get_hadron_mask => pythia8_get_hadron_mask
     procedure :: get_shower_particles => whizard_pythia8_get_shower_particles
     procedure :: get_hadron_particles => whizard_pythia8_get_hadron_particles
     procedure, private :: get_particles => whizard_pythia8_get_particles
     procedure, private :: import_pythia_particles => pythia8_import_pythia_particles
     procedure, private :: is_beam_particle => pythia8_is_beam_particle
     procedure, private :: get_parent_child_relation => pythia8_get_parent_child_relation
  end type pythia8_t


  interface
     function new_pythia8 (print_banner) bind(C) result (cptr)
       import
       type(c_ptr) :: cptr
       logical(c_bool), value, intent(in) :: print_banner
     end function new_pythia8
  end interface

  interface
     subroutine pythia8_delete (pythia) bind(C)
       import
       type(c_ptr), value :: pythia
     end subroutine pythia8_delete
  end interface

  interface
     function pythia8_set_lhaup_ptr (cptr, whizard_lha) bind(C) result (flag)
       import
       type(c_ptr), intent(in), value :: cptr
       type(c_ptr), intent(in), value :: whizard_lha
       logical(c_bool) :: flag
     end function pythia8_set_lhaup_ptr
  end interface

  interface
     function pythia8_set_rndm_engine_ptr (pythia, whizard_rndm) bind(C) result (flag)
       import
       type(c_ptr), intent(in), value :: pythia
       type(c_ptr), intent(in), value :: whizard_rndm
       logical(c_bool) :: flag
     end function pythia8_set_rndm_engine_ptr
  end interface
  interface
     function pythia8_read_file (cptr, filename, subrun) bind(C) result (flag)
       import
       type(c_ptr), intent(in), value :: cptr
       character(kind=c_char), dimension(*), intent(in) :: filename
       integer(c_int), intent(in), value :: subrun
       logical(c_bool) :: flag
     end function pythia8_read_file
  end interface

  interface
     function pythia8_read_string (cptr, str) bind(C) result (flag)
       import
       type(c_ptr), intent(in), value :: cptr
       character(kind=c_char), dimension(*), intent(in) :: str
       logical(c_bool) :: flag
     end function pythia8_read_string
  end interface

  interface
     function pythia8_init (cptr) bind(C) result (flag)
       import
       type(c_ptr), intent(in), value :: cptr
       logical(c_bool) :: flag
     end function pythia8_init
  end interface

  interface
     function pythia8_next (cptr) bind(C) result (flag)
       import
       type(c_ptr), intent(in), value :: cptr
       logical(c_bool) :: flag
     end function pythia8_next
  end interface

  interface
     subroutine pythia8_list_lha_event (cptr) bind(C)
       import
       type(c_ptr), intent(in), value :: cptr
     end subroutine pythia8_list_lha_event
  end interface

  interface
     subroutine pythia8_list_event (cptr) bind(C)
       import
       type(c_ptr), intent(in), value :: cptr
     end subroutine pythia8_list_event
  end interface

  interface
     function pythia8_get_event_size (cptr) bind(C) result(n)
       import
       type(c_ptr), intent(in), value :: cptr
       integer(c_int) :: n
     end function pythia8_get_event_size
  end interface

  interface
     function pythia8_get_single_event (cptr, index) bind(C) result (particle)
       import
       type(c_ptr), intent(in), value :: cptr
       integer(c_int), intent(in), value :: index
       type(lha_particle_t) :: particle
     end function pythia8_get_single_event
  end interface

  interface
     function pythia8_get_particle_status (cptr, index) bind(C) result (status)
       import
       type(c_ptr), intent(in), value :: cptr
       integer(c_int), intent(in), value :: index
       integer(c_int) :: status
     end function pythia8_get_particle_status
  end interface

  interface
     function pythia8_get_particle_id (cptr, index) bind(C) result (status)
       import
       type(c_ptr), intent(in), value :: cptr
       integer(c_int), intent(in), value :: index
       integer(c_int) :: status
     end function pythia8_get_particle_id
  end interface

  interface
     subroutine pythia8_get_particle_momentum (cptr, index, momentum) bind(C)
       import
       type(c_ptr), intent(in), value :: cptr
       integer(c_int), intent(in), value :: index
       real(c_double), dimension(*), intent(out) :: momentum
     end subroutine pythia8_get_particle_momentum
  end interface
  interface
     function pythia8_get_n_mothers (cptr, i_prt) bind(C) result (n_mothers)
       import
       type(c_ptr), intent(in), value :: cptr
       integer(c_int), value, intent(in) :: i_prt
       integer(c_int) :: n_mothers
     end function pythia8_get_n_mothers

     subroutine pythia8_get_mother_array (cptr, i_prt, n_mothers, mother) bind (C)
       import
       type(c_ptr), intent(in), value :: cptr
       integer(c_int), value :: i_prt
       integer(c_int), value :: n_mothers
       integer(c_int), dimension(*), intent(out) :: mother
     end subroutine pythia8_get_mother_array

     function pythia8_get_n_daughters (cptr, i_prt) bind(C) result (n_daughters)
       import
       type(c_ptr), intent(in), value :: cptr
       integer(c_int), value, intent(in) :: i_prt
       integer(c_int) :: n_daughters
     end function pythia8_get_n_daughters

     subroutine pythia8_get_daughter_array (cptr, i_prt, n_daughters, daughter) bind (C)
       import
       type(c_ptr), intent(in), value :: cptr
       integer(c_int), value :: i_prt
       integer(c_int), value :: n_daughters
       integer(c_int), dimension(*), intent(out) :: daughter
     end subroutine pythia8_get_daughter_array

     function pythia8_get_status_hepmc (cptr, i_prt) bind(C) result (status)
       import
       type(c_ptr), intent(in), value :: cptr
       integer(c_int), value, intent(in) :: i_prt
       integer(c_int) :: status
     end function pythia8_get_status_hepmc

     subroutine pythia8_get_decay_vertex (cptr, i_prt, time, vertex) bind(C)
       import
       type(c_ptr), intent(in), value :: cptr
       integer(c_int), value, intent(in) :: i_prt
       real(c_double), intent(out) :: time
       real(c_double), dimension(3), intent(out) :: vertex
     end subroutine pythia8_get_decay_vertex

     subroutine pythia8_get_production_vertex (cptr, i_prt, time, vertex) bind(C)
       import
       type(c_ptr), intent(in), value :: cptr
       integer(c_int), value, intent(in) :: i_prt
       real(c_double), intent(out) :: time
       real(c_double), dimension(3), intent(out) :: vertex
     end subroutine pythia8_get_production_vertex
  end interface

contains

  subroutine whizard_pythia8_init (pythia, verbose)
    class(pythia8_t), intent(out) :: pythia
    logical, intent(in), optional :: verbose
    logical(c_bool) :: verbose_opt
    verbose_opt = .false.
    if (present (verbose)) verbose_opt = verbose
    pythia%cptr = new_pythia8 (verbose_opt)
    if (.not. verbose_opt) &
         call pythia%read_string (var_str ("Print:quiet = on"))
  end subroutine whizard_pythia8_init

  subroutine whizard_pythia8_final (pythia)
    class(pythia8_t), intent(inout) :: pythia
    call pythia8_delete (pythia%cptr)
  end subroutine whizard_pythia8_final

  subroutine whizard_pythia8_set_lhaup_ptr (pythia, whizard_lha)
    class(pythia8_t), intent(inout) :: pythia
    type(whizard_lha_t), intent(in) :: whizard_lha
    logical(c_bool) :: flag
    flag = pythia8_set_lhaup_ptr (pythia%cptr, whizard_lha%get_ptr ())
  end subroutine whizard_pythia8_set_lhaup_ptr

  subroutine whizard_pythia8_import_rng (pythia, rng)
    class(pythia8_t), intent(inout), target :: pythia
    class(rng_t), allocatable, intent(in), target :: rng
    logical :: flag
    pythia%rndm%rng => rng
    flag = pythia8_set_rndm_engine_ptr (pythia%cptr, c_loc (pythia%rndm))
    if (.not. flag) then
       call msg_bug ("[whizard_pythia8_import_rng] Cannot export RNG to Pythia8.")
    end if
  end subroutine whizard_pythia8_import_rng

  function whizard_rndm_generate (whizard_rndm) bind(C) result (c_x)
    type(c_ptr), intent(in), value :: whizard_rndm
    real(c_double) :: c_x
    real(default) :: x
    type(whizard_rndm_t), pointer :: f_whizard_rndm
    call c_f_pointer (whizard_rndm, f_whizard_rndm)
    if (.not. associated (f_whizard_rndm)) then
       call msg_bug ("[whizard_rndm_generate] Cannot import pointer to RNG object from Pythia8.")
    end if
    call f_whizard_rndm%rng%generate (x)
    c_x = real (x, c_double)
  end function whizard_rndm_generate

  subroutine pythia8_set_rng_seed (pythia, r)
    class(pythia8_t), intent(inout) :: pythia
    real(default), intent(in) :: r
    real(default), parameter :: MAX_SEED = 900000000._default
    character(len=10) :: buffer; type(string_t) :: string
    write (buffer, "(I10)") floor (r * MAX_SEED)
    string = var_str ("Random:seed = " // buffer)
    call pythia%read_string (string)
  end subroutine pythia8_set_rng_seed

  subroutine whizard_pythia8_read_file (pythia, filename, subrun)
    class(pythia8_t), intent(inout) :: pythia
    type(string_t), intent(in) :: filename
    integer, intent(in), optional :: subrun
    character(len(filename) + 1, kind=c_char) :: c_filename
    integer(c_int) :: c_subrun
    logical(c_bool) :: flag
    c_filename = filename // c_null_char
    if (present (subrun)) then
       c_subrun = subrun
    else
       c_subrun = -1
    end if
    flag = pythia8_read_file (pythia%cptr, c_filename, c_subrun)
  end subroutine whizard_pythia8_read_file

  subroutine whizard_pythia8_read_string (pythia, str)
    class(pythia8_t), intent(in) :: pythia
    type(string_t), intent(in) :: str
    character(len(str) + 1, kind=c_char) :: c_str
    logical(c_bool) :: flag
    c_str = char (str) // c_null_char
    flag = pythia8_read_string (pythia%cptr, c_str)
  end subroutine whizard_pythia8_read_string

  subroutine whizard_pythia8_parse_and_set_config (pythia, config)
    class(pythia8_t), intent(in) :: pythia
    type(string_t), intent(in) :: config
    type(string_t) :: stream, line, token
    if (debug_on) call msg_debug (D_SHOWER, "whizard_pythia8_parse_and_set_config")
    if (len (config) == 0) return
    stream = config
    do while (len (stream) > 0)
       call split (stream, line, new_line("A"))
       if (debug_active (D_SHOWER)) &
            print *, "LINE: ", char(line), " | ", char(stream)
       if (index (line, ";") == 0) then
          call pythia%read_string (trim(line))
       else
          token = line
            do while (len (line) > 0)
               call split (line, token, ";")
               if (debug_active (D_SHOWER)) &
                    print *, "-> ", char(token), " | ", char(line)
               call pythia%read_string (trim(token))
            end do
         end if
      end do

  end subroutine whizard_pythia8_parse_and_set_config

  subroutine whizard_pythia8_init_pythia (pythia)
    class(pythia8_t), intent(in) :: pythia
    logical(c_bool) :: flag
    flag = pythia8_init (pythia%cptr)
    if (.not. flag)  then
       call msg_fatal ("[whizard_pythia8_init_pythia] Pythia8 initialisation failed.")
    end if
  end subroutine whizard_pythia8_init_pythia

  subroutine whizard_pythia8_next (pythia, flag)
    class(pythia8_t), intent(inout) :: pythia
    logical, intent(out), optional :: flag
    logical(c_bool) :: c_flag
    c_flag = pythia8_next (pythia%cptr)
    if (present (flag)) flag = c_flag
  end subroutine whizard_pythia8_next

  subroutine whizard_pythia8_list_lha_event (pythia)
    class(pythia8_t), intent(in) :: pythia
    call pythia8_list_lha_event (pythia%cptr)
  end subroutine whizard_pythia8_list_lha_event

  subroutine whizard_pythia8_list_event (pythia)
    class(pythia8_t), intent(in) :: pythia
    call pythia8_list_event (pythia%cptr)
  end subroutine whizard_pythia8_list_event

  function whizard_pythia8_get_event_size (pythia) result(n)
    class(pythia8_t), intent(in) :: pythia
    integer :: n
    integer(c_int) :: c_n
    c_n = pythia8_get_event_size (pythia%cptr)
    n = c_n - 1
  end function whizard_pythia8_get_event_size

  function whizard_pythia8_get_single_event (pythia, index) result (particle)
    class(pythia8_t), intent(in) :: pythia
    integer, intent(in) :: index
    type(lha_particle_t) :: particle
    integer(c_int) :: c_index
    c_index = index
    particle = pythia8_get_single_event (pythia%cptr, c_index)
  end function whizard_pythia8_get_single_event

  function whizard_pythia8_get_particle_status (pythia, index) result (status)
    class(pythia8_t), intent(in) :: pythia
    integer, intent(in) :: index
    integer :: status
    status = pythia8_get_particle_status (pythia%cptr, int(index, c_int))
  end function whizard_pythia8_get_particle_status

  function whizard_pythia8_get_particle_id (pythia, index) result (id)
    class(pythia8_t), intent(in) :: pythia
    integer, intent(in) :: index
    integer :: id
    id = pythia8_get_particle_id (pythia%cptr, int(index, c_int))
  end function whizard_pythia8_get_particle_id

  function whizard_pythia8_get_particle_momentum &
       (pythia, index) result (momentum)
    class(pythia8_t), intent(in) :: pythia
    integer, intent(in) :: index
    real(default), dimension(4) :: momentum
    real(c_double), dimension(4) :: c_momentum
    call pythia8_get_particle_momentum (pythia%cptr, index, c_momentum)
    momentum = real (c_momentum, kind=default)
  end function whizard_pythia8_get_particle_momentum

  subroutine whizard_pythia8_get_final_colored_ME_momenta &
         (pythia, momenta)
    class(pythia8_t), intent(in) :: pythia
    type(vector4_t), dimension(:), allocatable, intent(out) :: momenta
    logical, dimension(:), allocatable :: mask
    integer, parameter :: PYTHIA8_HARD_PROCESS_OUTGOING = 23
    integer :: i, j, n_particles, id, status
    if (debug_on) call msg_debug (D_TRANSFORMS, "whizard_pythia8_get_final_colored_ME_momenta")
    n_particles = pythia%get_event_size ()
    allocate (mask(n_particles), source=.false.)
    do i = 1, n_particles
       status = pythia%get_particle_status (i)
       id = pythia%get_particle_id (i)
       if (abs (status) == PYTHIA8_HARD_PROCESS_OUTGOING &
            .and. (id == 21 .or. abs(id) <= 6)) mask(i) = .true.
       ! Particle record is ordered. First beam, beam remnants, then hard process, second hard process and so on...
       ! if (abs(status) > 30) exit
    end do
    if (all (.not. mask)) return
    allocate (momenta(count(mask)))
    j = 1
    do i = 1, n_particles
       if (.not. mask(i)) cycle
       momenta(j) = pythia%get_particle_momentum(i)
       j = j + 1
    end do
  end subroutine whizard_pythia8_get_final_colored_ME_momenta

  subroutine pythia8_get_shower_mask (pythia, pset, mask, recover_beams)
    class(pythia8_t), intent(in) :: pythia
    type(particle_set_t), intent(in) :: pset
    logical, dimension(:), allocatable, intent(out) :: mask
    logical, intent(in) :: recover_beams
    integer :: skip_beams, i
    type(lha_particle_t) :: c_prt
    if (allocated (mask)) deallocate (mask)
    allocate (mask(pythia%get_event_size ()), source=.true.)
    if (.not. recover_beams) then
       skip_beams = 2
       mask(1:2) = .false.
    else
       skip_beams = 0
    end if
    do i = 1 + skip_beams, size(mask)
       c_prt = pythia%get_single_event (i)
       ! Search for unchanged entries
       mask(i) = (reverse_find_particle (c_prt, pset) == 0)
    end do
  end subroutine pythia8_get_shower_mask

  pure function reverse_find_particle (c_prt, particle_set) result (idx)
    type(lha_particle_t), intent(in) :: c_prt
    type(particle_set_t), intent(in) :: particle_set
    integer :: idx
    type(vector4_t) :: momentum
    momentum = real (c_prt%momentum, default)
    idx = particle_set%reverse_find_particle &
         (c_prt%id, momentum, tiny_10, tiny_07)
  end function reverse_find_particle

  subroutine pythia8_get_hadron_mask (pythia, pset, mask)
    class(pythia8_t), intent(in) :: pythia
    type(particle_set_t), intent(in) :: pset
    logical, dimension(:), allocatable, intent(out) :: mask
    integer :: i
    integer(c_int) :: c_i_prt
    type(lha_particle_t) :: c_prt
    type(vector4_t) :: momentum
    if (allocated (mask)) deallocate (mask)
    allocate (mask(pythia%get_event_size ()), source=.true.)
    do i = 1, size(mask)
       c_prt = pythia%get_single_event (i)
       ! Search for unchanged entries
       mask(i) = (reverse_find_particle (c_prt, pset) == 0)
    end do
  end subroutine pythia8_get_hadron_mask

  subroutine whizard_pythia8_get_shower_particles &
       (pythia, model, model_fallback, particle_set, helicity, recover_beams)
    class(pythia8_t), intent(in) :: pythia
    class(model_data_t), intent(in), target :: model, model_fallback
    type(particle_set_t), intent(inout) :: particle_set
    integer, intent(in), optional :: helicity
    logical, intent(in), optional :: recover_beams
    integer :: n_particles, n_old
    logical, dimension(:), allocatable :: mask
    type(particle_t), dimension(:), allocatable :: particle
    integer :: helicity_opt
    logical :: recover_beams_opt
    if (debug_on) call msg_debug (D_SHOWER, "whizard_pythia8_get_particle_set")
    recover_beams_opt = .false.; if (present (recover_beams)) &
         recover_beams_opt = recover_beams
    helicity_opt = PRT_UNPOLARIZED; if (present (helicity)) &
         helicity_opt = helicity
    call pythia%get_shower_mask (particle_set, mask, recover_beams)
    n_particles = pythia%get_event_size ()
    if (.not. recover_beams_opt) n_particles = n_particles - 2
    allocate (particle(n_particles))
    call pythia%get_particles (model, model_fallback, mask, particle, particle_set, &
        helicity_opt, recover_beams_opt)
  end subroutine whizard_pythia8_get_shower_particles

  subroutine whizard_pythia8_get_hadron_particles &
       (pythia, model, model_fallback, particle_set, helicity)
    class(pythia8_t), intent(in) :: pythia
    class(model_data_t), intent(in), target :: model, model_fallback
    type(particle_set_t), intent(inout) :: particle_set
    integer, intent(in), optional :: helicity
    integer :: n_particles
    logical, dimension(:), allocatable :: mask
    type(particle_t), dimension(:), allocatable :: particle
    integer, dimension(:), allocatable :: pythia_idx, whizard_idx
    integer :: helicity_opt
    if (debug_on) call msg_debug (D_TRANSFORMS, "whizard_pythia8_get_particle_set")
    helicity_opt = PRT_UNPOLARIZED; if (present (helicity)) &
         helicity_opt = helicity
    call pythia%get_hadron_mask (particle_set, mask)
    n_particles = pythia%get_event_size ()
    allocate (particle(n_particles))
    call pythia%get_particles (model, model_fallback, mask, particle, particle_set, &
        helicity_opt, recover_beams = .true.)
  end subroutine whizard_pythia8_get_hadron_particles

  subroutine whizard_pythia8_get_particles (&
     pythia, model, model_fallback, mask, particle, particle_set, &
     helicity, recover_beams)
    class(pythia8_t), intent(in) :: pythia
    class(model_data_t), intent(in), target :: model, model_fallback
    logical, dimension(:), intent(in) :: mask
    type(particle_t), dimension(:), allocatable, intent(inout) :: particle
    type(particle_set_t), intent(inout) :: particle_set
    integer, intent(in) :: helicity
    logical, intent(in) :: recover_beams
    integer, dimension(:), allocatable :: pythia_idx, whizard_idx
    if (debug_on) call msg_debug (D_SHOWER, "whizard_pythia8_get_particles")
    call pythia%import_pythia_particles (&
         model, model_fallback, mask, particle, particle_set, &
         pythia_idx, whizard_idx, helicity, recover_beams)
    call particle_set%replace (particle)
    call pythia%get_parent_child_relation (&
         pythia_idx, whizard_idx, particle_set, recover_beams)
    where ((particle_set%prt%status == PRT_OUTGOING .or. &
            particle_set%prt%status == PRT_VIRTUAL .or. &
            particle_set%prt%status == PRT_BEAM_REMNANT) .and. &
            particle_set%prt%has_children ()) &
            particle_set%prt%status = PRT_RESONANT
  end subroutine whizard_pythia8_get_particles

subroutine pythia8_import_pythia_particles (&
     pythia, model, model_fallback, mask, particle, particle_set, &
     pythia_idx, whizard_idx, helicity, recover_beams)
    class(pythia8_t), intent(in) :: pythia
    class(model_data_t), intent(in), target :: model, model_fallback
    logical, dimension(:), intent(in) :: mask
    type(particle_t), dimension(:), intent(inout) :: particle
    type(particle_set_t), intent(in) :: particle_set
    integer, dimension(:), allocatable, intent(out) :: pythia_idx, whizard_idx
    integer, intent(in) :: helicity
    logical, intent(in) :: recover_beams
    integer :: i_whizard, i_pythia, idx, skip_beams
    real(default) :: time
    real(default), dimension(3) :: vertex
    type(vector4_t) :: momentum
    type(lha_particle_t) :: c_prt
    allocate (whizard_idx(size (mask)), source = 0)
    allocate (pythia_idx(size (particle)), source = 0)
    i_whizard = 0;
    if (recover_beams) then
       skip_beams = 0
    else
       skip_beams = 2
    end if
    ADD_PARTICLE: do i_pythia = 1 + skip_beams, size(mask)
       idx = -1
       c_prt = pythia%get_single_event (i_pythia)
       ! Check on exisiting particle entry
       if (.not. mask(i_pythia)) then
          ! Retrieve particle from original particle_set
          whizard_idx(i_pythia) = reverse_find_particle (c_prt, particle_set)
          idx = reverse_find_particle (c_prt, particle_set)
          ! Skip entry completely
          if (idx > 0) then
             i_whizard = i_whizard + 1
             particle(i_whizard) = particle_set%get_particle (idx)
             whizard_idx(i_pythia) = i_whizard
             pythia_idx(i_whizard) = i_pythia
             if (debug2_active (D_SHOWER)) then
                print *, "Reverse search for particle ", i_pythia, " with PDG: ", c_prt%id
                print *, "Momentum: ", c_prt%momentum
                print *, "Found: ", whizard_idx(i_pythia)
             end if
             cycle ADD_PARTICLE
          end if
          ! Fallthrough: We could not retrieve the particle from our set, retrieve it from PYTHIA8.
       end if
       ! idx is exactly zero iff the reverse particle search failed.
       if (mask(i_pythia) .or. idx == 0) then
          i_whizard = i_whizard + 1
          call get_particle_status (i_pythia, c_prt, particle(i_whizard))
          call fill_particle (model, model_fallback, helicity, c_prt, particle(i_whizard))
          ! call get_particle_color (i_color, dangling_color, c_prt%status, particle(i_whizard))
          call get_particle_vertex (i_pythia, particle(i_whizard))
          whizard_idx(i_pythia) = i_whizard
          pythia_idx(i_whizard) = i_pythia
       end if
       if (debug2_active (D_SHOWER)) then
          print *, "Shower: ", mask(i_pythia)
          print *, "i_pythia: ", i_pythia, " -> i_whizard: ", whizard_idx(i_pythia)
       end if
    end do ADD_PARTICLE
  contains
    subroutine get_particle_status (i_pythia, c_particle, particle)
      integer, intent(in) :: i_pythia
      type(lha_particle_t), intent(in) :: c_particle
      type(particle_t), intent(inout) :: particle
      integer :: whizard_status
      select case (pythia8_get_status_hepmc &
           (pythia%cptr, int(i_pythia, c_int)))
      case(1); whizard_status = (PRT_OUTGOING)
      case(2); whizard_status = (PRT_RESONANT)
      case(4); whizard_status = (PRT_BEAM)
      case default;
         if (c_particle%status < 0) &
              whizard_status = PRT_VIRTUAL
      end select
      if (debug2_active (D_SHOWER) .or. debug2_active (D_TRANSFORMS)) then
         write (*, "(1X,A,1X,I0)") "Particle's status code:", i_pythia
         write (*, "(1X,3(A,1X,I0,1X))") "HEPMC:", pythia8_get_status_hepmc (pythia%cptr, int(i_pythia, c_int)), &
              "PYTHIA:", c_particle%status, &
              "WHIZARD:", whizard_status
      end if
      call particle%set_status (whizard_status)
    end subroutine get_particle_status

    subroutine fill_particle &
         (model_in, model_fallback, polarization, c_prt, prt)
      type(lha_particle_t), intent(in) :: c_prt
      class(model_data_t), intent(in), target :: model_in, model_fallback
      type(particle_t), intent(inout) :: prt
      integer, intent(in) :: polarization
      integer :: whizard_status, hmax
      class(model_data_t), pointer :: model
      type(flavor_t) :: flv
      type(color_t) :: color
      type(helicity_t) :: hel
      type(vector4_t) :: p
      integer :: col, acol
      call find_model (model, c_prt%id, model_in, model_fallback)
      call flv%init (c_prt%id, model)
      col = max (c_prt%color(1), 0)
      acol = max (c_prt%color(2), 0)
      call color%init_col_acl (col, acol)
      if (flv%is_beam_remnant ()) &
           call prt%set_status (PRT_BEAM_REMNANT)
      call prt%set_flavor (flv); call prt%set_color (color)
      p = real (c_prt%momentum, kind=default)
      call prt%set_momentum (p, real (c_prt%mass**2, default))
      select case (polarization)
      case (PRT_DEFINITE_HELICITY)
         if (abs (c_prt%spin) <= 1.) then
            hmax = flv%get_spin_type () / 2.
            call hel%init (sign (hmax, nint (c_prt%spin)))
            call prt%set_helicity (hel)
         end if
      case (PRT_GENERIC_POLARIZATION)
         call msg_fatal ("[whizard_pythia8_get_particle_set]" // &
              "generic polarization with Pythia8 not defined.")
      case (PRT_UNPOLARIZED)
      case default
         call msg_bug ("[whizard_pythia8_get_particle]" // &
              "Helicity handling is undefined.")
      end select
      if (.not. vanishes (real (c_prt%tau, kind=default))) &
           call prt%set_lifetime (real (c_prt%tau, kind=default))
    end subroutine fill_particle

    subroutine get_particle_vertex (i_pythia, particle)
      integer, intent(in) :: i_pythia
      type(particle_t), intent(inout) :: particle
      real(c_double) :: time
      real(c_double), dimension(3) :: vertex
      type(vector4_t) :: vtx4
      call pythia8_get_production_vertex (pythia%cptr, i_pythia, time, vertex)
      vtx4 = vector4_moving (real (time, kind=default), &
           vector3_moving (real (vertex, kind=default)))
      if (vtx4 /= vector4_null) call particle%set_vertex (vtx4)
    end subroutine get_particle_vertex
  end subroutine pythia8_import_pythia_particles

  function pythia8_is_beam_particle (pythia, i_pythia) result (flag)
    class(pythia8_t), intent(in) :: pythia
    integer, intent(in) :: i_pythia
    logical :: flag
    integer(c_int) :: c_i_pythia
    integer, parameter :: HEPMC_BEAM_PRT = 4
    c_i_pythia = int(i_pythia, c_int)
    flag = HEPMC_BEAM_PRT == pythia8_get_status_hepmc (pythia%cptr, c_i_pythia)
  end function pythia8_is_beam_particle

  subroutine pythia8_get_parent_child_relation (&
       pythia, pythia_idx, whizard_idx, particle_set, recover_beams)
    class(pythia8_t), intent(in) :: pythia
    integer, dimension(:), intent(in) :: pythia_idx, whizard_idx
    type(particle_set_t), intent(inout) :: particle_set
    logical, intent(in) :: recover_beams
    integer(c_int) :: c_n_parents, c_n_children, c_i_pythia
    integer, dimension(:), allocatable :: parent, child
    integer :: i_pythia, i, skip_beams
    if (debug_on) call msg_debug (D_SHOWER, "pythia8_get_parent_child_relation")
    skip_beams = 0; if (recover_beams) skip_beams = 2
    do i_pythia = 1 + skip_beams, size(whizard_idx)
       if (whizard_idx(i_pythia) == 0) cycle
       c_i_pythia = int(i_pythia, c_int)
       c_n_parents = pythia8_get_n_mothers (pythia%cptr, c_i_pythia)
       c_n_children = pythia8_get_n_daughters (pythia%cptr, c_i_pythia)
       allocate (parent(c_n_parents), child(c_n_children))
       parent = 0; child = 0
       if (c_n_parents > 0) then
          call pythia8_get_mother_array (pythia%cptr, c_i_pythia, c_n_parents, parent)
          if (count (parent > 0) > 0) then
             if (debug2_active (D_SHOWER) .or. debug2_active (D_TRANSFORMS)) then
                write (*, "(1X,A,1X,I0)") "Particle's parents ", whizard_idx(i_pythia)
                do i = 1, c_n_parents
                   if (parent(i) > 0) &
                        write (*, "(1X,I0,1X,'(',I0,')',1X)", advance="no") parent(i), whizard_idx(parent(i))
                end do
                write (*, *)
             end if
             call particle_set%prt(whizard_idx(i_pythia))%set_parents (&
                  whizard_idx(pack(parent, parent > 0)))
          end if
       end if
       if (c_n_children > 0) then
          call pythia8_get_daughter_array (pythia%cptr, c_i_pythia, c_n_children, child)
          if (count (child > 0) > 0) then
             if (debug2_active (D_SHOWER) .or. debug2_active (D_TRANSFORMS)) then
                write (*, "(1X,A,1X,I0)") "Particle's children ", whizard_idx(i_pythia)
                do i = 1, c_n_children
                   if (child(i) > 0) &
                        write (*, "(1X,I0,1X,'(',I0,')',1X)", advance="no") child(i), whizard_idx(child(i))
                end do
                write (*, *)
             end if
             call particle_set%prt(whizard_idx(i_pythia))%set_children (&
                  whizard_idx(pack (child, child > 0)))
          end if
       end if
       deallocate (parent, child)
    end do
  end subroutine pythia8_get_parent_child_relation


end module pythia8
