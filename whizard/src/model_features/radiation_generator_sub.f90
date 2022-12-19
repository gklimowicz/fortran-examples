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

submodule (radiation_generator) radiation_generator_s

  use io_units
  use diagnostics
  use physics_defs, only: PHOTON, GLUON
  use string_utils, only: split_string, string_contains_word
  use flavors

  implicit none

contains

  module subroutine pdg_states_init (states)
    class(pdg_states_t), intent(inout) :: states
    nullify (states%next)
  end subroutine pdg_states_init

  module subroutine pdg_states_add (states, pdg)
    class(pdg_states_t), intent(inout), target :: states
    type(pdg_array_t), dimension(:), intent(in) :: pdg
    type(pdg_states_t), pointer :: current_state
    select type (states)
    type is (pdg_states_t)
      current_state => states
      do
        if (associated (current_state%next)) then
          current_state => current_state%next
        else
          allocate (current_state%next)
          nullify(current_state%next%next)
          current_state%pdg = pdg
          exit
        end if
      end do
    end select
  end subroutine pdg_states_add

  module function pdg_states_get_n_states (states) result (n)
    class(pdg_states_t), intent(in), target :: states
    integer :: n
    type(pdg_states_t), pointer :: current_state
    n = 0
    select type(states)
    type is (pdg_states_t)
      current_state => states
      do
        if (associated (current_state%next)) then
          n = n+1
          current_state => current_state%next
        else
          exit
        end if
      end do
    end select
  end function pdg_states_get_n_states

  module subroutine prt_queue_null (queue)
    class(prt_queue_t), intent(out) :: queue
    queue%next => null ()
    queue%previous => null ()
    queue%front => null ()
    queue%current_prt => null ()
    queue%back => null ()
    queue%n_lists = 0
    if (allocated (queue%prt_string))  deallocate (queue%prt_string)
  end subroutine prt_queue_null

  module subroutine prt_queue_append (queue, prt_string)
    class(prt_queue_t), intent(inout) :: queue
    type(string_t), intent(in), dimension(:) :: prt_string
    type(prt_queue_t), pointer :: new_element => null ()
    type(prt_queue_t), pointer :: current_back => null ()
    allocate (new_element)
    allocate (new_element%prt_string(size (prt_string)))
    new_element%prt_string = prt_string
    if (associated (queue%back)) then
       current_back => queue%back
       current_back%next => new_element
       new_element%previous => current_back
       queue%back => new_element
    else
       !!! Initial entry
       queue%front => new_element
       queue%back => queue%front
       queue%current_prt => queue%front
    end if
    queue%n_lists = queue%n_lists + 1
  end subroutine prt_queue_append

  module subroutine prt_queue_get (queue, prt_string)
    class(prt_queue_t), intent(inout) :: queue
    type(string_t), dimension(:), allocatable, intent(out) :: prt_string
    if (associated (queue%current_prt)) then
       prt_string = queue%current_prt%prt_string
       if (associated (queue%current_prt%next)) &
          queue%current_prt => queue%current_prt%next
    else
       prt_string = " "
    end if
  end subroutine prt_queue_get

  module subroutine prt_queue_get_last (queue, prt_string)
    class(prt_queue_t), intent(in) :: queue
    type(string_t), dimension(:), allocatable, intent(out) :: prt_string
    if (associated (queue%back)) then
       allocate (prt_string(size (queue%back%prt_string)))
       prt_string = queue%back%prt_string
    else
       prt_string = " "
    end if
  end subroutine prt_queue_get_last

  module subroutine prt_queue_reset (queue)
    class(prt_queue_t), intent(inout) :: queue
    queue%current_prt => queue%front
  end subroutine prt_queue_reset

  module function prt_queue_check_for_same_prt_strings (queue) result (val)
    class(prt_queue_t), intent(inout) :: queue
    logical :: val
    type(string_t), dimension(:), allocatable :: prt_string
    integer, dimension(:,:), allocatable :: i_particle
    integer :: n_d, n_dbar, n_u, n_ubar, n_s, n_sbar, n_gl, n_e, n_ep, n_mu, n_mup, n_A
    integer :: i, j
    call queue%reset ()
    allocate (i_particle (queue%n_lists, 12))
    do i = 1, queue%n_lists
       call queue%get (prt_string)
       n_d = count_particle (prt_string, 1)
       n_dbar = count_particle (prt_string, -1)
       n_u = count_particle (prt_string, 2)
       n_ubar = count_particle (prt_string, -2)
       n_s = count_particle (prt_string, 3)
       n_sbar = count_particle (prt_string, -3)
       n_gl = count_particle (prt_string, 21)
       n_e = count_particle (prt_string, 11)
       n_ep = count_particle (prt_string, -11)
       n_mu = count_particle (prt_string, 13)
       n_mup = count_particle (prt_string, -13)
       n_A = count_particle (prt_string, 22)
       i_particle (i, 1) = n_d
       i_particle (i, 2) = n_dbar
       i_particle (i, 3) = n_u
       i_particle (i, 4) = n_ubar
       i_particle (i, 5) = n_s
       i_particle (i, 6) = n_sbar
       i_particle (i, 7) = n_gl
       i_particle (i, 8) = n_e
       i_particle (i, 9) = n_ep
       i_particle (i, 10) = n_mu
       i_particle (i, 11) = n_mup
       i_particle (i, 12) = n_A
    end do
    val = .false.
    do i = 1, queue%n_lists
       do j = 1, queue%n_lists
          if  (i == j) cycle
          val = val .or. all (i_particle (i,:) == i_particle(j,:))
       end do
    end do
  contains
    function count_particle (prt_string, pdg) result (n)
      type(string_t), dimension(:), intent(in) :: prt_string
      integer, intent(in) :: pdg
      integer :: n
      integer :: i
      type(string_t) :: prt_ref
      n = 0
      select case (pdg)
      case (1)
         prt_ref = "d"
      case (-1)
         prt_ref = "dbar"
      case (2)
         prt_ref = "u"
      case (-2)
         prt_ref = "ubar"
      case (3)
         prt_ref = "s"
      case (-3)
         prt_ref = "sbar"
      case (21)
         prt_ref = "gl"
      case (11)
         prt_ref = "e-"
      case (-11)
         prt_ref = "e+"
      case (13)
         prt_ref = "mu-"
      case (-13)
         prt_ref = "mu+"
      case (22)
         prt_ref = "A"
      end select
      do i = 1, size (prt_string)
         if (prt_string(i) == prt_ref) n = n+1
      end do
    end function count_particle

  end function prt_queue_check_for_same_prt_strings

  module function prt_queue_contains (queue, prt_string) result (val)
    class(prt_queue_t), intent(in) :: queue
    type(string_t), intent(in), dimension(:) :: prt_string
    logical :: val
    type(prt_queue_t), pointer :: current => null()
    if (associated (queue%front)) then
       current => queue%front
    else
       call msg_fatal ("Trying to access empty particle queue")
    end if
    val = .false.
    do
       if (size (current%prt_string) == size (prt_string)) then
          if (all (current%prt_string == prt_string)) then
             val = .true.
             exit
          end if
       end if
       if (associated (current%next)) then
          current => current%next
       else
          exit
       end if
    end do
  end function prt_queue_contains

  module subroutine prt_queue_write (queue, unit)
    class(prt_queue_t), intent(in) :: queue
    integer, optional :: unit
    type(prt_queue_t), pointer :: current => null ()
    integer :: i, j, u
    u = given_output_unit (unit)
    if (associated (queue%front)) then
       current => queue%front
    else
       write (u, "(A)") "[Particle queue is empty]"
       return
    end if
    j = 1
    do
       write (u, "(I2,A,1X)", advance = 'no') j , ":"
       do i = 1, size (current%prt_string)
          write (u, "(A,1X)", advance = 'no') char (current%prt_string(i))
       end do
       write (u, "(A)")
       if (associated (current%next)) then
          current => current%next
          j = j+1
       else
          exit
       end if
    end do
  end subroutine prt_queue_write

  subroutine sort_prt (prt, model)
    type(string_t), dimension(:), intent(inout) :: prt
    class(model_data_t), intent(in), target :: model
    type(pdg_array_t), dimension(:), allocatable :: pdg
    type(flavor_t) :: flv
    integer :: i
    call create_pdg_array (prt, model, pdg)
    call sort_pdg (pdg)
    do i = 1, size (pdg)
       call flv%init (pdg(i)%get(), model)
       prt(i) = flv%get_name ()
    end do
  end subroutine sort_prt

  subroutine sort_pdg (pdg)
    type(pdg_array_t), dimension(:), intent(inout) :: pdg
    integer, dimension(:), allocatable :: i_pdg
    integer :: i
    allocate (i_pdg (size (pdg)))
    do i = 1, size (pdg)
       i_pdg(i) = pdg(i)%get ()
    end do
    i_pdg = sort_abs (i_pdg)
    do i = 1, size (pdg)
       call pdg(i)%set (1, i_pdg(i))
    end do
  end subroutine sort_pdg

  subroutine create_pdg_array (prt, model, pdg)
    type (string_t), dimension(:), intent(in) :: prt
    class (model_data_t), intent(in), target :: model
    type(pdg_array_t), dimension(:), allocatable, intent(out) :: pdg
    type(flavor_t) :: flv
    integer :: i
    allocate (pdg (size (prt)))
    do i = 1, size (prt)
       call flv%init (prt(i), model)
       pdg(i) = flv%get_pdg ()
    end do
  end subroutine create_pdg_array

  module subroutine reshuffle_list_write (rlist)
    class(reshuffle_list_t), intent(in) :: rlist
    type(reshuffle_list_t), pointer :: current => null ()
    integer :: i
    print *, 'Content of reshuffling list: '
    if (associated (rlist%next)) then
       current => rlist%next
       i = 1
       do
         print *, 'i: ', i, 'list: ', current%ii
         i = i + 1
         if (associated (current%next)) then
            current => current%next
         else
            exit
         end if
       end do
    else
       print *, '[EMPTY]'
    end if
  end subroutine reshuffle_list_write

  module subroutine reshuffle_list_append (rlist, ii)
    class(reshuffle_list_t), intent(inout) :: rlist
    integer, dimension(:), allocatable, intent(in) :: ii
    type(reshuffle_list_t), pointer :: current
    if (associated (rlist%next)) then
       current => rlist%next
       do
          if (associated (current%next)) then
             current => current%next
          else
             allocate (current%next)
             allocate (current%next%ii (size (ii)))
             current%next%ii = ii
             exit
          end if
       end do
    else
       allocate (rlist%next)
       allocate (rlist%next%ii (size (ii)))
       rlist%next%ii = ii
    end if
  end subroutine reshuffle_list_append

  elemental module function reshuffle_list_is_empty (rlist) result (is_empty)
    logical :: is_empty
    class(reshuffle_list_t), intent(in) :: rlist
    is_empty = .not. associated (rlist%next)
  end function reshuffle_list_is_empty

  module function reshuffle_list_get (rlist, index) result (ii)
    integer, dimension(:), allocatable :: ii
    class(reshuffle_list_t), intent(inout) :: rlist
    integer, intent(in) :: index
    type(reshuffle_list_t), pointer :: current => null ()
    integer :: i
    current => rlist%next
    do i = 1, index - 1
       if (associated (current%next)) then
          current => current%next
       else
          call msg_fatal ("Index exceeds size of reshuffling list")
       end if
    end do
    allocate (ii (size (current%ii)))
    ii = current%ii
  end function reshuffle_list_get

  module subroutine reshuffle_list_reset (rlist)
    class(reshuffle_list_t), intent(inout) :: rlist
    rlist%next => null ()
  end subroutine reshuffle_list_reset

  module subroutine radiation_generator_init_pdg_list &
       (generator, pl_in, pl_out, pl_excluded_gauge_splittings, qcd, qed)
    class(radiation_generator_t), intent(inout) :: generator
    type(pdg_list_t), intent(in) :: pl_in, pl_out
    type(pdg_list_t), intent(in) :: pl_excluded_gauge_splittings
    logical, intent(in), optional :: qcd, qed
    if (present (qcd))  generator%qcd_enabled = qcd
    if (present (qed))  generator%qed_enabled = qed
    generator%pl_in = pl_in
    generator%pl_out = pl_out
    generator%pl_excluded_gauge_splittings = pl_excluded_gauge_splittings
    generator%is_gluon = pl_in%search_for_particle (GLUON)
    generator%fs_gluon = pl_out%search_for_particle (GLUON)
    generator%is_photon = pl_in%search_for_particle (PHOTON)
    generator%fs_photon = pl_out%search_for_particle (PHOTON)
    generator%mass_sum = 0._default
    call generator%pdg_raw%init ()
  end subroutine radiation_generator_init_pdg_list

  module subroutine radiation_generator_init_pdg_array &
       (generator, pdg_in, pdg_out, pdg_excluded_gauge_splittings, qcd, qed)
    class(radiation_generator_t), intent(inout) :: generator
    type(pdg_array_t), intent(in), dimension(:) :: pdg_in, pdg_out
    type(pdg_array_t), intent(in), dimension(:) :: pdg_excluded_gauge_splittings
    logical, intent(in), optional :: qcd, qed
    type(pdg_list_t) :: pl_in, pl_out
    type(pdg_list_t) :: pl_excluded_gauge_splittings
    integer :: i
    call pl_in%init(size (pdg_in))
    call pl_out%init(size (pdg_out))
    do i = 1, size (pdg_in)
       call pl_in%set (i, pdg_in(i))
    end do
    do i = 1, size (pdg_out)
       call pl_out%set (i, pdg_out(i))
    end do
    call pl_excluded_gauge_splittings%init(size (pdg_excluded_gauge_splittings))
    do i = 1, size (pdg_excluded_gauge_splittings)
       call pl_excluded_gauge_splittings%set &
            (i, pdg_excluded_gauge_splittings(i))
    end do
    call generator%init (pl_in, pl_out, pl_excluded_gauge_splittings, qcd, qed)
  end subroutine radiation_generator_init_pdg_array

  module subroutine radiation_generator_set_initial_state_emissions (generator)
     class(radiation_generator_t), intent(inout) :: generator
     generator%only_final_state = .false.
  end subroutine radiation_generator_set_initial_state_emissions

  module subroutine radiation_generator_setup_if_table (generator, model)
    class(radiation_generator_t), intent(inout) :: generator
    class(model_data_t), intent(in), target :: model
    type(pdg_list_t), dimension(:), allocatable :: pl_in, pl_out

    allocate (pl_in(1), pl_out(1))

    pl_in(1) = generator%pl_in
    pl_out(1) = generator%pl_out

    call generator%if_table%init &
         (model, pl_in, pl_out, generator%constraints)
  end subroutine radiation_generator_setup_if_table

  module subroutine radiation_generator_reset_particle_content_pdg_list (generator, pl)
    class(radiation_generator_t), intent(inout) :: generator
    type(pdg_list_t), intent(in) :: pl
    generator%pl_out = pl
    generator%fs_gluon = pl%search_for_particle (GLUON)
    generator%fs_photon = pl%search_for_particle (PHOTON)
  end subroutine radiation_generator_reset_particle_content_pdg_list

  module subroutine radiation_generator_reset_particle_content_pdg_array (generator, pdg)
    class(radiation_generator_t), intent(inout) :: generator
    type(pdg_array_t), intent(in), dimension(:) :: pdg
    type(pdg_list_t) :: pl
    integer :: i
    call pl%init (size (pdg))
    do i = 1, size (pdg)
       call pl%set (i, pdg(i))
    end do
    call generator%reset_particle_content (pl)
  end subroutine radiation_generator_reset_particle_content_pdg_array

  module subroutine radiation_generator_reset_reshuffle_list (generator)
    class(radiation_generator_t), intent(inout) :: generator
    call generator%reshuffle_list%reset ()
  end subroutine radiation_generator_reset_reshuffle_list

  module subroutine radiation_generator_set_n (generator, n_in, n_out, n_loops)
    class(radiation_generator_t), intent(inout) :: generator
    integer, intent(in) :: n_in, n_out, n_loops
    generator%n_tot = n_in + n_out + 1
    generator%n_in = n_in
    generator%n_out = n_out
    generator%n_loops = n_loops
  end subroutine radiation_generator_set_n

  module subroutine radiation_generator_set_constraints &
       (generator, set_n_loop, set_mass_sum, &
        set_selected_particles, set_required_particles)
    class(radiation_generator_t), intent(inout), target :: generator
    logical, intent(in) :: set_n_loop
    logical, intent(in) :: set_mass_sum
    logical, intent(in) :: set_selected_particles
    logical, intent(in) :: set_required_particles
    logical :: set_no_photon_induced = .true.
    integer :: i, j, n, n_constraints
    type(pdg_list_t) :: pl_req, pl_insert
    type(pdg_list_t) :: pl_antiparticles
    type(pdg_array_t) :: pdg_gluon, pdg_photon
    type(pdg_array_t) :: pdg_add, pdg_tmp
    integer :: last_index
    integer :: n_new_particles, n_skip
    integer, dimension(:), allocatable :: i_skip
    integer :: n_nlo_correction_types

    n_nlo_correction_types = count ([generator%qcd_enabled, generator%qed_enabled])
    if (generator%is_photon) set_no_photon_induced = .false.

    allocate (i_skip (generator%n_tot))
    i_skip = -1

    n_constraints = 2 + count([set_n_loop, set_mass_sum, &
         set_selected_particles, set_required_particles, set_no_photon_induced])
    associate (constraints => generator%constraints)
      n = 1
      call constraints%init (n_constraints)
      call constraints%set (n, constrain_n_tot (generator%n_tot))
      n = 2
      call constraints%set (n, constrain_couplings (generator%qcd_enabled, &
           generator%qed_enabled, n_nlo_correction_types))
      n = n + 1
      if (set_no_photon_induced) then
         call constraints%set (n, constrain_photon_induced_processes (generator%n_in))
         n = n + 1
      end if
      if (set_n_loop) then
         call constraints%set (n, constrain_n_loop(generator%n_loops))
         n = n + 1
      end if
      if (set_mass_sum) then
         call constraints%set (n, constrain_mass_sum(generator%mass_sum))
         n = n + 1
      end if
      if (set_required_particles) then
         if (generator%fs_gluon .or. generator%fs_photon) then
            do i = 1, generator%n_out
               pdg_tmp = generator%pl_out%get(i)
               if (pdg_tmp%search_for_particle (GLUON) &
                    .or. pdg_tmp%search_for_particle (PHOTON)) then
                  i_skip(i) = i
               end if
            end do

            n_skip = count (i_skip > 0)
            call pl_req%init (generator%n_out-n_skip)
         else
            call pl_req%init (generator%n_out)
         end if
         j = 1
         do i = 1, generator%n_out
            if (any (i == i_skip)) cycle
            call pl_req%set (j, generator%pl_out%get(i))
            j = j + 1
         end do
         call constraints%set (n, constrain_require (pl_req))
         n = n + 1
      end if
      if (set_selected_particles) then
         if (generator%only_final_state ) then
            call pl_insert%init (generator%n_out + n_nlo_correction_types)
            do i = 1, generator%n_out
               call pl_insert%set(i, generator%pl_out%get(i))
            end do
            last_index = generator%n_out + 1
         else
            call generator%pl_in%create_antiparticles (pl_antiparticles, n_new_particles)
            call pl_insert%init (generator%n_tot + n_new_particles &
                 + n_nlo_correction_types)
            do i = 1, generator%n_in
               call pl_insert%set(i, generator%pl_in%get(i))
            end do
            do i = 1, generator%n_out
               j = i + generator%n_in
               call pl_insert%set(j, generator%pl_out%get(i))
            end do
            do i = 1, n_new_particles
               j = i + generator%n_in + generator%n_out
               call pl_insert%set(j, pl_antiparticles%get(i))
            end do
            last_index = generator%n_tot + n_new_particles + 1
         end if
         pdg_gluon = GLUON; pdg_photon = PHOTON
         if (generator%qcd_enabled) then
            pdg_add = pdg_gluon
            call pl_insert%set (last_index, pdg_add)
            last_index = last_index + 1
         end if
         if (generator%qed_enabled) then
            pdg_add = pdg_photon
            call pl_insert%set (last_index, pdg_add)
         end if
         call constraints%set (n, constrain_splittings (pl_insert, &
              generator%pl_excluded_gauge_splittings))
      end if
    end associate
  end subroutine radiation_generator_set_constraints

  module subroutine radiation_generator_find_splittings (generator)
    class(radiation_generator_t), intent(inout) :: generator
    integer :: i
    type(pdg_array_t), dimension(:), allocatable :: pdg_in, pdg_out, pdg_tmp
    integer, dimension(:), allocatable :: reshuffle_list

    call generator%pl_in%create_pdg_array (pdg_in)
    call generator%pl_out%create_pdg_array (pdg_out)

    associate (if_table => generator%if_table)
       call if_table%radiate (generator%constraints, do_not_check_regular = .true.)

       do i = 1, if_table%get_length ()
          call if_table%get_pdg_out (i, pdg_tmp)
          if (size (pdg_tmp) == generator%n_tot) then
             call pdg_reshuffle (pdg_out, pdg_tmp, reshuffle_list)
             call generator%reshuffle_list%append (reshuffle_list)
          end if
       end do
    end associate

  contains

    subroutine pdg_reshuffle (pdg_born, pdg_real, list)
      type(pdg_array_t), intent(in), dimension(:) :: pdg_born, pdg_real
      integer, intent(out), dimension(:), allocatable :: list
      type(pdg_sorter_t), dimension(:), allocatable :: sort_born
      type(pdg_sorter_t), dimension(:), allocatable :: sort_real
      integer :: i_min, n_in, n_born, n_real
      integer :: ib, ir

      n_in = generator%n_in
      n_born = size (pdg_born)
      n_real = size (pdg_real)
      allocate (list (n_real - n_in))
      allocate (sort_born (n_born))
      allocate (sort_real (n_real - n_in))
      sort_born%pdg = pdg_born%get ()
      sort_real%pdg = pdg_real(n_in + 1 : n_real)%get()
      do ib = 1, n_born
         if (any (sort_born(ib)%pdg == sort_real%pdg)) &
            call associate_born_indices (sort_born(ib), sort_real, ib, n_real)
      end do
      i_min = maxval (sort_real%associated_born) + 1
      do ir = 1, n_real - n_in
         if (sort_real(ir)%associated_born == 0) then
            sort_real(ir)%associated_born = i_min
            i_min = i_min + 1
         end if
      end do
      list = sort_real%associated_born
    end subroutine pdg_reshuffle

    subroutine associate_born_indices (sort_born, sort_real, ib, n_real)
      type(pdg_sorter_t), intent(in) :: sort_born
      type(pdg_sorter_t), intent(inout), dimension(:) :: sort_real
      integer, intent(in) :: ib, n_real
      integer :: ir
      do ir = 1, n_real - generator%n_in
         if (sort_born%pdg == sort_real(ir)%pdg &
            .and..not. sort_real(ir)%checked) then
            sort_real(ir)%associated_born = ib
            sort_real(ir)%checked = .true.
            exit
        end if
      end do
    end subroutine associate_born_indices
  end subroutine radiation_generator_find_splittings

  module subroutine radiation_generator_generate_real_particle_strings &
       (generator, prt_tot_in, prt_tot_out)
    class(radiation_generator_t), intent(inout) :: generator
    type(string_t), intent(out), dimension(:), allocatable :: prt_tot_in, prt_tot_out
    type(prt_array_t), dimension(:), allocatable :: prt_in, prt_out
    type(prt_array_t), dimension(:), allocatable :: prt_out0, prt_in0
    type(pdg_array_t), dimension(:), allocatable :: pdg_tmp, pdg_out, pdg_in
    type(pdg_list_t), dimension(:), allocatable :: pl_in, pl_out
    type(prt_array_t) :: prt_out0_tmp, prt_in0_tmp
    integer :: i, j
    integer, dimension(:), allocatable :: reshuffle_list_local
    type(reshuffle_list_t) :: reshuffle_list
    integer :: flv
    type(string_t), dimension(:), allocatable :: buf
    integer :: i_buf

    flv = 0
    allocate (prt_in0(0), prt_out0(0))
    associate (if_table => generator%if_table)
       do i = 1, if_table%get_length ()
          call if_table%get_pdg_out (i, pdg_tmp)
          if (size (pdg_tmp) == generator%n_tot) then
             call if_table%get_particle_string (i, &
                  prt_in0_tmp%prt, prt_out0_tmp%prt)
             prt_in0 = [prt_in0, prt_in0_tmp]
             prt_out0 = [prt_out0, prt_out0_tmp]
             flv = flv + 1
          end if
       end do
    end associate

    allocate (prt_in(size (prt_in0)), prt_out(size (prt_out0)))
    do i = 1, flv
       allocate (prt_in(i)%prt (generator%n_in))
       allocate (prt_out(i)%prt (generator%n_tot - generator%n_in))
    end do
    allocate (prt_tot_in (generator%n_in))
    allocate (prt_tot_out (generator%n_tot - generator%n_in))
    allocate (buf (generator%n_tot))
    buf = ""

    do j = 1, flv
       do i = 1, generator%n_in
          prt_in(j)%prt(i) = prt_in0(j)%prt(i)
          call fill_buffer (buf(i), prt_in0(j)%prt(i))
       end do
    end do
    prt_tot_in = buf(1 : generator%n_in)

    do j = 1, flv
       allocate (reshuffle_list_local (size (generator%reshuffle_list%get(j))))
       reshuffle_list_local = generator%reshuffle_list%get(j)
       do i = 1, size (reshuffle_list_local)
          prt_out(j)%prt(reshuffle_list_local(i)) = prt_out0(j)%prt(i)
          i_buf = reshuffle_list_local(i) + generator%n_in
          call fill_buffer (buf(i_buf), &
               prt_out(j)%prt(reshuffle_list_local(i)))
       end do
       !!! Need to deallocate here because in the next iteration the reshuffling
       !!! list can have a different size
       deallocate (reshuffle_list_local)
    end do
    prt_tot_out = buf(generator%n_in + 1 : generator%n_tot)
    if (debug2_active (D_CORE)) then
       print *, 'Generated initial state: '
       do i = 1, size (prt_tot_in)
          print *, char (prt_tot_in(i))
       end do
       print *, 'Generated final state: '
       do i = 1, size (prt_tot_out)
          print *, char (prt_tot_out(i))
       end do
    end if


  contains

    subroutine fill_buffer (buffer, particle)
      type(string_t), intent(inout) :: buffer
      type(string_t), intent(in) :: particle
      logical :: particle_present
      if (len (buffer) > 0) then
         particle_present = check_for_substring (char(buffer), particle)
         if (.not. particle_present) buffer = buffer // ":" // particle
      else
         buffer = buffer // particle
      end if
    end subroutine fill_buffer

    function check_for_substring (buffer, substring) result (exist)
      character(len=*), intent(in) :: buffer
      type(string_t), intent(in) :: substring
      character(len=50) :: buffer_internal
      logical :: exist
      integer :: i_first, i_last
      exist = .false.
      i_first = 1; i_last = 1
      do
         if (buffer(i_last:i_last) == ":") then
            buffer_internal = buffer (i_first : i_last - 1)
            if (buffer_internal == char (substring)) then
               exist = .true.
               exit
            end if
            i_first = i_last + 1; i_last = i_first + 1
            if (i_last > len(buffer)) exit
         else if (i_last == len(buffer)) then
            buffer_internal = buffer (i_first : i_last)
            exist = buffer_internal == char (substring)
            exit
         else
            i_last = i_last + 1
            if (i_last > len(buffer)) exit
         end if
      end do
    end function check_for_substring
  end subroutine radiation_generator_generate_real_particle_strings

  module function radiation_generator_contains_emissions (generator) result (has_em)
    logical :: has_em
    class(radiation_generator_t), intent(in) :: generator
    has_em = .not. generator%reshuffle_list%is_empty ()
  end function radiation_generator_contains_emissions

  module subroutine radiation_generator_generate (generator, prt_in, prt_out)
    class(radiation_generator_t), intent(inout) :: generator
    type(string_t), intent(out), dimension(:), allocatable :: prt_in, prt_out
    call generator%find_splittings ()
    call generator%generate_real_particle_strings (prt_in, prt_out)
  end subroutine radiation_generator_generate

  module subroutine radiation_generator_generate_multiple &
       (generator, max_multiplicity, model)
    class(radiation_generator_t), intent(inout) :: generator
    integer, intent(in) :: max_multiplicity
    class(model_data_t), intent(in), target :: model
    if (max_multiplicity <= generator%n_out) &
         call msg_fatal ("GKS states: Multiplicity is not large enough!")
    call generator%first_emission (model)
    call generator%reset_reshuffle_list ()
    if (max_multiplicity - generator%n_out > 1) &
         call generator%append_emissions (max_multiplicity, model)
  end subroutine radiation_generator_generate_multiple

  module subroutine radiation_generator_first_emission (generator, model)
    class(radiation_generator_t), intent(inout) :: generator
    class(model_data_t), intent(in), target :: model
    type(string_t), dimension(:), allocatable :: prt_in, prt_out
    call generator%setup_if_table (model)
    call generator%generate (prt_in, prt_out)
    call generator%prt_queue%null ()
    call generator%prt_queue%append (prt_out)
  end subroutine radiation_generator_first_emission

  module subroutine radiation_generator_append_emissions &
       (generator, max_multiplicity, model)
    class(radiation_generator_t), intent(inout) :: generator
    integer, intent(in) :: max_multiplicity
    class(model_data_t), intent(in), target :: model
    type(string_t), dimension(:), allocatable :: prt_fetched
    type(string_t), dimension(:), allocatable :: prt_in
    type(string_t), dimension(:), allocatable :: prt_out
    type(pdg_array_t), dimension(:), allocatable :: pdg_new_out
    integer :: current_multiplicity, i, j, n_longest_length
    type(prt_table_t), dimension(:), allocatable :: prt_table_out
    do
       call generator%prt_queue%get (prt_fetched)
       current_multiplicity = size (prt_fetched)
       if (current_multiplicity == max_multiplicity) exit
       call create_pdg_array (prt_fetched, model, &
            pdg_new_out)
       call generator%reset_particle_content (pdg_new_out)
       call generator%set_n (2, current_multiplicity, 0)
       call generator%set_constraints (.false., .false., .true., .true.)
       call generator%setup_if_table (model)
       call generator%generate (prt_in, prt_out)
       n_longest_length = get_length_of_longest_tuple (prt_out)
       call separate_particles (prt_out, prt_table_out)
       do i = 1, n_longest_length
          if (.not. any (prt_table_out(i)%prt == " ")) then
             call sort_prt (prt_table_out(i)%prt, model)
             if (.not. generator%prt_queue%contains (prt_table_out(i)%prt)) then
                call generator%prt_queue%append (prt_table_out(i)%prt)
             end if
          end if
       end do
       call generator%reset_reshuffle_list ()
    end do

  contains
    subroutine separate_particles (prt, prt_table)
      type(string_t), intent(in), dimension(:) :: prt
      type(string_t), dimension(:), allocatable :: prt_tmp
      type(prt_table_t), intent(out), dimension(:), allocatable :: prt_table
      integer :: i, j
      logical, dimension(:), allocatable :: tuples_occured
      allocate (prt_table (n_longest_length))
      do i = 1, n_longest_length
         allocate (prt_table(i)%prt (size (prt)))
      end do
      allocate (tuples_occured (size (prt)))
      do j = 1, size (prt)
         call split_string (prt(j), var_str (":"), prt_tmp)
         do i = 1, n_longest_length
            if (i <= size (prt_tmp)) then
               prt_table(i)%prt(j) = prt_tmp(i)
            else
               prt_table(i)%prt(j) = " "
            end if
         end do
         if (n_longest_length > 1) &
              tuples_occured(j) = prt_table(1)%prt(j) /= " " &
              .and. prt_table(2)%prt(j) /= " "
      end do
      if (any (tuples_occured)) then
         do j = 1, size (tuples_occured)
            if (.not. tuples_occured(j)) then
               do i = 2, n_longest_length
                  prt_table(i)%prt(j) = prt_table(1)%prt(j)
               end do
            end if
         end do
      end if
    end subroutine separate_particles

  function get_length_of_longest_tuple (prt) result (longest_length)
      type(string_t), intent(in), dimension(:) :: prt
      integer :: longest_length, i
      type(prt_table_t), dimension(:), allocatable :: prt_table
      allocate (prt_table (size (prt)))
      longest_length = 0
      do i = 1, size (prt)
         call split_string (prt(i), var_str (":"), prt_table(i)%prt)
         if (size (prt_table(i)%prt) > longest_length) &
              longest_length = size (prt_table(i)%prt)
      end do
    end function get_length_of_longest_tuple
  end subroutine radiation_generator_append_emissions
  module subroutine radiation_generator_reset_queue (generator)
    class(radiation_generator_t), intent(inout) :: generator
    call generator%prt_queue%reset ()
  end subroutine radiation_generator_reset_queue

  module function radiation_generator_get_n_gks_states (generator) result (n)
    class(radiation_generator_t), intent(in) :: generator
    integer :: n
    n = generator%prt_queue%n_lists
  end function radiation_generator_get_n_gks_states

  module function radiation_generator_get_next_state (generator) result (prt_string)
    class(radiation_generator_t), intent(inout) :: generator
    type(string_t), dimension(:), allocatable :: prt_string
    call generator%prt_queue%get (prt_string)
  end function radiation_generator_get_next_state

  module subroutine radiation_generator_get_emitter_indices (generator, indices)
    class(radiation_generator_t), intent(in) :: generator
    integer, dimension(:), allocatable, intent(out) :: indices
    type(pdg_array_t), dimension(:), allocatable :: pdg_in, pdg_out
    integer, dimension(:), allocatable :: flv_in, flv_out
    integer, dimension(:), allocatable :: emitters
    integer :: i, j
    integer :: n_in, n_out
    
    call generator%pl_in%create_pdg_array (pdg_in)
    call generator%pl_out%create_pdg_array (pdg_out)
    
    n_in = size (pdg_in); n_out = size (pdg_out)
    allocate (flv_in (n_in), flv_out (n_out))
    forall (i=1:n_in) flv_in(i) = pdg_in(i)%get()
    forall (i=1:n_out) flv_out(i) = pdg_out(i)%get()
    
    call generator%if_table%get_emitters (generator%constraints, emitters)
    allocate (indices (size (emitters)))
    
    j = 1
    do i = 1, n_in + n_out
       if (i <= n_in) then
          if (any (flv_in(i) == emitters)) then
             indices (j) = i
             j = j + 1
          end if
       else
          if (any (flv_out(i-n_in) == emitters)) then
             indices (j) = i
             j = j + 1
          end if
       end if
    end do
  end subroutine radiation_generator_get_emitter_indices

  module function radiation_generator_get_raw_states (generator) result (raw_states)
    class(radiation_generator_t), intent(in), target :: generator
    integer, dimension(:,:), allocatable :: raw_states
    type(pdg_states_t), pointer :: state
    integer :: n_states, n_particles
    integer :: i_state
    integer :: j
    state => generator%pdg_raw
    n_states = generator%pdg_raw%get_n_states ()
    n_particles = size (generator%pdg_raw%pdg)
    allocate (raw_states (n_particles, n_states))
    do i_state = 1, n_states
      do j = 1, n_particles
        raw_states (j, i_state) = state%pdg(j)%get ()
      end do
        state => state%next
    end do
  end function radiation_generator_get_raw_states

  module subroutine radiation_generator_save_born_raw (generator, pdg_in, pdg_out)
    class(radiation_generator_t), intent(inout) :: generator
    type(pdg_array_t), dimension(:), allocatable, intent(in) :: pdg_in, pdg_out
    generator%pdg_in_born = pdg_in
    generator%pdg_out_born = pdg_out
  end subroutine radiation_generator_save_born_raw
  module function radiation_generator_get_born_raw (generator) result (flv_born)
    class(radiation_generator_t), intent(in) :: generator
    integer, dimension(:,:), allocatable :: flv_born
    integer :: i_part, n_particles
    n_particles = size (generator%pdg_in_born) + size (generator%pdg_out_born)
    allocate (flv_born (n_particles, 1))
    flv_born(1,1) = generator%pdg_in_born(1)%get ()
    flv_born(2,1) = generator%pdg_in_born(2)%get ()
    do i_part = 3, n_particles
      flv_born(i_part, 1) = generator%pdg_out_born(i_part-2)%get ()
    end do
  end function radiation_generator_get_born_raw


end submodule radiation_generator_s

