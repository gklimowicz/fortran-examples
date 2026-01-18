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

module whizard_lha
  use kinds, only: default

  use, intrinsic :: iso_c_binding
  use diagnostics
  use format_utils, only: write_separator
  use lorentz
  use io_units, only: given_output_unit
  use polarizations
  use particles
  use subevents, only: PRT_BEAM, PRT_INCOMING, PRT_OUTGOING, &
       PRT_UNDEFINED, PRT_VIRTUAL, PRT_RESONANT, PRT_BEAM_REMNANT


  implicit none
  private

  public :: lha_particle_t
  public :: whizard_lha_t



  type, bind(C) :: lha_particle_t
     integer(c_int) :: id, status
     integer(c_int), dimension(2) :: mother
     integer(c_int), dimension(2) :: color
     real(c_double), dimension(4) :: momentum
     real(c_double) :: mass, tau, spin
  end type lha_particle_t

  type :: whizard_lha_t
     private
     type(c_ptr) :: cptr
     logical :: new_event = .false.
   contains
     procedure :: init => whizard_lha_init
     procedure :: final => whizard_lha_final
     procedure :: get_ptr => whizard_lha_get_ptr
     procedure :: set_init => whizard_lha_set_init
     procedure :: set_process_parameters => whizard_lha_set_process_parameters
     procedure :: list_init => whizard_lha_list_init
     procedure :: list_event => whizard_lha_list_event
     procedure :: set_event_process => whizard_lha_set_event_process
     procedure :: set_event => whizard_lha_set_event
  end type whizard_lha_t


  interface lha_particle_write
     module procedure lha_particle_write_single, lha_particle_write_array
  end interface lha_particle_write
  interface
     function new_whizard_lha () bind(C) result (cptr)
       import
       type(c_ptr) :: cptr
     end function new_whizard_lha
  end interface

    interface
     subroutine lhaup_whizard_delete (cptr) bind(C)
       import
       ! Attribute value cannot have intent(inout).
       type(c_ptr), value :: cptr
     end subroutine lhaup_whizard_delete
  end interface

  interface
     function lhaup_whizard_set_init (cptr, beam_pdg, beam_energy, n_processes, unweighted, negative_weights) bind(C) result (flag)
       import
       type(c_ptr), value :: cptr
       integer(c_int), dimension(2), intent(in) :: beam_pdg
       real(c_double), dimension(2), intent(in) :: beam_energy
       integer(c_int), intent(in), value :: n_processes
       logical(c_bool), intent(in), value :: unweighted, negative_weights
       logical(c_bool) :: flag
     end function lhaup_whizard_set_init
  end interface

  interface
     function lhaup_whizard_set_process_parameters &
          (cptr, process_id, cross_section, error, max_weight) &
          bind(C) result (flag)
       import
       type(c_ptr), value :: cptr
       integer(c_int), intent(in), value :: process_id
       real(c_double), intent(in), value :: cross_section, error, max_weight
       logical(c_bool) :: flag
     end function lhaup_whizard_set_process_parameters
  end interface

  interface
     subroutine lhaup_whizard_list_init (cptr) bind(C)
       import
       type(c_ptr), value :: cptr
     end subroutine lhaup_whizard_list_init
  end interface

  interface
     subroutine lhaup_whizard_list_event (cptr) bind(C)
       import
       type(c_ptr), value :: cptr
     end subroutine lhaup_whizard_list_event
  end interface

  interface
     subroutine lhaup_whizard_set_event_process &
          (cptr, process_id, scale, alpha_qcd, alpha_qed, weight) bind(C)
       import
       type(c_ptr), value :: cptr
       integer(c_int), intent(in), value :: process_id
       real(c_double), intent(in), value :: scale, alpha_qcd, alpha_qed, weight
     end subroutine lhaup_whizard_set_event_process
  end interface

  interface
     function lhaup_whizard_set_event (cptr, process_id, n_particles, &
          particle_set) bind(C) result (flag)
       import
       type(c_ptr), value :: cptr
       integer(c_int), intent(in), value :: process_id
       integer(c_int), intent(in), value :: n_particles
       ! IMPORTANT NOTE: Assumed-size array has to be defined by *.
       type(lha_particle_t), dimension(*), intent(in) :: particle_set
       logical(c_bool) :: flag
     end function lhaup_whizard_set_event
  end interface


contains

  subroutine lha_particle_write_single (particle, unit)
    type(lha_particle_t) :: particle
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(I9,1X)", advance="no") particle%id
    write (u, "(I4,1X)", advance="no") particle%status
    write (u, "(2(I5,1X))", advance="no") particle%mother(1), particle%mother(2)
    write (u, "(2(I5,1X))", advance="no") particle%color(1), particle%color(2)
    write (u, "(5(F11.3,1X))", advance="no") particle%momentum(2), particle%momentum(3), &
         particle%momentum(4), particle%momentum(1), particle%mass
    write (u, "(F8.3,1X)", advance="no") particle%tau
    write (u, "(F8.3,1X)", advance="no") particle%tau
    write (u, "(A)")
  end subroutine lha_particle_write_single

  subroutine lha_particle_write_array (particle_set, unit)
    type(lha_particle_t), dimension(:), intent(in) :: particle_set
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(1X,A)") "LHA Particle set:"
    call write_separator (u)
    write (u, "((A4,1X),(A9,1X),(A4,1X),(1X,A10,1X),(1X,A10,1X),5(A11,1X),2(A8,1X))") &
         "No", "ID", "Stat", "Mothers", "Colours", &
         "P(1)", "P(2)", "P(3)", "E", "M", "Tau", "Spin"
    if (size (particle_set) == 0) then
       write (u, "(3X,A)") "[empty]"
    else
       do i = 1, size(particle_set)
          write (u, "(I4,1X)", advance="no") i
          call lha_particle_write_single (particle_set(i), unit)
       end do
    end if
    call write_separator (u)
  end subroutine lha_particle_write_array

  subroutine whizard_lha_init (whizard_lha)
    class(whizard_lha_t), intent(out) :: whizard_lha
    whizard_lha%cptr = new_whizard_lha ()
  end subroutine whizard_lha_init

  subroutine whizard_lha_final (whizard_lha)
    class(whizard_lha_t), intent(inout)  :: whizard_lha
    call lhaup_whizard_delete (whizard_lha%cptr)
  end subroutine whizard_lha_final

  function whizard_lha_get_ptr (whizard_lha) result (cptr)
    class(whizard_lha_t), intent(in) :: whizard_lha
    type(c_ptr) :: cptr
    cptr = whizard_lha%cptr
  end function whizard_lha_get_ptr

  subroutine whizard_lha_set_init (whizard_lha, beam_pdg, beam_energy, n_processes, unweighted, negative_weights)
    class(whizard_lha_t), intent(inout) :: whizard_lha
    integer, dimension(2), intent(in) :: beam_pdg
    real(default), dimension(2), intent(in) :: beam_energy
    integer, intent(in) :: n_processes
    logical, intent(in) :: unweighted
    logical, intent(in) :: negative_weights
    logical(c_bool) :: flag
    integer(c_int) :: c_n_processes
    integer(c_int), dimension(2) :: c_beam_pdg
    real(c_double), dimension(2) :: c_beam_energy
    logical(c_bool) :: c_unweighted, c_negative_weights
    c_beam_pdg = int (beam_pdg, c_int)
    c_beam_energy = real (beam_energy, c_double)
    c_n_processes = int (n_processes, c_int)
    c_unweighted = unweighted 
    c_negative_weights = negative_weights
    flag = lhaup_whizard_set_init (whizard_lha%cptr, c_beam_pdg, &
         c_beam_energy, c_n_processes, c_unweighted, c_negative_weights)
    if (.not. flag)  then
       call msg_fatal ("[whizard_lha_set_init] could not " // & 
            "initialize the LHAUpWhizard interface.")
    end if
  end subroutine whizard_lha_set_init

  ! get this directly from event_sample_data_t
  subroutine whizard_lha_set_process_parameters (whizard_lha, process_id, cross_section, error, max_weight)
    class(whizard_lha_t), intent(inout) :: whizard_lha
    integer, intent(in) :: process_id
    real(default), intent(in), optional :: cross_section, error, max_weight
    real(default), parameter :: pb_per_fb = 1.e-3_default
    integer(c_int) :: c_process_id
    real(c_double) :: c_cross_section, c_error, c_max_weight
    logical(c_bool) :: flag
    c_process_id = int (process_id, c_int)
    if (present (cross_section)) then
       c_cross_section = real (cross_section * pb_per_fb, c_double)
    else
       c_cross_section = 0._c_double
    end if
    if (present (error)) then
       c_error = real (error * pb_per_fb, c_double)
    else
       c_error = 0._c_double
    end if
    if (present (max_weight)) then
       c_max_weight = real (max_weight, c_double)
    else
       c_max_weight = 0._c_double
    end if
    flag = lhaup_whizard_set_process_parameters (whizard_lha%cptr, &
         c_process_id, c_cross_section, c_error, c_max_weight)
    if (.not. flag) then
       call msg_fatal ("[whizard_lha_add_process] could not add a process.")
    end if
  end subroutine whizard_lha_set_process_parameters

  subroutine whizard_lha_list_init (whizard_lha)
    class(whizard_lha_t), intent(in) :: whizard_lha
    call lhaup_whizard_list_init (whizard_lha%cptr)
  end subroutine whizard_lha_list_init

  subroutine whizard_lha_list_event (whizard_lha)
    class(whizard_lha_t), intent(in) :: whizard_lha
    call lhaup_whizard_list_event (whizard_lha%cptr)
  end subroutine whizard_lha_list_event

  subroutine whizard_lha_set_event_process &
       (whizard_lha, process_id, scale, alpha_qcd, alpha_qed, weight)
    class(whizard_lha_t), intent(inout) :: whizard_lha
    integer, intent(in) :: process_id
    real(default), intent(in) :: scale, alpha_qcd, alpha_qed, weight
    integer(c_int) :: c_process_id
    real(c_double) :: c_scale, c_alpha_qcd, c_alpha_qed, c_weight
    c_scale = real (scale, c_double)
    c_alpha_qcd = real (alpha_qcd, c_double)
    c_alpha_qed = real (alpha_qed, c_double)
    c_weight = real (weight, c_double)
    c_process_id = int (process_id, c_int)
    call lhaup_whizard_set_event_process (whizard_lha%cptr, &
         c_process_id, c_scale, c_alpha_qcd, c_alpha_qed, c_weight)
    whizard_lha%new_event = .true.
  end subroutine whizard_lha_set_event_process

  subroutine whizard_lha_set_event (whizard_lha, process_id, particle_set,&
       keep_beams, keep_remnants, polarization)
    class(whizard_lha_t), intent(inout) :: whizard_lha
    integer, intent(in) :: process_id
    type(particle_set_t), intent(in) :: particle_set
    logical, intent(in), optional :: keep_beams, keep_remnants, polarization
    type(particle_set_t) :: pset
    logical :: kr, pol
    type(lha_particle_t), dimension(:), allocatable :: c_particle_set
    integer(c_int) :: c_process_id, c_n_particles
    logical(c_bool) :: flag
    kr = .true.; if (present (keep_remnants))  kr = keep_remnants
    pol = .true.; if (present (polarization))  pol = polarization
    if (.not. whizard_lha%new_event) then
       call msg_bug ("[whizard_lha_set_event] new event was not prepared.")
    end if
    call particle_set%filter_particles (pset, real_parents = .true., &
         keep_beams = keep_beams, keep_virtuals = .false.)
    if  (debug_active (D_SHOWER) .or. debug_active(D_TRANSFORMS)) then
       print *, "After particle_set%filter: pset"
       call pset%write (summary = .true., compressed = .true.)
    end if
    allocate (c_particle_set (pset%get_n_tot ()))
    call fill_c_particle_set (pset, c_particle_set, kr, pol)
    if (debug_active (D_SHOWER) .or. debug_active (D_TRANSFORMS)) &
         call lha_particle_write (c_particle_set)
    c_n_particles = pset%get_n_tot (); c_process_id = process_id
    flag = lhaup_whizard_set_event (whizard_lha%cptr, c_process_id, c_n_particles, c_particle_set)
    whizard_lha%new_event = .false.
  contains
    subroutine fill_c_particle_set (particle_set, c_particle_set, keep_remnants, polarization)
      type(particle_set_t), intent(in) :: particle_set
      type(lha_particle_t), dimension(:), intent(out) :: c_particle_set
      logical, intent(in) :: keep_remnants, polarization
      integer :: i, status
      integer, dimension(:), allocatable :: parent
      integer, dimension(2) :: color
      type(vector4_t) :: p
      do i = 1, particle_set%get_n_tot  ()
         associate (c_prt => c_particle_set, prt => particle_set%prt(i))
           c_prt(i)%id = prt%get_pdg ()
           status = prt%get_status ()
           if (keep_remnants .and. status == PRT_BEAM_REMNANT &
                .and. prt%get_n_children () == 0) then
              status = PRT_OUTGOING
           end if
           select case (status)
           case (PRT_BEAM);         c_prt(i)%status = -9
           case (PRT_INCOMING);     c_prt(i)%status = -1
           case (PRT_OUTGOING);     c_prt(i)%status =  1
           case (PRT_RESONANT);     c_prt(i)%status =  2
           case (PRT_VIRTUAL);      c_prt(i)%status =  3
           case default;            c_prt(i)%status =  0
           end select
           parent = prt%get_parents ()
           select case (size (parent))
           case (0)
              c_prt(i)%mother(1) = 0; c_prt(i)%mother(2) = 0
           case (1)
              c_prt(i)%mother(1) = parent(1); c_prt(i)%mother(2) = 0
           case (2)
              c_prt(i)%mother(1) = parent(1); c_prt(i)%mother(2) = parent(2)
           case default
              call msg_bug("[fill_c_particle_set] Too many parents. &
                   &Please contact the WHIZARD developers.")
           end select
           color = prt%get_color ()
           where (color > 0)
              c_prt(i)%color = 500 + color
           elsewhere
              c_prt(i)%color = 0
           end where
           p = prt%get_momentum ()
           c_prt(i)%momentum = p%p
           c_prt(i)%mass = invariant_mass(p)
           c_prt(i)%tau = prt%get_lifetime ()
           c_prt(i)%spin = 9
           if (polarization) then
              select case (prt%get_polarization_status ())
              case (PRT_GENERIC_POLARIZATION)
                 if (prt%get_n_parents () == 1) then
                    parent = prt%get_parents ()
                    c_prt(i)%spin = polarization_to_spin &
                         (prt%get_momentum (), prt%get_polarization (), &
                         particle_set%prt(parent(1))%get_momentum ())
                 end if
              end select
           end if
         end associate
      end do
    end subroutine fill_c_particle_set

    real(default) function polarization_to_spin (p, pol, p_mother) result (spin)
      type(vector4_t), intent(in) :: p
      type(polarization_t), intent(in) :: pol
      type(vector4_t), intent(in) :: p_mother
      type(vector3_t) :: s3, p3
      type(vector4_t) :: s4
      ! TODO sbrass move the conversion of polarization to spin to a better place (with documentation)
      s3 = vector3_moving (pol%get_axis ())
      p3 = space_part (p)
      s4 = rotation_to_2nd (3, p3) * vector4_moving (0._default, s3)
      spin = enclosed_angle_ct (s4, p_mother)
    end function polarization_to_spin

  end subroutine whizard_lha_set_event


end module whizard_lha
