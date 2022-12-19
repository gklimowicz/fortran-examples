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

submodule (isr_epa_handler) isr_epa_handler_s

  use diagnostics, only: msg_fatal
  use diagnostics, only: msg_bug
  use io_units
  use format_defs, only: FMT_12, FMT_19
  use format_utils, only: write_separator
  use format_utils, only: pac_fmt
  use physics_defs, only: PHOTON
  use recoil_kinematics, only: initial_transformation
  use recoil_kinematics, only: generate_recoil
  use recoil_kinematics, only: recoil_transformation

  implicit none

contains

  function rad_mode_string (mode) result (string)
    type(string_t) :: string
    integer, intent(in) :: mode
    select case (mode)
    case (BEAM_RAD_NONE);  string = "---"
    case (BEAM_RAD_ISR);   string = "ISR"
    case (BEAM_RAD_EPA);   string = "EPA"
    case default;  string = "???"
    end select
  end function rad_mode_string

  module function evt_isr_epa_get_mode_string (evt) result (string)
    type(string_t) :: string
    class(evt_isr_epa_t), intent(in) :: evt
    select case (evt%mode)
    case (ISR_TRIVIAL_COLLINEAR)
       string = "trivial, collinear"
    case (ISR_PAIR_RECOIL)
       string = "pair recoil"
    case default
       string = "[undefined]"
    end select
  end function evt_isr_epa_get_mode_string

  module subroutine evt_isr_epa_set_mode_string (evt, string)
    class(evt_isr_epa_t), intent(inout) :: evt
    type(string_t), intent(in) :: string
    select case (char (string))
    case ("trivial")
       evt%mode = ISR_TRIVIAL_COLLINEAR
    case ("recoil")
       evt%mode = ISR_PAIR_RECOIL
    case default
       call msg_fatal ("ISR handler: mode '" // char (string) &
            // "' is undefined")
    end select
  end subroutine evt_isr_epa_set_mode_string

  module subroutine evt_isr_epa_write_name (evt, unit)
    class(evt_isr_epa_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Event transform: ISR/EPA handler"
  end subroutine evt_isr_epa_write_name

  module subroutine evt_isr_epa_write_mode (evt, unit)
    class(evt_isr_epa_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A,1x,I0,':',1x,A)")  "Insertion mode =", evt%mode, &
         char (evt%get_mode_string ())
  end subroutine evt_isr_epa_write_mode

  module subroutine evt_isr_epa_write_input (evt, unit, testflag)
    class(evt_isr_epa_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    character(len=7) :: fmt
    integer :: u
    u = given_output_unit (unit)
    call pac_fmt (fmt, FMT_19, FMT_12, testflag)
    if (evt%isr_active) then
       write (u, "(3x,A,1x," // fmt // ")") "ISR: Q_max =", evt%isr_q_max
       write (u, "(3x,A,1x," // fmt // ")") "     m     =", evt%isr_mass
       write (u, "(3x,A,1x,L1)") "     keep m=", evt%isr_keep_mass
    else
       write (u, "(3x,A)") "ISR: [inactive]"
    end if
    if (evt%epa_active) then
       write (u, "(3x,A,1x," // fmt // ")") "EPA: Q_max =", evt%epa_q_max
       write (u, "(3x,A,1x," // fmt // ")") "     m     =", evt%epa_mass
    else
       write (u, "(3x,A)") "EPA: [inactive]"
    end if
  end subroutine evt_isr_epa_write_input

  module subroutine evt_isr_epa_write_data (evt, unit, testflag)
    class(evt_isr_epa_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    character(len=7), parameter :: FMTL_19 = "A3,16x"
    character(len=7), parameter :: FMTL_12 = "A3,9x"
    character(len=7) :: fmt, fmtl
    integer :: u
    u = given_output_unit (unit)
    call pac_fmt (fmt, FMT_19, FMT_12, testflag)
    call pac_fmt (fmtl, FMTL_19, FMTL_12, testflag)
    select case (evt%mode)
    case (ISR_PAIR_RECOIL)
       write (u, "(1x,A)")  "Event:"
       write (u, "(3x,A,2(1x," // fmtl // "))") &
            "mode  = ", &
            char (rad_mode_string (evt%rad_mode(1))), &
            char (rad_mode_string (evt%rad_mode(2)))
       write (u, "(3x,A,2(1x," // fmt // "))") "Q_max =", evt%q_max
       write (u, "(3x,A,2(1x," // fmt // "))") "m     =", evt%m
       write (u, "(3x,A,2(1x," // fmt // "))") "x     =", evt%xc
       write (u, "(3x,A,2(1x," // fmt // "))") "xb    =", evt%xcb
       write (u, "(3x,A,1x," // fmt // ")")    "sqrts =", evt%sqrts
       call write_separator (u)
       write (u, "(A)")  "Lorentz boost (partons before radiation &
            &c.m. -> lab) ="
       call evt%lti%write (u, testflag)
       write (u, "(A)")  "Lorentz transformation (collinear partons &
            &-> partons with recoil in c.m.) ="
       call evt%lto%write (u, testflag)
       write (u, "(A)")  "Combined transformation (partons &
            &-> partons with recoil in lab frame) ="
       call evt%lt%write (u, testflag)
    end select
  end subroutine evt_isr_epa_write_data

  module subroutine evt_isr_epa_write &
       (evt, unit, verbose, more_verbose, testflag)
    class(evt_isr_epa_t), intent(in) :: evt
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose, more_verbose, testflag
    logical :: show_mass
    integer :: u, i
    u = given_output_unit (unit)
    if (present (testflag)) then
       show_mass = .not. testflag
    else
       show_mass = .true.
    end if
    call write_separator (u, 2)
    call evt%write_name (u)
    call write_separator (u, 2)
    call evt%write_mode (u)
    call evt%write_input (u, testflag=testflag)
    call evt%write_data (u, testflag=testflag)
    call write_separator (u)
    call evt%base_write (u, testflag = testflag, show_set = .false.)
    if (all (evt%i_beam > 0)) then
       call write_separator (u)
       write (u, "(A,2(1x,I0))")  "Partons before radiation:", evt%i_beam
       do i = 1, 2
          call evt%beam(i)%write (u, testflag=testflag)
       end do
       call write_separator (u)
       write (u, "(A)")  "... boosted to c.m.:"
       do i = 1, 2
          call evt%pi(i)%write (u, show_mass=show_mass, testflag=testflag)
       end do
    end if
    if (all (evt%i_radiated > 0)) then
       call write_separator (u)
       write (u, "(A,2(1x,I0))")  "Radiated particles, collinear:", &
            evt%i_radiated
       do i = 1, 2
          call evt%radiated(i)%write (u, testflag=testflag)
       end do
       call write_separator (u)
       write (u, "(A)")  "... boosted to c.m.:"
       do i = 1, 2
          call evt%ki(i)%write (u, show_mass=show_mass, testflag=testflag)
       end do
       call write_separator (u)
       write (u, "(A)")  "... with kT:"
       do i = 1, 2
          call evt%km(i)%write (u, show_mass=show_mass, testflag=testflag)
       end do
    end if
    if (all (evt%i_parton > 0)) then
       call write_separator (u)
       write (u, "(A,2(1x,I0))")  "Partons after radiation, collinear:", &
            evt%i_parton
       do i = 1, 2
          call evt%parton(i)%write (u, testflag=testflag)
       end do
       call write_separator (u)
       write (u, "(A)")  "... boosted to c.m.:"
       do i = 1, 2
          call evt%qi(i)%write (u, show_mass=show_mass, testflag=testflag)
       end do
       call write_separator (u)
       write (u, "(A)")  "... with qT, off-shell:"
       do i = 1, 2
          call evt%qm(i)%write (u, show_mass=show_mass, testflag=testflag)
       end do
       call write_separator (u)
       write (u, "(A)")  "... projected on-shell:"
       do i = 1, 2
          call evt%qo(i)%write (u, show_mass=show_mass, testflag=testflag)
       end do
       call write_separator (u)
    end if
    if (evt%particle_set_exists)  &
         call evt%particle_set%write &
         (u, summary = .true., compressed = .true., testflag = testflag)
    call write_separator (u)
  end subroutine evt_isr_epa_write

  module subroutine evt_isr_epa_import_rng (evt, rng)
    class(evt_isr_epa_t), intent(inout) :: evt
    class(rng_t), allocatable, intent(inout) :: rng
    call move_alloc (from = rng, to = evt%rng)
  end subroutine evt_isr_epa_import_rng

  module subroutine evt_isr_epa_set_data_isr (evt, sqrts, q_max, m, keep_mass)
    class(evt_isr_epa_t), intent(inout) :: evt
    real(default), intent(in) :: sqrts
    real(default), intent(in) :: q_max
    real(default), intent(in) :: m
    logical, intent(in) :: keep_mass
    if (sqrts <= 0) then
       call msg_fatal ("ISR handler: sqrts value must be positive")
    end if
    if (q_max <= 0 .or. q_max > sqrts) then
       evt%isr_q_max = sqrts
    else
       evt%isr_q_max = q_max
    end if
    if (m > 0) then
       evt%isr_mass = m
    else
       call msg_fatal ("ISR handler: ISR_mass value must be positive")
    end if
    evt%isr_active = .true.
    evt%isr_keep_mass = keep_mass
  end subroutine evt_isr_epa_set_data_isr

  module subroutine evt_isr_epa_set_data_epa (evt, sqrts, q_max, m)
    class(evt_isr_epa_t), intent(inout) :: evt
    real(default), intent(in) :: sqrts
    real(default), intent(in) :: q_max
    real(default), intent(in) :: m
    if (sqrts <= 0) then
       call msg_fatal ("EPA handler: sqrts value must be positive")
    end if
    if (q_max <= 0 .or. q_max > sqrts) then
       evt%epa_q_max = sqrts
    else
       evt%epa_q_max = q_max
    end if
    if (m > 0) then
       evt%epa_mass = m
    else
       call msg_fatal ("EPA handler: EPA_mass value must be positive")
    end if
    evt%epa_active = .true.
  end subroutine evt_isr_epa_set_data_epa

  module subroutine identify_radiated (evt)
    class(evt_isr_epa_t), intent(inout) :: evt
    integer :: i, k
    k = 2
    FIND_LAST_RADIATED: do i = evt%particle_set%get_n_tot (), 1, -1
       associate (prt => evt%particle_set%prt(i))
         if (prt%is_beam_remnant ()) then
            evt%i_radiated(k) = i
            evt%radiated(k) = prt
            k = k - 1
            if (k == 0)  exit FIND_LAST_RADIATED
         end if
       end associate
    end do FIND_LAST_RADIATED
    if (k /= 0)  call err_count
  contains
    subroutine err_count
      call evt%particle_set%write ()
      call msg_fatal ("ISR/EPA handler: &
           &event does not contain two radiated particles")
    end subroutine err_count
  end subroutine identify_radiated

  module subroutine identify_partons (evt)
    class(evt_isr_epa_t), intent(inout) :: evt
    integer, dimension(:), allocatable :: parent, child
    integer :: i, j
    if (all (evt%i_radiated > 0)) then
       do i = 1, 2
          parent = evt%radiated(i)%get_parents ()
          if (size (parent) /= 1)  call err_mismatch
          evt%i_beam(i) = parent(1)
          evt%beam(i) = evt%particle_set%prt(parent(1))
          associate (prt => evt%beam(i))
            child = prt%get_children ()
            if (size (child) /= 2)  call err_mismatch
            do j = 1, 2
               if (child(j) /= evt%i_radiated(i)) then
                  evt%i_parton(i) = child(j)
                  evt%parton(i) = evt%particle_set%prt(child(j))
               end if
            end do
          end associate
       end do
    end if
  contains
    subroutine err_mismatch
      call evt%particle_set%write ()
      call msg_bug ("ISR/EPA handler: mismatch in parent-child relations")
    end subroutine err_mismatch
  end subroutine identify_partons

  module subroutine evt_isr_epa_check_radiation (evt)
    class(evt_isr_epa_t), intent(inout) :: evt
    type(flavor_t) :: flv
    integer :: i
    do i = 1, 2
       flv = evt%radiated(i)%get_flv ()
       if (flv%get_pdg () == PHOTON) then
          if (evt%isr_active) then
             evt%rad_mode(i) = BEAM_RAD_ISR
          else
             call err_isr_init
          end if
       else
          flv = evt%parton(i)%get_flv ()
          if (flv%get_pdg () == PHOTON) then
             if (evt%epa_active) then
                evt%rad_mode(i) = BEAM_RAD_EPA
             else
                call err_epa_init
             end if
          else
             call err_no_photon
          end if
       end if
    end do
  contains
    subroutine err_isr_init
      call evt%particle_set%write ()
      call msg_fatal ("ISR/EPA handler: &
           &event contains radiated photon, but ISR is not initialized")
    end subroutine err_isr_init
    subroutine err_epa_init
      call evt%particle_set%write ()
      call msg_fatal ("ISR/EPA handler: &
           &event contains incoming photon, but EPA is not initialized")
    end subroutine err_epa_init
    subroutine err_no_photon
      call evt%particle_set%write ()
      call msg_fatal ("ISR/EPA handler: &
           &event does not appear to be ISR or EPA - missing photon")
    end subroutine err_no_photon
  end subroutine evt_isr_epa_check_radiation

  module subroutine evt_isr_epa_set_recoil_parameters (evt)
    class(evt_isr_epa_t), intent(inout) :: evt
    integer :: i
    do i = 1, 2
       select case (evt%rad_mode(i))
       case (BEAM_RAD_ISR)
          evt%q_max(i) = evt%isr_q_max
          evt%m(i) = evt%isr_mass
       case (BEAM_RAD_EPA)
          evt%q_max(i) = evt%epa_q_max
          evt%m(i) = evt%epa_mass
       end select
    end do
  end subroutine evt_isr_epa_set_recoil_parameters

  module subroutine boost_to_cm (evt)
    class(evt_isr_epa_t), intent(inout) :: evt
    type(vector4_t), dimension(2) :: p
    type(vector4_t), dimension(2) :: k
    type(vector4_t), dimension(2) :: q
    logical :: ok
    p = evt%beam%get_momentum ()
    k = evt%radiated%get_momentum ()
    q = evt%parton%get_momentum ()
    call initial_transformation (p, evt%sqrts, evt%lti, ok)
    if (.not. ok)  call err_non_collinear
    evt%pi = inverse (evt%lti) * p
    evt%ki = inverse (evt%lti) * k
    evt%qi = inverse (evt%lti) * q
  contains
    subroutine err_non_collinear
      call evt%particle_set%write ()
      call msg_fatal ("ISR/EPA handler: &
           &partons before radiation are not collinear")
    end subroutine err_non_collinear
  end subroutine boost_to_cm

  module subroutine infer_x (evt)
    class(evt_isr_epa_t), intent(inout) :: evt
    real(default) :: E_parent, E_radiated, E_parton
    integer :: i
    if (all (evt%i_radiated > 0)) then
       do i = 1, 2
          E_parent = energy (evt%pi(i))
          E_radiated = energy (evt%ki(i))
          E_parton = energy (evt%qi(i))
          if (E_parent > 0) then
             evt%xc(i) = E_parton / E_parent
             evt%xcb(i)= E_radiated / E_parent
          else
             call err_energy
          end if
       end do
    end if
  contains
    subroutine err_energy
      call evt%particle_set%write ()
      call msg_bug ("ISR/EPA handler: non-positive energy in splitting")
    end subroutine err_energy
  end subroutine infer_x

  module subroutine evt_generate_recoil (evt, ok)
    class(evt_isr_epa_t), intent(inout) :: evt
    logical, intent(out) :: ok
    real(default), dimension(4) :: r
    real(default), dimension(2) :: m, mo
    integer :: i
    call evt%rng%generate (r)
    m = 0
    mo = 0
    do i = 1, 2
       select case (evt%rad_mode(i))
       case (BEAM_RAD_ISR)
          m(i) = evt%m(i)
          if (evt%isr_keep_mass)  mo(i) = m(i)
       case (BEAM_RAD_EPA)
          m(i) = evt%xc(i) * evt%m(i)
       end select
    end do
    call generate_recoil (evt%sqrts, evt%q_max, m, mo, evt%xc, evt%xcb, r, &
         evt%km, evt%qm, evt%qo, ok)
  end subroutine evt_generate_recoil

  module subroutine replace_radiated (evt)
    class(evt_isr_epa_t), intent(inout) :: evt
    integer :: i
    do i = 1, 2
       associate (prt => evt%particle_set%prt(evt%i_radiated(i)))
         call prt%set_momentum (evt%lti * evt%km(i))
       end associate
    end do
  end subroutine replace_radiated

  module subroutine replace_partons (evt)
    class(evt_isr_epa_t), intent(inout) :: evt
    integer :: i
    do i = 1, 2
       associate (prt => evt%particle_set%prt(evt%i_parton(i)))
         call prt%set_momentum (evt%lti * evt%qo(i))
       end associate
    end do
  end subroutine replace_partons

  module subroutine evt_transform_outgoing (evt)
    class(evt_isr_epa_t), intent(inout) :: evt
    logical, dimension(:), allocatable :: mask
    call recoil_transformation (evt%sqrts, evt%xc, evt%qo, evt%lto)
    evt%lt = evt%lti * evt%lto * inverse (evt%lti)
    allocate (mask (evt%particle_set%get_n_tot ()), source=.false.)
    call transform_children (evt%i_parton(1))
  contains
    recursive subroutine transform_children (i)
      integer, intent(in) :: i
      integer :: j, n_child, c
      integer, dimension(:), allocatable :: child
      child = evt%particle_set%prt(i)%get_children ()
      do j = 1, size (child)
         c = child(j)
         if (.not. mask(c)) then
            associate (prt => evt%particle_set%prt(c))
              call prt%set_momentum (evt%lt * prt%get_momentum ())
              mask(c) = .true.
              call transform_children (c)
            end associate
         end if
      end do
    end subroutine transform_children
  end subroutine evt_transform_outgoing

  module subroutine evt_isr_epa_generate_weighted (evt, probability)
    class(evt_isr_epa_t), intent(inout) :: evt
    real(default), intent(inout) :: probability
    logical :: valid
    call evt%particle_set%final ()
    evt%particle_set = evt%previous%particle_set
    evt%particle_set_exists = .true.
    select case (evt%mode)
    case (ISR_TRIVIAL_COLLINEAR)
       probability = 1
       valid = .true.
    case (ISR_PAIR_RECOIL)
       call evt%identify_radiated ()
       call evt%identify_partons ()
       call evt%check_radiation ()
       call evt%set_recoil_parameters ()
       call evt%boost_to_cm ()
       call evt%infer_x ()
       call evt%generate_recoil (valid)
       if (valid) then
          probability = 1
       else
          probability = 0
       end if
    case default
       call msg_bug ("ISR/EPA handler: generate weighted: unsupported mode")
    end select
    evt%particle_set_exists = .false.
  end subroutine evt_isr_epa_generate_weighted

  module subroutine evt_isr_epa_make_particle_set &
       (evt, factorization_mode, keep_correlations, r)
    class(evt_isr_epa_t), intent(inout) :: evt
    integer, intent(in) :: factorization_mode
    logical, intent(in) :: keep_correlations
    real(default), dimension(:), intent(in), optional :: r
    select case (evt%mode)
    case (ISR_TRIVIAL_COLLINEAR)
    case (ISR_PAIR_RECOIL)
       call evt%replace_radiated ()
       call evt%replace_partons ()
       call evt%transform_outgoing ()
    case default
       call msg_bug ("ISR/EPA handler: make particle set: unsupported mode")
    end select
    evt%particle_set_exists = .true.
  end subroutine evt_isr_epa_make_particle_set

  module subroutine evt_isr_epa_prepare_new_event (evt, i_mci, i_term)
    class(evt_isr_epa_t), intent(inout) :: evt
    integer, intent(in) :: i_mci, i_term
    call evt%reset ()
  end subroutine evt_isr_epa_prepare_new_event


end submodule isr_epa_handler_s

