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

submodule (beams) beams_s

  use io_units
  use format_defs, only: FMT_19
  use numeric_utils
  use diagnostics
  use md5

  implicit none

contains

  subroutine beam_data_init (beam_data, n)
    type(beam_data_t), intent(out) :: beam_data
    integer, intent(in) :: n
    beam_data%n = n
    allocate (beam_data%flv (n))
    allocate (beam_data%mass (n))
    allocate (beam_data%pmatrix (n))
    allocate (beam_data%p_cm (n))
    allocate (beam_data%p (n))
    beam_data%initialized = .true.
  end subroutine beam_data_init

  module subroutine beam_data_final (beam_data)
    class(beam_data_t), intent(inout) :: beam_data
    beam_data%initialized = .false.
  end subroutine beam_data_final

  module subroutine beam_data_write (beam_data, unit, verbose, write_md5sum)
    class(beam_data_t), intent(in) :: beam_data
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose, write_md5sum
    integer :: prt_name_len
    logical :: verb, write_md5
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    verb = .false.;  if (present (verbose))  verb = verbose
    write_md5 = verb;  if (present (write_md5sum)) write_md5 = write_md5sum
    if (.not. beam_data%initialized) then
       write (u, "(1x,A)") "Beam data: [undefined]"
       return
    end if
    prt_name_len = maxval (len (beam_data%flv%get_name ()))
    select case (beam_data%n)
    case (1)
       write (u, "(1x,A)") "Beam data (decay):"
       if (verb) then
          call write_prt (1)
          call beam_data%pmatrix(1)%write (u)
          write (u, *) "R.f. momentum:"
          call vector4_write (beam_data%p_cm(1), u)
          write (u, *) "Lab momentum:"
          call vector4_write (beam_data%p(1), u)
       else
          call write_prt (1)
       end if
    case (2)
       write (u, "(1x,A)") "Beam data (collision):"
       if (verb) then
          call write_prt (1)
          call beam_data%pmatrix(1)%write (u)
          call write_prt (2)
          call beam_data%pmatrix(2)%write (u)
          call write_sqrts
          write (u, *) "C.m. momenta:"
          call vector4_write (beam_data%p_cm(1), u)
          call vector4_write (beam_data%p_cm(2), u)
          write (u, *) "Lab momenta:"
          call vector4_write (beam_data%p(1), u)
          call vector4_write (beam_data%p(2), u)
       else
          call write_prt (1)
          call write_prt (2)
          call write_sqrts
       end if
    end select
    if (allocated (beam_data%L_cm_to_lab)) then
       if (verb) then
          call lorentz_transformation_write (beam_data%L_cm_to_lab, u)
       else
          write (u, "(1x,A)")  "Beam structure: lab and c.m. frame differ"
       end if
    end if
    if (write_md5) then
       write (u, *) "MD5 sum: ", beam_data%md5sum
    end if
  contains
    subroutine write_sqrts
      character(80) :: sqrts_str
      write (sqrts_str, "(" // FMT_19 // ")")  beam_data%sqrts
      write (u, "(3x,A)")  "sqrts = " // trim (adjustl (sqrts_str)) // " GeV"
    end subroutine write_sqrts
    subroutine write_prt (i)
      integer, intent(in) :: i
      character(80) :: name_str, mass_str
      write (name_str, "(A)")  char (beam_data%flv(i)%get_name ())
      write (mass_str, "(ES13.7)")  beam_data%mass(i)
      write (u, "(3x,A)", advance="no") &
           name_str(:prt_name_len) // "  (mass = " &
           // trim (adjustl (mass_str)) // " GeV)"
      if (beam_data%pmatrix(i)%is_polarized ()) then
         write (u, "(2x,A)")  "polarized"
      else
         write (u, *)
      end if
    end subroutine write_prt
  end subroutine beam_data_write

  module function beam_data_are_valid (beam_data) result (flag)
    class(beam_data_t), intent(in) :: beam_data
    logical :: flag
    flag = beam_data%initialized
  end function beam_data_are_valid

  module subroutine beam_data_check_scattering (beam_data, sqrts)
    class(beam_data_t), intent(in) :: beam_data
    real(default), intent(in), optional :: sqrts
    if (beam_data_are_valid (beam_data)) then
       if (present (sqrts)) then
          if (.not. nearly_equal (sqrts, beam_data%sqrts)) then
             call msg_error ("Current setting of sqrts is inconsistent " &
                  // "with beam setup (ignored).")
          end if
       end if
    else
       call msg_bug ("Beam setup: invalid beam data")
    end if
  end subroutine beam_data_check_scattering

  module function beam_data_get_n_in (beam_data) result (n_in)
    class(beam_data_t), intent(in) :: beam_data
    integer :: n_in
    n_in = beam_data%n
  end function beam_data_get_n_in

  module function beam_data_get_flavor (beam_data) result (flv)
    class(beam_data_t), intent(in) :: beam_data
    type(flavor_t), dimension(:), allocatable :: flv
    allocate (flv (beam_data%n))
    flv = beam_data%flv
  end function beam_data_get_flavor

  module function beam_data_get_energy (beam_data) result (e)
    class(beam_data_t), intent(in) :: beam_data
    real(default), dimension(:), allocatable :: e
    integer :: i
    allocate (e (beam_data%n))
    if (beam_data%initialized) then
       do i = 1, beam_data%n
          e(i) = energy (beam_data%p(i))
       end do
    else
       e = 0
    end if
  end function beam_data_get_energy

  module function beam_data_get_sqrts (beam_data) result (sqrts)
    class(beam_data_t), intent(in) :: beam_data
    real(default) :: sqrts
    sqrts = beam_data%sqrts
  end function beam_data_get_sqrts

  module function beam_data_get_polarization (beam_data) result (pol)
    class(beam_data_t), intent(in) :: beam_data
    real(default), dimension(beam_data%n) :: pol
    pol = beam_data%pmatrix%get_simple_pol ()
  end function beam_data_get_polarization

  module function beam_data_get_helicity_state_matrix &
       (beam_data) result (state_hel)
    type(state_matrix_t) :: state_hel
    class(beam_data_t), intent(in) :: beam_data
    type(polarization_t), dimension(:), allocatable :: pol
    integer :: i
    allocate (pol (beam_data%n))
    do i = 1, beam_data%n
       call pol(i)%init_pmatrix (beam_data%pmatrix(i))
    end do
    call combine_polarization_states (pol, state_hel)
  end function beam_data_get_helicity_state_matrix

  module function beam_data_is_initialized (beam_data) result (initialized)
    logical :: initialized
    class(beam_data_t), intent(in) :: beam_data
    initialized = any (beam_data%pmatrix%exists ())
  end function beam_data_is_initialized

  module function beam_data_get_md5sum &
       (beam_data, sqrts) result (md5sum_beams)
    class(beam_data_t), intent(in) :: beam_data
    real(default), intent(in) :: sqrts
    character(32) :: md5sum_beams
    character(80) :: buffer
    if (beam_data%md5sum /= "") then
       md5sum_beams = beam_data%md5sum
    else
       write (buffer, *)  sqrts
       md5sum_beams = md5sum (buffer)
    end if
  end function beam_data_get_md5sum

  module subroutine beam_data_init_structure &
       (beam_data, structure, sqrts, model, decay_rest_frame)
    class(beam_data_t), intent(out) :: beam_data
    type(beam_structure_t), intent(in) :: structure
    integer :: n_beam
    real(default), intent(in) :: sqrts
    class(model_data_t), intent(in), target :: model
    logical, intent(in), optional :: decay_rest_frame
    type(flavor_t), dimension(:), allocatable :: flv
    n_beam = structure%get_n_beam ()
    allocate (flv (n_beam))
    call flv%init (structure%get_prt (), model)
    if (structure%asymmetric ()) then
       if (structure%polarized ()) then
          call beam_data%init_momenta (structure%get_momenta (), flv, &
               structure%get_smatrix (), structure%get_pol_f ())
       else
          call beam_data%init_momenta (structure%get_momenta (), flv)
       end if
    else
       select case (n_beam)
       case (1)
          if (structure%polarized ()) then
             call beam_data%init_decay (flv, &
                  structure%get_smatrix (), structure%get_pol_f (), &
                  rest_frame = decay_rest_frame)
          else
             call beam_data%init_decay (flv, &
                  rest_frame = decay_rest_frame)
          end if
       case (2)
          if (structure%polarized ()) then
             call beam_data%init_sqrts (sqrts, flv, &
                  structure%get_smatrix (), structure%get_pol_f ())
          else
             call beam_data%init_sqrts (sqrts, flv)
          end if
       case default
          call msg_bug ("Beam data: invalid beam structure object")
       end select
    end if
  end subroutine beam_data_init_structure

  module subroutine beam_data_init_sqrts &
       (beam_data, sqrts, flv, smatrix, pol_f)
    class(beam_data_t), intent(out) :: beam_data
    real(default), intent(in) :: sqrts
    type(flavor_t), dimension(:), intent(in) :: flv
    type(smatrix_t), dimension(:), intent(in), optional :: smatrix
    real(default), dimension(:), intent(in), optional :: pol_f
    real(default), dimension(size(flv)) :: E, p
    call beam_data_init (beam_data, size (flv))
    beam_data%sqrts = sqrts
    beam_data%lab_is_cm = .true.
    select case (beam_data%n)
    case (1)
       E = sqrts;  p = 0
       beam_data%p_cm = vector4_moving (E, p, 3)
       beam_data%p = beam_data%p_cm
    case (2)
       beam_data%p_cm = colliding_momenta (sqrts, flv%get_mass ())
       beam_data%p = colliding_momenta (sqrts, flv%get_mass ())
    end select
    call beam_data_finish_initialization (beam_data, flv, smatrix, pol_f)
  end subroutine beam_data_init_sqrts

  module subroutine beam_data_init_momenta &
       (beam_data, p3, flv, smatrix, pol_f)
    class(beam_data_t), intent(out) :: beam_data
    type(vector3_t), dimension(:), intent(in) :: p3
    type(flavor_t), dimension(:), intent(in) :: flv
    type(smatrix_t), dimension(:), intent(in), optional :: smatrix
    real(default), dimension(:), intent(in), optional :: pol_f
    type(vector4_t) :: p0
    type(vector4_t), dimension(:), allocatable :: p, p_cm_rot
    real(default), dimension(size(p3)) :: e
    real(default), dimension(size(flv)) :: m
    type(lorentz_transformation_t) :: L_boost, L_rot
    call beam_data_init (beam_data, size (flv))
    m = flv%get_mass ()
    e = sqrt (p3 ** 2 + m ** 2)
    allocate (p (beam_data%n))
    p = vector4_moving (e, p3)
    p0 = sum (p)
    beam_data%p = p
    beam_data%lab_is_cm = .false.
    beam_data%sqrts = p0 ** 1
    L_boost = boost (p0, beam_data%sqrts)
    allocate (p_cm_rot (beam_data%n))
    p_cm_rot = inverse (L_boost) * p
    allocate (beam_data%L_cm_to_lab)
    select case (beam_data%n)
    case (1)
       beam_data%L_cm_to_lab = L_boost
       beam_data%p_cm = vector4_at_rest (beam_data%sqrts)
    case (2)
       L_rot = rotation_to_2nd (3, space_part (p_cm_rot(1)))
       beam_data%L_cm_to_lab = L_boost * L_rot
       beam_data%p_cm = &
            colliding_momenta (beam_data%sqrts, flv%get_mass ())
    end select
    call beam_data_finish_initialization (beam_data, flv, smatrix, pol_f)
  end subroutine beam_data_init_momenta

  subroutine beam_data_finish_initialization (beam_data, flv, smatrix, pol_f)
    type(beam_data_t), intent(inout) :: beam_data
    type(flavor_t), dimension(:), intent(in) :: flv
    type(smatrix_t), dimension(:), intent(in), optional :: smatrix
    real(default), dimension(:), intent(in), optional :: pol_f
    integer :: i
    do i = 1, beam_data%n
       beam_data%flv(i) = flv(i)
       beam_data%mass(i) = flv(i)%get_mass ()
       if (present (smatrix)) then
          if (size (smatrix) /= beam_data%n) &
               call msg_fatal ("Beam data: &
               &polarization density array has wrong dimension")
          beam_data%pmatrix(i) = smatrix(i)
          if (present (pol_f)) then
             if (size (pol_f) /= size (smatrix)) &
                  call msg_fatal ("Beam data: &
                  &polarization fraction array has wrong dimension")
             call beam_data%pmatrix(i)%normalize (flv(i), pol_f(i))
          else
             call beam_data%pmatrix(i)%normalize (flv(i), 1._default)
          end if
       else
          call beam_data%pmatrix(i)%init (2, 0)
          call beam_data%pmatrix(i)%normalize (flv(i), 0._default)
       end if
    end do
    call beam_data%compute_md5sum ()
  end subroutine beam_data_finish_initialization

  module subroutine beam_data_compute_md5sum (beam_data)
    class(beam_data_t), intent(inout) :: beam_data
    integer :: unit
    unit = free_unit ()
    open (unit = unit, status = "scratch", action = "readwrite")
    call beam_data%write (unit, write_md5sum = .false., &
       verbose = .true.)
    rewind (unit)
    beam_data%md5sum = md5sum (unit)
    close (unit)
  end subroutine beam_data_compute_md5sum

  module subroutine beam_data_init_decay &
       (beam_data, flv, smatrix, pol_f, rest_frame)
    class(beam_data_t), intent(out) :: beam_data
    type(flavor_t), dimension(1), intent(in) :: flv
    type(smatrix_t), dimension(1), intent(in), optional :: smatrix
    real(default), dimension(:), intent(in), optional :: pol_f
    logical, intent(in), optional :: rest_frame
    real(default), dimension(1) :: m
    m = flv%get_mass ()
    if (present (smatrix)) then
       call beam_data%init_sqrts (m(1), flv, smatrix, pol_f)
    else
       call beam_data%init_sqrts (m(1), flv, smatrix, pol_f)
    end if
    if (present (rest_frame))  beam_data%lab_is_cm = rest_frame
  end subroutine beam_data_init_decay

  module subroutine beam_init (beam, beam_data)
    type(beam_t), intent(out) :: beam
    type(beam_data_t), intent(in), target :: beam_data
    logical, dimension(beam_data%n) :: polarized, diagonal
    type(quantum_numbers_mask_t), dimension(beam_data%n) :: mask, mask_d
    type(state_matrix_t), target :: state_hel, state_fc, state_tmp
    type(state_iterator_t) :: it_hel, it_tmp
    type(quantum_numbers_t), dimension(:), allocatable :: qn
    complex(default) :: value
    real(default), parameter :: tolerance = 100 * epsilon (1._default)
    polarized = beam_data%pmatrix%is_polarized ()
    diagonal = beam_data%pmatrix%is_diagonal ()
    mask = quantum_numbers_mask (.false., .false., &
         mask_h = .not. polarized, &
         mask_hd = diagonal)
    mask_d = quantum_numbers_mask (.false., .false., .false., &
         mask_hd = polarized .and. diagonal)
    call beam%int%basic_init &
         (0, 0, beam_data%n, mask = mask, store_values = .true.)
    state_hel = beam_data%get_helicity_state_matrix ()
    allocate (qn (beam_data%n))
    call qn%init (beam_data%flv, color_from_flavor (beam_data%flv, 1))
    call state_fc%init ()
    call state_fc%add_state (qn)
    call merge_state_matrices (state_hel, state_fc, state_tmp)
    call it_hel%init (state_hel)
    call it_tmp%init (state_tmp)
    do while (it_hel%is_valid ())
       qn = it_tmp%get_quantum_numbers ()
       value = it_hel%get_matrix_element ()
       if (any (qn%are_redundant (mask_d))) then
          ! skip off-diagonal elements for diagonal polarization
       else if (abs (value) <= tolerance) then
          ! skip zero entries
       else
          call beam%int%add_state (qn, value = value)
       end if
       call it_hel%advance ()
       call it_tmp%advance ()
    end do
    call beam%int%freeze ()
    call beam%int%set_momenta (beam_data%p, outgoing = .true.)
    call state_hel%final ()
    call state_fc%final ()
    call state_tmp%final ()
  end subroutine beam_init

  module subroutine beam_final (beam)
    type(beam_t), intent(inout) :: beam
    call beam%int%final ()
  end subroutine beam_final

  module subroutine beam_write &
       (beam, unit, verbose, show_momentum_sum, show_mass, col_verbose)
    type(beam_t), intent(in) :: beam
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose, show_momentum_sum, show_mass
    logical, intent(in), optional :: col_verbose
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    select case (beam%int%get_n_out ())
    case (1);  write (u, *) "Decaying particle:"
    case (2);  write (u, *) "Colliding beams:"
    end select
    call beam%int%basic_write &
         (unit, verbose = verbose, show_momentum_sum = &
          show_momentum_sum, show_mass = show_mass, &
          col_verbose = col_verbose)
  end subroutine beam_write

  module subroutine beam_assign (beam_out, beam_in)
    type(beam_t), intent(out) :: beam_out
    type(beam_t), intent(in) :: beam_in
    beam_out%int = beam_in%int
  end subroutine beam_assign

  module subroutine interaction_set_source_link_beam (int, i, beam1, i1)
    type(interaction_t), intent(inout) :: int
    type(beam_t), intent(in), target :: beam1
    integer, intent(in) :: i, i1
    call int%set_source_link (i, beam1%int, i1)
  end subroutine interaction_set_source_link_beam

  module function beam_get_int_ptr (beam) result (int)
    type(interaction_t), pointer :: int
    type(beam_t), intent(in), target :: beam
    int => beam%int
  end function beam_get_int_ptr

  module subroutine beam_set_momenta (beam, p)
    type(beam_t), intent(inout) :: beam
    type(vector4_t), dimension(:), intent(in) :: p
    call beam%int%set_momenta (p)
  end subroutine beam_set_momenta


end submodule beams_s

