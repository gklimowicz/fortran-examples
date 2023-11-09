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

submodule (sf_circe1) sf_circe1_s

  use io_units
  use format_defs, only: FMT_17, FMT_19
  use diagnostics
  use physics_defs, only: ELECTRON, PHOTON
  use lorentz
  use colors
  use quantum_numbers
  use state_matrices

  implicit none

contains

  module subroutine circe1_data_init &
       (data, model, pdg_in, sqrts, eps, out_photon, &
        ver, rev, acc, chat, with_radiation)
    class(circe1_data_t), intent(out) :: data
    class(model_data_t), intent(in), target :: model
    type(pdg_array_t), dimension(2), intent(in) :: pdg_in
    real(default), intent(in) :: sqrts
    real(default), intent(in) :: eps
    logical, dimension(2), intent(in) :: out_photon
    character(*), intent(in) :: acc
    integer, intent(in) :: ver, rev, chat
    logical, intent(in) :: with_radiation
    data%model => model
    if (any (pdg_in%get_length () /= 1)) then
       call msg_fatal ("CIRCE1: incoming beam particles must be unique")
    end if
    call data%flv_in(1)%init (pdg_in(1)%get (1), model)
    call data%flv_in(2)%init (pdg_in(2)%get (1), model)
    data%pdg_in = data%flv_in%get_pdg ()
    data%m_in = data%flv_in%get_mass ()
    data%sqrts = sqrts
    data%eps = eps
    data%photon = out_photon
    data%ver = ver
    data%rev = rev
    data%acc = acc
    data%chat = chat
    data%with_radiation = with_radiation
    call data%check ()
    call circex (0.d0, 0.d0, dble (data%sqrts), &
         data%acc, data%ver, data%rev, data%chat)
  end subroutine circe1_data_init

  module subroutine circe1_data_set_generator_mode (data, rng_factory)
    class(circe1_data_t), intent(inout) :: data
    class(rng_factory_t), intent(inout), allocatable :: rng_factory
    data%generate = .true.
    call move_alloc (from = rng_factory, to = data%rng_factory)
  end subroutine circe1_data_set_generator_mode

  module subroutine circe1_data_check (data)
    class(circe1_data_t), intent(in) :: data
    type(flavor_t) :: flv_electron, flv_photon
    call flv_electron%init (ELECTRON, data%model)
    call flv_photon%init (PHOTON, data%model)
    if (.not. flv_electron%is_defined () &
         .or. .not. flv_photon%is_defined ()) then
       call msg_fatal ("CIRCE1: model must contain photon and electron")
    end if
    if (any (abs (data%pdg_in) /= ELECTRON) &
         .or. (data%pdg_in(1) /= - data%pdg_in(2))) then
       call msg_fatal ("CIRCE1: applicable only for e+e- or e-e+ collisions")
    end if
    if (data%eps <= 0) then
       call msg_error ("CIRCE1: circe1_eps = 0: integration will &
            &miss x=1 peak")
    end if
  end subroutine circe1_data_check

  module subroutine circe1_data_write (data, unit, verbose)
    class(circe1_data_t), intent(in) :: data
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u
    logical :: verb
    verb = .false.;  if (present (verbose))  verb = verbose
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(1x,A)") "CIRCE1 data:"
    write (u, "(3x,A,2(1x,A))") "prt_in   =", &
         char (data%flv_in(1)%get_name ()), &
         char (data%flv_in(2)%get_name ())
    write (u, "(3x,A,2(1x,L1))")  "photon   =", data%photon
    write (u, "(3x,A,L1)")        "generate = ", data%generate
    write (u, "(3x,A,2(1x," // FMT_19 // "))") "m_in     =", data%m_in
    write (u, "(3x,A," // FMT_19 // ")") "sqrts    = ", data%sqrts
    write (u, "(3x,A," // FMT_19 // ")") "eps      = ", data%eps
    write (u, "(3x,A,I0)") "ver      = ", data%ver
    write (u, "(3x,A,I0)") "rev      = ", data%rev
    write (u, "(3x,A,A)")  "acc      = ", data%acc
    write (u, "(3x,A,I0)") "chat     = ", data%chat
    write (u, "(3x,A,L1)") "with rad.= ", data%with_radiation
    if (data%generate) then
       if (verb) then
          call data%rng_factory%write (u)
       end if
    end if
  end subroutine circe1_data_write

  module function circe1_data_is_generator (data) result (flag)
    class(circe1_data_t), intent(in) :: data
    logical :: flag
    flag = data%generate
  end function circe1_data_is_generator

  module function circe1_data_get_n_par (data) result (n)
    class(circe1_data_t), intent(in) :: data
    integer :: n
    n = 2
  end function circe1_data_get_n_par

  module subroutine circe1_data_get_pdg_out (data, pdg_out)
    class(circe1_data_t), intent(in) :: data
    type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    integer :: i, n
    n = 2
    do i = 1, n
       if (data%photon(i)) then
          pdg_out(i) = PHOTON
       else
          pdg_out(i) = data%pdg_in(i)
       end if
    end do
  end subroutine circe1_data_get_pdg_out

  module function circe1_data_get_pdg_int (data) result (pdg)
    class(circe1_data_t), intent(in) :: data
    integer, dimension(2) :: pdg
    integer :: i
    do i = 1, 2
       if (data%photon(i)) then
          pdg(i) = PHOTON
       else
          pdg(i) = data%pdg_in(i)
       end if
    end do
  end function circe1_data_get_pdg_int

  module function circe1_data_get_beam_file (data) result (file)
    class(circe1_data_t), intent(in) :: data
    type(string_t) :: file
    file = "CIRCE1: " // data%acc
  end function circe1_data_get_beam_file

  module subroutine rng_obj_generate (rng_obj, u)
    class(rng_obj_t), intent(inout) :: rng_obj
    real(double), intent(out) :: u
    real(default) :: x
    call rng_obj%rng%generate (x)
    u = x
  end subroutine rng_obj_generate

  module function circe1_type_string (object) result (string)
    class(circe1_t), intent(in) :: object
    type(string_t) :: string
    if (associated (object%data)) then
       string = "CIRCE1: beamstrahlung"
    else
       string = "CIRCE1: [undefined]"
    end if
  end function circe1_type_string

  module subroutine circe1_write (object, unit, testflag)
    class(circe1_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    if (associated (object%data)) then
       call object%data%write (u)
       if (object%data%generate)  call object%rng_obj%rng%write (u)
       if (object%status >= SF_DONE_KINEMATICS) then
          write (u, "(3x,A,2(1x," // FMT_17 // "))")  "x =", object%x
          write (u, "(3x,A,2(1x," // FMT_17 // "))")  "xb=", object%xb
          if (object%status >= SF_FAILED_EVALUATION) then
             write (u, "(3x,A,1x," // FMT_17 // ")")  "f =", object%f
          end if
       end if
       call object%base_write (u, testflag)
    else
       write (u, "(1x,A)")  "CIRCE1 data: [undefined]"
    end if
  end subroutine circe1_write

  module subroutine circe1_init (sf_int, data)
    class(circe1_t), intent(out) :: sf_int
    class(sf_data_t), intent(in), target :: data
    logical, dimension(6) :: mask_h
    type(quantum_numbers_mask_t), dimension(6) :: mask
    integer, dimension(6) :: hel_lock
    type(polarization_t), target :: pol1, pol2
    type(quantum_numbers_t), dimension(1) :: qn_fc1, qn_fc2
    type(flavor_t) :: flv_photon
    type(color_t) :: col0
    real(default), dimension(2) :: mi2, mr2, mo2
    type(quantum_numbers_t) :: qn_hel1, qn_hel2, qn_photon, qn1, qn2
    type(quantum_numbers_t), dimension(6) :: qn
    type(polarization_iterator_t) :: it_hel1, it_hel2
    hel_lock = 0
    mask_h = .false.
    select type (data)
    type is (circe1_data_t)
       mi2 = data%m_in**2
       if (data%with_radiation) then
          if (data%photon(1)) then
             hel_lock(1) = 3;  hel_lock(3) = 1;  mask_h(5) = .true.
             mr2(1) = mi2(1)
             mo2(1) = 0._default
          else
             hel_lock(1) = 5;  hel_lock(5) = 1;  mask_h(3) = .true.
             mr2(1) = 0._default
             mo2(1) = mi2(1)
          end if
          if (data%photon(2)) then
             hel_lock(2) = 4;  hel_lock(4) = 2;  mask_h(6) = .true.
             mr2(2) = mi2(2)
             mo2(2) = 0._default
          else
             hel_lock(2) = 6;  hel_lock(6) = 2;  mask_h(4) = .true.
             mr2(2) = 0._default
             mo2(2) = mi2(2)
          end if
          mask = quantum_numbers_mask (.false., .false., mask_h)
          call sf_int%base_init (mask, mi2, mr2, mo2, &
               hel_lock = hel_lock)
          sf_int%data => data
          call flv_photon%init (PHOTON, data%model)
          call col0%init ()
          call qn_photon%init (flv_photon, col0)
          call pol1%init_generic (data%flv_in(1))
          call qn_fc1(1)%init (flv = data%flv_in(1), col = col0)
          call pol2%init_generic (data%flv_in(2))
          call qn_fc2(1)%init (flv = data%flv_in(2), col = col0)
          call it_hel1%init (pol1)

          do while (it_hel1%is_valid ())
             qn_hel1 = it_hel1%get_quantum_numbers ()
             qn1 = qn_hel1 .merge. qn_fc1(1)
             qn(1) = qn1
             if (data%photon(1)) then
                qn(3) = qn1;  qn(5) = qn_photon
             else
                qn(3) = qn_photon;  qn(5) = qn1
             end if
             call it_hel2%init (pol2)
             do while (it_hel2%is_valid ())
                qn_hel2 = it_hel2%get_quantum_numbers ()
                qn2 = qn_hel2 .merge. qn_fc2(1)
                qn(2) = qn2
                if (data%photon(2)) then
                   qn(4) = qn2;  qn(6) = qn_photon
                else
                   qn(4) = qn_photon;  qn(6) = qn2
                end if
                call qn(3:4)%tag_radiated ()
                call sf_int%add_state (qn)
                call it_hel2%advance ()
             end do
             call it_hel1%advance ()
          end do
!           call pol1%final ()
!           call pol2%final ()
          call sf_int%freeze ()
          call sf_int%set_incoming ([1,2])
          call sf_int%set_radiated ([3,4])
          call sf_int%set_outgoing ([5,6])
       else
          if (data%photon(1)) then
             mask_h(3) = .true.
             mo2(1) = 0._default
          else
             hel_lock(1) = 3;  hel_lock(3) = 1
             mo2(1) = mi2(1)
          end if
          if (data%photon(2)) then
             mask_h(4) = .true.
             mo2(2) = 0._default
          else
             hel_lock(2) = 4;  hel_lock(4) = 2
             mo2(2) = mi2(2)
          end if
          mask = quantum_numbers_mask (.false., .false., mask_h)
          call sf_int%base_init (mask(1:4), mi2, [real(default) :: ], mo2, &
               hel_lock = hel_lock(1:4))
          sf_int%data => data
          call flv_photon%init (PHOTON, data%model)
          call col0%init ()
          call qn_photon%init (flv_photon, col0)
          call pol1%init_generic (data%flv_in(1))
          call qn_fc1(1)%init (flv = data%flv_in(1), col = col0)
          call pol2%init_generic (data%flv_in(2))
          call qn_fc2(1)%init (flv = data%flv_in(2), col = col0)
          call it_hel1%init (pol1)

          do while (it_hel1%is_valid ())
             qn_hel1 = it_hel1%get_quantum_numbers ()
             qn1 = qn_hel1 .merge. qn_fc1(1)
             qn(1) = qn1
             if (data%photon(1)) then
                qn(3) = qn_photon
             else
                qn(3) = qn1
             end if
             call it_hel2%init (pol2)
             do while (it_hel2%is_valid ())
                qn_hel2 = it_hel2%get_quantum_numbers ()
                qn2 = qn_hel2 .merge. qn_fc2(1)
                qn(2) = qn2
                if (data%photon(2)) then
                   qn(4) = qn_photon
                else
                   qn(4) = qn2
                end if
                call sf_int%add_state (qn(1:4))
                call it_hel2%advance ()
             end do
             call it_hel1%advance ()
          end do
!           call pol1%final ()
!           call pol2%final ()
          call sf_int%freeze ()
          call sf_int%set_incoming ([1,2])
          call sf_int%set_outgoing ([3,4])
       end if
       sf_int%status = SF_INITIAL
    end select
    if (sf_int%data%generate) then
       call sf_int%data%rng_factory%make (sf_int%rng_obj%rng)
    end if
  end subroutine circe1_init

  module function circe1_is_generator (sf_int) result (flag)
    class(circe1_t), intent(in) :: sf_int
    logical :: flag
    flag = sf_int%data%is_generator ()
  end function circe1_is_generator

  module subroutine circe1_generate_free (sf_int, r, rb,  x_free)
    class(circe1_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(inout) :: x_free

    if (sf_int%data%generate) then
       call circe_generate (r, sf_int%data%get_pdg_int (), sf_int%rng_obj)
       rb = 1 - r
       x_free = x_free * product (r)
    else
       r = 0
       rb= 1
    end if
  end subroutine circe1_generate_free

  subroutine circe_generate (x, pdg, rng_obj)
    real(default), dimension(2), intent(out) :: x
    integer, dimension(2), intent(in) :: pdg
    class(rng_obj_t), intent(inout) :: rng_obj
    real(double) :: xc1, xc2
    select case (abs (pdg(1)))
    case (ELECTRON)
       select case (abs (pdg(2)))
       case (ELECTRON)
          call gircee (xc1, xc2, rng_obj = rng_obj)
       case (PHOTON)
          call girceg (xc1, xc2, rng_obj = rng_obj)
       end select
    case (PHOTON)
       select case (abs (pdg(2)))
       case (ELECTRON)
          call girceg (xc2, xc1, rng_obj = rng_obj)
       case (PHOTON)
          call gircgg (xc1, xc2, rng_obj = rng_obj)
       end select
    end select
    x = [xc1, xc2]
  end subroutine circe_generate

  module subroutine circe1_complete_kinematics &
       (sf_int, x, xb, f, r, rb, map)
    class(circe1_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: x
    real(default), dimension(:), intent(out) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: r
    real(default), dimension(:), intent(in) :: rb
    logical, intent(in) :: map
    x = r
    xb = rb
    sf_int%x = x
    sf_int%xb= xb
    f = 1
    if (sf_int%data%with_radiation) then
       call sf_int%split_momenta (x, xb)
    else
       call sf_int%reduce_momenta (x)
    end if
    select case (sf_int%status)
    case (SF_FAILED_KINEMATICS);  f = 0
    end select
  end subroutine circe1_complete_kinematics

  module subroutine circe1_inverse_kinematics &
       (sf_int, x, xb, f, r, rb, map, set_momenta)
    class(circe1_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(in) :: x
    real(default), dimension(:), intent(in) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: r
    real(default), dimension(:), intent(out) :: rb
    logical, intent(in) :: map
    logical, intent(in), optional :: set_momenta
    logical :: set_mom
    set_mom = .false.;  if (present (set_momenta))  set_mom = set_momenta
    r = x
    rb = xb
    sf_int%x = x
    sf_int%xb= xb
    f = 1
    if (set_mom) then
       call sf_int%split_momenta (x, xb)
       select case (sf_int%status)
       case (SF_FAILED_KINEMATICS);  f = 0
       end select
    end if
  end subroutine circe1_inverse_kinematics

  module subroutine circe1_apply (sf_int, scale, negative_sf, rescale, i_sub)
    class(circe1_t), intent(inout) :: sf_int
    real(default), intent(in) :: scale
    logical, intent(in), optional :: negative_sf
    class(sf_rescale_t), intent(in), optional :: rescale
    integer, intent(in), optional :: i_sub
    real(default), dimension(2) :: xb
    real(double), dimension(2) :: xc
    real(double), parameter :: one = 1
    associate (data => sf_int%data)
      xc = sf_int%x
      xb = sf_int%xb
      if (data%generate) then
         sf_int%f = 1
      else
         sf_int%f = 0
         if (all (sf_int%continuum)) then
            sf_int%f = circe (xc(1), xc(2), data%pdg_in(1), data%pdg_in(2))
         end if
         if (sf_int%continuum(2) .and. sf_int%peak(1)) then
            sf_int%f = sf_int%f &
                 + circe (one, xc(2), data%pdg_in(1), data%pdg_in(2)) &
                 * peak (xb(1), data%eps)
         end if
         if (sf_int%continuum(1) .and. sf_int%peak(2)) then
            sf_int%f = sf_int%f &
                 + circe (xc(1), one, data%pdg_in(1), data%pdg_in(2)) &
                 * peak (xb(2), data%eps)
         end if
         if (all (sf_int%peak)) then
            sf_int%f = sf_int%f &
                 + circe (one, one, data%pdg_in(1), data%pdg_in(2)) &
                 * peak (xb(1), data%eps) * peak (xb(2), data%eps)
         end if
      end if
    end associate
    call sf_int%set_matrix_element (cmplx (sf_int%f, kind=default))
    sf_int%status = SF_EVALUATED
  end subroutine circe1_apply

  function peak (x, eps) result (f)
    real(default), intent(in) :: x, eps
    real(default) :: f
    f = exp (-x / eps) / eps
  end function peak


end submodule sf_circe1_s

