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

submodule (sf_circe2) sf_circe2_s

  use io_units
  use format_defs, only: FMT_19
  use numeric_utils
  use diagnostics
  use physics_defs, only: PHOTON, ELECTRON, MUON
  use lorentz
  use colors
  use helicities
  use quantum_numbers
  use state_matrices

  implicit none

contains

  module subroutine circe2_data_init (data, os_data, model, pdg_in, &
       sqrts, polarized, beam_pol, file, design)
    class(circe2_data_t), intent(out) :: data
    type(os_data_t), intent(in) :: os_data
    class(model_data_t), intent(in), target :: model
    type(pdg_array_t), dimension(2), intent(in) :: pdg_in
    real(default), intent(in) :: sqrts
    logical, intent(in) :: polarized, beam_pol
    type(string_t), intent(in) :: file, design
    integer :: h
    data%model => model
    if (any (pdg_in%get_length () /= 1)) then
       call msg_fatal ("CIRCE2: incoming beam particles must be unique")
    end if
    call data%flv_in(1)%init (pdg_in(1)%get (1), model)
    call data%flv_in(2)%init (pdg_in(2)%get (1), model)
    data%pdg_in = data%flv_in%get_pdg ()
    data%sqrts = sqrts
    data%polarized = polarized
    data%beams_polarized = beam_pol
    data%filename = file
    data%design = design
    call data%check_file (os_data)
    call circe2_load (circe2_global_state, trim (char(data%file)), &
            trim (char(data%design)), data%sqrts, data%error)
    call data%check ()
    data%lumi = circe2_luminosity (circe2_global_state, data%pdg_in, [0, 0])
    if (vanishes (data%lumi)) then
       call msg_fatal ("CIRCE2: luminosity vanishes for specified beams.")
    end if
    if (data%polarized) then
       do h = 1, 4
          data%lumi_hel_frac(h) = &
               circe2_luminosity (circe2_global_state, data%pdg_in, &
                            [data%h1(h), data%h2(h)]) &
               / data%lumi
       end do
    end if
  end subroutine circe2_data_init

  module subroutine circe2_data_set_generator_mode (data, rng_factory)
    class(circe2_data_t), intent(inout) :: data
    class(rng_factory_t), intent(inout), allocatable :: rng_factory
    call move_alloc (from = rng_factory, to = data%rng_factory)
  end subroutine circe2_data_set_generator_mode

  module subroutine circe2_check_file (data, os_data)
    class(circe2_data_t), intent(inout) :: data
    type(os_data_t), intent(in) :: os_data
    logical :: exist
    type(string_t) :: file
    file = data%filename
    if (file == "") &
         call msg_fatal ("CIRCE2: $circe2_file is not set")
    inquire (file = char (file), exist = exist)
    if (exist) then
       data%file = file
    else
       file = os_data%whizard_circe2path // "/" // data%filename
       inquire (file = char (file), exist = exist)
       if (exist) then
          data%file = file
       else
          call msg_fatal ("CIRCE2: data file '" // char (data%filename) &
               // "' not found")
       end if
    end if
  end subroutine circe2_check_file

  module subroutine circe2_data_check (data)
    class(circe2_data_t), intent(in) :: data
    type(flavor_t) :: flv_photon, flv_electron, flv_muon
    call flv_photon%init (PHOTON, data%model)
    if (.not. flv_photon%is_defined ()) then
       call msg_fatal ("CIRCE2: model must contain photon")
    end if
    if (any (abs (data%pdg_in) /= PHOTON .and. abs (data%pdg_in) /= &
         ELECTRON .and. abs (data%pdg_in) /= MUON)) then
       call msg_fatal ("CIRCE2: applicable only for e+e-, mu+mu- or " // &
            "photon collisions")
    end if
    if (any (abs (data%pdg_in) == ELECTRON)) then
       call flv_electron%init (ELECTRON, data%model)
       if (.not. flv_electron%is_defined ()) then
          call msg_fatal ("CIRCE2: model must contain electron")
       end if
    end if
    if (any (abs (data%pdg_in) == MUON)) then
       call flv_muon%init (MUON, data%model)
       if (.not. flv_muon%is_defined ()) then
          call msg_fatal ("CIRCE2: model must contain muon")
       end if
    end if
    select case (data%error)
    case (-1)
       call msg_fatal ("CIRCE2: data file not found.")
    case (-2)
       call msg_fatal ("CIRCE2: beam setup does not match data file.")
    case (-3)
       call msg_fatal ("CIRCE2: invalid format of data file.")
    case (-4)
       call msg_fatal ("CIRCE2: data file too large.")
    end select
  end subroutine circe2_data_check

  module subroutine circe2_data_write (data, unit, verbose)
    class(circe2_data_t), intent(in) :: data
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: u, h
    logical :: verb
    verb = .false.;  if (present (verbose))  verb = verbose
    u = given_output_unit (unit)
    write (u, "(1x,A)") "CIRCE2 data:"
    write (u, "(3x,A,A)")       "file   = ", char(data%filename)
    write (u, "(3x,A,A)")       "design = ", char(data%design)
    write (u, "(3x,A," // FMT_19 // ")") "sqrts  = ", data%sqrts
    write (u, "(3x,A,A,A,A)")   "prt_in = ", &
         char (data%flv_in(1)%get_name ()), &
         ", ", char (data%flv_in(2)%get_name ())
    write (u, "(3x,A,L1)")      "polarized  = ", data%polarized
    write (u, "(3x,A,L1)")      "beams pol. = ", data%beams_polarized
    write (u, "(3x,A," // FMT_19 // ")") "luminosity = ", data%lumi
    if (data%polarized) then
       do h = 1, 4
          write (u, "(6x,'(',I2,1x,I2,')',1x,'=',1x)", advance="no") &
               data%h1(h), data%h2(h)
          write (u, "(6x, " // FMT_19 // ")")  data%lumi_hel_frac(h)
       end do
    end if
    if (verb) then
       call data%rng_factory%write (u)
    end if
  end subroutine circe2_data_write

  module function circe2_data_is_generator (data) result (flag)
    class(circe2_data_t), intent(in) :: data
    logical :: flag
    flag = .true.
  end function circe2_data_is_generator

  module function circe2_data_get_n_par (data) result (n)
    class(circe2_data_t), intent(in) :: data
    integer :: n
    n = 2
  end function circe2_data_get_n_par

  module subroutine circe2_data_get_pdg_out (data, pdg_out)
    class(circe2_data_t), intent(in) :: data
    type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    integer :: i, n
    n = 2
    do i = 1, n
       pdg_out(i) = data%pdg_in(i)
    end do
  end subroutine circe2_data_get_pdg_out

  module function circe2_data_get_beam_file (data) result (file)
    class(circe2_data_t), intent(in) :: data
    type(string_t) :: file
    file = "CIRCE2: " // data%filename
  end function circe2_data_get_beam_file

  module subroutine rng_obj_generate (rng_obj, u)
    class(rng_obj_t), intent(inout) :: rng_obj
    real(default), intent(out) :: u
    real(default) :: x
    call rng_obj%rng%generate (x)
    u = x
  end subroutine rng_obj_generate

  module function circe2_type_string (object) result (string)
    class(circe2_t), intent(in) :: object
    type(string_t) :: string
    if (associated (object%data)) then
       string = "CIRCE2: " // object%data%design
    else
       string = "CIRCE2: [undefined]"
    end if
  end function circe2_type_string

  module subroutine circe2_write (object, unit, testflag)
    class(circe2_t), intent(in) :: object
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: testflag
    integer :: u
    u = given_output_unit (unit)
    if (associated (object%data)) then
       call object%data%write (u)
       call object%base_write (u, testflag)
    else
       write (u, "(1x,A)")  "CIRCE2 data: [undefined]"
    end if
  end subroutine circe2_write

  module subroutine circe2_init (sf_int, data)
    class(circe2_t), intent(out) :: sf_int
    class(sf_data_t), intent(in), target :: data
    logical, dimension(4) :: mask_h
    real(default), dimension(2) :: m2_array
    real(default), dimension(0) :: null_array
    type(quantum_numbers_mask_t), dimension(4) :: mask
    type(quantum_numbers_t), dimension(4) :: qn
    type(helicity_t) :: hel
    type(color_t) :: col0
    integer :: h
    select type (data)
    type is (circe2_data_t)
       if (data%polarized .and. data%beams_polarized) then
          call msg_fatal ("CIRCE2: Beam polarization can't be set &
               &for polarized data file")
       else if (data%beams_polarized) then
          call msg_warning ("CIRCE2: User-defined beam polarization set &
               &for unpolarized CIRCE2 data file")
       end if
       mask_h(1:2) = .not. data%beams_polarized
       mask_h(3:4) = .not. (data%polarized .or. data%beams_polarized)
       mask = quantum_numbers_mask (.false., .false., mask_h)
       m2_array(:) = (data%flv_in(:)%get_mass ())**2
       call sf_int%base_init (mask, m2_array, null_array, m2_array)
       sf_int%data => data
       if (data%polarized) then
          if (vanishes (sum (data%lumi_hel_frac)) .or. &
               any (data%lumi_hel_frac < 0)) then
             call msg_fatal ("CIRCE2: Helicity-dependent lumi " &
                  // "fractions all vanish or",  &
                  [var_str ("are negative: Please inspect the " &
                  // "CIRCE2 file or "), &
                   var_str ("switch off the polarized" // &
                  " option for CIRCE2.")])
          else
             call sf_int%selector%init (data%lumi_hel_frac)
          end if
       end if
       call col0%init ()
       if (data%beams_polarized) then
          do h = 1, 4
             call hel%init (data%h1(h))
             call qn(1)%init &
                  (flv = data%flv_in(1), col = col0, hel = hel)
             call qn(3)%init &
                  (flv = data%flv_in(1), col = col0, hel = hel)
             call hel%init (data%h2(h))
             call qn(2)%init &
                  (flv = data%flv_in(2), col = col0, hel = hel)
             call qn(4)%init &
                  (flv = data%flv_in(2), col = col0, hel = hel)
             call sf_int%add_state (qn)
          end do
       else if (data%polarized) then
          call qn(1)%init (flv = data%flv_in(1), col = col0)
          call qn(2)%init (flv = data%flv_in(2), col = col0)
          do h = 1, 4
             call hel%init (data%h1(h))
             call qn(3)%init &
                  (flv = data%flv_in(1), col = col0, hel = hel)
             call hel%init (data%h2(h))
             call qn(4)%init &
                  (flv = data%flv_in(2), col = col0, hel = hel)
             call sf_int%add_state (qn)
          end do
       else
          call qn(1)%init (flv = data%flv_in(1), col = col0)
          call qn(2)%init (flv = data%flv_in(2), col = col0)
          call qn(3)%init (flv = data%flv_in(1), col = col0)
          call qn(4)%init (flv = data%flv_in(2), col = col0)
          call sf_int%add_state (qn)
       end if
       call sf_int%freeze ()
       call sf_int%set_incoming ([1,2])
       call sf_int%set_outgoing ([3,4])
       call sf_int%data%rng_factory%make (sf_int%rng_obj%rng)
       sf_int%status = SF_INITIAL
    end select
  end subroutine circe2_init

  module function circe2_is_generator (sf_int) result (flag)
    class(circe2_t), intent(in) :: sf_int
    logical :: flag
    flag = sf_int%data%is_generator ()
  end function circe2_is_generator

  module subroutine circe2_generate_whizard_free (sf_int, r, rb, x_free)
    class(circe2_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: r, rb
    real(default), intent(inout) :: x_free
    integer :: h_sel
    if (sf_int%data%polarized) then
       call sf_int%selector%generate (sf_int%rng_obj%rng, h_sel)
    else
       h_sel = 0
    end if
    sf_int%h_sel = h_sel
    call circe2_generate_whizard (r, sf_int%data%pdg_in, &
         [sf_int%data%h1(h_sel), sf_int%data%h2(h_sel)], &
         sf_int%rng_obj)
    rb = 1 - r
    x_free = x_free * product (r)
  end subroutine circe2_generate_whizard_free

  module subroutine circe2_generate_whizard (x, pdg, hel, rng_obj)
    real(default), dimension(2), intent(out) :: x
    integer, dimension(2), intent(in) :: pdg
    integer, dimension(2), intent(in) :: hel
    class(rng_obj_t), intent(inout) :: rng_obj
    call circe2_generate (circe2_global_state, rng_obj, x, pdg, hel)
  end subroutine circe2_generate_whizard

  module subroutine circe2_complete_kinematics (sf_int, x, xb, f, r, rb, map)
    class(circe2_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(out) :: x
    real(default), dimension(:), intent(out) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(in) :: r
    real(default), dimension(:), intent(in) :: rb
    logical, intent(in) :: map
    if (map) then
       call msg_fatal ("CIRCE2: map flag not supported")
    else
       x = r
       xb= rb
       f = 1
    end if
    call sf_int%reduce_momenta (x)
  end subroutine circe2_complete_kinematics

  module subroutine circe2_inverse_kinematics &
       (sf_int, x, xb, f, r, rb, map, set_momenta)
    class(circe2_t), intent(inout) :: sf_int
    real(default), dimension(:), intent(in) :: x
    real(default), dimension(:), intent(in) :: xb
    real(default), intent(out) :: f
    real(default), dimension(:), intent(out) :: r
    real(default), dimension(:), intent(out) :: rb
    logical, intent(in) :: map
    logical, intent(in), optional :: set_momenta
    logical :: set_mom
    set_mom = .false.;  if (present (set_momenta))  set_mom = set_momenta
    if (map) then
       call msg_fatal ("CIRCE2: map flag not supported")
    else
       r = x
       rb= xb
       f = 1
    end if
    if (set_mom) then
       call sf_int%reduce_momenta (x)
    end if
  end subroutine circe2_inverse_kinematics

  module subroutine circe2_apply (sf_int, scale, negative_sf, rescale, i_sub)
    class(circe2_t), intent(inout) :: sf_int
    real(default), intent(in) :: scale
    logical, intent(in), optional :: negative_sf
    class(sf_rescale_t), intent(in), optional :: rescale
    integer, intent(in), optional :: i_sub
    complex(default) :: f
    associate (data => sf_int%data)
      f = 1
      if (data%beams_polarized) then
         call sf_int%set_matrix_element (f)
      else if (data%polarized) then
         call sf_int%set_matrix_element (sf_int%h_sel, f)
      else
         call sf_int%set_matrix_element (1, f)
      end if
    end associate
    sf_int%status = SF_EVALUATED
  end subroutine circe2_apply


end submodule sf_circe2_s

