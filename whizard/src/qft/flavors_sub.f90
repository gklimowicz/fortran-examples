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

submodule (flavors) flavors_s

  use io_units
  use diagnostics
  use physics_defs, only: INVALID
  use physics_defs, only: HADRON_REMNANT
  use physics_defs, only: HADRON_REMNANT_SINGLET
  use physics_defs, only: HADRON_REMNANT_TRIPLET
  use physics_defs, only: HADRON_REMNANT_OCTET

  implicit none

contains

  elemental module subroutine flavor_init_empty (flv)
    class(flavor_t), intent(inout) :: flv
    flv%f = UNDEFINED
    flv%hard_process = .false.
    flv%radiated = .false.
    flv%field_data => null ()
  end subroutine flavor_init_empty

  elemental module subroutine flavor_init (flv, f)
    class(flavor_t), intent(inout) :: flv
    integer, intent(in) :: f
    flv%f = f
    flv%hard_process = .false.
    flv%radiated = .false.
    flv%field_data => null ()
  end subroutine flavor_init

  impure elemental module subroutine flavor_init_field_data (flv, field_data)
    class(flavor_t), intent(inout) :: flv
    type(field_data_t), intent(in), target :: field_data
    flv%f = field_data%get_pdg ()
    flv%hard_process = .false.
    flv%radiated = .false.
    flv%field_data => field_data
  end subroutine flavor_init_field_data

  impure elemental module subroutine flavor_init_model (flv, f, model)
    class(flavor_t), intent(inout) :: flv
    integer, intent(in) :: f
    class(model_data_t), intent(in), target :: model
    flv%f = f
    flv%hard_process = .false.
    flv%radiated = .false.
    flv%field_data => model%get_field_ptr (f, check=.true.)
  end subroutine flavor_init_model

  impure elemental module subroutine flavor_init_model_alt (flv, f, model, alt_model)
    class(flavor_t), intent(inout) :: flv
    integer, intent(in) :: f
    class(model_data_t), intent(in), target :: model, alt_model
    flv%f = f
    flv%hard_process = .false.
    flv%radiated = .false.
    flv%field_data => model%get_field_ptr (f, check=.false.)
    if (.not. associated (flv%field_data)) then
       flv%field_data => alt_model%get_field_ptr (f, check=.false.)
       if (.not. associated (flv%field_data)) then
          write (msg_buffer, "(A,1x,I0,1x,A,1x,A,1x,A,1x,A)") &
               "Particle with code", f, &
               "found neither in model", char (model%get_name ()), &
               "nor in model", char (alt_model%get_name ())
          call msg_fatal ()
       end if
    end if
  end subroutine flavor_init_model_alt

  impure elemental module subroutine flavor_init_name_model (flv, name, model)
    class(flavor_t), intent(inout) :: flv
    type(string_t), intent(in) :: name
    class(model_data_t), intent(in), target :: model
    flv%f = model%get_pdg (name)
    flv%hard_process = .false.
    flv%radiated = .false.
    flv%field_data => model%get_field_ptr (name, check=.true.)
  end subroutine flavor_init_name_model

  elemental module subroutine flavor_tag_radiated (flv)
    class(flavor_t), intent(inout) :: flv
    flv%radiated = .true.
  end subroutine flavor_tag_radiated

  elemental module subroutine flavor_tag_hard_process (flv, hard)
    class(flavor_t), intent(inout) :: flv
    logical, intent(in), optional :: hard
    if (present (hard)) then
       flv%hard_process = hard
    else
       flv%hard_process = .true.
    end if
  end subroutine flavor_tag_hard_process

  elemental module subroutine flavor_undefine (flv)
    class(flavor_t), intent(inout) :: flv
    flv%f = UNDEFINED
    flv%field_data => null ()
  end subroutine flavor_undefine

  module subroutine flavor_write (flv, unit)
    class(flavor_t), intent(in) :: flv
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit);  if (u < 0)  return
    if (associated (flv%field_data)) then
       write (u, "(A)", advance="no")  "f("
    else
       write (u, "(A)", advance="no")  "p("
    end if
    write (u, "(I0)", advance="no")  flv%f
    if (flv%radiated) then
       write (u, "('*')", advance="no")
    end if
    if (msg_level (D_FLAVOR) >= DEBUG) then
       if (flv%hard_process) then
          write (u, "('#')", advance="no")
       end if
    end if
    write (u, "(A)", advance="no")  ")"
  end subroutine flavor_write

  module subroutine flavor_write_array (flv, unit)
    type(flavor_t), intent(in), dimension(:) :: flv
    integer, intent(in), optional :: unit
    integer :: u, i_flv
    u = given_output_unit (unit); if (u < 0) return
    do i_flv = 1, size (flv)
       call flv(i_flv)%write (u)
       if (i_flv /= size (flv)) write (u,"(A)", advance = "no") " / "
    end do
    write (u,"(A)")
  end subroutine flavor_write_array

  module subroutine flavor_write_raw (flv, u)
    class(flavor_t), intent(in) :: flv
    integer, intent(in) :: u
    write (u) flv%f
    write (u) flv%radiated
    write (u) flv%hard_process
  end subroutine flavor_write_raw

  module subroutine flavor_read_raw (flv, u, iostat)
    class(flavor_t), intent(out) :: flv
    integer, intent(in) :: u
    integer, intent(out), optional :: iostat
    read (u, iostat=iostat) flv%f
    if (present (iostat)) then
       if (iostat /= 0)  return
    end if
    read (u, iostat=iostat)  flv%radiated
    read (u, iostat=iostat)  flv%hard_process
  end subroutine flavor_read_raw

  impure elemental module subroutine flavor_set_model_single (flv, model)
    class(flavor_t), intent(inout) :: flv
    class(model_data_t), intent(in), target :: model
    if (flv%f /= UNDEFINED) &
         flv%field_data => model%get_field_ptr (flv%f)
  end subroutine flavor_set_model_single

  elemental module function flavor_is_defined (flv) result (defined)
    class(flavor_t), intent(in) :: flv
    logical :: defined
    defined = flv%f /= UNDEFINED
  end function flavor_is_defined

  elemental module function flavor_is_valid (flv) result (valid)
    class(flavor_t), intent(in) :: flv
    logical :: valid
    valid = flv%f /= INVALID
  end function flavor_is_valid

  elemental module function flavor_is_associated (flv) result (flag)
    class(flavor_t), intent(in) :: flv
    logical :: flag
    flag = associated (flv%field_data)
  end function flavor_is_associated

  elemental module function flavor_is_radiated (flv) result (flag)
    class(flavor_t), intent(in) :: flv
    logical :: flag
    flag = flv%radiated
  end function flavor_is_radiated

  elemental module function flavor_is_hard_process (flv) result (flag)
    class(flavor_t), intent(in) :: flv
    logical :: flag
    flag = flv%hard_process
  end function flavor_is_hard_process

  elemental module function flavor_get_pdg (flv) result (f)
    integer :: f
    class(flavor_t), intent(in) :: flv
    f = flv%f
  end function flavor_get_pdg

  elemental module function flavor_get_pdg_anti (flv) result (f)
    integer :: f
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       if (flv%field_data%has_antiparticle ()) then
          f = -flv%f
       else
          f = flv%f
       end if
    else
       f = 0
    end if
  end function flavor_get_pdg_anti

  elemental module function flavor_get_pdg_abs (flv) result (f)
    integer :: f
    class(flavor_t), intent(in) :: flv
    f = abs (flv%f)
  end function flavor_get_pdg_abs

  elemental module function flavor_is_visible (flv) result (flag)
    logical :: flag
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       flag = flv%field_data%is_visible ()
    else
       flag = .false.
    end if
  end function flavor_is_visible

  elemental module function flavor_is_parton (flv) result (flag)
    logical :: flag
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       flag = flv%field_data%is_parton ()
    else
       flag = .false.
    end if
  end function flavor_is_parton

  elemental module function flavor_is_beam_remnant (flv) result (flag)
    logical :: flag
    class(flavor_t), intent(in) :: flv
    select case (abs (flv%f))
    case (HADRON_REMNANT, &
         HADRON_REMNANT_SINGLET, HADRON_REMNANT_TRIPLET, HADRON_REMNANT_OCTET)
       flag = .true.
    case default
       flag = .false.
    end select
  end function flavor_is_beam_remnant

  elemental module function flavor_is_gauge (flv) result (flag)
    logical :: flag
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       flag = flv%field_data%is_gauge ()
    else
       flag = .false.
    end if
  end function flavor_is_gauge

  elemental module function flavor_is_left_handed (flv) result (flag)
    logical :: flag
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       if (flv%f > 0) then
          flag = flv%field_data%is_left_handed ()
       else
          flag = flv%field_data%is_right_handed ()
       end if
    else
       flag = .false.
    end if
  end function flavor_is_left_handed

  elemental module function flavor_is_right_handed (flv) result (flag)
    logical :: flag
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       if (flv%f > 0) then
          flag = flv%field_data%is_right_handed ()
       else
          flag = flv%field_data%is_left_handed ()
       end if
    else
       flag = .false.
    end if
  end function flavor_is_right_handed

  elemental module function flavor_is_antiparticle (flv) result (flag)
    logical :: flag
    class(flavor_t), intent(in) :: flv
    flag = flv%f < 0
  end function flavor_is_antiparticle

  elemental module function flavor_has_antiparticle (flv) result (flag)
    logical :: flag
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       flag = flv%field_data%has_antiparticle ()
    else
       flag = .false.
    end if
  end function flavor_has_antiparticle

  elemental module function flavor_is_stable (flv) result (flag)
    logical :: flag
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       flag = flv%field_data%is_stable (anti = flv%f < 0)
    else
       flag = .true.
    end if
  end function flavor_is_stable

  module subroutine flavor_get_decays (flv, decay)
    class(flavor_t), intent(in) :: flv
    type(string_t), dimension(:), intent(out), allocatable :: decay
    logical :: anti
    anti = flv%f < 0
    if (.not. flv%field_data%is_stable (anti)) then
       call flv%field_data%get_decays (decay, anti)
    end if
  end subroutine flavor_get_decays

  elemental module function flavor_decays_isotropically (flv) result (flag)
    logical :: flag
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       flag = flv%field_data%decays_isotropically (anti = flv%f < 0)
    else
       flag = .true.
    end if
  end function flavor_decays_isotropically

  elemental module function flavor_decays_diagonal (flv) result (flag)
    logical :: flag
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       flag = flv%field_data%decays_diagonal (anti = flv%f < 0)
    else
       flag = .true.
    end if
  end function flavor_decays_diagonal

  elemental module function flavor_has_decay_helicity (flv) result (flag)
    logical :: flag
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       flag = flv%field_data%has_decay_helicity (anti = flv%f < 0)
    else
       flag = .false.
    end if
  end function flavor_has_decay_helicity

  elemental module function flavor_get_decay_helicity (flv) result (hel)
    integer :: hel
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       hel = flv%field_data%decay_helicity (anti = flv%f < 0)
    else
       hel = 0
    end if
  end function flavor_get_decay_helicity

  elemental module function flavor_is_polarized (flv) result (flag)
    logical :: flag
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       flag = flv%field_data%is_polarized (anti = flv%f < 0)
    else
       flag = .false.
    end if
  end function flavor_is_polarized

  elemental module function flavor_get_name (flv) result (name)
    type(string_t) :: name
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       name = flv%field_data%get_name (flv%f < 0)
    else
       name = "?"
    end if
  end function flavor_get_name

  elemental module function flavor_get_tex_name (flv) result (name)
    type(string_t) :: name
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       name = flv%field_data%get_tex_name (flv%f < 0)
    else
       name = "?"
    end if
  end function flavor_get_tex_name

  elemental module function flavor_get_spin_type (flv) result (type)
    integer :: type
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       type = flv%field_data%get_spin_type ()
    else
       type = 1
    end if
  end function flavor_get_spin_type

  elemental module function flavor_get_multiplicity (flv) result (type)
    integer :: type
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       type = flv%field_data%get_multiplicity ()
    else
       type = 1
    end if
  end function flavor_get_multiplicity

  elemental module function flavor_get_isospin_type (flv) result (type)
    integer :: type
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       type = flv%field_data%get_isospin_type ()
    else
       type = 1
    end if
  end function flavor_get_isospin_type

  elemental module function flavor_get_charge_type (flv) result (type)
    integer :: type
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       type = flv%field_data%get_charge_type ()
    else
       type = 1
    end if
  end function flavor_get_charge_type

  elemental module function flavor_get_color_type (flv) result (type)
    integer :: type
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       if (flavor_is_antiparticle (flv)) then
          type = - flv%field_data%get_color_type ()
       else
          type = flv%field_data%get_color_type ()
       end if
       select case (type)
       case (-1,-8);  type = abs (type)
       end select
    else
       type = 1
    end if
  end function flavor_get_color_type

  elemental module function flavor_get_charge (flv) result (charge)
    real(default) :: charge
    class(flavor_t), intent(in) :: flv
    integer :: charge_type
    if (associated (flv%field_data)) then
       charge_type = flv%get_charge_type ()
       if (charge_type == 0 .or. charge_type == 1) then
          charge = 0
       else
          if (flavor_is_antiparticle (flv)) then
             charge = - flv%field_data%get_charge ()
          else
             charge = flv%field_data%get_charge ()
          end if
       end if
    else
       charge = 0
    end if
  end function flavor_get_charge

  elemental module function flavor_get_mass (flv) result (mass)
    real(default) :: mass
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       mass = flv%field_data%get_mass ()
    else
       mass = 0
    end if
  end function flavor_get_mass

  elemental module function flavor_get_width (flv) result (width)
    real(default) :: width
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       width = flv%field_data%get_width ()
    else
       width = 0
    end if
  end function flavor_get_width

  elemental module function flavor_get_isospin (flv) result (isospin)
    real(default) :: isospin
    class(flavor_t), intent(in) :: flv
    if (associated (flv%field_data)) then
       if (flavor_is_antiparticle (flv)) then
          isospin = - flv%field_data%get_isospin ()
       else
          isospin = flv%field_data%get_isospin ()
       end if
    else
       isospin = 0
    end if
  end function flavor_get_isospin

  elemental module function flavor_match (flv1, flv2) result (eq)
    logical :: eq
    class(flavor_t), intent(in) :: flv1, flv2
    if (flv1%f /= UNDEFINED .and. flv2%f /= UNDEFINED) then
       eq = flv1%f == flv2%f
    else
       eq = .true.
    end if
  end function flavor_match

  elemental module function flavor_eq (flv1, flv2) result (eq)
    logical :: eq
    class(flavor_t), intent(in) :: flv1, flv2
    if (flv1%f /= UNDEFINED .and. flv2%f /= UNDEFINED) then
       eq = flv1%f == flv2%f
    else if (flv1%f == UNDEFINED .and. flv2%f == UNDEFINED) then
       eq = .true.
    else
       eq = .false.
    end if
  end function flavor_eq

  elemental module function flavor_neq (flv1, flv2) result (neq)
    logical :: neq
    class(flavor_t), intent(in) :: flv1, flv2
    if (flv1%f /= UNDEFINED .and. flv2%f /= UNDEFINED) then
       neq = flv1%f /= flv2%f
    else if (flv1%f == UNDEFINED .and. flv2%f == UNDEFINED) then
       neq = .false.
    else
       neq = .true.
    end if
  end function flavor_neq

  module function merge_flavors0 (flv1, flv2) result (flv)
    type(flavor_t) :: flv
    type(flavor_t), intent(in) :: flv1, flv2
    if (flavor_is_defined (flv1) .and. flavor_is_defined (flv2)) then
       if (flv1 == flv2) then
          flv = flv1
       else
          flv%f = INVALID
       end if
    else if (flavor_is_defined (flv1)) then
       flv = flv1
    else if (flavor_is_defined (flv2)) then
       flv = flv2
    end if
  end function merge_flavors0

  module function merge_flavors1 (flv1, flv2) result (flv)
    type(flavor_t), dimension(:), intent(in) :: flv1, flv2
    type(flavor_t), dimension(size(flv1)) :: flv
    integer :: i
    do i = 1, size (flv1)
       flv(i) = flv1(i) .merge. flv2(i)
    end do
  end function merge_flavors1

  module function color_from_flavor0 (flv, c_seed, reverse) result (col)
    type(color_t) :: col
    type(flavor_t), intent(in) :: flv
    integer, intent(in), optional :: c_seed
    logical, intent(in), optional :: reverse
    integer, save :: c = 1
    logical :: rev
    if (present (c_seed))  c = c_seed
    rev = .false.;  if (present (reverse)) rev = reverse
    select case (flavor_get_color_type (flv))
    case (1)
       call col%init ()
    case (3)
       call col%init ([c]);  c = c + 1
    case (-3)
       call col%init ([-c]);  c = c + 1
    case (8)
       if (rev) then
          call col%init ([c+1, -c]);  c = c + 2
       else
          call col%init ([c, -(c+1)]);  c = c + 2
       end if
    end select
  end function color_from_flavor0

  module function color_from_flavor1 (flv, c_seed, reverse) result (col)
    type(flavor_t), dimension(:), intent(in) :: flv
    integer, intent(in), optional :: c_seed
    logical, intent(in), optional :: reverse
    type(color_t), dimension(size(flv)) :: col
    integer :: i
    col(1) = color_from_flavor0 (flv(1), c_seed, reverse)
    do i = 2, size (flv)
       col(i) = color_from_flavor0 (flv(i), reverse=reverse)
    end do
  end function color_from_flavor1

  module function flavor_anti (flv) result (aflv)
    type(flavor_t) :: aflv
    class(flavor_t), intent(in) :: flv
    if (flavor_has_antiparticle (flv)) then
       aflv%f = - flv%f
    else
       aflv%f = flv%f
    end if
    aflv%field_data => flv%field_data
  end function flavor_anti


end submodule flavors_s

