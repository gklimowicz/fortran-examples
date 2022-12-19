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

submodule (model_data) model_data_s

  use format_defs, only: FMT_19
  use io_units
  use diagnostics
  use md5
  use hashes, only: hash

  implicit none

  integer, parameter :: VERTEX_TABLE_SCALE_FACTOR = 60

contains

  module subroutine par_write (par, unit)
    class(modelpar_data_t), intent(in) :: par
    integer, intent(in), optional :: unit
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A,1x,A)", advance="no")  char (par%name), "= "
    select type (par)
    type is (modelpar_real_t)
       write (u, "(" // FMT_19 // ")", advance="no")  par%value
    type is (modelpar_complex_t)
       write (u, "(" // FMT_19 // ",1x,'+',1x," // FMT_19 // ",1x,'I')", &
            advance="no")  par%value
    end select
  end subroutine par_write

  module subroutine par_show (par, l, u)
    class(modelpar_data_t), intent(in) :: par
    integer, intent(in) :: l, u
    character(len=l) :: buffer
    buffer = par%name
    select type (par)
    type is (modelpar_real_t)
       write (u, "(4x,A,1x,'=',1x," // FMT_19 // ")")  buffer, par%value
    type is (modelpar_complex_t)
       write (u, "(4x,A,1x,'=',1x," // FMT_19 // ",1x,'+',1x," &
            // FMT_19 // ",1x,'I')")  buffer, par%value
    end select
  end subroutine par_show

  module subroutine modelpar_data_init_real (par, name, value)
    class(modelpar_data_t), intent(out) :: par
    type(string_t), intent(in) :: name
    real(default), intent(in) :: value
    par%name = name
    par = value
  end subroutine modelpar_data_init_real

  module subroutine modelpar_data_init_complex (par, name, value)
    class(modelpar_data_t), intent(out) :: par
    type(string_t), intent(in) :: name
    complex(default), intent(in) :: value
    par%name = name
    par = value
  end subroutine modelpar_data_init_complex

  elemental module subroutine modelpar_data_set_real (par, value)
    class(modelpar_data_t), intent(inout) :: par
    real(default), intent(in) :: value
    select type (par)
    type is (modelpar_real_t)
       par%value = value
    type is (modelpar_complex_t)
       par%value = value
    end select
  end subroutine modelpar_data_set_real

  elemental module subroutine modelpar_data_set_complex (par, value)
    class(modelpar_data_t), intent(inout) :: par
    complex(default), intent(in) :: value
    select type (par)
    type is (modelpar_real_t)
       par%value = value
    type is (modelpar_complex_t)
       par%value = value
    end select
  end subroutine modelpar_data_set_complex

  module function modelpar_data_get_name (par) result (name)
    class(modelpar_data_t), intent(in) :: par
    type(string_t) :: name
    name = par%name
  end function modelpar_data_get_name

  elemental module function modelpar_data_get_real (par) result (value)
    class(modelpar_data_t), intent(in), target :: par
    real(default) :: value
    select type (par)
    type is (modelpar_real_t)
       value = par%value
    type is (modelpar_complex_t)
       value = par%value
    end select
  end function modelpar_data_get_real

  elemental module function modelpar_data_get_complex (par) result (value)
    class(modelpar_data_t), intent(in), target :: par
    complex(default) :: value
    select type (par)
    type is (modelpar_real_t)
       value = par%value
    type is (modelpar_complex_t)
       value = par%value
    end select
  end function modelpar_data_get_complex

  module function modelpar_data_get_real_ptr (par) result (ptr)
    class(modelpar_data_t), intent(in), target :: par
    real(default), pointer :: ptr
    select type (par)
    type is (modelpar_real_t)
       ptr => par%value
    class default
       ptr => null ()
    end select
  end function modelpar_data_get_real_ptr

  module function modelpar_data_get_complex_ptr (par) result (ptr)
    class(modelpar_data_t), intent(in), target :: par
    complex(default), pointer :: ptr
    select type (par)
    type is (modelpar_complex_t)
       ptr => par%value
    class default
       ptr => null ()
    end select
  end function modelpar_data_get_complex_ptr

  module subroutine field_data_init (prt, longname, pdg)
    class(field_data_t), intent(out) :: prt
    type(string_t), intent(in) :: longname
    integer, intent(in) :: pdg
    prt%longname = longname
    prt%pdg = pdg
    prt%tex_name = ""
    prt%tex_anti = ""
  end subroutine field_data_init

  module subroutine field_data_copy_from (prt, prt_src)
    class(field_data_t), intent(inout) :: prt
    class(field_data_t), intent(in) :: prt_src
    prt%visible = prt_src%visible
    prt%parton = prt_src%parton
    prt%gauge = prt_src%gauge
    prt%left_handed = prt_src%left_handed
    prt%right_handed = prt_src%right_handed
    prt%p_is_stable =             prt_src%p_is_stable
    prt%p_decays_isotropically =  prt_src%p_decays_isotropically
    prt%p_decays_diagonal =       prt_src%p_decays_diagonal
    prt%p_has_decay_helicity =    prt_src%p_has_decay_helicity
    prt%p_decay_helicity =        prt_src%p_decay_helicity
    prt%p_decays_diagonal =       prt_src%p_decays_diagonal
    prt%a_is_stable =             prt_src%a_is_stable
    prt%a_decays_isotropically =  prt_src%a_decays_isotropically
    prt%a_decays_diagonal =       prt_src%a_decays_diagonal
    prt%a_has_decay_helicity =    prt_src%a_has_decay_helicity
    prt%a_decay_helicity =        prt_src%a_decay_helicity
    prt%p_polarized =             prt_src%p_polarized
    prt%a_polarized =             prt_src%a_polarized
    prt%spin_type = prt_src%spin_type
    prt%isospin_type = prt_src%isospin_type
    prt%charge_type = prt_src%charge_type
    prt%color_type = prt_src%color_type
    prt%has_anti = prt_src%has_anti
    if (allocated (prt_src%name)) then
       if (allocated (prt%name))  deallocate (prt%name)
       allocate (prt%name (size (prt_src%name)), source = prt_src%name)
    end if
    if (allocated (prt_src%anti)) then
       if (allocated (prt%anti))  deallocate (prt%anti)
       allocate (prt%anti (size (prt_src%anti)), source = prt_src%anti)
    end if
    prt%tex_name = prt_src%tex_name
    prt%tex_anti = prt_src%tex_anti
    if (allocated (prt_src%p_decay)) then
       if (allocated (prt%p_decay))  deallocate (prt%p_decay)
       allocate (prt%p_decay (size (prt_src%p_decay)), source = prt_src%p_decay)
    end if
    if (allocated (prt_src%a_decay)) then
       if (allocated (prt%a_decay))  deallocate (prt%a_decay)
       allocate (prt%a_decay (size (prt_src%a_decay)), source = prt_src%a_decay)
    end if
  end subroutine field_data_copy_from

  module subroutine field_data_set (prt, &
       is_visible, is_parton, is_gauge, is_left_handed, is_right_handed, &
       p_is_stable, p_decays_isotropically, p_decays_diagonal, &
       p_decay_helicity, &
       a_is_stable, a_decays_isotropically, a_decays_diagonal, &
       a_decay_helicity, &
       p_polarized, a_polarized, &
       name, anti, tex_name, tex_anti, &
       spin_type, isospin_type, charge_type, color_type, &
       mass_data, width_data, &
       p_decay, a_decay)
    class(field_data_t), intent(inout) :: prt
    logical, intent(in), optional :: is_visible, is_parton, is_gauge
    logical, intent(in), optional :: is_left_handed, is_right_handed
    logical, intent(in), optional :: p_is_stable
    logical, intent(in), optional :: p_decays_isotropically, p_decays_diagonal
    integer, intent(in), optional :: p_decay_helicity
    logical, intent(in), optional :: a_is_stable
    logical, intent(in), optional :: a_decays_isotropically, a_decays_diagonal
    integer, intent(in), optional :: a_decay_helicity
    logical, intent(in), optional :: p_polarized, a_polarized
    type(string_t), dimension(:), intent(in), optional :: name, anti
    type(string_t), intent(in), optional :: tex_name, tex_anti
    integer, intent(in), optional :: spin_type, isospin_type
    integer, intent(in), optional :: charge_type, color_type
    class(modelpar_data_t), intent(in), pointer, optional :: mass_data, width_data
    type(string_t), dimension(:), intent(in), optional :: p_decay, a_decay
    if (present (is_visible))  prt%visible = is_visible
    if (present (is_parton))  prt%parton = is_parton
    if (present (is_gauge))  prt%gauge = is_gauge
    if (present (is_left_handed))  prt%left_handed = is_left_handed
    if (present (is_right_handed))  prt%right_handed = is_right_handed
    if (present (p_is_stable))  prt%p_is_stable = p_is_stable
    if (present (p_decays_isotropically)) &
          prt%p_decays_isotropically = p_decays_isotropically
    if (present (p_decays_diagonal)) &
          prt%p_decays_diagonal = p_decays_diagonal
    if (present (p_decay_helicity)) then
       prt%p_has_decay_helicity = .true.
       prt%p_decay_helicity = p_decay_helicity
    end if
    if (present (a_is_stable))  prt%a_is_stable = a_is_stable
    if (present (a_decays_isotropically)) &
          prt%a_decays_isotropically = a_decays_isotropically
    if (present (a_decays_diagonal)) &
          prt%a_decays_diagonal = a_decays_diagonal
    if (present (a_decay_helicity)) then
       prt%a_has_decay_helicity = .true.
       prt%a_decay_helicity = a_decay_helicity
    end if
    if (present (p_polarized)) prt%p_polarized = p_polarized
    if (present (a_polarized)) prt%a_polarized = a_polarized
    if (present (name)) then
       if (allocated (prt%name))  deallocate (prt%name)
       allocate (prt%name (size (name)), source = name)
    end if
    if (present (anti)) then
       if (allocated (prt%anti))  deallocate (prt%anti)
       allocate (prt%anti (size (anti)), source = anti)
       prt%has_anti = .true.
    end if
    if (present (tex_name))  prt%tex_name = tex_name
    if (present (tex_anti))  prt%tex_anti = tex_anti
    if (present (spin_type))  prt%spin_type = spin_type
    if (present (isospin_type))  prt%isospin_type = isospin_type
    if (present (charge_type))  prt%charge_type = charge_type
    if (present (color_type))  prt%color_type = color_type
    if (present (mass_data)) then
       prt%mass_data => mass_data
       if (associated (mass_data)) then
          prt%mass_val => mass_data%get_real_ptr ()
       else
          prt%mass_val => null ()
       end if
    end if
    if (present (width_data)) then
       prt%width_data => width_data
       if (associated (width_data)) then
          prt%width_val => width_data%get_real_ptr ()
       else
          prt%width_val => null ()
       end if
    end if
    if (present (spin_type) .or. present (mass_data)) then
       call prt%set_multiplicity ()
    end if
    if (present (p_decay)) then
       if (allocated (prt%p_decay))  deallocate (prt%p_decay)
       if (size (p_decay) > 0) &
            allocate (prt%p_decay (size (p_decay)), source = p_decay)
    end if
    if (present (a_decay)) then
       if (allocated (prt%a_decay))  deallocate (prt%a_decay)
       if (size (a_decay) > 0) &
            allocate (prt%a_decay (size (a_decay)), source = a_decay)
    end if
  end subroutine field_data_set

  module subroutine field_data_set_multiplicity (prt)
    class(field_data_t), intent(inout) :: prt
    if (prt%spin_type /= SCALAR) then
       if (associated (prt%mass_data)) then
          prt%multiplicity = prt%spin_type
       else if (prt%left_handed .or. prt%right_handed) then
          prt%multiplicity = 1
       else
          prt%multiplicity = 2
       end if
    end if
  end subroutine field_data_set_multiplicity

  module subroutine field_data_set_mass (prt, mass)
    class(field_data_t), intent(inout) :: prt
    real(default), intent(in) :: mass
    if (associated (prt%mass_val))  prt%mass_val = mass
  end subroutine field_data_set_mass

  module subroutine field_data_set_width (prt, width)
    class(field_data_t), intent(inout) :: prt
    real(default), intent(in) :: width
    if (associated (prt%width_val))  prt%width_val = width
  end subroutine field_data_set_width

  elemental module subroutine field_data_freeze (prt)
    class(field_data_t), intent(inout) :: prt
    if (.not. allocated (prt%name))  allocate (prt%name (0))
    if (.not. allocated (prt%anti))  allocate (prt%anti (0))
  end subroutine field_data_freeze

  module subroutine field_data_write (prt, unit)
    class(field_data_t), intent(in) :: prt
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    write (u, "(3x,A,1x,A)", advance="no") "particle", char (prt%longname)
    write (u, "(1x,I0)", advance="no") prt%pdg
    if (.not. prt%visible) write (u, "(2x,A)", advance="no") "invisible"
    if (prt%parton)  write (u, "(2x,A)", advance="no") "parton"
    if (prt%gauge)  write (u, "(2x,A)", advance="no") "gauge"
    if (prt%left_handed)  write (u, "(2x,A)", advance="no") "left"
    if (prt%right_handed)  write (u, "(2x,A)", advance="no") "right"
    write (u, *)
    write (u, "(5x,A)", advance="no") "name"
    if (allocated (prt%name)) then
       do i = 1, size (prt%name)
          write (u, "(1x,A)", advance="no")  '"' // char (prt%name(i)) // '"'
       end do
       write (u, *)
       if (prt%has_anti) then
          write (u, "(5x,A)", advance="no") "anti"
          do i = 1, size (prt%anti)
             write (u, "(1x,A)", advance="no")  '"' // char (prt%anti(i)) // '"'
          end do
          write (u, *)
       end if
       if (prt%tex_name /= "") then
          write (u, "(5x,A)")  &
               "tex_name " // '"' // char (prt%tex_name) // '"'
       end if
       if (prt%has_anti .and. prt%tex_anti /= "") then
          write (u, "(5x,A)")  &
               "tex_anti " // '"' // char (prt%tex_anti) // '"'
       end if
    else
       write (u, "(A)")  "???"
    end if
    write (u, "(5x,A)", advance="no") "spin "
    select case (mod (prt%spin_type - 1, 2))
    case (0);  write (u, "(I0)", advance="no") (prt%spin_type-1) / 2
    case default;  write (u, "(I0,A)", advance="no") prt%spin_type-1, "/2"
    end select
    ! write (u, "(2x,A,I1,A)") "! [multiplicity = ", prt%multiplicity, "]"
    if (abs (prt%isospin_type) /= 1) then
       write (u, "(2x,A)", advance="no") "isospin "
       select case (mod (abs (prt%isospin_type) - 1, 2))
       case (0);  write (u, "(I0)", advance="no") &
            sign (abs (prt%isospin_type) - 1, prt%isospin_type) / 2
       case default;  write (u, "(I0,A)", advance="no") &
            sign (abs (prt%isospin_type) - 1, prt%isospin_type), "/2"
       end select
    end if
    if (abs (prt%charge_type) /= 1) then
       write (u, "(2x,A)", advance="no") "charge "
       select case (mod (abs (prt%charge_type) - 1, 3))
       case (0);  write (u, "(I0)", advance="no") &
            sign (abs (prt%charge_type) - 1, prt%charge_type) / 3
       case default;  write (u, "(I0,A)", advance="no") &
            sign (abs (prt%charge_type) - 1, prt%charge_type), "/3"
       end select
    end if
    if (prt%color_type /= 1) then
       write (u, "(2x,A,I0)", advance="no") "color ", prt%color_type
    end if
    write (u, *)
    if (associated (prt%mass_data)) then
       write (u, "(5x,A)", advance="no") &
            "mass " // char (prt%mass_data%get_name ())
       if (associated (prt%width_data)) then
          write (u, "(2x,A)") &
               "width " // char (prt%width_data%get_name ())
       else
          write (u, *)
       end if
    end if
    call prt%write_decays (u)
  end subroutine field_data_write

  module subroutine field_data_write_decays (prt, unit)
    class(field_data_t), intent(in) :: prt
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    if (.not. prt%p_is_stable) then
       if (allocated (prt%p_decay)) then
          write (u, "(5x,A)", advance="no") "p_decay"
          do i = 1, size (prt%p_decay)
             write (u, "(1x,A)", advance="no")  char (prt%p_decay(i))
          end do
          if (prt%p_decays_isotropically) then
             write (u, "(1x,A)", advance="no")  "isotropic"
          else if (prt%p_decays_diagonal) then
             write (u, "(1x,A)", advance="no")  "diagonal"
          else if (prt%p_has_decay_helicity) then
             write (u, "(1x,A,I0)", advance="no")  "helicity = ", &
                  prt%p_decay_helicity
          end if
          write (u, *)
       end if
    else if (prt%p_polarized) then
       write (u, "(5x,A)")  "p_polarized"
    end if
    if (.not. prt%a_is_stable) then
       if (allocated (prt%a_decay)) then
          write (u, "(5x,A)", advance="no") "a_decay"
          do i = 1, size (prt%a_decay)
             write (u, "(1x,A)", advance="no")  char (prt%a_decay(i))
          end do
          if (prt%a_decays_isotropically) then
             write (u, "(1x,A)", advance="no")  "isotropic"
          else if (prt%a_decays_diagonal) then
             write (u, "(1x,A)", advance="no")  "diagonal"
          else if (prt%a_has_decay_helicity) then
             write (u, "(1x,A,I0)", advance="no")  "helicity = ", &
                  prt%a_decay_helicity
          end if
          write (u, *)
       end if
    else if (prt%a_polarized) then
       write (u, "(5x,A)")  "a_polarized"
    end if
  end subroutine field_data_write_decays

  module subroutine field_data_show (prt, l, u)
    class(field_data_t), intent(in) :: prt
    integer, intent(in) :: l, u
    character(len=l) :: buffer
    integer :: i
    type(string_t), dimension(:), allocatable :: decay
    buffer = prt%get_name (.false.)
    write (u, "(4x,A,1x,I8)", advance="no")  buffer, &
         prt%get_pdg ()
    if (prt%is_polarized ()) then
       write (u, "(3x,A)")  "polarized"
    else if (.not. prt%is_stable ()) then
       write (u, "(3x,A)", advance="no")  "decays:"
       call prt%get_decays (decay)
       do i = 1, size (decay)
          write (u, "(1x,A)", advance="no")  char (decay(i))
       end do
       write (u, *)
    else
       write (u, *)
    end if
    if (prt%has_antiparticle ()) then
       buffer = prt%get_name (.true.)
       write (u, "(4x,A,1x,I8)", advance="no")  buffer, &
            prt%get_pdg_anti ()
       if (prt%is_polarized (.true.)) then
          write (u, "(3x,A)")  "polarized"
       else if (.not. prt%is_stable (.true.)) then
          write (u, "(3x,A)", advance="no")  "decays:"
          call prt%get_decays (decay, .true.)
          do i = 1, size (decay)
             write (u, "(1x,A)", advance="no")  char (decay(i))
          end do
          write (u, *)
       else
          write (u, *)
       end if
    end if
  end subroutine field_data_show

  elemental module function field_data_get_pdg (prt) result (pdg)
    integer :: pdg
    class(field_data_t), intent(in) :: prt
    pdg = prt%pdg
  end function field_data_get_pdg

  elemental module function field_data_get_pdg_anti (prt) result (pdg)
    integer :: pdg
    class(field_data_t), intent(in) :: prt
    if (prt%has_anti) then
       pdg = - prt%pdg
    else
       pdg = prt%pdg
    end if
  end function field_data_get_pdg_anti

  elemental module function field_data_is_visible (prt) result (flag)
    logical :: flag
    class(field_data_t), intent(in) :: prt
    flag = prt%visible
  end function field_data_is_visible

  elemental module function field_data_is_parton (prt) result (flag)
    logical :: flag
    class(field_data_t), intent(in) :: prt
    flag = prt%parton
  end function field_data_is_parton

  elemental module function field_data_is_gauge (prt) result (flag)
    logical :: flag
    class(field_data_t), intent(in) :: prt
    flag = prt%gauge
  end function field_data_is_gauge

  elemental module function field_data_is_left_handed (prt) result (flag)
    logical :: flag
    class(field_data_t), intent(in) :: prt
    flag = prt%left_handed
  end function field_data_is_left_handed

  elemental module function field_data_is_right_handed (prt) result (flag)
    logical :: flag
    class(field_data_t), intent(in) :: prt
    flag = prt%right_handed
  end function field_data_is_right_handed

  elemental module function field_data_has_antiparticle (prt) result (flag)
    logical :: flag
    class(field_data_t), intent(in) :: prt
    flag = prt%has_anti
  end function field_data_has_antiparticle

  elemental module function field_data_is_stable (prt, anti) result (flag)
    logical :: flag
    class(field_data_t), intent(in) :: prt
    logical, intent(in), optional :: anti
    if (present (anti)) then
       if (anti) then
          flag = prt%a_is_stable
       else
          flag = prt%p_is_stable
       end if
    else
       flag = prt%p_is_stable
    end if
  end function field_data_is_stable

  module subroutine field_data_get_decays (prt, decay, anti)
    class(field_data_t), intent(in) :: prt
    type(string_t), dimension(:), intent(out), allocatable :: decay
    logical, intent(in), optional :: anti
    if (present (anti)) then
       if (anti) then
          allocate (decay (size (prt%a_decay)), source = prt%a_decay)
       else
          allocate (decay (size (prt%p_decay)), source = prt%p_decay)
       end if
    else
       allocate (decay (size (prt%p_decay)), source = prt%p_decay)
    end if
  end subroutine field_data_get_decays

  elemental module function field_data_decays_isotropically &
       (prt, anti) result (flag)
    logical :: flag
    class(field_data_t), intent(in) :: prt
    logical, intent(in), optional :: anti
    if (present (anti)) then
       if (anti) then
          flag = prt%a_decays_isotropically
       else
          flag = prt%p_decays_isotropically
       end if
    else
       flag = prt%p_decays_isotropically
    end if
  end function field_data_decays_isotropically

  elemental module function field_data_decays_diagonal &
       (prt, anti) result (flag)
    logical :: flag
    class(field_data_t), intent(in) :: prt
    logical, intent(in), optional :: anti
    if (present (anti)) then
       if (anti) then
          flag = prt%a_decays_diagonal
       else
          flag = prt%p_decays_diagonal
       end if
    else
       flag = prt%p_decays_diagonal
    end if
  end function field_data_decays_diagonal

  elemental module function field_data_has_decay_helicity &
       (prt, anti) result (flag)
    logical :: flag
    class(field_data_t), intent(in) :: prt
    logical, intent(in), optional :: anti
    if (present (anti)) then
       if (anti) then
          flag = prt%a_has_decay_helicity
       else
          flag = prt%p_has_decay_helicity
       end if
    else
       flag = prt%p_has_decay_helicity
    end if
  end function field_data_has_decay_helicity

  elemental module function field_data_decay_helicity &
       (prt, anti) result (hel)
    integer :: hel
    class(field_data_t), intent(in) :: prt
    logical, intent(in), optional :: anti
    if (present (anti)) then
       if (anti) then
          hel = prt%a_decay_helicity
       else
          hel = prt%p_decay_helicity
       end if
    else
       hel = prt%p_decay_helicity
    end if
  end function field_data_decay_helicity

  elemental module function field_data_is_polarized (prt, anti) result (flag)
    logical :: flag
    class(field_data_t), intent(in) :: prt
    logical, intent(in), optional :: anti
    logical :: a
    if (present (anti)) then
       a = anti
    else
       a = .false.
    end if
    if (a) then
       flag = prt%a_polarized
    else
       flag = prt%p_polarized
    end if
  end function field_data_is_polarized

  pure module function field_data_get_longname (prt) result (name)
    type(string_t) :: name
    class(field_data_t), intent(in) :: prt
    name = prt%longname
  end function field_data_get_longname

  pure module function field_data_get_name (prt, is_antiparticle) result (name)
    type(string_t) :: name
    class(field_data_t), intent(in) :: prt
    logical, intent(in) :: is_antiparticle
    name = prt%longname
    if (is_antiparticle) then
       if (prt%has_anti) then
          if (allocated (prt%anti)) then
             if (size(prt%anti) > 0) name = prt%anti(1)
          end if
       else
          if (allocated (prt%name)) then
             if (size (prt%name) > 0) name = prt%name(1)
          end if
       end if
    else
       if (allocated (prt%name)) then
          if (size (prt%name) > 0) name = prt%name(1)
       end if
    end if
  end function field_data_get_name

  module subroutine field_data_get_name_array (prt, is_antiparticle, name)
    class(field_data_t), intent(in) :: prt
    logical, intent(in) :: is_antiparticle
    type(string_t), dimension(:), allocatable, intent(inout) :: name
    if (allocated (name))  deallocate (name)
    if (is_antiparticle) then
       if (prt%has_anti) then
          allocate (name (size (prt%anti)))
          name = prt%anti
       else
          allocate (name (0))
       end if
    else
       allocate (name (size (prt%name)))
       name = prt%name
    end if
  end subroutine field_data_get_name_array

  elemental module function field_data_get_tex_name &
       (prt, is_antiparticle) result (name)
    type(string_t) :: name
    class(field_data_t), intent(in) :: prt
    logical, intent(in) :: is_antiparticle
    if (is_antiparticle) then
       if (prt%has_anti) then
          name = prt%tex_anti
       else
          name = prt%tex_name
       end if
    else
       name = prt%tex_name
    end if
    if (name == "")  name = prt%get_name (is_antiparticle)
  end function field_data_get_tex_name

  module function field_data_matches_name &
       (field, name, is_antiparticle) result (flag)
    class(field_data_t), intent(in) :: field
    type(string_t), intent(in) :: name
    logical, intent(in) :: is_antiparticle
    logical :: flag
    if (is_antiparticle) then
       if (field%has_anti) then
          flag = any (name == field%anti)
       else
          flag = .false.
       end if
    else
       flag = name == field%longname .or. any (name == field%name)
    end if
  end function field_data_matches_name

  elemental module function field_data_get_spin_type (prt) result (type)
    integer :: type
    class(field_data_t), intent(in) :: prt
    type = prt%spin_type
  end function field_data_get_spin_type

  elemental module function field_data_get_multiplicity (prt) result (type)
    integer :: type
    class(field_data_t), intent(in) :: prt
    type = prt%multiplicity
  end function field_data_get_multiplicity

  elemental module function field_data_get_isospin_type (prt) result (type)
    integer :: type
    class(field_data_t), intent(in) :: prt
    type = prt%isospin_type
  end function field_data_get_isospin_type

  elemental module function field_data_get_charge_type (prt) result (type)
    integer :: type
    class(field_data_t), intent(in) :: prt
    type = prt%charge_type
  end function field_data_get_charge_type

  elemental module function field_data_get_color_type (prt) result (type)
    integer :: type
    class(field_data_t), intent(in) :: prt
    type = prt%color_type
  end function field_data_get_color_type

  elemental module function field_data_get_charge (prt) result (charge)
    real(default) :: charge
    class(field_data_t), intent(in) :: prt
    if (prt%charge_type /= 0) then
       charge = real (sign ((abs(prt%charge_type) - 1), &
             prt%charge_type), default) / 3
    else
       charge = 0
    end if
  end function field_data_get_charge

  elemental module function field_data_get_isospin (prt) result (isospin)
    real(default) :: isospin
    class(field_data_t), intent(in) :: prt
    if (prt%isospin_type /= 0) then
       isospin = real (sign (abs(prt%isospin_type) - 1, &
            prt%isospin_type), default) / 2
    else
       isospin = 0
    end if
  end function field_data_get_isospin

  elemental module function field_data_get_mass (prt) result (mass)
    real(default) :: mass
    class(field_data_t), intent(in) :: prt
    if (associated (prt%mass_val)) then
       mass = abs (prt%mass_val)
    else
       mass = 0
    end if
  end function field_data_get_mass

  elemental module function field_data_get_mass_sign (prt) result (sgn)
    integer :: sgn
    class(field_data_t), intent(in) :: prt
    if (associated (prt%mass_val)) then
       sgn = sign (1._default, prt%mass_val)
    else
       sgn = 0
    end if
  end function field_data_get_mass_sign

  elemental module function field_data_get_width (prt) result (width)
    real(default) :: width
    class(field_data_t), intent(in) :: prt
    if (associated (prt%width_val)) then
       width = prt%width_val
    else
       width = 0
    end if
  end function field_data_get_width

  module subroutine find_model (model, PDG, model_A, model_B)
    class(model_data_t), pointer, intent(out) :: model
    integer, intent(in) :: PDG
    class(model_data_t), intent(in), target :: model_A, model_B
    character(len=10) :: buffer
    if (model_A%test_field (PDG)) then
       model => model_A
    else if (model_B%test_field (PDG)) then
       model => model_B
    else
       call model_A%write ()
       call model_B%write ()
       write (buffer, "(I10)") PDG
       call msg_fatal ("Parton " // buffer // &
            " not found in the given model files")
    end if
  end subroutine find_model

  module subroutine vertex_write (vtx, unit)
    class(vertex_t), intent(in) :: vtx
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(3x,A)", advance="no")  "vertex"
    do i = 1, size (vtx%prt)
       if (associated (vtx%prt(i)%p)) then
          write (u, "(1x,A)", advance="no") &
               '"' // char (vtx%prt(i)%p%get_name (vtx%pdg(i) < 0)) &
                   // '"'
       else
          write (u, "(1x,I7)", advance="no") vtx%pdg(i)
       end if
    end do
    write (u, *)
  end subroutine vertex_write

  module subroutine vertex_init (vtx, pdg, model)
    class(vertex_t), intent(out) :: vtx
    integer, dimension(:), intent(in) :: pdg
    type(model_data_t), intent(in), target, optional :: model
    integer :: i
    allocate (vtx%pdg (size (pdg)))
    allocate (vtx%prt (size (pdg)))
    vtx%trilinear = size (pdg) == 3
    vtx%pdg = pdg
    if (present (model)) then
       do i = 1, size (pdg)
          vtx%prt(i)%p => model%get_field_ptr (pdg(i))
       end do
    end if
  end subroutine vertex_init

  module subroutine vertex_copy_from (vtx, old_vtx, new_model)
    class(vertex_t), intent(out) :: vtx
    class(vertex_t), intent(in) :: old_vtx
    type(model_data_t), intent(in), target, optional :: new_model
    call vtx%init (old_vtx%pdg, new_model)
  end subroutine vertex_copy_from

  module subroutine vertex_get_match (vtx, pdg1, pdg2)
    class(vertex_t), intent(in) :: vtx
    integer, intent(in) :: pdg1
    integer, dimension(:), allocatable, intent(out) :: pdg2
    integer :: i, j
    do i = 1, size (vtx%pdg)
       if (vtx%pdg(i) == pdg1) then
          allocate (pdg2 (size (vtx%pdg) - 1))
          do j = 1, i-1
             pdg2(j) = anti (j)
          end do
          do j = i, size (pdg2)
             pdg2(j) = anti (j+1)
          end do
          exit
       end if
    end do
  contains
    function anti (i) result (pdg)
      integer, intent(in) :: i
      integer :: pdg
      if (vtx%prt(i)%p%has_antiparticle ()) then
         pdg = - vtx%pdg(i)
      else
         pdg = vtx%pdg(i)
      end if
    end function anti
  end subroutine vertex_get_match

  module subroutine vertex_iterator_init (it, model, pdg, save_pdg_index)
    class(vertex_iterator_t), intent(out) :: it
    class(model_data_t), intent(in), target :: model
    integer, dimension(:), intent(in) :: pdg
    logical, intent(in) :: save_pdg_index
    it%model => model
    allocate (it%pdg (size (pdg)), source = pdg)
    it%save_pdg_index = save_pdg_index
  end subroutine vertex_iterator_init

  module subroutine vertex_iterator_get_next_match (it, pdg_match)
    class(vertex_iterator_t), intent(inout) :: it
    integer, dimension(:), allocatable, intent(out) :: pdg_match
    integer :: i, j
    do i = it%vertex_index + 1, size (it%model%vtx)
       do j = it%pdg_index + 1, size (it%pdg)
          call it%model%vtx(i)%get_match (it%pdg(j), pdg_match)
          if (it%save_pdg_index) then
             if (allocated (pdg_match) .and. j < size (it%pdg)) then
                it%pdg_index = j
                return
             else if (allocated (pdg_match) .and. j == size (it%pdg)) then
                it%vertex_index = i
                it%pdg_index = 0
                return
             end if
          else if (allocated (pdg_match)) then
             it%vertex_index = i
             return
          end if
       end do
    end do
    it%vertex_index = 0
    it%pdg_index = 0
  end subroutine vertex_iterator_get_next_match

  function vertex_table_size (n_vtx) result (n)
    integer(i32) :: n
    integer, intent(in) :: n_vtx
    integer :: i, s
    s = VERTEX_TABLE_SCALE_FACTOR * n_vtx
    n = 1
    do i = 1, 31
       n = ishft (n, 1)
       s = ishft (s,-1)
       if (s == 0)  exit
    end do
  end function vertex_table_size

  function hash2 (pdg1, pdg2)
    integer(i32) :: hash2
    integer, intent(in) :: pdg1, pdg2
    integer(i8), dimension(1) :: mold
    hash2 = hash (transfer ([pdg1, pdg2], mold))
  end function hash2

  module subroutine vertex_table_write (vt, unit)
    class(vertex_table_t), intent(in) :: vt
    integer, intent(in), optional :: unit
    integer :: u, i
    character(9) :: size_pdg3
    u = given_output_unit (unit)
    write (u, "(A)") "vertex hash table:"
    write (u, "(A,I7)") "  size = ", size (vt%entry)
    write (u, "(A,I7)") "  used = ", count (vt%entry%n /= 0)
    write (u, "(A,I7)") "  coll = ", vt%n_collisions
    do i = lbound (vt%entry, 1), ubound (vt%entry, 1)
       if (vt%entry(i)%n /= 0) then
          write (size_pdg3, "(I7)") size (vt%entry(i)%pdg3)
          write (u, "(A,1x,I7,1x,A,2(1x,I7),A," // &
               size_pdg3 // "(1x,I7))")  &
               "  ", i, ":", vt%entry(i)%pdg1, &
               vt%entry(i)%pdg2, "->", vt%entry(i)%pdg3
       end if
    end do
  end subroutine vertex_table_write

  module subroutine vertex_table_init (vt, prt, vtx)
    class(vertex_table_t), intent(out) :: vt
    type(field_data_t), dimension(:), intent(in) :: prt
    type(vertex_t), dimension(:), intent(in) :: vtx
    integer :: n_vtx, vt_size, i, p1, p2, p3
    integer, dimension(3) :: p
    n_vtx = size (vtx)
    vt_size = vertex_table_size (count (vtx%trilinear))
    vt%mask = vt_size - 1
    allocate (vt%entry (0:vt_size-1))
    do i = 1, n_vtx
       if (vtx(i)%trilinear) then
          p = vtx(i)%pdg
          p1 = p(1);  p2 = p(2)
          call create (hash2 (p1, p2))
          if (p(2) /= p(3)) then
             p2 = p(3)
             call create (hash2 (p1, p2))
          end if
          if (p(1) /= p(2)) then
             p1 = p(2);  p2 = p(1)
             call create (hash2 (p1, p2))
             if (p(1) /= p(3)) then
                p2 = p(3)
                call create (hash2 (p1, p2))
             end if
          end if
          if (p(1) /= p(3)) then
             p1 = p(3);  p2 = p(1)
             call create (hash2 (p1, p2))
             if (p(1) /= p(2)) then
                p2 = p(2)
                call create (hash2 (p1, p2))
             end if
          end if
       end if
    end do
    do i = 0, vt_size - 1
       allocate (vt%entry(i)%pdg3 (vt%entry(i)%n))
    end do
    vt%entry%n = 0
    do i = 1, n_vtx
       if (vtx(i)%trilinear) then
          p = vtx(i)%pdg
          p1 = p(1);  p2 = p(2);  p3 = p(3)
          call register (hash2 (p1, p2))
          if (p(2) /= p(3)) then
             p2 = p(3);  p3 = p(2)
             call register (hash2 (p1, p2))
          end if
          if (p(1) /= p(2)) then
             p1 = p(2);  p2 = p(1);  p3 = p(3)
             call register (hash2 (p1, p2))
             if (p(1) /= p(3)) then
                p2 = p(3);  p3 = p(1)
                call register (hash2 (p1, p2))
             end if
          end if
          if (p(1) /= p(3)) then
             p1 = p(3);  p2 = p(1);  p3 = p(2)
             call register (hash2 (p1, p2))
             if (p(1) /= p(2)) then
                p2 = p(2);  p3 = p(1)
                call register (hash2 (p1, p2))
             end if
          end if
       end if
    end do
  contains
    recursive subroutine create (hashval)
      integer(i32), intent(in) :: hashval
      integer :: h
      h = iand (hashval, vt%mask)
      if (vt%entry(h)%n == 0) then
         vt%entry(h)%pdg1 = p1
         vt%entry(h)%pdg2 = p2
         vt%entry(h)%n = 1
      else if (vt%entry(h)%pdg1 == p1 .and. vt%entry(h)%pdg2 == p2) then
         vt%entry(h)%n = vt%entry(h)%n + 1
      else
         vt%n_collisions = vt%n_collisions + 1
         call create (hashval + 1)
      end if
    end subroutine create
    recursive subroutine register (hashval)
      integer(i32), intent(in) :: hashval
      integer :: h
      h = iand (hashval, vt%mask)
      if (vt%entry(h)%pdg1 == p1 .and. vt%entry(h)%pdg2 == p2) then
         vt%entry(h)%n = vt%entry(h)%n + 1
         vt%entry(h)%pdg3(vt%entry(h)%n) = p3
      else
         call register (hashval + 1)
      end if
    end subroutine register
  end subroutine vertex_table_init

  module subroutine vertex_table_match (vt, pdg1, pdg2, pdg3)
    class(vertex_table_t), intent(in) :: vt
    integer, intent(in) :: pdg1, pdg2
    integer, dimension(:), allocatable, intent(out) :: pdg3
    call match (hash2 (pdg1, pdg2))
  contains
    recursive subroutine match (hashval)
      integer(i32), intent(in) :: hashval
      integer :: h
      h = iand (hashval, vt%mask)
      if (vt%entry(h)%n == 0) then
         allocate (pdg3 (0))
      else if (vt%entry(h)%pdg1 == pdg1 .and. vt%entry(h)%pdg2 == pdg2) then
         allocate (pdg3 (size (vt%entry(h)%pdg3)))
         pdg3 = vt%entry(h)%pdg3
      else
         call match (hashval + 1)
      end if
    end subroutine match
  end subroutine vertex_table_match

  module function vertex_table_check (vt, pdg1, pdg2, pdg3) result (flag)
    class(vertex_table_t), intent(in) :: vt
    integer, intent(in) :: pdg1, pdg2, pdg3
    logical :: flag
    flag = check (hash2 (pdg1, pdg2))
  contains
    recursive function check (hashval) result (flag)
      integer(i32), intent(in) :: hashval
      integer :: h
      logical :: flag
      h = iand (hashval, vt%mask)
      if (vt%entry(h)%n == 0) then
         flag = .false.
      else if (vt%entry(h)%pdg1 == pdg1 .and. vt%entry(h)%pdg2 == pdg2) then
         flag = any (vt%entry(h)%pdg3 == pdg3)
      else
         flag = check (hashval + 1)
      end if
    end function check
  end function vertex_table_check

  module subroutine model_data_final (model)
    class(model_data_t), intent(inout) :: model
    if (associated (model%par_real)) then
       deallocate (model%par_real)
    end if
    if (associated (model%par_complex)) then
       deallocate (model%par_complex)
    end if
  end subroutine model_data_final

  module subroutine model_data_write (model, unit, verbose, &
       show_md5sum, show_variables, show_parameters, &
       show_particles, show_vertices, show_scheme)
    class(model_data_t), intent(in) :: model
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    logical, intent(in), optional :: show_md5sum
    logical, intent(in), optional :: show_variables
    logical, intent(in), optional :: show_parameters
    logical, intent(in), optional :: show_particles
    logical, intent(in), optional :: show_vertices
    logical, intent(in), optional :: show_scheme
    logical :: show_sch, show_par, show_prt, show_vtx
    integer :: u, i
    u = given_output_unit (unit)
    show_sch = .false.;  if (present (show_scheme)) &
         show_sch = show_scheme
    show_par = .true.;  if (present (show_parameters)) &
         show_par = show_parameters
    show_prt = .true.;  if (present (show_particles)) &
         show_prt = show_particles
    show_vtx = .true.;  if (present (show_vertices)) &
         show_vtx = show_vertices
    if (show_sch) then
       write (u, "(3x,A,1X,I0)")  "scheme =", model%scheme
    end if
    if (show_par) then
       do i = 1, size (model%par_real)
          call model%par_real(i)%write (u)
          write (u, "(A)")
       end do
       do i = 1, size (model%par_complex)
          call model%par_complex(i)%write (u)
          write (u, "(A)")
       end do
    end if
    if (show_prt) then
       write (u, "(A)")
       call model%write_fields (u)
    end if
    if (show_vtx) then
       write (u, "(A)")
       call model%write_vertices (u, verbose)
    end if
  end subroutine model_data_write

  module subroutine model_data_init (model, name, &
       n_par_real, n_par_complex, n_field, n_vtx)
    class(model_data_t), intent(out) :: model
    type(string_t), intent(in) :: name
    integer, intent(in) :: n_par_real, n_par_complex
    integer, intent(in) :: n_field
    integer, intent(in) :: n_vtx
    model%name = name
    allocate (model%par_real (n_par_real))
    allocate (model%par_complex (n_par_complex))
    allocate (model%field (n_field))
    allocate (model%vtx (n_vtx))
  end subroutine model_data_init

  module subroutine model_data_set_scheme_num (model, scheme)
    class(model_data_t), intent(inout) :: model
    integer, intent(in) :: scheme
    model%scheme = scheme
  end subroutine model_data_set_scheme_num

  module subroutine model_data_freeze_fields (model)
    class(model_data_t), intent(inout) :: model
    call model%field%freeze ()
  end subroutine model_data_freeze_fields

  module subroutine model_data_copy (model, src)
    class(model_data_t), intent(inout), target :: model
    class(model_data_t), intent(in), target :: src
    class(modelpar_data_t), pointer :: data, src_data
    integer :: i
    model%scheme = src%scheme
    model%par_real = src%par_real
    model%par_complex = src%par_complex
    do i = 1, size (src%field)
       associate (field => model%field(i), src_field => src%field(i))
         call field%init (src_field%get_longname (), src_field%get_pdg ())
         call field%copy_from (src_field)
         src_data => src_field%mass_data
         if (associated (src_data)) then
            data => model%get_par_data_ptr (src_data%get_name ())
            call field%set (mass_data = data)
         end if
         src_data => src_field%width_data
         if (associated (src_data)) then
            data => model%get_par_data_ptr (src_data%get_name ())
            call field%set (width_data = data)
         end if
         call field%set_multiplicity ()
       end associate
    end do
    do i = 1, size (src%vtx)
       call model%vtx(i)%copy_from (src%vtx(i), model)
    end do
    call model%freeze_vertices ()
  end subroutine model_data_copy

  module function model_data_get_name (model) result (name)
    class(model_data_t), intent(in) :: model
    type(string_t) :: name
    name = model%name
  end function model_data_get_name

  module function model_data_get_scheme_num (model) result (scheme)
    class(model_data_t), intent(in) :: model
    integer :: scheme
    scheme = model%scheme
  end function model_data_get_scheme_num

  module function model_data_get_parameters_md5sum (model) result (par_md5sum)
    character(32) :: par_md5sum
    class(model_data_t), intent(in) :: model
    real(default), dimension(:), allocatable :: par
    type(field_data_t), pointer :: field
    integer :: unit, i
    allocate (par (model%get_n_real ()))
    call model%real_parameters_to_array (par)
    unit = free_unit ()
    open (unit, status="scratch", action="readwrite")
    if (model%scheme /= 0)  write (unit, "(I0)")  model%scheme
    write (unit, "(" // FMT_19 // ")")  par
    do i = 1, model%get_n_field ()
       field => model%get_field_ptr_by_index (i)
       if (.not. field%is_stable (.false.) .or. .not. field%is_stable (.true.) &
            .or. field%is_polarized (.false.) .or. field%is_polarized (.true.))&
            then
          write (unit, "(3x,A)") char (field%get_longname ())
          call field%write_decays (unit)
       end if
    end do
    rewind (unit)
    par_md5sum = md5sum (unit)
    close (unit)
  end function model_data_get_parameters_md5sum

  module function model_data_get_md5sum (model) result (md5sum)
    class(model_data_t), intent(in) :: model
    character(32) :: md5sum
    md5sum = model%get_parameters_md5sum ()
  end function model_data_get_md5sum

  module subroutine model_data_init_par_real (model, i, name, value)
    class(model_data_t), intent(inout) :: model
    integer, intent(in) :: i
    type(string_t), intent(in) :: name
    real(default), intent(in) :: value
    call model%par_real(i)%init (name, value)
  end subroutine model_data_init_par_real

  module subroutine model_data_init_par_complex (model, i, name, value)
    class(model_data_t), intent(inout) :: model
    integer, intent(in) :: i
    type(string_t), intent(in) :: name
    complex(default), intent(in) :: value
    call model%par_complex(i)%init (name, value)
  end subroutine model_data_init_par_complex

  module function model_data_get_n_real (model) result (n)
    class(model_data_t), intent(in) :: model
    integer :: n
    n = size (model%par_real)
  end function model_data_get_n_real

  module function model_data_get_n_complex (model) result (n)
    class(model_data_t), intent(in) :: model
    integer :: n
    n = size (model%par_complex)
  end function model_data_get_n_complex

  module subroutine model_data_real_par_to_array (model, array)
    class(model_data_t), intent(in) :: model
    real(default), dimension(:), intent(inout) :: array
    array = model%par_real%get_real ()
  end subroutine model_data_real_par_to_array

  module subroutine model_data_complex_par_to_array (model, array)
    class(model_data_t), intent(in) :: model
    complex(default), dimension(:), intent(inout) :: array
    array = model%par_complex%get_complex ()
  end subroutine model_data_complex_par_to_array

  module subroutine model_data_real_par_from_array (model, array)
    class(model_data_t), intent(inout) :: model
    real(default), dimension(:), intent(in) :: array
    model%par_real = array
  end subroutine model_data_real_par_from_array

  module subroutine model_data_complex_par_from_array (model, array)
    class(model_data_t), intent(inout) :: model
    complex(default), dimension(:), intent(in) :: array
    model%par_complex = array
  end subroutine model_data_complex_par_from_array

  module subroutine model_data_real_par_to_c_array (model, array)
    class(model_data_t), intent(in) :: model
    real(c_default_float), dimension(:), intent(inout) :: array
    array = model%par_real%get_real ()
  end subroutine model_data_real_par_to_c_array

  module subroutine model_data_real_par_from_c_array (model, array)
    class(model_data_t), intent(inout) :: model
    real(c_default_float), dimension(:), intent(in) :: array
    model%par_real = real (array, default)
  end subroutine model_data_real_par_from_c_array

  module function model_data_get_par_real_ptr_index (model, i) result (ptr)
    class(model_data_t), intent(inout) :: model
    integer, intent(in) :: i
    class(modelpar_data_t), pointer :: ptr
    ptr => model%par_real(i)
  end function model_data_get_par_real_ptr_index

  module function model_data_get_par_complex_ptr_index (model, i) result (ptr)
    class(model_data_t), intent(inout) :: model
    integer, intent(in) :: i
    class(modelpar_data_t), pointer :: ptr
    ptr => model%par_complex(i)
  end function model_data_get_par_complex_ptr_index

  module function model_data_get_par_data_ptr_name (model, name) result (ptr)
    class(model_data_t), intent(in) :: model
    type(string_t), intent(in) :: name
    class(modelpar_data_t), pointer :: ptr
    integer :: i
    do i = 1, size (model%par_real)
       if (model%par_real(i)%name == name) then
          ptr => model%par_real(i)
          return
       end if
    end do
    do i = 1, size (model%par_complex)
       if (model%par_complex(i)%name == name) then
          ptr => model%par_complex(i)
          return
       end if
    end do
    ptr => null ()
  end function model_data_get_par_data_ptr_name

  module function model_data_get_par_real_value (model, name) result (value)
    class(model_data_t), intent(in) :: model
    type(string_t), intent(in) :: name
    class(modelpar_data_t), pointer :: par
    real(default) :: value
    par => model%get_par_data_ptr (name)
    value = par%get_real ()
  end function model_data_get_par_real_value

  module function model_data_get_par_complex_value (model, name) result (value)
    class(model_data_t), intent(in) :: model
    type(string_t), intent(in) :: name
    class(modelpar_data_t), pointer :: par
    complex(default) :: value
    par => model%get_par_data_ptr (name)
    value = par%get_complex ()
  end function model_data_get_par_complex_value

  module subroutine model_data_set_par_real (model, name, value)
    class(model_data_t), intent(inout) :: model
    type(string_t), intent(in) :: name
    real(default), intent(in) :: value
    class(modelpar_data_t), pointer :: par
    par => model%get_par_data_ptr (name)
    par = value
  end subroutine model_data_set_par_real

  module subroutine model_data_set_par_complex (model, name, value)
    class(model_data_t), intent(inout) :: model
    type(string_t), intent(in) :: name
    complex(default), intent(in) :: value
    class(modelpar_data_t), pointer :: par
    par => model%get_par_data_ptr (name)
    par = value
  end subroutine model_data_set_par_complex

  module subroutine model_data_write_fields (model, unit)
    class(model_data_t), intent(in) :: model
    integer, intent(in), optional :: unit
    integer :: i
    do i = 1, size (model%field)
       call model%field(i)%write (unit)
    end do
  end subroutine model_data_write_fields

  module function model_data_get_n_field (model) result (n)
    class(model_data_t), intent(in) :: model
    integer :: n
    n = size (model%field)
  end function model_data_get_n_field

  module function model_data_get_field_pdg_index (model, i) result (pdg)
    class(model_data_t), intent(in) :: model
    integer, intent(in) :: i
    integer :: pdg
    pdg = model%field(i)%get_pdg ()
  end function model_data_get_field_pdg_index

  module function model_data_get_field_pdg_name &
       (model, name, check) result (pdg)
    class(model_data_t), intent(in) :: model
    type(string_t), intent(in) :: name
    logical, intent(in), optional :: check
    integer :: pdg
    integer :: i
    do i = 1, size (model%field)
       associate (field => model%field(i))
         if (field%matches_name (name, .false.)) then
            pdg = field%get_pdg ()
            return
         else if (field%matches_name (name, .true.)) then
            pdg = - field%get_pdg ()
            return
         end if
       end associate
    end do
    pdg = 0
    call model%field_error (check, name)
  end function model_data_get_field_pdg_name

  module subroutine model_data_get_all_pdg (model, pdg)
    class(model_data_t), intent(in) :: model
    integer, dimension(:), allocatable, intent(inout) :: pdg
    integer :: n0, n1, i, k
    n0 = size (model%field)
    n1 = n0 + count (model%field%has_antiparticle ())
    allocate (pdg (n1))
    pdg(1:n0) = model%field%get_pdg ()
    k = n0
    do i = 1, size (model%field)
       associate (field => model%field(i))
         if (field%has_antiparticle ()) then
            k = k + 1
            pdg(k) = - field%get_pdg ()
         end if
       end associate
    end do
  end subroutine model_data_get_all_pdg

  module function model_data_get_field_array_ptr (model) result (ptr)
    class(model_data_t), intent(in), target :: model
    type(field_data_t), dimension(:), pointer :: ptr
    ptr => model%field
  end function model_data_get_field_array_ptr

  module function model_data_get_field_ptr_name &
       (model, name, check) result (ptr)
    class(model_data_t), intent(in), target :: model
    type(string_t), intent(in) :: name
    logical, intent(in), optional :: check
    type(field_data_t), pointer :: ptr
    integer :: i
    do i = 1, size (model%field)
       if (model%field(i)%matches_name (name, .false.)) then
          ptr => model%field(i)
          return
       else if (model%field(i)%matches_name (name, .true.)) then
          ptr => model%field(i)
          return
       end if
    end do
    ptr => null ()
    call model%field_error (check, name)
  end function model_data_get_field_ptr_name

  module function model_data_get_field_ptr_pdg (model, pdg, check) result (ptr)
    class(model_data_t), intent(in), target :: model
    integer, intent(in) :: pdg
    logical, intent(in), optional :: check
    type(field_data_t), pointer :: ptr
    integer :: i, pdg_abs
    if (pdg == 0) then
       ptr => null ()
       return
    end if
    pdg_abs = abs (pdg)
    do i = 1, size (model%field)
       if (abs(model%field(i)%get_pdg ()) == pdg_abs) then
          ptr => model%field(i)
          return
       end if
    end do
    ptr => null ()
    call model%field_error (check, pdg=pdg)
  end function model_data_get_field_ptr_pdg

  module function model_data_get_field_ptr_index (model, i) result (ptr)
    class(model_data_t), intent(in), target :: model
    integer, intent(in) :: i
    type(field_data_t), pointer :: ptr
    ptr => model%field(i)
  end function model_data_get_field_ptr_index

  module function model_data_test_field_pdg (model, pdg, check) result (exist)
    class(model_data_t), intent(in), target :: model
    integer, intent(in) :: pdg
    logical, intent(in), optional :: check
    logical :: exist
    exist = associated (model%get_field_ptr (pdg, check))
  end function model_data_test_field_pdg

  module subroutine model_data_field_error (model, check, name, pdg)
    class(model_data_t), intent(in) :: model
    logical, intent(in), optional :: check
    type(string_t), intent(in), optional :: name
    integer, intent(in), optional :: pdg
    if (present (check)) then
       if (check) then
          if (present (name)) then
             write (msg_buffer, "(A,1x,A,1x,A,1x,A)") &
                  "No particle with name", char (name), &
                  "is contained in model", char (model%name)
          else if (present (pdg)) then
             write (msg_buffer, "(A,1x,I0,1x,A,1x,A)") &
                  "No particle with PDG code", pdg, &
                  "is contained in model", char (model%name)
          else
             write (msg_buffer, "(A,1x,A,1x,A)") &
                  "Particle missing", &
                  "in model", char (model%name)
          end if
          call msg_fatal ()
       end if
    end if
  end subroutine model_data_field_error

  module subroutine model_data_set_field_mass_pdg (model, pdg, value)
    class(model_data_t), intent(inout) :: model
    integer, intent(in) :: pdg
    real(default), intent(in) :: value
    type(field_data_t), pointer :: field
    field => model%get_field_ptr (pdg, check = .true.)
    call field%set_mass (value)
  end subroutine model_data_set_field_mass_pdg

  module subroutine model_data_set_field_width_pdg (model, pdg, value)
    class(model_data_t), intent(inout) :: model
    integer, intent(in) :: pdg
    real(default), intent(in) :: value
    type(field_data_t), pointer :: field
    field => model%get_field_ptr (pdg, check = .true.)
    call field%set_width (value)
  end subroutine model_data_set_field_width_pdg

  module subroutine model_data_set_unstable &
       (model, pdg, decay, isotropic, diagonal, decay_helicity)
    class(model_data_t), intent(inout), target :: model
    integer, intent(in) :: pdg
    type(string_t), dimension(:), intent(in) :: decay
    logical, intent(in), optional :: isotropic, diagonal
    integer, intent(in), optional :: decay_helicity
    type(field_data_t), pointer :: field
    field => model%get_field_ptr (pdg)
    if (pdg > 0) then
       call field%set ( &
            p_is_stable = .false., p_decay = decay, &
            p_decays_isotropically = isotropic, &
            p_decays_diagonal = diagonal, &
            p_decay_helicity = decay_helicity)
    else
       call field%set ( &
            a_is_stable = .false., a_decay = decay, &
            a_decays_isotropically = isotropic, &
            a_decays_diagonal = diagonal, &
            a_decay_helicity = decay_helicity)
    end if
  end subroutine model_data_set_unstable

  module subroutine model_data_set_stable (model, pdg)
    class(model_data_t), intent(inout), target :: model
    integer, intent(in) :: pdg
    type(field_data_t), pointer :: field
    field => model%get_field_ptr (pdg)
    if (pdg > 0) then
       call field%set (p_is_stable = .true.)
    else
       call field%set (a_is_stable = .true.)
    end if
  end subroutine model_data_set_stable

  module subroutine model_data_set_polarized (model, pdg)
    class(model_data_t), intent(inout), target :: model
    integer, intent(in) :: pdg
    type(field_data_t), pointer :: field
    field => model%get_field_ptr (pdg)
    if (pdg > 0) then
       call field%set (p_polarized = .true.)
    else
       call field%set (a_polarized = .true.)
    end if
  end subroutine model_data_set_polarized

  module subroutine model_data_set_unpolarized (model, pdg)
    class(model_data_t), intent(inout), target :: model
    integer, intent(in) :: pdg
    type(field_data_t), pointer :: field
    field => model%get_field_ptr (pdg)
    if (pdg > 0) then
       call field%set (p_polarized = .false.)
    else
       call field%set (a_polarized = .false.)
    end if
  end subroutine model_data_set_unpolarized

  module subroutine model_clear_unstable (model)
    class(model_data_t), intent(inout), target :: model
    integer :: i
    type(field_data_t), pointer :: field
    do i = 1, model%get_n_field ()
       field => model%get_field_ptr_by_index (i)
       call field%set (p_is_stable = .true.)
       if (field%has_antiparticle ()) then
          call field%set (a_is_stable = .true.)
       end if
    end do
  end subroutine model_clear_unstable

  module subroutine model_clear_polarized (model)
    class(model_data_t), intent(inout), target :: model
    integer :: i
    type(field_data_t), pointer :: field
    do i = 1, model%get_n_field ()
       field => model%get_field_ptr_by_index (i)
       call field%set (p_polarized = .false.)
       if (field%has_antiparticle ()) then
          call field%set (a_polarized = .false.)
       end if
    end do
  end subroutine model_clear_polarized

  module subroutine model_data_write_vertices (model, unit, verbose)
    class(model_data_t), intent(in) :: model
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: verbose
    integer :: i, u
    u = given_output_unit (unit)
    do i = 1, size (model%vtx)
       call vertex_write (model%vtx(i), unit)
    end do
    if (present (verbose)) then
       if (verbose) then
          write (u, *)
          call vertex_table_write (model%vt, unit)
       end if
    end if
  end subroutine model_data_write_vertices

  module subroutine model_data_set_vertex_pdg (model, i, pdg)
    class(model_data_t), intent(inout), target :: model
    integer, intent(in) :: i
    integer, dimension(:), intent(in) :: pdg
    call vertex_init (model%vtx(i), pdg, model)
  end subroutine model_data_set_vertex_pdg

  module subroutine model_data_set_vertex_names (model, i, name)
    class(model_data_t), intent(inout), target :: model
    integer, intent(in) :: i
    type(string_t), dimension(:), intent(in) :: name
    integer, dimension(size(name)) :: pdg
    integer :: j
    do j = 1, size (name)
       pdg(j) = model%get_pdg (name(j))
    end do
    call model%set_vertex (i, pdg)
  end subroutine model_data_set_vertex_names

  module subroutine model_data_freeze_vertices (model)
    class(model_data_t), intent(inout) :: model
    call model%vt%init (model%field, model%vtx)
  end subroutine model_data_freeze_vertices

  module function model_data_get_n_vtx (model) result (n)
    class(model_data_t), intent(in) :: model
    integer :: n
    n = size (model%vtx)
  end function model_data_get_n_vtx

  module subroutine model_data_match_vertex (model, pdg1, pdg2, pdg3)
    class(model_data_t), intent(in) :: model
    integer, intent(in) :: pdg1, pdg2
    integer, dimension(:), allocatable, intent(out) :: pdg3
    call model%vt%match (pdg1, pdg2, pdg3)
  end subroutine model_data_match_vertex

  module function model_data_check_vertex &
       (model, pdg1, pdg2, pdg3) result (flag)
    logical :: flag
    class(model_data_t), intent(in) :: model
    integer, intent(in) :: pdg1, pdg2, pdg3
    flag = model%vt%check (pdg1, pdg2, pdg3)
  end function model_data_check_vertex

  module subroutine model_data_init_test (model)
    class(model_data_t), intent(out) :: model
    type(field_data_t), pointer :: field
    integer, parameter :: n_real = 4
    integer, parameter :: n_field = 2
    integer, parameter :: n_vertex = 2
    integer :: i
    call model%init (var_str ("Test"), &
         n_real, 0, n_field, n_vertex)
    i = 0
    i = i + 1
    call model%init_par (i, var_str ("gy"), 1._default)
    i = i + 1
    call model%init_par (i, var_str ("ms"), 125._default)
    i = i + 1
    call model%init_par (i, var_str ("ff"), 1.5_default)
    i = i + 1
    call model%init_par (i, var_str ("mf"), 1.5_default * 125._default)
    i = 0
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("SCALAR"), 25)
    call field%set (spin_type=1)
    call field%set (mass_data=model%get_par_real_ptr (2))
    call field%set (name = [var_str ("s")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("FERMION"), 6)
    call field%set (spin_type=2)
    call field%set (mass_data=model%get_par_real_ptr (4))
    call field%set (name = [var_str ("f")], anti = [var_str ("fbar")])
    call model%freeze_fields ()
    i = 0
    i = i + 1
    call model%set_vertex (i, [var_str ("fbar"), var_str ("f"), var_str ("s")])
    i = i + 1
    call model%set_vertex (i, [var_str ("s"), var_str ("s"), var_str ("s")])
    call model%freeze_vertices ()
  end subroutine model_data_init_test

  module subroutine model_data_init_qed_test (model)
    class(model_data_t), intent(out) :: model
    type(field_data_t), pointer :: field
    integer, parameter :: n_real = 1
    integer, parameter :: n_field = 2
    integer :: i
    call model%init (var_str ("QED_test"), &
         n_real, 0, n_field, 0)
    i = 0
    i = i + 1
    call model%init_par (i, var_str ("me"), 0.000510997_default)
    i = 0
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("E_LEPTON"), 11)
    call field%set (spin_type=2, charge_type=-4)
    call field%set (mass_data=model%get_par_real_ptr (1))
    call field%set (name = [var_str ("e-")], anti = [var_str ("e+")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("PHOTON"), 22)
    call field%set (spin_type=3)
    call field%set (name = [var_str ("A")])
    call model%freeze_fields ()
    call model%freeze_vertices ()
  end subroutine model_data_init_qed_test

  module subroutine model_data_init_sm_test (model)
    class(model_data_t), intent(out) :: model
    type(field_data_t), pointer :: field
    integer, parameter :: n_real = 11
    integer, parameter :: n_field = 19
    integer, parameter :: n_vtx = 9
    integer :: i
    call model%init (var_str ("SM_test"), &
         n_real, 0, n_field, n_vtx)
    i = 0
    i = i + 1
    call model%init_par (i, var_str ("mZ"), 91.1882_default)
    i = i + 1
    call model%init_par (i, var_str ("mW"), 80.419_default)
    i = i + 1
    call model%init_par (i, var_str ("me"), 0.000510997_default)
    i = i + 1
    call model%init_par (i, var_str ("mmu"), 0.105658389_default)
    i = i + 1
    call model%init_par (i, var_str ("mb"), 4.2_default)
    i = i + 1
    call model%init_par (i, var_str ("mtop"), 173.1_default)
    i = i + 1
    call model%init_par (i, var_str ("wZ"), 2.443_default)
    i = i + 1
    call model%init_par (i, var_str ("wW"), 2.049_default)
    i = i + 1
    call model%init_par (i, var_str ("ee"), 0.3079561542961_default)
    i = i + 1
    call model%init_par (i, var_str ("cw"), 8.819013863636E-01_default)
    i = i + 1
    call model%init_par (i, var_str ("sw"), 4.714339240339E-01_default)
    i = 0
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("D_QUARK"), 1)
    call field%set (spin_type=2, color_type=3, charge_type=-2, isospin_type=-2)
    call field%set (name = [var_str ("d")], anti = [var_str ("dbar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("U_QUARK"), 2)
    call field%set (spin_type=2, color_type=3, charge_type=3, isospin_type=2)
    call field%set (name = [var_str ("u")], anti = [var_str ("ubar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("S_QUARK"), 3)
    call field%set (spin_type=2, color_type=3, charge_type=-2, isospin_type=-2)
    call field%set (name = [var_str ("s")], anti = [var_str ("sbar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("C_QUARK"), 4)
    call field%set (spin_type=2, color_type=3, charge_type=3, isospin_type=2)
    call field%set (name = [var_str ("c")], anti = [var_str ("cbar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("B_QUARK"), 5)
    call field%set (spin_type=2, color_type=3, charge_type=-2, isospin_type=-2)
    call field%set (mass_data=model%get_par_real_ptr (5))
    call field%set (name = [var_str ("b")], anti = [var_str ("bbar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("T_QUARK"), 6)
    call field%set (spin_type=2, color_type=3, charge_type=3, isospin_type=2)
    call field%set (mass_data=model%get_par_real_ptr (6))
    call field%set (name = [var_str ("t")], anti = [var_str ("tbar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("E_LEPTON"), 11)
    call field%set (spin_type=2)
    call field%set (mass_data=model%get_par_real_ptr (3))
    call field%set (name = [var_str ("e-")], anti = [var_str ("e+")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("E_NEUTRINO"), 12)
    call field%set (spin_type=2, is_left_handed=.true.)
    call field%set (name = [var_str ("nue")], anti = [var_str ("nuebar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("MU_LEPTON"), 13)
    call field%set (spin_type=2)
    call field%set (mass_data=model%get_par_real_ptr (4))
    call field%set (name = [var_str ("mu-")], anti = [var_str ("mu+")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("MU_NEUTRINO"), 14)
    call field%set (spin_type=2, is_left_handed=.true.)
    call field%set (name = [var_str ("numu")], anti = [var_str ("numubar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("GLUON"), 21)
    call field%set (spin_type=3, color_type=8)
    call field%set (name = [var_str ("gl")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("PHOTON"), 22)
    call field%set (spin_type=3)
    call field%set (name = [var_str ("A")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("Z_BOSON"), 23)
    call field%set (spin_type=3)
    call field%set (mass_data=model%get_par_real_ptr (1))
    call field%set (width_data=model%get_par_real_ptr (7))
    call field%set (name = [var_str ("Z")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("W_BOSON"), 24)
    call field%set (spin_type=3)
    call field%set (mass_data=model%get_par_real_ptr (2))
    call field%set (width_data=model%get_par_real_ptr (8))
    call field%set (name = [var_str ("W+")], anti = [var_str ("W-")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("HIGGS"), 25)
    call field%set (spin_type=1)
!    call field%set (mass_data=model%get_par_real_ptr (2))
!    call field%set (width_data=model%get_par_real_ptr (8))
    call field%set (name = [var_str ("H")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("PROTON"), 2212)
    call field%set (spin_type=2)
    call field%set (name = [var_str ("p")], anti = [var_str ("pbar")])
!    call field%set (mass_data=model%get_par_real_ptr (12))
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("HADRON_REMNANT_SINGLET"), 91)
    call field%set (color_type=1)
    call field%set (name = [var_str ("hr1")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("HADRON_REMNANT_TRIPLET"), 92)
    call field%set (color_type=3)
    call field%set (name = [var_str ("hr3")], anti = [var_str ("hr3bar")])
    i = i + 1
    field => model%get_field_ptr_by_index (i)
    call field%init (var_str ("HADRON_REMNANT_OCTET"), 93)
    call field%set (color_type=8)
    call field%set (name = [var_str ("hr8")])
    call model%freeze_fields ()
    i = 0
    i = i + 1
    call model%set_vertex (i, [var_str ("dbar"), var_str ("d"), var_str ("A")])
    i = i + 1
    call model%set_vertex (i, [var_str ("ubar"), var_str ("u"), var_str ("A")])
    i = i + 1
    call model%set_vertex (i, [var_str ("gl"), var_str ("gl"), var_str ("gl")])
    i = i + 1
    call model%set_vertex (i, [var_str ("dbar"), var_str ("d"), var_str ("gl")])
    i = i + 1
    call model%set_vertex (i, [var_str ("ubar"), var_str ("u"), var_str ("gl")])
    i = i + 1
    call model%set_vertex (i, [var_str ("dbar"), var_str ("d"), var_str ("Z")])
    i = i + 1
    call model%set_vertex (i, [var_str ("ubar"), var_str ("u"), var_str ("Z")])
    i = i + 1
    call model%set_vertex (i, [var_str ("ubar"), var_str ("d"), var_str ("W+")])
    i = i + 1
    call model%set_vertex (i, [var_str ("dbar"), var_str ("u"), var_str ("W-")])
    call model%freeze_vertices ()
  end subroutine model_data_init_sm_test


end submodule model_data_s

