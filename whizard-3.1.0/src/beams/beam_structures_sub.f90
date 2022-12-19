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

submodule (beam_structures) beam_structures_s

  use io_units
  use format_defs, only: FMT_19
  use diagnostics

  implicit none

contains

  module function beam_structure_entry_to_string (object) result (string)
    class(beam_structure_entry_t), intent(in) :: object
    type(string_t) :: string
    if (object%is_valid) then
       string = object%name
    else
       string = "none"
    end if
  end function beam_structure_entry_to_string

  module subroutine beam_structure_final_sf (object)
    class(beam_structure_t), intent(inout) :: object
    if (allocated (object%prt))  deallocate (object%prt)
    if (allocated (object%record))  deallocate (object%record)
    object%n_beam = 0
  end subroutine beam_structure_final_sf

  module subroutine beam_structure_write (object, unit)
    class(beam_structure_t), intent(in) :: object
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(1x,A,A)")  "Beam structure: ", char (object%to_string ())
    if (allocated (object%smatrix)) then
       do i = 1, size (object%smatrix)
          write (u, "(3x,A,I0,A)") "polarization (beam ", i, "):"
          call object%smatrix(i)%write (u, indent=2)
       end do
    end if
    if (allocated (object%pol_f)) then
       write (u, "(3x,A,F10.7,:,',',F10.7)")  "polarization degree =", &
            object%pol_f
    end if
    if (allocated (object%p)) then
       write (u, "(3x,A," // FMT_19 // ",:,','," // FMT_19 // &
            ")")  "momentum =", object%p
    end if
    if (allocated (object%theta)) then
       write (u, "(3x,A," // FMT_19 // ",:,','," // FMT_19 // &
            ")")  "angle th =", object%theta
    end if
    if (allocated (object%phi)) then
       write (u, "(3x,A," // FMT_19 // ",:,','," // FMT_19 // &
            ")")  "angle ph =", object%phi
    end if
  end subroutine beam_structure_write

  module function beam_structure_to_string (object, sf_only) result (string)
    class(beam_structure_t), intent(in) :: object
    logical, intent(in), optional :: sf_only
    type(string_t) :: string
    integer :: i, j
    logical :: with_beams
    with_beams = .true.;  if (present (sf_only))  with_beams = .not. sf_only
    select case (object%n_beam)
    case (1)
       if (with_beams) then
          string = object%prt(1)
       else
          string = ""
       end if
    case (2)
       if (with_beams) then
          string = object%prt(1) // ", " // object%prt(2)
       else
          string = ""
       end if
       if (allocated (object%record)) then
          if (size (object%record) > 0) then
             if (with_beams)  string = string // " => "
             do i = 1, size (object%record)
                if (i > 1)  string = string // " => "
                do j = 1, size (object%record(i)%entry)
                   if (j > 1)  string = string // ", "
                   string = string // object%record(i)%entry(j)%to_string ()
                end do
             end do
          end if
       end if
    case default
       string = "[any particles]"
    end select
  end function beam_structure_to_string

  module subroutine beam_structure_init_sf (beam_structure, prt, dim_array)
    class(beam_structure_t), intent(inout) :: beam_structure
    type(string_t), dimension(:), intent(in) :: prt
    integer, dimension(:), intent(in), optional :: dim_array
    integer :: i
    call beam_structure%final_sf ()
    beam_structure%n_beam = size (prt)
    allocate (beam_structure%prt (size (prt)))
    beam_structure%prt = prt
    if (present (dim_array)) then
       allocate (beam_structure%record (size (dim_array)))
       do i = 1, size (dim_array)
          allocate (beam_structure%record(i)%entry (dim_array(i)))
       end do
    else
       allocate (beam_structure%record (0))
    end if
  end subroutine beam_structure_init_sf

  module subroutine beam_structure_set_sf (beam_structure, i, j, name)
    class(beam_structure_t), intent(inout) :: beam_structure
    integer, intent(in) :: i, j
    type(string_t), intent(in) :: name
    associate (entry => beam_structure%record(i)%entry(j))
      entry%name = name
      entry%is_valid = .true.
    end associate
  end subroutine beam_structure_set_sf

  module subroutine beam_structure_expand (beam_structure, strfun_mode)
    class(beam_structure_t), intent(inout) :: beam_structure
    procedure(strfun_mode_fun) :: strfun_mode
    type(beam_structure_record_t), dimension(:), allocatable :: new
    integer :: n_record, i, j
    if (.not. allocated (beam_structure%record))  return
    do i = 1, size (beam_structure%record)
       associate (entry => beam_structure%record(i)%entry)
         do j = 1, size (entry)
            select case (strfun_mode (entry(j)%name))
            case (0);  entry(j)%is_valid = .false.
            end select
         end do
       end associate
    end do
    n_record = 0
    do i = 1, size (beam_structure%record)
       associate (entry => beam_structure%record(i)%entry)
         select case (size (entry))
         case (1)
            if (entry(1)%is_valid) then
               select case (strfun_mode (entry(1)%name))
               case (1);  n_record = n_record + 2
               case (2);  n_record = n_record + 1
               end select
            end if
         case (2)
            do j = 1, 2
               if (entry(j)%is_valid) then
                  select case (strfun_mode (entry(j)%name))
                  case (1);  n_record = n_record + 1
                  case (2)
                     call beam_structure%write ()
                     call msg_fatal ("Pair spectrum used as &
                          &single-particle structure function")
                  end select
               end if
            end do
         end select
       end associate
    end do
    allocate (new (n_record))
    n_record = 0
    do i = 1, size (beam_structure%record)
       associate (entry => beam_structure%record(i)%entry)
         select case (size (entry))
         case (1)
            if (entry(1)%is_valid) then
               select case (strfun_mode (entry(1)%name))
               case (1)
                  n_record = n_record + 1
                  allocate (new(n_record)%entry (2))
                  new(n_record)%entry(1) = entry(1)
                  n_record = n_record + 1
                  allocate (new(n_record)%entry (2))
                  new(n_record)%entry(2) = entry(1)
               case (2)
                  n_record = n_record + 1
                  allocate (new(n_record)%entry (1))
                  new(n_record)%entry(1) = entry(1)
               end select
            end if
         case (2)
            do j = 1, 2
               if (entry(j)%is_valid) then
                  n_record = n_record + 1
                  allocate (new(n_record)%entry (2))
                  new(n_record)%entry(j) = entry(j)
               end if
            end do
         end select
       end associate
    end do
    call move_alloc (from = new, to = beam_structure%record)
  end subroutine beam_structure_expand

  module subroutine beam_structure_final_pol (beam_structure)
    class(beam_structure_t), intent(inout) :: beam_structure
    if (allocated (beam_structure%smatrix))  deallocate (beam_structure%smatrix)
    if (allocated (beam_structure%pol_f))  deallocate (beam_structure%pol_f)
  end subroutine beam_structure_final_pol

  module subroutine beam_structure_init_pol (beam_structure, n)
    class(beam_structure_t), intent(inout) :: beam_structure
    integer, intent(in) :: n
    if (allocated (beam_structure%smatrix))  deallocate (beam_structure%smatrix)
    allocate (beam_structure%smatrix (n))
    if (.not. allocated (beam_structure%pol_f)) &
         allocate (beam_structure%pol_f (n), source = 1._default)
  end subroutine beam_structure_init_pol

  elemental module function beam_structure_has_polarized_beams &
       (beam_structure) result (pol)
    logical :: pol
    class(beam_structure_t), intent(in) :: beam_structure
    if (allocated (beam_structure%pol_f)) then
       pol = any (beam_structure%pol_f /= 0)
    else
       pol = .false.
    end if
  end function beam_structure_has_polarized_beams

  module subroutine beam_structure_set_smatrix (beam_structure, i, smatrix)
    class(beam_structure_t), intent(inout) :: beam_structure
    integer, intent(in) :: i
    type(smatrix_t), intent(in) :: smatrix
    beam_structure%smatrix(i) = smatrix
  end subroutine beam_structure_set_smatrix

  module subroutine beam_structure_init_smatrix (beam_structure, i, n_entry)
    class(beam_structure_t), intent(inout) :: beam_structure
    integer, intent(in) :: i
    integer, intent(in) :: n_entry
    call beam_structure%smatrix(i)%init (2, n_entry)
  end subroutine beam_structure_init_smatrix

  module subroutine beam_structure_set_sentry &
       (beam_structure, i, i_entry, index, value)
    class(beam_structure_t), intent(inout) :: beam_structure
    integer, intent(in) :: i
    integer, intent(in) :: i_entry
    integer, dimension(:), intent(in) :: index
    complex(default), intent(in) :: value
    call beam_structure%smatrix(i)%set_entry (i_entry, index, value)
  end subroutine beam_structure_set_sentry

  module subroutine beam_structure_set_pol_f (beam_structure, f)
    class(beam_structure_t), intent(inout) :: beam_structure
    real(default), dimension(:), intent(in) :: f
    if (allocated (beam_structure%pol_f))  deallocate (beam_structure%pol_f)
    allocate (beam_structure%pol_f (size (f)), source = f)
  end subroutine beam_structure_set_pol_f

  module subroutine beam_structure_final_mom (beam_structure)
    class(beam_structure_t), intent(inout) :: beam_structure
    if (allocated (beam_structure%p))  deallocate (beam_structure%p)
    if (allocated (beam_structure%theta))  deallocate (beam_structure%theta)
    if (allocated (beam_structure%phi))  deallocate (beam_structure%phi)
  end subroutine beam_structure_final_mom

  module subroutine beam_structure_set_momentum (beam_structure, p)
    class(beam_structure_t), intent(inout) :: beam_structure
    real(default), dimension(:), intent(in) :: p
    if (allocated (beam_structure%p))  deallocate (beam_structure%p)
    allocate (beam_structure%p (size (p)), source = p)
  end subroutine beam_structure_set_momentum

  module subroutine beam_structure_set_theta (beam_structure, theta)
    class(beam_structure_t), intent(inout) :: beam_structure
    real(default), dimension(:), intent(in) :: theta
    if (allocated (beam_structure%theta))  deallocate (beam_structure%theta)
    allocate (beam_structure%theta (size (theta)), source = theta)
  end subroutine beam_structure_set_theta

  module subroutine beam_structure_set_phi (beam_structure, phi)
    class(beam_structure_t), intent(inout) :: beam_structure
    real(default), dimension(:), intent(in) :: phi
    if (allocated (beam_structure%phi))  deallocate (beam_structure%phi)
    allocate (beam_structure%phi (size (phi)), source = phi)
  end subroutine beam_structure_set_phi

  module function beam_structure_is_set (beam_structure) result (flag)
    class(beam_structure_t), intent(in) :: beam_structure
    logical :: flag
    flag = beam_structure%n_beam > 0 .or. beam_structure%asymmetric ()
  end function beam_structure_is_set

  module function beam_structure_get_n_beam (beam_structure) result (n)
    class(beam_structure_t), intent(in) :: beam_structure
    integer :: n
    n = beam_structure%n_beam
  end function beam_structure_get_n_beam

  module function beam_structure_get_prt (beam_structure) result (prt)
    class(beam_structure_t), intent(in) :: beam_structure
    type(string_t), dimension(:), allocatable :: prt
    allocate (prt (size (beam_structure%prt)))
    prt = beam_structure%prt
  end function beam_structure_get_prt

  module function beam_structure_get_n_record (beam_structure) result (n)
    class(beam_structure_t), intent(in) :: beam_structure
    integer :: n
    if (allocated (beam_structure%record)) then
       n = size (beam_structure%record)
    else
       n = 0
    end if
  end function beam_structure_get_n_record

  module function beam_structure_get_i_entry &
       (beam_structure, i) result (i_entry)
    class(beam_structure_t), intent(in) :: beam_structure
    integer, intent(in) :: i
    integer, dimension(:), allocatable :: i_entry
    associate (record => beam_structure%record(i))
      select case (size (record%entry))
      case (1)
         if (record%entry(1)%is_valid) then
            allocate (i_entry (2), source = [1, 2])
         else
            allocate (i_entry (0))
         end if
      case (2)
         if (all (record%entry%is_valid)) then
            allocate (i_entry (2), source = [1, 2])
         else if (record%entry(1)%is_valid) then
            allocate (i_entry (1), source = [1])
         else if (record%entry(2)%is_valid) then
            allocate (i_entry (1), source = [2])
         else
            allocate (i_entry (0))
         end if
      end select
    end associate
  end function beam_structure_get_i_entry

  module function beam_structure_get_name (beam_structure, i) result (name)
    type(string_t) :: name
    class(beam_structure_t), intent(in) :: beam_structure
    integer, intent(in) :: i
    associate (record => beam_structure%record(i))
      if (record%entry(1)%is_valid) then
         name = record%entry(1)%name
      else if (size (record%entry) == 2) then
         name = record%entry(2)%name
      end if
    end associate
  end function beam_structure_get_name

  module function beam_structure_has_pdf (beam_structure) result (has_pdf)
    logical :: has_pdf
    class(beam_structure_t), intent(in) :: beam_structure
    integer :: i
    type(string_t) :: name
    has_pdf = .false.
    do i = 1, beam_structure%get_n_record ()
       name = beam_structure%get_name (i)
       has_pdf = has_pdf .or. name == var_str ("pdf_builtin") .or. name == var_str ("lhapdf")
    end do
  end function beam_structure_has_pdf

  module function beam_structure_contains (beam_structure, name) result (flag)
    class(beam_structure_t), intent(in) :: beam_structure
    character(*), intent(in) :: name
    logical :: flag
    integer :: i, j
    flag = .false.
    if (allocated (beam_structure%record)) then
       do i = 1, size (beam_structure%record)
          do j = 1, size (beam_structure%record(i)%entry)
             flag = beam_structure%record(i)%entry(j)%name == name
             if (flag)  return
          end do
       end do
    end if
  end function beam_structure_contains

  module function beam_structure_polarized (beam_structure) result (flag)
    class(beam_structure_t), intent(in) :: beam_structure
    logical :: flag
    flag = allocated (beam_structure%smatrix)
  end function beam_structure_polarized

  module function beam_structure_get_smatrix (beam_structure) result (smatrix)
    class(beam_structure_t), intent(in) :: beam_structure
    type(smatrix_t), dimension(:), allocatable :: smatrix
    allocate (smatrix (size (beam_structure%smatrix)), &
         source = beam_structure%smatrix)
  end function beam_structure_get_smatrix

  module function beam_structure_get_pol_f (beam_structure) result (pol_f)
    class(beam_structure_t), intent(in) :: beam_structure
    real(default), dimension(:), allocatable :: pol_f
    allocate (pol_f (size (beam_structure%pol_f)), &
         source = beam_structure%pol_f)
  end function beam_structure_get_pol_f

  module function beam_structure_asymmetric (beam_structure) result (flag)
    class(beam_structure_t), intent(in) :: beam_structure
    logical :: flag
    flag = allocated (beam_structure%p) &
         .or. allocated (beam_structure%theta) &
         .or. allocated (beam_structure%phi)
  end function beam_structure_asymmetric

  module function beam_structure_get_momenta (beam_structure) result (p)
    class(beam_structure_t), intent(in) :: beam_structure
    type(vector3_t), dimension(:), allocatable :: p
    real(default), dimension(:), allocatable :: theta, phi
    integer :: n, i
    if (allocated (beam_structure%p)) then
       n = size (beam_structure%p)
       if (allocated (beam_structure%theta)) then
          if (size (beam_structure%theta) == n) then
             allocate (theta (n), source = beam_structure%theta)
          else
             call msg_fatal ("Beam structure: mismatch in momentum vs. &
                  &angle theta specification")
          end if
       else
          allocate (theta (n), source = 0._default)
       end if
       if (allocated (beam_structure%phi)) then
          if (size (beam_structure%phi) == n) then
             allocate (phi (n), source = beam_structure%phi)
          else
             call msg_fatal ("Beam structure: mismatch in momentum vs. &
                  &angle phi specification")
          end if
       else
          allocate (phi (n), source = 0._default)
       end if
       allocate (p (n))
       do i = 1, n
          p(i) = beam_structure%p(i) * vector3_moving ([ &
               sin (theta(i)) * cos (phi(i)), &
               sin (theta(i)) * sin (phi(i)), &
               cos (theta(i))])
       end do
       if (n == 2)  p(2) = - p(2)
    else
       call msg_fatal ("Beam structure: angle theta/phi specified but &
            &momentum/a p undefined")
    end if
  end function beam_structure_get_momenta

  module subroutine beam_structure_check_against_n_in &
       (beam_structure, n_in, applies)
    class(beam_structure_t), intent(in) :: beam_structure
    integer, intent(in) :: n_in
    logical, intent(out) :: applies
    if (beam_structure%is_set ()) then
       if (n_in == beam_structure%get_n_beam ()) then
          applies = .true.
       else if (beam_structure%get_n_beam () == 0) then
          call msg_fatal &
               ("Asymmetric beams: missing beam particle specification")
          applies = .false.
       else
          call msg_fatal &
               ("Mismatch of process and beam setup (scattering/decay)")
          applies = .false.
       end if
    else
       applies = .false.
    end if
  end subroutine beam_structure_check_against_n_in


end submodule beam_structures_s

