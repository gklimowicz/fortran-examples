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

submodule (polarizations) polarizations_s

  use io_units
  use format_defs, only: FMT_19
  use diagnostics
  use helicities

  implicit none

contains

  module subroutine polarization_init (pol, spin_type, multiplicity, &
       anti, left_handed, right_handed)
    class(polarization_t), intent(out) :: pol
    integer, intent(in) :: spin_type
    integer, intent(in) :: multiplicity
    logical, intent(in) :: anti
    logical, intent(in) :: left_handed
    logical, intent(in) :: right_handed
    pol%spin_type = spin_type
    pol%multiplicity = multiplicity
    pol%anti = anti
    select case (pol%multiplicity)
    case (1)
       if (left_handed) then
          pol%chirality = -1
       else if (right_handed) then
          pol%chirality = 1
       end if
    end select
    select case (pol%chirality)
    case (0)
       call pol%bv%init_unpolarized (spin_type)
    end select
  end subroutine polarization_init

  module subroutine polarization_init_flv (pol, flv)
    class(polarization_t), intent(out) :: pol
    type(flavor_t), intent(in) :: flv
    call pol%init ( &
         spin_type = flv%get_spin_type (), &
         multiplicity = flv%get_multiplicity (), &
         anti = flv%is_antiparticle (), &
         left_handed = flv%is_left_handed (), &
         right_handed = flv%is_right_handed ())
  end subroutine polarization_init_flv

  module subroutine polarization_init_generic (pol, spin_type, multiplicity, &
       anti, left_handed, right_handed)
    class(polarization_t), intent(out) :: pol
    integer, intent(in) :: spin_type
    integer, intent(in) :: multiplicity
    logical, intent(in) :: anti
    logical, intent(in) :: left_handed
    logical, intent(in) :: right_handed
    call pol%init (spin_type, multiplicity, &
         anti, left_handed, right_handed)
    select case (pol%chirality)
    case (0)
       if (pol%multiplicity == pol%bv%get_n_states ()) then
          call pol%bv%init (spin_type)
       else
          call pol%bv%init_max_weight (spin_type)
       end if
    end select
  end subroutine polarization_init_generic

  module subroutine polarization_init_generic_flv (pol, flv)
    class(polarization_t), intent(out) :: pol
    type(flavor_t), intent(in) :: flv
    call pol%init_generic ( &
         spin_type = flv%get_spin_type (), &
         multiplicity = flv%get_multiplicity (), &
         anti = flv%is_antiparticle (), &
         left_handed = flv%is_left_handed (), &
         right_handed = flv%is_right_handed ())
  end subroutine polarization_init_generic_flv

  module subroutine polarization_write (pol, unit, state_matrix, all_states, tolerance)
    class(polarization_t), intent(in) :: pol
    integer, intent(in), optional :: unit
    logical, intent(in), optional :: state_matrix, all_states
    real(default), intent(in), optional :: tolerance
    logical :: state_m
    type(state_matrix_t) :: state
    real(default), dimension(:), allocatable :: a
    integer :: u, i
    u = given_output_unit (unit);  if (u < 0)  return
    state_m = .false.;  if (present (state_matrix))  state_m = state_matrix
    if (pol%anti) then
       write (u, "(1x,A,I1,A,I1,A,L1,A)")  &
            "Polarization: [spin_type = ", pol%spin_type, &
            ", mult = ", pol%multiplicity, ", anti = ", pol%anti, "]"
    else
       write (u, "(1x,A,I1,A,I1,A)")  &
            "Polarization: [spin_type = ", pol%spin_type, &
            ", mult = ", pol%multiplicity, "]"
    end if
    if (state_m) then
       call pol%to_state (state, all_states, tolerance)
       call state%write (unit=unit)
       call state%final ()
    else if (pol%chirality == 1) then
       write (u, "(1x,A)")  "chirality = +"
    else if (pol%chirality == -1) then
       write (u, "(1x,A)")  "chirality = -"
    else if (pol%bv%is_polarized ()) then
       call pol%bv%to_array (a)
       do i = 1, size (a)
          write (u, "(1x,I2,':',1x,F10.7)")  i, a(i)
       end do
    else
       write (u, "(1x,A)")  "[unpolarized]"
    end if
  end subroutine polarization_write

  module subroutine polarization_write_raw (pol, u)
    class(polarization_t), intent(in) :: pol
    integer, intent(in) :: u
    write (u) pol%spin_type
    write (u) pol%multiplicity
    write (u) pol%chirality
    write (u) pol%anti
    call pol%bv%write_raw (u)
  end subroutine polarization_write_raw

  module subroutine polarization_read_raw (pol, u, iostat)
    class(polarization_t), intent(out) :: pol
    integer, intent(in) :: u
    integer, intent(out), optional :: iostat
    read (u, iostat=iostat) pol%spin_type
    read (u, iostat=iostat) pol%multiplicity
    read (u, iostat=iostat) pol%chirality
    read (u, iostat=iostat) pol%anti
    call pol%bv%read_raw (u, iostat)
  end subroutine polarization_read_raw

  module function polarization_is_polarized (pol) result (polarized)
    class(polarization_t), intent(in) :: pol
    logical :: polarized
    polarized = pol%chirality /= 0 .or. pol%bv%is_polarized ()
  end function polarization_is_polarized

  module function polarization_is_diagonal (pol) result (diagonal)
    class(polarization_t), intent(in) :: pol
    logical :: diagonal
    select case (pol%chirality)
    case (0)
       diagonal = pol%bv%is_diagonal ()
    case default
       diagonal = .true.
    end select
  end function polarization_is_diagonal

  module subroutine polarization_init_state_matrix (pol, state)
    class(polarization_t), intent(out) :: pol
    type(state_matrix_t), intent(in), target :: state
    type(state_iterator_t) :: it
    type(flavor_t) :: flv
    type(helicity_t) :: hel
    integer :: d, h1, h2, i, j
    complex(default), dimension(:,:), allocatable :: r
    complex(default) :: me
    real(default) :: trace
    call it%init (state)
    flv = it%get_flavor (1)
    hel = it%get_helicity (1)
    if (hel%is_defined ()) then
       call pol%init_generic (flv)
       select case (pol%chirality)
       case (0)
          trace = 0
          d = pol%bv%get_n_states ()
          allocate (r (d, d), source = (0._default, 0._default))
          do while (it%is_valid ())
             hel = it%get_helicity (1)
             call hel%get_indices (h1, h2)
             i = pol%bv%hel_index (h1)
             j = pol%bv%hel_index (h2)
             me = it%get_matrix_element ()
             r(i,j) = me
             if (i == j)  trace = trace + real (me)
             call it%advance ()
          end do
          if (trace /= 0)  call pol%bv%set (r / trace)
       end select
    else
       call pol%init (flv)
    end if
  end subroutine polarization_init_state_matrix

  module subroutine polarization_to_state_matrix (pol, state, all_states, tolerance)
    class(polarization_t), intent(in), target :: pol
    type(state_matrix_t), intent(out) :: state
    logical, intent(in), optional :: all_states
    real(default), intent(in), optional :: tolerance
    type(polarization_iterator_t) :: it
    type(quantum_numbers_t), dimension(1) :: qn
    complex(default) :: value
    call it%init (pol, all_states, tolerance)
    call state%init (store_values = .true.)
    do while (it%is_valid ())
       value = it%get_value ()
       qn(1) = it%get_quantum_numbers ()
       call state%add_state (qn, value = value)
       call it%advance ()
    end do
    call state%freeze ()
  end subroutine polarization_to_state_matrix

  module subroutine polarization_init_unpolarized (pol, flv)
    class(polarization_t), intent(out) :: pol
    type(flavor_t), intent(in) :: flv
    call pol%init (flv)
  end subroutine polarization_init_unpolarized

  module subroutine polarization_init_circular (pol, flv, f)
    class(polarization_t), intent(out) :: pol
    type(flavor_t), intent(in) :: flv
    real(default), intent(in) :: f
    call pol%init (flv)
    select case (pol%chirality)
    case (0)
       call pol%bv%init_vector (pol%spin_type, &
            [0._default, 0._default, f])
    end select
  end subroutine polarization_init_circular

  module subroutine polarization_init_transversal (pol, flv, phi, f)
    class(polarization_t), intent(out) :: pol
    type(flavor_t), intent(in) :: flv
    real(default), intent(in) :: phi, f
    call pol%init (flv)
    select case (pol%chirality)
    case (0)
       if (pol%anti) then
          call pol%bv%init_vector (pol%spin_type, &
               [f * cos (phi), f * sin (phi), 0._default])
       else
          call pol%bv%init_vector (pol%spin_type, &
               [f * cos (phi),-f * sin (phi), 0._default])
       end if
    end select
  end subroutine polarization_init_transversal

  module subroutine polarization_init_axis (pol, flv, alpha)
    class(polarization_t), intent(out) :: pol
    type(flavor_t), intent(in) :: flv
    real(default), dimension(3), intent(in) :: alpha
    call pol%init (flv)
    select case (pol%chirality)
    case (0)
       if (pol%anti) then
          call pol%bv%init_vector (pol%spin_type, &
               [alpha(1), alpha(2), alpha(3)])
       else
          call pol%bv%init_vector (pol%spin_type, &
               [alpha(1),-alpha(2), alpha(3)])
       end if
    end select
  end subroutine polarization_init_axis

  module subroutine polarization_init_angles (pol, flv, r, theta, phi)
    class(polarization_t), intent(out) :: pol
    type(flavor_t), intent(in) :: flv
    real(default), intent(in) :: r, theta, phi
    real(default), dimension(3) :: alpha
    real(default), parameter :: eps = 10 * epsilon (1._default)

    alpha(1) = r * sin (theta) * cos (phi)
    alpha(2) = r * sin (theta) * sin (phi)
    alpha(3) = r * cos (theta)
    where (abs (alpha) < eps)  alpha = 0
    call pol%init_axis (flv, alpha)
  end subroutine polarization_init_angles

  module subroutine polarization_init_longitudinal (pol, flv, f)
    class(polarization_t), intent(out) :: pol
    type(flavor_t), intent(in) :: flv
    real(default), intent(in) :: f
    real(default), dimension(:), allocatable :: rd
    integer :: s, d
    s = flv%get_spin_type ()
    select case (s)
    case (VECTOR, TENSOR)
       call pol%init_generic (flv)
       if (pol%bv%is_polarized ()) then
          d = pol%bv%get_n_states ()
          allocate (rd (d), source = 0._default)
          rd(pol%bv%hel_index (0)) = f
          call pol%bv%set (rd)
       end if
    case default
       call pol%init_unpolarized (flv)
    end select
  end subroutine polarization_init_longitudinal

  module subroutine polarization_init_diagonal (pol, flv, rd)
    class(polarization_t), intent(out) :: pol
    type(flavor_t), intent(in) :: flv
    real(default), dimension(:), intent(in) :: rd
    real(default) :: trace
    call pol%init_generic (flv)
    if (pol%bv%is_polarized ()) then
       trace = sum (rd)
       if (trace /= 0)  call pol%bv%set (rd / trace)
    end if
  end subroutine polarization_init_diagonal

  module subroutine combine_polarization_states (pol, state)
    type(polarization_t), dimension(:), intent(in), target :: pol
    type(state_matrix_t), intent(out) :: state
    type(state_matrix_t), dimension(size(pol)), target :: pol_state
    integer :: i
    do i = 1, size (pol)
       call pol(i)%to_state (pol_state(i))
    end do
    call outer_multiply (pol_state, state)
    do i = 1, size (pol)
       call pol_state(i)%final ()
    end do
  end subroutine combine_polarization_states

  module function polarization_get_axis (pol) result (alpha)
    class(polarization_t), intent(in), target :: pol
    real(default), dimension(3) :: alpha
    select case (pol%chirality)
    case (0)
       call pol%bv%to_vector (alpha)
       if (.not. pol%anti)  alpha(2) = - alpha(2)
    case (-1)
       alpha = [0._default, 0._default, -1._default]
    case (1)
       alpha = [0._default, 0._default, 1._default]
    end select
  end function polarization_get_axis

  module subroutine polarization_to_angles (pol, r, theta, phi)
    class(polarization_t), intent(in) :: pol
    real(default), intent(out) :: r, theta, phi
    real(default), dimension(3) :: alpha
    real(default) :: norm, r12
    alpha = pol%get_axis ()
    norm = sum (alpha**2)
    r = sqrt (norm)
    if (norm > 0) then
       r12 = sqrt (alpha(1)**2 + alpha(2)**2)
       theta = atan2 (r12, alpha(3))
       if (any (alpha(1:2) /= 0)) then
          phi = atan2 (alpha(2), alpha(1))
       else
          phi = 0
       end if
    else
       theta = 0
       phi = 0
    end if
  end subroutine polarization_to_angles

  module subroutine polarization_iterator_write (it, unit)
    class(polarization_iterator_t), intent(in) :: it
    integer, intent(in), optional :: unit
    integer :: u, i
    u = given_output_unit (unit)
    write (u, "(1X,A)")  "Polarization iterator:"
    write (u, "(3X,A,L1)")  "assigned = ", associated (it%pol)
    write (u, "(3X,A,L1)")  "valid    = ", it%valid
    if (it%valid) then
       write (u, "(3X,A,2(1X,I2))")  "i, j     = ", it%i, it%j
       write (u, "(3X,A,2(1X,I2))")  "h1, h2   = ", it%h1, it%h2
       write (u, "(3X,A)", advance="no")  "value    = "
       write (u, *)  it%value
       if (allocated (it%r)) then
          do i = 1, size (it%r, 2)
             write (u, *)  it%r(i,:)
          end do
       end if
    end if
  end subroutine polarization_iterator_write

  module subroutine polarization_iterator_init (it, pol, all_states, tolerance)
    class(polarization_iterator_t), intent(out) :: it
    type(polarization_t), intent(in), target :: pol
    logical, intent(in), optional :: all_states
    real(default), intent(in), optional :: tolerance
    integer :: d
    logical :: only_max_weight
    it%pol => pol
    if (present (all_states)) then
       if (.not. all_states) then
          if (present (tolerance)) then
             it%tolerance = tolerance
          else
             it%tolerance = 0
          end if
       end if
    end if
    select case (pol%chirality)
    case (0)
       d = pol%bv%get_n_states ()
       only_max_weight = pol%multiplicity < d
       it%polarized = pol%bv%is_polarized ()
       if (it%polarized) then
          it%i = d
          it%j = it%i
          it%h1 = pol%bv%hel_value (it%i)
          it%h2 = it%h1
          call pol%bv%to_matrix (it%r, only_max_weight)
          it%value = it%r(it%i, it%j)
       else
          it%value = 1._default / d
       end if
       it%valid = .true.
    case (1,-1)
       it%polarized = .true.
       select case (pol%spin_type)
       case (SPINOR)
          it%h1 = pol%chirality
       case (VECTORSPINOR)
          it%h1 = 2 * pol%chirality
       end select
       it%h2 = it%h1
       it%valid = .true.
    end select
    if (it%valid .and. abs (it%value) <= it%tolerance)  call it%advance ()
  end subroutine polarization_iterator_init

  recursive module subroutine polarization_iterator_advance (it)
    class(polarization_iterator_t), intent(inout) :: it
    if (it%valid) then
       select case (it%pol%chirality)
       case (0)
          if (it%polarized) then
             if (it%j > 1) then
                it%j = it%j - 1
                it%h2 = it%pol%bv%hel_value (it%j)
                it%value = it%r(it%i, it%j)
             else if (it%i > 1) then
                it%j = it%pol%bv%get_n_states ()
                it%h2 = it%pol%bv%hel_value (it%j)
                it%i = it%i - 1
                it%h1 = it%pol%bv%hel_value (it%i)
                it%value = it%r(it%i, it%j)
             else
                it%valid = .false.
             end if
          else
             it%valid = .false.
          end if
       case default
          it%valid = .false.
       end select
       if (it%valid .and. abs (it%value) <= it%tolerance)  call it%advance ()
    end if
  end subroutine polarization_iterator_advance

  module function polarization_iterator_is_valid (it) result (is_valid)
    logical :: is_valid
    class(polarization_iterator_t), intent(in) :: it
    is_valid = it%valid
  end function polarization_iterator_is_valid

  module function polarization_iterator_get_value (it) result (value)
    complex(default) :: value
    class(polarization_iterator_t), intent(in) :: it
    if (it%valid) then
       value = it%value
    else
       value = 0
    end if
  end function polarization_iterator_get_value

  module function polarization_iterator_get_quantum_numbers (it) result (qn)
    class(polarization_iterator_t), intent(in) :: it
    type(helicity_t) :: hel
    type(quantum_numbers_t) :: qn
    if (it%polarized) then
       call hel%init (it%h2, it%h1)
    end if
    call qn%init (hel)
  end function polarization_iterator_get_quantum_numbers

  module subroutine smatrix_write (object, unit, indent)
    class(smatrix_t), intent(in) :: object
    integer, intent(in), optional :: unit, indent
    integer :: u, i, ind
    u = given_output_unit (unit)
    ind = 0;  if (present (indent))  ind = indent
    if (allocated (object%value)) then
       if (size (object%value) > 0) then
          do i = 1, object%n_entry
             write (u, "(1x,A,'@(')", advance="no")  repeat ("  ", ind)
             write (u, "(SP,9999(I2.1,':',1x))", advance="no") &
                  object%index(:,i)
             write (u, "('('," // FMT_19 // ",','," // FMT_19 // &
                  ",'))')")  object%value(i)
          end do
       else
          write (u, "(1x,A)", advance="no")  repeat ("  ", ind)
          write (u, "(A)")  "[empty matrix]"
       end if
    else
       write (u, "(1x,A)", advance="no")  repeat ("  ", ind)
       write (u, "(A)")  "[undefined matrix]"
    end if
  end subroutine smatrix_write

  module subroutine smatrix_init (smatrix, dim, n_entry)
    class(smatrix_t), intent(out) :: smatrix
    integer, intent(in) :: dim
    integer, intent(in) :: n_entry
    smatrix%dim = dim
    smatrix%n_entry = n_entry
    allocate (smatrix%index (dim, n_entry))
    allocate (smatrix%value (n_entry))
  end subroutine smatrix_init

  module subroutine smatrix_set_entry (smatrix, i, index, value)
    class(smatrix_t), intent(inout) :: smatrix
    integer, intent(in) :: i
    integer, dimension(:), intent(in) :: index
    complex(default), intent(in) :: value
    smatrix%index(:,i) = index
    smatrix%value(i) = value
  end subroutine smatrix_set_entry

  elemental module function smatrix_exists (smatrix) result (exist)
    logical :: exist
    class(smatrix_t), intent(in) :: smatrix
    exist = .not. all (smatrix%value == 0)
  end function smatrix_exists

  module subroutine pmatrix_write (object, unit, indent)
    class(pmatrix_t), intent(in) :: object
    integer, intent(in), optional :: unit, indent
    integer :: u
    u = given_output_unit (unit)
    write (u, "(1x,A)")  "Polarization: spin density matrix"
    write (u, "(3x,A,I0)")  "spin type     = ", object%spin_type
    write (u, "(3x,A,I0)")  "multiplicity  = ", object%multiplicity
    write (u, "(3x,A,L1)")  "massive       = ", object%massive
    write (u, "(3x,A,I0)")  "chirality     = ", object%chirality
    write (u, "(3x,A,F10.7)")  "pol.degree    =", object%degree
    write (u, "(3x,A,L1)")  "pure state    = ", object%pure
    call object%smatrix_t%write (u, 1)
  end subroutine pmatrix_write

  module subroutine pmatrix_assign_from_smatrix (pmatrix, smatrix)
    class(pmatrix_t), intent(out) :: pmatrix
    type(smatrix_t), intent(in) :: smatrix
    pmatrix%smatrix_t = smatrix
  end subroutine pmatrix_assign_from_smatrix

  module subroutine pmatrix_normalize (pmatrix, flv, degree, tolerance)
    class(pmatrix_t), intent(inout) :: pmatrix
    type(flavor_t), intent(in) :: flv
    real(default), intent(in), optional :: degree
    real(default), intent(in), optional :: tolerance
    integer :: i, hmax
    logical :: fermion, ok
    real(default) :: trace, trace_sq
    real(default) :: tol
    tol = 0;  if (present (tolerance))  tol = tolerance
    pmatrix%spin_type = flv%get_spin_type ()
    pmatrix%massive = flv%get_mass () /= 0
    if (.not. pmatrix%massive) then
       if (flv%is_left_handed ()) then
          pmatrix%chirality = -1
       else if (flv%is_right_handed ()) then
          pmatrix%chirality = +1
       end if
    end if
    if (pmatrix%spin_type == SCALAR) then
       pmatrix%multiplicity = 1
    else if (pmatrix%massive) then
       pmatrix%multiplicity = pmatrix%spin_type
    else if (pmatrix%chirality == 0) then
       pmatrix%multiplicity = 2
    else
       pmatrix%multiplicity = 1
    end if
    if (present (degree)) then
       if (degree < 0 .or. degree > 1) &
            call msg_error ("polarization degree must be between 0 and 1")
       pmatrix%degree = degree
    end if
    if (size (pmatrix%index, 1) /= 2)  call error ("wrong array rank")
    fermion = mod (pmatrix%spin_type, 2) == 0
    hmax = pmatrix%spin_type / 2
    if (pmatrix%n_entry > 0) then
       if (fermion) then
          if (pmatrix%massive) then
             ok = all (pmatrix%index /= 0) &
                  .and. all (abs (pmatrix%index) <= hmax)
          else if (pmatrix%chirality == -1) then
             ok = all (pmatrix%index == -hmax)
          else if (pmatrix%chirality == +1) then
             ok = all (pmatrix%index == +hmax)
          else
             ok = all (abs (pmatrix%index) == hmax)
          end if
       else
          if (pmatrix%massive) then
             ok = all (abs (pmatrix%index) <= hmax)
          else
             ok = all (abs (pmatrix%index) == hmax)
          end if
       end if
       if (.not. ok)  call error ("illegal index value")
    else
       pmatrix%degree = 0
       pmatrix%pure = pmatrix%multiplicity == 1
       return
    end if
    trace = 0
    do i = 1, pmatrix%n_entry
       associate (index => pmatrix%index(:,i), value => pmatrix%value(i))
         if (index(1) == index(2)) then
            if (abs (aimag (value)) > tol)  call error ("diagonal must be real")
            value = real (value, kind=default)
            trace = trace + value

         else if (any (pmatrix%index(1,:) == index(2) &
              .and.    pmatrix%index(2,:) == index(1))) then
            call error ("redundant off-diagonal entry")
         else if (index(2) < index (1)) then
            index = index([2,1])
            value = conjg (value)
         end if
       end associate
    end do
    if (abs (trace) <= tol)  call error ("trace must not vanish")
    trace = real (trace, kind=default)
    pmatrix%value = pmatrix%value / trace * pmatrix%degree
    trace_sq = (1 - pmatrix%degree ** 2) / pmatrix%multiplicity
    do i = 1, pmatrix%n_entry
       associate (index => pmatrix%index(:,i), value => pmatrix%value(i))
         if (index(1) == index(2)) then
            trace_sq = trace_sq + abs (value) ** 2
         else
            trace_sq = trace_sq + 2 * abs (value) ** 2
         end if
       end associate
    end do
    if (pmatrix%multiplicity == 1) then
       pmatrix%pure = .true.
    else if (abs (trace_sq - 1) <= tol) then
       pmatrix%pure = .true.
    else if (trace_sq - 1 > tol .or. trace_sq < -tol) then
       print *, "Trace of matrix square = ", trace_sq
       call error ("not permissible as density matrix")
    end if
  contains
    subroutine error (msg)
      character(*), intent(in) :: msg
      call pmatrix%write ()
      call msg_fatal ("Spin density matrix: " // msg)
    end subroutine error
  end subroutine pmatrix_normalize

  elemental module function pmatrix_is_polarized (pmatrix) result (flag)
    class(pmatrix_t), intent(in) :: pmatrix
    logical :: flag
    flag = pmatrix%degree > 0
  end function pmatrix_is_polarized

  elemental module function pmatrix_is_diagonal (pmatrix) result (flag)
    class(pmatrix_t), intent(in) :: pmatrix
    logical :: flag
    flag = all (pmatrix%index(1,:) == pmatrix%index(2,:))
  end function pmatrix_is_diagonal

  elemental module function pmatrix_get_simple_pol (pmatrix) result (pol)
    class(pmatrix_t), intent(in) :: pmatrix
    real(default) :: pol
    if (pmatrix%is_polarized ()) then
       select case (size (pmatrix%value))
       case (0)
          pol = 0
       case (1)
          pol = pmatrix%index (1,1) * pmatrix%degree
       case (2)
          pol = 42
       end select
    else
       pol = 0
    end if
  end function pmatrix_get_simple_pol

  module subroutine polarization_init_pmatrix (pol, pmatrix)
    class(polarization_t), intent(out) :: pol
    type(pmatrix_t), intent(in) :: pmatrix
    integer :: d, i, j, k, h1, h2
    complex(default), dimension(:,:), allocatable :: r
    call pol%init_generic ( &
         spin_type = pmatrix%spin_type, &
         multiplicity = pmatrix%multiplicity, &
         anti = .false., &                     !!! SUFFICIENT?
         left_handed = pmatrix%chirality < 0, &
         right_handed = pmatrix%chirality > 0)
    if (pol%bv%is_polarized ()) then
       d = pol%bv%get_n_states ()
       allocate (r (d, d), source = (0._default, 0._default))
       if (d == pmatrix%multiplicity) then
          do i = 1, d
             r(i,i) = (1 - pmatrix%degree) / d
          end do
       else if (d > pmatrix%multiplicity) then
          r(1,1) = (1 - pmatrix%degree) / 2
          r(d,d) = r(1,1)
       end if
       do k = 1, size (pmatrix%value)
          h1 = pmatrix%index(1,k)
          h2 = pmatrix%index(2,k)
          i = pol%bv%hel_index (h1)
          j = pol%bv%hel_index (h2)
          r(i,j) = r(i,j) + pmatrix%value(k)
          r(j,i) = conjg (r(i,j))
       end do
       call pol%bv%set (r)
    end if
  end subroutine polarization_init_pmatrix


end submodule polarizations_s

