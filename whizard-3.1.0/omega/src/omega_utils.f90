!  omegalib.nw --
!
!  Copyright (C) 1999-2022 by
!      Wolfgang Kilian <kilian@physik.uni-siegen.de>
!      Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!      Juergen Reuter <juergen.reuter@desy.de>
!      with contributions from                                                                                                                                    
!      Fabian Bach <fabian.bach@t-online.de>                                                                                                                 
!      Bijan Chokoufe Nejad <bijan.chokoufe@desy.de>                                                                                                              
!      Christian Speckner <cnspeckn@googlemail.com>     
!
!  WHIZARD is free software; you can redistribute it and/or modify it
!  under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2, or (at your option)
!  any later version.
!
!  WHIZARD is distributed in the hope that it will be useful, but
!  WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module omega_utils
  use kinds
  use omega_vectors
  use omega_polarizations
  implicit none
  private
  public :: omega_update_helicity_selection
  public :: omega_report_helicity_selection
  public :: omega_ward_warn, omega_ward_panic
  public :: omega_slavnov_warn, omega_slavnov_panic
  public :: omega_check_arguments_warn, omega_check_arguments_panic
  public :: omega_check_helicities_warn, omega_check_helicities_panic
  private :: omega_check_helicity
  public :: omega_check_momenta_warn, omega_check_momenta_panic
  private :: check_momentum_conservation, check_mass_shell
  integer, parameter, private :: MOMENTUM_TOLERANCE = 10000
  integer, parameter, private :: ON_SHELL_TOLERANCE = 1000000
  integer, parameter, public :: omega_utils_2010_01_A = 0
contains
  pure subroutine omega_update_helicity_selection &
               (count, amp, max_abs, sum_abs, mask, threshold, cutoff, mask_dirty)
    integer, intent(inout) :: count
    complex(kind=default), dimension(:,:,:), intent(in) :: amp
    real(kind=default), dimension(:), intent(inout) :: max_abs
    real(kind=default), intent(inout) :: sum_abs
    logical, dimension(:), intent(inout) :: mask
    real(kind=default), intent(in) :: threshold
    integer, intent(in) :: cutoff
    logical, intent(out) :: mask_dirty
    integer :: h
    real(kind=default) :: avg
    mask_dirty = .false.
    if (threshold > 0) then
       count = count + 1
       if (count <= cutoff) then
          forall (h = lbound (amp, 3) : ubound (amp, 3))
             max_abs(h) = max (max_abs(h), maxval (abs (amp(:,:,h))))
          end forall
          sum_abs = sum_abs + sum (abs (amp))
          if (count == cutoff) then
             avg = sum_abs / size (amp) / cutoff
             mask = max_abs >= threshold * epsilon (avg) * avg
             mask_dirty = .true.
          end if
       end if
    end if
  end subroutine omega_update_helicity_selection
  subroutine omega_report_helicity_selection (mask, spin_states, threshold, unit)
    logical, dimension(:), intent(in) :: mask
    integer, dimension(:,:), intent(in) :: spin_states
    real(kind=default), intent(in) :: threshold
    integer, intent(in), optional :: unit
    integer :: u
    integer :: h, i
    if (present(unit)) then
       u = unit
    else
       u = 6
    end if
    if (u >= 0) then
       write (unit = u, &
            fmt = "('| ','Contributing Helicity Combinations: ', I5, ' of ', I5)") &
            count (mask), size (mask)
       write (unit = u, &
            fmt = "('| ','Threshold: amp / avg > ', E9.2, ' = ', E9.2, ' * epsilon()')") &
            threshold * epsilon (threshold), threshold
       i = 0
       do h = 1, size (mask)
          if (mask(h)) then
             i = i + 1
             write (unit = u, fmt = "('| ',I4,': ',20I4)") i, spin_states (:, h)
          end if
       end do
    end if
  end subroutine omega_report_helicity_selection
  subroutine omega_ward_warn (name, m, k, e)
    character(len=*), intent(in) :: name
    real(kind=default), intent(in) :: m
    type(momentum), intent(in) :: k
    type(vector), intent(in) :: e
    type(vector) :: ek
    real(kind=default) :: abs_eke, abs_ek_abs_e
    ek = eps (m, k, 4)
    abs_eke = abs (ek * e)
    abs_ek_abs_e = abs (ek) * abs (e)
    print *, name, ":", abs_eke / abs_ek_abs_e, abs (ek), abs (e)
    if (abs_eke > 1000 * epsilon (abs_ek_abs_e)) then
       print *, "O'Mega: warning: non-transverse vector field: ", &
            name, ":", abs_eke / abs_ek_abs_e, abs (e)
    end if
  end subroutine omega_ward_warn
  subroutine omega_ward_panic (name, m, k, e)
    character(len=*), intent(in) :: name
    real(kind=default), intent(in) :: m
    type(momentum), intent(in) :: k
    type(vector), intent(in) :: e
    type(vector) :: ek
    real(kind=default) :: abs_eke, abs_ek_abs_e
    ek = eps (m, k, 4)
    abs_eke = abs (ek * e)
    abs_ek_abs_e = abs (ek) * abs (e)
    if (abs_eke > 1000 * epsilon (abs_ek_abs_e)) then
       print *, "O'Mega: panic: non-transverse vector field: ", &
            name, ":", abs_eke / abs_ek_abs_e, abs (e)
       stop
    end if
  end subroutine omega_ward_panic
  subroutine omega_slavnov_warn (name, m, k, e, phi)
    character(len=*), intent(in) :: name
    real(kind=default), intent(in) :: m
    type(momentum), intent(in) :: k
    type(vector), intent(in) :: e
    complex(kind=default), intent(in) :: phi
    type(vector) :: ek
    real(kind=default) :: abs_eke, abs_ek_abs_e
    ek = eps (m, k, 4)
    abs_eke = abs (ek * e - phi)
    abs_ek_abs_e = abs (ek) * abs (e)
    print *, name, ":", abs_eke / abs_ek_abs_e, abs (ek), abs (e)
    if (abs_eke > 1000 * epsilon (abs_ek_abs_e)) then
       print *, "O'Mega: warning: non-transverse vector field: ", &
            name, ":", abs_eke / abs_ek_abs_e, abs (e)
    end if
  end subroutine omega_slavnov_warn
  subroutine omega_slavnov_panic (name, m, k, e, phi)
    character(len=*), intent(in) :: name
    real(kind=default), intent(in) :: m
    type(momentum), intent(in) :: k
    type(vector), intent(in) :: e
    complex(kind=default), intent(in) :: phi
    type(vector) :: ek
    real(kind=default) :: abs_eke, abs_ek_abs_e
    ek = eps (m, k, 4)
    abs_eke = abs (ek * e - phi)
    abs_ek_abs_e = abs (ek) * abs (e)
    if (abs_eke > 1000 * epsilon (abs_ek_abs_e)) then
       print *, "O'Mega: panic: non-transverse vector field: ", &
            name, ":", abs_eke / abs_ek_abs_e, abs (e)
       stop
    end if
  end subroutine omega_slavnov_panic
  subroutine omega_check_arguments_warn (n, k)
    integer, intent(in) :: n
    real(kind=default), dimension(0:,:), intent(in) :: k
    integer :: i
    i = size(k,dim=1)
    if (i /= 4) then
       print *, "O'Mega: warning: wrong # of dimensions:", i
    end if
    i = size(k,dim=2)
    if (i /= n) then
       print *, "O'Mega: warning: wrong # of momenta:", i, &
            ", expected", n
    end if
  end subroutine omega_check_arguments_warn
  subroutine omega_check_arguments_panic (n, k)
    integer, intent(in) :: n
    real(kind=default), dimension(0:,:), intent(in) :: k
    logical :: error
    integer :: i
    error = .false.
    i = size(k,dim=1)
    if (i /= n) then
       print *, "O'Mega: warning: wrong # of dimensions:", i
       error = .true.
    end if
    i = size(k,dim=2)
    if (i /= n) then
       print *, "O'Mega: warning: wrong # of momenta:", i, &
            ", expected", n
       error = .true.
    end if
    if (error) then
       stop
    end if
  end subroutine omega_check_arguments_panic
  function omega_check_helicity (m, smax, s) result (error)
    real(kind=default), intent(in) :: m
    integer, intent(in) :: smax, s
    logical :: error
    select case (smax)
    case (0)
       error = (s /= 0)
    case (1)
       error = (abs (s) /= 1)
    case (2)
       if (m == 0.0_default) then
          error = .not. (abs (s) == 1 .or. abs (s) == 4)
       else
          error = .not. (abs (s) <= 1 .or. abs (s) == 4)
       end if
    case (4)
       error = .true.
    case default
       error = .true.
    end select
  end function omega_check_helicity
  subroutine omega_check_helicities_warn (m, smax, s)
    real(kind=default), dimension(:), intent(in) :: m
    integer, dimension(:), intent(in) :: smax, s
    integer :: i
    do i = 1, size (m)
       if (omega_check_helicity (m(i), smax(i), s(i))) then
          print *, "O'Mega: warning: invalid helicity", s(i)
       end if
    end do
  end subroutine omega_check_helicities_warn
  subroutine omega_check_helicities_panic (m, smax, s)
    real(kind=default), dimension(:), intent(in) :: m
    integer, dimension(:), intent(in) :: smax, s
    logical :: error
    logical :: error1
    integer :: i
    error = .false.
    do i = 1, size (m)
       error1 = omega_check_helicity (m(i), smax(i), s(i))
       if (error1) then
          print *, "O'Mega: panic: invalid helicity", s(i)
          error = .true.
       end if
    end do
    if (error) then
       stop
    end if
  end subroutine omega_check_helicities_panic
  function check_momentum_conservation (k) result (error)
    real(kind=default), dimension(0:,:), intent(in) :: k
    logical :: error
    error = any (abs (sum (k(:,3:), dim = 2) - k(:,1) - k(:,2)) > &
         MOMENTUM_TOLERANCE * epsilon (maxval (abs (k), dim = 2)))
    if (error) then
       print *, sum (k(:,3:), dim = 2) - k(:,1) - k(:,2)
       print *, MOMENTUM_TOLERANCE * epsilon (maxval (abs (k), dim = 2)), &
            maxval (abs (k), dim = 2)
    end if
  end function check_momentum_conservation
  function check_mass_shell (m, k) result (error)
    real(kind=default), intent(in) :: m
    real(kind=default), dimension(0:), intent(in) :: k
    real(kind=default) :: e2
    logical :: error
    e2 = k(1)**2 + k(2)**2 + k(3)**2 + m**2
    error = abs (k(0)**2 - e2) > ON_SHELL_TOLERANCE * epsilon (max (k(0)**2, e2))
    if (error) then
       print *, k(0)**2 - e2
       print *, ON_SHELL_TOLERANCE * epsilon (max (k(0)**2, e2)), max (k(0)**2, e2)
    end if
  end function check_mass_shell
  subroutine omega_check_momenta_warn (m, k)
    real(kind=default), dimension(:), intent(in) :: m
    real(kind=default), dimension(0:,:), intent(in) :: k
    integer :: i
    if (check_momentum_conservation (k)) then
       print *, "O'Mega: warning: momentum not conserved"
    end if
    do i = 1, size(m)
       if (check_mass_shell (m(i), k(:,i))) then
          print *, "O'Mega: warning: particle #", i, "not on-shell"
       end if
    end do
  end subroutine omega_check_momenta_warn
  subroutine omega_check_momenta_panic (m, k)
    real(kind=default), dimension(:), intent(in) :: m
    real(kind=default), dimension(0:,:), intent(in) :: k
    logical :: error
    logical :: error1
    integer :: i
    error = check_momentum_conservation (k)
    if (error) then
       print *, "O'Mega: panic: momentum not conserved"
    end if
    do i = 1, size(m)
       error1 = check_mass_shell (m(i), k(0:,i))
       if (error1) then
          print *, "O'Mega: panic: particle #", i, "not on-shell"
          error = .true.
       end if
    end do
    if (error) then
       stop
    end if
  end subroutine omega_check_momenta_panic
end module omega_utils
