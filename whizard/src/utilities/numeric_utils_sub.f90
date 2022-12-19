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

submodule (numeric_utils) numeric_utils_s

  use string_utils
  use constants
  use format_defs

  implicit none

contains

  module subroutine int_workspace_init (work, a, b, limit)
    class(int_workspace_t), intent(out) :: work
    real(default), intent(in) :: a, b
    integer, intent(in) :: limit
    work%limit = limit
    allocate (work%alist (limit), work%blist (limit), &
         work%rlist (limit), work%elist(limit), &
         work%order(limit), work%level(limit))
    work%alist(1) = a
    work%blist(1) = b
  end subroutine int_workspace_init

  module subroutine int_workspace_set_initial (work, res, err)
    class(int_workspace_t), intent(inout) :: work
    real(default), intent(in) :: res, err
    work%size = 1
    work%order(1) = 1
    work%rlist(1) = res
    work%elist(1) = err
  end subroutine int_workspace_set_initial

  module subroutine int_workspace_update (work, a1, b1, area1, error1, &
       a2, b2, area2, error2)
    class(int_workspace_t), intent(inout) :: work
    real(default), intent(in) :: a1, b1, area1, error1, &
         a2, b2, area2, error2
    integer :: i_max, i_new, new_level
    i_max = work%i
    i_new = work%size + 1
    new_level = work%level(i_max) + 1
    ! append the newly-created intervals to the list
    if (error2 > error1) then
       ! work%blist(maxerr) is already b2
       work%alist(i_max) = a2
       work%rlist(i_max) = area2
       work%elist(i_max) = error2
       work%level(i_max) = new_level
       work%alist(i_new) = a1
       work%blist(i_new) = b1
       work%rlist(i_new) = area1
       work%elist(i_new) = error1
       work%level(i_new) = new_level
    else
       ! work%alist(maxerr) is already a1
       work%blist(i_max) = b1
       work%rlist(i_max) = area1
       work%elist(i_max) = error1
       work%level(i_max) = new_level
       work%alist(i_new) = a2
       work%blist(i_new) = b2
       work%rlist(i_new) = area2
       work%elist(i_new) = error2
       work%level(i_new) = new_level
    end if
    work%size = work%size + 1
    if (new_level > work%maximum_level) &
         work%maximum_level = new_level
    call work%sort ()
  end subroutine int_workspace_update

  module subroutine int_workspace_sort (work)
    class(int_workspace_t), intent(inout) :: work
    integer :: last, limit, i, k, top, i_nrmax, i_maxerr
    real(default) :: errmax, errmin
    last = work%size
    limit = work%limit
    i_nrmax = work%nrmax
    i_maxerr = work%order(i_nrmax)
    ! Check whether the list contains more than two error estimates
    if (last < 3) then
       work%order(1) = 1
       work%order(2) = 2
       work%i = i_maxerr
       return
    end if
    errmax = work%elist(i_maxerr)
    ! This part of the routine is only executed if, due to a difficult
    ! integrand, subdivision increased the error estimate. In the normal
    ! case the insert procedure should start after the nrmax-th largest
    ! error estimate.
    DESCEND_NRMAX: do while (i_nrmax > 1) 
       if (errmax > work%elist(work%order(i_nrmax - 1))) then
          work%order(i_nrmax) = work%order(i_nrmax - 1)
          i_nrmax = i_nrmax - 1
       else
          exit DESCEND_NRMAX
       end if
    end do DESCEND_NRMAX
    ! Compute the number of elements in the list to be maintained in
    ! descending order. This number depends on the number of
    ! subdivisions still allowed.
    if (last < (limit/2 + 2)) then
       top = last
    else
       top = limit - last + 1
    end if
    ! Insert errmax by traversing the list top-down, starting
    ! comparison from the element elist(order(i_nrmax+1)).
    i = i_nrmax + 1  
    ! The order of the tests in the following line is important to
    ! prevent a segmentation fault
    ASCEND_TOP: do while (i < top)
       if (errmax < work%elist(work%order(i))) then
          work%order(i-1) = work%order(i)
          i = i + 1
       else
          exit ASCEND_TOP
       end if
    end do ASCEND_TOP
    work%order(i-1) = i_maxerr
    ! Insert errmin by traversing the list bottom-up
    errmin = work%elist(last)
    k = top - 1
    DESCEND_K: do while (k > i-2)
       if (errmin >= work%elist(work%order(k))) then
          work%order(k+1) = work%order(k)
          k = k - 1
       else
          exit DESCEND_K
       end if
    end do DESCEND_K
    work%order(k+1) = last
    ! Set i_max and e_max
    i_maxerr = work%order(i_nrmax)
    work%i = i_maxerr
    work%nrmax = i_nrmax
  end subroutine int_workspace_sort

  module subroutine assert (unit, ok, description, exit_on_fail)
    integer, intent(in) :: unit
    logical, intent(in) :: ok
    character(*), intent(in), optional :: description
    logical, intent(in), optional :: exit_on_fail
    logical :: ef
    ef = .false.;  if (present (exit_on_fail)) ef = exit_on_fail
    if (.not. ok) then
       if (present(description)) then
          write (unit, "(A)") "* FAIL: " // description
       else
          write (unit, "(A)") "* FAIL: Assertion error"
       end if
       if (ef)  stop 1
    end if
  end subroutine assert

  module subroutine assert_equal_integer (unit, lhs, rhs, description, exit_on_fail)
    integer, intent(in) :: unit
    integer, intent(in) :: lhs, rhs
    character(*), intent(in), optional :: description
    logical, intent(in), optional :: exit_on_fail
    type(string_t) :: desc
    logical :: ok
    ok = lhs == rhs
    desc = '';  if (present (description)) desc = var_str(description) // ": "
    call assert (unit, ok, char(desc // str (lhs) // " /= " // str (rhs)), exit_on_fail)
  end subroutine assert_equal_integer

  module subroutine assert_equal_integers (unit, lhs, rhs, description, exit_on_fail)
    integer, intent(in) :: unit
    integer, dimension(:), intent(in) :: lhs, rhs
    character(*), intent(in), optional :: description
    logical, intent(in), optional :: exit_on_fail
    type(string_t) :: desc
    logical :: ok
    ok = all(lhs == rhs)
    desc = '';  if (present (description)) desc = var_str(description) // ": "
    call assert (unit, ok, char(desc // str (lhs) // " /= " // str (rhs)), exit_on_fail)
  end subroutine assert_equal_integers

  module subroutine assert_equal_real (unit, lhs, rhs, description, &
                                abs_smallness, rel_smallness, exit_on_fail)
    integer, intent(in) :: unit
    real(default), intent(in) :: lhs, rhs
    character(*), intent(in), optional :: description
    real(default), intent(in), optional :: abs_smallness, rel_smallness
    logical, intent(in), optional :: exit_on_fail
    type(string_t) :: desc
    logical :: ok
    ok = nearly_equal (lhs, rhs, abs_smallness, rel_smallness)
    desc = '';  if (present (description)) desc = var_str(description) // ": "
    call assert (unit, ok, char(desc // str (lhs) // " /= " // str (rhs)), exit_on_fail)
  end subroutine assert_equal_real

  module subroutine assert_equal_reals (unit, lhs, rhs, description, &
                                abs_smallness, rel_smallness, exit_on_fail)
    integer, intent(in) :: unit
    real(default), dimension(:), intent(in) :: lhs, rhs
    character(*), intent(in), optional :: description
    real(default), intent(in), optional :: abs_smallness, rel_smallness
    logical, intent(in), optional :: exit_on_fail
    type(string_t) :: desc
    logical :: ok
    ok = all(nearly_equal (lhs, rhs, abs_smallness, rel_smallness))
    desc = '';  if (present (description)) desc = var_str(description) // ": "
    call assert (unit, ok, char(desc // str (lhs) // " /= " // str (rhs)), exit_on_fail)
  end subroutine assert_equal_reals

  module subroutine assert_equal_complex (unit, lhs, rhs, description, &
                                abs_smallness, rel_smallness, exit_on_fail)
    integer, intent(in) :: unit
    complex(default), intent(in) :: lhs, rhs
    character(*), intent(in), optional :: description
    real(default), intent(in), optional :: abs_smallness, rel_smallness
    logical, intent(in), optional :: exit_on_fail
    type(string_t) :: desc
    logical :: ok
    ok = nearly_equal (real(lhs), real(rhs), abs_smallness, rel_smallness) &
         .and. nearly_equal (aimag(lhs), aimag(rhs), abs_smallness, rel_smallness)
    desc = '';  if (present (description)) desc = var_str(description) // ": "
    call assert (unit, ok, char(desc // str (lhs) // " /= " // str (rhs)), exit_on_fail)
  end subroutine assert_equal_complex

  module subroutine assert_equal_complexs (unit, lhs, rhs, description, &
                                abs_smallness, rel_smallness, exit_on_fail)
    integer, intent(in) :: unit
    complex(default), dimension(:), intent(in) :: lhs, rhs
    character(*), intent(in), optional :: description
    real(default), intent(in), optional :: abs_smallness, rel_smallness
    logical, intent(in), optional :: exit_on_fail
    type(string_t) :: desc
    logical :: ok
    ok = all (nearly_equal (real(lhs), real(rhs), abs_smallness, rel_smallness)) &
         .and. all (nearly_equal (aimag(lhs), aimag(rhs), abs_smallness, rel_smallness))
    desc = '';  if (present (description)) desc = var_str(description) // ": "
    call assert (unit, ok, char(desc // str (lhs) // " /= " // str (rhs)), exit_on_fail)
  end subroutine assert_equal_complexs

  elemental function ieee_is_nan (x) result (yorn)
    logical :: yorn
    real(default), intent(in) :: x
    yorn = (x /= x)
  end function ieee_is_nan

  elemental module function nearly_equal_real &
       (a, b, abs_smallness, rel_smallness) result (r)
    logical :: r
    real(default), intent(in) :: a, b
    real(default), intent(in), optional :: abs_smallness, rel_smallness
    real(default) :: abs_a, abs_b, diff, abs_small, rel_small
    abs_a = abs (a)
    abs_b = abs (b)
    diff = abs (a - b)
    ! shortcut, handles infinities and nans
    if (a == b) then
       r = .true.
       return
    else if (ieee_is_nan (a) .or. ieee_is_nan (b) .or. ieee_is_nan (diff)) then
       r = .false.
       return
    end if
    abs_small = tiny_13; if (present (abs_smallness)) abs_small = abs_smallness
    rel_small = tiny_10; if (present (rel_smallness)) rel_small = rel_smallness
    if (abs_a < abs_small .and. abs_b < abs_small) then
       r = diff < abs_small
    else
       r = diff / max (abs_a, abs_b) < rel_small
    end if
  end function nearly_equal_real

  elemental module function nearly_equal_complex &
       (a, b, abs_smallness, rel_smallness) result (r)
    logical :: r
    complex(default), intent(in) :: a, b
    real(default), intent(in), optional :: abs_smallness, rel_smallness
    r = nearly_equal_real (real (a), real (b), abs_smallness, rel_smallness) .and. &
        nearly_equal_real (aimag (a), aimag(b), abs_smallness, rel_smallness)
  end function nearly_equal_complex

  elemental module function vanishes_real &
       (x, abs_smallness, rel_smallness) result (r)
    logical :: r
    real(default), intent(in) :: x
    real(default), intent(in), optional :: abs_smallness, rel_smallness
    r = nearly_equal (x, zero, abs_smallness, rel_smallness)
  end function vanishes_real

  elemental module function vanishes_complex &
       (x, abs_smallness, rel_smallness) result (r)
    logical :: r
    complex(default), intent(in) :: x
    real(default), intent(in), optional :: abs_smallness, rel_smallness
    r = vanishes_real (abs (x), abs_smallness, rel_smallness)
  end function vanishes_complex

  pure module function expanded_amp2 (amp_tree, amp_blob) result (amp2)
    real(default) :: amp2
    complex(default), dimension(:), intent(in) :: amp_tree, amp_blob
    amp2 = sum (amp_tree * conjg (amp_tree) + &
                amp_tree * conjg (amp_blob) + &
                amp_blob * conjg (amp_tree))
  end function expanded_amp2

  elemental module function abs2 (c) result (c2)
    real(default) :: c2
    complex(default), intent(in) :: c
    c2 = real (c * conjg(c))
  end function abs2

  module function remove_array_element_logical &
       (array, index) result (array_reduced)
    logical, intent(in), dimension(:) :: array
    integer, intent(in) :: index
    logical, dimension(:), allocatable :: array_reduced
    integer :: i
    allocate (array_reduced(0))
    do i = 1, size (array)
       if (i /= index) then
          array_reduced = [array_reduced, [array(i)]]
       end if
    end do
  end function remove_array_element_logical

  module function remove_duplicates_from_int_array &
       (array) result (array_unique)
    integer, intent(in), dimension(:) :: array
    integer, dimension(:), allocatable :: array_unique
    integer :: i
    allocate (array_unique(0))
    do i = 1, size (array)
       if (any (array_unique == array(i))) cycle
       array_unique = [array_unique, [array(i)]]
    end do
  end function remove_duplicates_from_int_array

  module subroutine extend_integer_array (list, incr, initial_value)
    integer, intent(inout), dimension(:), allocatable :: list
    integer, intent(in) :: incr
    integer, intent(in), optional :: initial_value
    integer, dimension(:), allocatable :: list_store
    integer :: n, ini
    ini = 0; if (present (initial_value)) ini = initial_value
    n = size (list)
    allocate (list_store (n))
    list_store = list
    deallocate (list)
    allocate (list (n+incr))
    list(1:n) = list_store
    list(1+n : n+incr) = ini
    deallocate (list_store)
  end subroutine extend_integer_array

  module subroutine crop_integer_array (list, i_crop)
    integer, intent(inout), dimension(:), allocatable :: list
    integer, intent(in) :: i_crop
    integer, dimension(:), allocatable :: list_store
    allocate (list_store (i_crop))
    list_store = list(1:i_crop)
    deallocate (list)
    allocate (list (i_crop))
    list = list_store
    deallocate (list_store)
  end subroutine crop_integer_array

  module function log_prec (x, xb) result (lx)
    real(default), intent(in) :: x, xb
    real(default) :: a1, a2, a3, lx
    a1 = xb
    a2 = a1 * xb / two
    a3 = a2 * xb * two / three
    if (abs (a3) < epsilon (a3)) then
       lx = - a1 - a2 - a3
    else
       lx = log (x)
    end if
  end function log_prec
  
  module subroutine split_integer_array (list1, list2)
    integer, intent(inout), dimension(:), allocatable :: list1, list2
    integer, dimension(:), allocatable :: list_store
    allocate (list_store (size (list1) - size (list2)))
    list2 = list1(:size (list2))
    list_store = list1 (size (list2) + 1:)
    deallocate (list1)
    allocate (list1 (size (list_store)))
    list1 = list_store
    deallocate (list_store)
  end subroutine split_integer_array

  module subroutine split_real_array (list1, list2)
    real(default), intent(inout), dimension(:), allocatable :: list1, list2
    real(default), dimension(:), allocatable :: list_store
    allocate (list_store (size (list1) - size (list2)))
    list2 = list1(:size (list2))
    list_store = list1 (size (list2) + 1:)
    deallocate (list1)
    allocate (list1 (size (list_store)))
    list1 = list_store
    deallocate (list_store)
  end subroutine split_real_array

  module function d1mach (i) result (d1)
    integer, intent(in) :: i
    real(default) :: b, x
    real(default) :: d1
    !***begin prologue  d1mach
    !***purpose  return floating point machine dependent constants.
    !***library   slatec
    !***category  r1
    !***type      single precision (d1mach-s, d1mach-d)
    !***keywords  machine constants
    !***author  fox, p. a., (bell labs)
    !           hall, a. d., (bell labs)
    !           schryer, n. l., (bell labs)
    !***description
    !
    !   d1mach can be used to obtain machine-dependent parameters for the
    !   local machine environment.  it is a function subprogram with one
    !   (input) argument, and can be referenced as follows:
    !
    !        a = d1mach(i)
    !
    !   where i=1,...,5.  the (output) value of a above is determined by
    !   the (input) value of i.  the results for various values of i are
    !   discussed below.
    !
    !   d1mach(1) = b**(emin-1), the smallest positive magnitude.
    !   d1mach(2) = b**emax*(1 - b**(-t)), the largest magnitude.
    !   d1mach(3) = b**(-t), the smallest relative spacing.
    !   d1mach(4) = b**(1-t), the largest relative spacing.
    !   d1mach(5) = log10(b)
    !
    !   assume single precision numbers are represented in the t-digit,
    !   base-b form
    !
    !              sign (b**e)*( (x(1)/b) + ... + (x(t)/b**t) )
    !
    !   where 0 .le. x(i) .lt. b for i=1,...,t, 0 .lt. x(1), and
    !   emin .le. e .le. emax.
    !
    !   the values of b, t, emin and emax are provided in i1mach as
    !   follows:
    !   i1mach(10) = b, the base.
    !   i1mach(11) = t, the number of base-b digits.
    !   i1mach(12) = emin, the smallest exponent e.
    !   i1mach(13) = emax, the largest exponent e.
    !
    !
    !***references  p. a. fox, a. d. hall and n. l. schryer, framework for
    !                 a portable library, acm transactions on mathematical
    !                 software 4, 2 (june 1978), pp. 177-188.
    !***routines called  xemsgr
    !***revision history  (yymmdd)
    !   790101  date written
    !   960329  modified for fortran 90 (be after suggestions by ehg)
    !***end prologue  d1mach
    !
    x = 1.0_default
    b = radix(x)
    select case (i)
    case (1)
       d1 = b**(minexponent(x)-1) ! the smallest positive magnitude.
    case (2)
       d1 = huge(X)               ! the largest magnitude.
    case (3)
       d1 = b**(-digits(x))       ! the smallest relative spacing.
    case (4)
       d1 = b**(1-digits(x))      ! the largest relative spacing.
    case (5)
       d1 = log10(b)
    case default
       d1 = b**(minexponent(x)-1)
    end select
  end function d1mach

  module subroutine dqk41 (f, a, b, result, abserr, resabs, resasc)
    !c***begin prologue  dqk41
    !c***date written   800101   (yymmdd)
    !c***revision date  830518   (yymmdd)
    !c***category no.  h2a1a2
    !c***keywords  41-point gauss-kronrod rules
    !c***author  piessens,robert,appl. math. & progr. div. - k.u.leuven
    !c           de doncker,elise,appl. math. & progr. div. - k.u.leuven
    !c***purpose  to compute i = integral of f over (a,b), with error
    !c                           estimate
    !c                       j = integral of abs(f) over (a,b)
    !c***description
    !c
    !c           integration rules
    !c           standard fortran subroutine
    !c           double precision version
    !c
    !c           parameters
    !c            on entry
    !c              f      - double precision
    !c                       function subprogram defining the integrand
    !c                       function f(x). the actual name for f needs to be
    !c                       declared e x t e r n a l in the calling program.
    !c
    !c              a      - double precision
    !c                       lower limit of integration
    !c
    !c              b      - double precision
    !c                       upper limit of integration
    !c
    !c            on return
    !c              result - double precision
    !c                       approximation to the integral i
    !c                       result is computed by applying the 41-point
    !c                       gauss-kronrod rule (resk) obtained by optimal
    !c                       addition of abscissae to the 20-point gauss
    !c                       rule (resg).
    !c
    !c              abserr - double precision
    !c                       estimate of the modulus of the absolute error,
    !c                       which should not exceed abs(i-result)
    !c
    !c              resabs - double precision
    !c                       approximation to the integral j
    !c
    !c              resasc - double precision
    !c                       approximation to the integal of abs(f-i/(b-a))
    !c                       over (a,b)
    !c
    !c***references  (none)
    !c***routines called  d1mach
    !c***end prologue  dqk41
    procedure(g_func) :: f
    real(default), intent(in) :: a, b
    real(default), intent(out) :: result, abserr, resabs, resasc
    real(default) :: absc, centr, dhlgth, dmax1, dmin1, &
         epmach, fc, fsum, fval1, fval2, hlgth, &
         resg, resk, reskh, uflow
    real(default), dimension(20) :: fv1, fv2
    real(default), dimension(21) :: xgk, wgk
    real(default), dimension(10) :: wg
    integer :: j, jtw, jtwm1
    !           the abscissae and weights are given for the interval (-1,1).
    !           because of symmetry only the positive abscissae and their
    !           corresponding weights are given.
    !
    !           xgk    - abscissae of the 41-point gauss-kronrod rule
    !                    xgk(2), xgk(4), ...  abscissae of the 20-point
    !                    gauss rule
    !                    xgk(1), xgk(3), ...  abscissae which are optimally
    !                    added to the 20-point gauss rule
    !
    !           wgk    - weights of the 41-point gauss-kronrod rule
    !
    !           wg     - weights of the 20-point gauss rule
    !
    !
    ! gauss quadrature weights and kronron quadrature abscissae and weights
    ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    ! bell labs, nov. 1981.
    !
    data wg  (  1) / 0.017614007139152118311861962351853_default /
    data wg  (  2) / 0.040601429800386941331039952274932_default /
    data wg  (  3) / 0.062672048334109063569506535187042_default /
    data wg  (  4) / 0.083276741576704748724758143222046_default /
    data wg  (  5) / 0.101930119817240435036750135480350_default /
    data wg  (  6) / 0.118194531961518417312377377711382_default /
    data wg  (  7) / 0.131688638449176626898494499748163_default /
    data wg  (  8) / 0.142096109318382051329298325067165_default /
    data wg  (  9) / 0.149172986472603746787828737001969_default /
    data wg  ( 10) / 0.152753387130725850698084331955098_default /

    data xgk (  1) / 0.998859031588277663838315576545863_default /
    data xgk (  2) / 0.993128599185094924786122388471320_default /
    data xgk (  3) / 0.981507877450250259193342994720217_default /
    data xgk (  4) / 0.963971927277913791267666131197277_default /
    data xgk (  5) / 0.940822633831754753519982722212443_default /
    data xgk (  6) / 0.912234428251325905867752441203298_default /
    data xgk (  7) / 0.878276811252281976077442995113078_default /
    data xgk (  8) / 0.839116971822218823394529061701521_default /
    data xgk (  9) / 0.795041428837551198350638833272788_default /
    data xgk ( 10) / 0.746331906460150792614305070355642_default /
    data xgk ( 11) / 0.693237656334751384805490711845932_default /
    data xgk ( 12) / 0.636053680726515025452836696226286_default /
    data xgk ( 13) / 0.575140446819710315342946036586425_default /
    data xgk ( 14) / 0.510867001950827098004364050955251_default /
    data xgk ( 15) / 0.443593175238725103199992213492640_default /
    data xgk ( 16) / 0.373706088715419560672548177024927_default /
    data xgk ( 17) / 0.301627868114913004320555356858592_default /
    data xgk ( 18) / 0.227785851141645078080496195368575_default /
    data xgk ( 19) / 0.152605465240922675505220241022678_default /
    data xgk ( 20) / 0.076526521133497333754640409398838_default /
    data xgk ( 21) / 0.000000000000000000000000000000000_default /

    data wgk (  1) / 0.003073583718520531501218293246031_default /
    data wgk (  2) / 0.008600269855642942198661787950102_default /
    data wgk (  3) / 0.014626169256971252983787960308868_default /
    data wgk (  4) / 0.020388373461266523598010231432755_default /
    data wgk (  5) / 0.025882133604951158834505067096153_default /
    data wgk (  6) / 0.031287306777032798958543119323801_default /
    data wgk (  7) / 0.036600169758200798030557240707211_default /
    data wgk (  8) / 0.041668873327973686263788305936895_default /
    data wgk (  9) / 0.046434821867497674720231880926108_default /
    data wgk ( 10) / 0.050944573923728691932707670050345_default /
    data wgk ( 11) / 0.055195105348285994744832372419777_default /
    data wgk ( 12) / 0.059111400880639572374967220648594_default /
    data wgk ( 13) / 0.062653237554781168025870122174255_default /
    data wgk ( 14) / 0.065834597133618422111563556969398_default /
    data wgk ( 15) / 0.068648672928521619345623411885368_default /
    data wgk ( 16) / 0.071054423553444068305790361723210_default /
    data wgk ( 17) / 0.073030690332786667495189417658913_default /
    data wgk ( 18) / 0.074582875400499188986581418362488_default /
    data wgk ( 19) / 0.075704497684556674659542775376617_default /
    data wgk ( 20) / 0.076377867672080736705502835038061_default /
    data wgk ( 21) / 0.076600711917999656445049901530102_default /
    !
    !
    !           list of major variables
    !           -----------------------
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc   - abscissa
    !           fval*  - function value
    !           resg   - result of the 20-point gauss formula
    !           resk   - result of the 41-point kronrod formula
    !           reskh  - approximation to mean value of f over (a,b), i.e.
    !                    to i/(b-a)
    !
    !           machine dependent constants
    !           ---------------------------
    !
    !           epmach is the largest relative spacing.
    !           uflow is the smallest positive magnitude.
    !
    !***first executable statement  dqk41
    epmach = d1mach(4)
    uflow = d1mach(1)

    centr = 0.5d+00*(a+b)
    hlgth = 0.5d+00*(b-a)
    dhlgth = abs(hlgth)
    !
    ! compute the 41-point gauss-kronrod approximation to
    ! the integral, and estimate the absolute error.
    !
    resg = 0.0d+00
    fc = f(centr)
    resk = wgk(21)*fc
    resabs = abs(resk)
    do j = 1, 10
       jtw = j*2
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do
    do j = 1,10
       jtwm1 = j*2-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
     end do
     reskh = resk*0.5d+00
     resasc = wgk(21)*abs(fc-reskh)
     do j = 1, 20
        resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
     end do
     result = resk*hlgth
     resabs = resabs*dhlgth
     resasc = resasc*dhlgth
     abserr = abs((resk-resg)*hlgth)
     if (resasc.ne.0.0d+00.and.abserr.ne.0.d+00) &
          abserr = resasc*dmin1(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
     if (resabs.gt.uflow/(0.5d+02*epmach)) abserr = dmax1 &
          ((epmach*0.5d+02)*resabs,abserr)
   end subroutine dqk41

  module subroutine dqk61 (f, a, b, result, abserr, resabs, resasc)
    !c***begin prologue  dqk61
    !c***purpose  to compute i = integral of f over (a,b) with error
    !c                           estimate
    !c                       j = integral of abs(f) over (a,b)
    !c***library   slatec (quadpack)
    !c***category  h2a1a2
    !c***type      double precision (qk61-s, dqk61-d)
    !c***keywords  61-point gauss-kronrod rules, quadpack, quadrature
    !c***author  piessens, robert
    !c             applied mathematics and programming division
    !c             k. u. leuven
    !c           de doncker, elise
    !c             applied mathematics and programming division
    !c             k. u. leuven
    !c***description
    !c
    !c        integration rule
    !c        standard fortran subroutine
    !c        double precision version
    !c
    !c
    !c        parameters
    !c         on entry
    !c           f      - double precision
    !c                    function subprogram defining the integrand
    !c                    function f(x). the actual name for f needs to be
    !c                    declared e x t e r n a l in the calling program.
    !c
    !c           a      - double precision
    !c                    lower limit of integration
    !c
    !c           b      - double precision
    !c                    upper limit of integration
    !c
    !c         on return
    !c           result - double precision
    !c                    approximation to the integral i
    !c                    result is computed by applying the 61-point
    !c                    kronrod rule (resk) obtained by optimal addition of
    !c                    abscissae to the 30-point gauss rule (resg).
    !c
    !c           abserr - double precision
    !c                    estimate of the modulus of the absolute error,
    !c                    which should equal or exceed abs(i-result)
    !c
    !c           resabs - double precision
    !c                    approximation to the integral j
    !c
    !c           resasc - double precision
    !c                    approximation to the integral of abs(f-i/(b-a))
    !c
    !c***references  (none)
    !c***routines called  d1mach
    !c***revision history  (yymmdd)
    !c   800101  date written
    !c   890531  changed all specific intrinsics to generic.  (wrb)
    !c   890531  revision date from version 3.2
    !c   891214  prologue converted to version 4.0 format.  (bab)
    !c***end prologue  dqk61
    !c
    procedure(g_func) :: f
    real(default), intent(in) :: a, b
    real(default), intent(out) :: result, abserr, resabs, resasc
    real(default) ::  absc, centr, dhlgth, &
         epmach, fc, fsum, fval1, fval2, hlgth, &
         resg, resk, reskh, uflow
    real(default), dimension(30) :: fv1, fv2
    real(default), dimension(31) :: xgk, wgk
    real(default), dimension(15) :: wg
    integer :: j, jtw, jtwm1
    !
    !           the abscissae and weights are given for the
    !           interval (-1,1). because of symmetry only the positive
    !           abscissae and their corresponding weights are given.
    !
    !           xgk   - abscissae of the 61-point kronrod rule
    !                   xgk(2), xgk(4)  ... abscissae of the 30-point
    !                   gauss rule
    !                   xgk(1), xgk(3)  ... optimally added abscissae
    !                   to the 30-point gauss rule
    !
    !           wgk   - weights of the 61-point kronrod rule
    !
    !           wg    - weights of the 30-point gauss rule
    !
    !
    ! gauss quadrature weights and kronrod quadrature abscissae and weights
    ! as evaluated with 80 decimal digit arithmetic by l. w. fullerton,
    ! bell labs, nov. 1981.
    !
    data wg  (  1) / 0.007968192496166605615465883474674_default /
    data wg  (  2) / 0.018466468311090959142302131912047_default /
    data wg  (  3) / 0.028784707883323369349719179611292_default /
    data wg  (  4) / 0.038799192569627049596801936446348_default /
    data wg  (  5) / 0.048402672830594052902938140422808_default /
    data wg  (  6) / 0.057493156217619066481721689402056_default /
    data wg  (  7) / 0.065974229882180495128128515115962_default /
    data wg  (  8) / 0.073755974737705206268243850022191_default /
    data wg  (  9) / 0.080755895229420215354694938460530_default /
    data wg  ( 10) / 0.086899787201082979802387530715126_default /
    data wg  ( 11) / 0.092122522237786128717632707087619_default /
    data wg  ( 12) / 0.096368737174644259639468626351810_default /
    data wg  ( 13) / 0.099593420586795267062780282103569_default /
    data wg  ( 14) / 0.101762389748405504596428952168554_default /
    data wg  ( 15) / 0.102852652893558840341285636705415_default /

    data xgk (  1) / 0.999484410050490637571325895705811_default /
    data xgk (  2) / 0.996893484074649540271630050918695_default /
    data xgk (  3) / 0.991630996870404594858628366109486_default /
    data xgk (  4) / 0.983668123279747209970032581605663_default /
    data xgk (  5) / 0.973116322501126268374693868423707_default /
    data xgk (  6) / 0.960021864968307512216871025581798_default /
    data xgk (  7) / 0.944374444748559979415831324037439_default /
    data xgk (  8) / 0.926200047429274325879324277080474_default /
    data xgk (  9) / 0.905573307699907798546522558925958_default /
    data xgk ( 10) / 0.882560535792052681543116462530226_default /
    data xgk ( 11) / 0.857205233546061098958658510658944_default /
    data xgk ( 12) / 0.829565762382768397442898119732502_default /
    data xgk ( 13) / 0.799727835821839083013668942322683_default /
    data xgk ( 14) / 0.767777432104826194917977340974503_default /
    data xgk ( 15) / 0.733790062453226804726171131369528_default /
    data xgk ( 16) / 0.697850494793315796932292388026640_default /
    data xgk ( 17) / 0.660061064126626961370053668149271_default /
    data xgk ( 18) / 0.620526182989242861140477556431189_default /
    data xgk ( 19) / 0.579345235826361691756024932172540_default /
    data xgk ( 20) / 0.536624148142019899264169793311073_default /
    data xgk ( 21) / 0.492480467861778574993693061207709_default /
    data xgk ( 22) / 0.447033769538089176780609900322854_default /
    data xgk ( 23) / 0.400401254830394392535476211542661_default /
    data xgk ( 24) / 0.352704725530878113471037207089374_default /
    data xgk ( 25) / 0.304073202273625077372677107199257_default /
    data xgk ( 26) / 0.254636926167889846439805129817805_default /
    data xgk ( 27) / 0.204525116682309891438957671002025_default /
    data xgk ( 28) / 0.153869913608583546963794672743256_default /
    data xgk ( 29) / 0.102806937966737030147096751318001_default /
    data xgk ( 30) / 0.051471842555317695833025213166723_default /
    data xgk ( 31) / 0.000000000000000000000000000000000_default /
    
    data wgk (  1) / 0.001389013698677007624551591226760_default /
    data wgk (  2) / 0.003890461127099884051267201844516_default /
    data wgk (  3) / 0.006630703915931292173319826369750_default /
    data wgk (  4) / 0.009273279659517763428441146892024_default /
    data wgk (  5) / 0.011823015253496341742232898853251_default /
    data wgk (  6) / 0.014369729507045804812451432443580_default /
    data wgk (  7) / 0.016920889189053272627572289420322_default /
    data wgk (  8) / 0.019414141193942381173408951050128_default /
    data wgk (  9) / 0.021828035821609192297167485738339_default /
    data wgk ( 10) / 0.024191162078080601365686370725232_default /
    data wgk ( 11) / 0.026509954882333101610601709335075_default /
    data wgk ( 12) / 0.028754048765041292843978785354334_default /
    data wgk ( 13) / 0.030907257562387762472884252943092_default /
    data wgk ( 14) / 0.032981447057483726031814191016854_default /
    data wgk ( 15) / 0.034979338028060024137499670731468_default /
    data wgk ( 16) / 0.036882364651821229223911065617136_default /
    data wgk ( 17) / 0.038678945624727592950348651532281_default /
    data wgk ( 18) / 0.040374538951535959111995279752468_default /
    data wgk ( 19) / 0.041969810215164246147147541285970_default /
    data wgk ( 20) / 0.043452539701356069316831728117073_default /
    data wgk ( 21) / 0.044814800133162663192355551616723_default /
    data wgk ( 22) / 0.046059238271006988116271735559374_default /
    data wgk ( 23) / 0.047185546569299153945261478181099_default /
    data wgk ( 24) / 0.048185861757087129140779492298305_default /
    data wgk ( 25) / 0.049055434555029778887528165367238_default /
    data wgk ( 26) / 0.049795683427074206357811569379942_default /
    data wgk ( 27) / 0.050405921402782346840893085653585_default /
    data wgk ( 28) / 0.050881795898749606492297473049805_default /
    data wgk ( 29) / 0.051221547849258772170656282604944_default /
    data wgk ( 30) / 0.051426128537459025933862879215781_default /
    data wgk ( 31) / 0.051494729429451567558340433647099_default /
    !
    !           list of major variables
    !           -----------------------
    !
    !           centr  - mid point of the interval
    !           hlgth  - half-length of the interval
    !           absc  - abscissa
    !           fval*  - function value
    !           resg   - result of the 30-point gauss rule
    !           resk   - result of the 61-point kronrod rule
    !           reskh  - approximation to the mean value of f
    !                    over (a,b), i.e. to i/(b-a)
    !
    !           machine dependent constants
    !           ---------------------------
    !
    !           epmach is the largest relative spacing.
    !           uflow is the smallest positive magnitude.
    !
    !***first executable statement  dqk61
    epmach = d1mach(4)
    uflow = d1mach(1)

    centr = 0.5d+00*(b+a)
    hlgth = 0.5d+00*(b-a)
    dhlgth = abs(hlgth)
    !
    ! compute the 61-point kronrod approximation to the
    ! integral, and estimate the absolute error.
    !
    resg = 0.0d+00
    fc = f(centr)
    resk = wgk(31)*fc
    resabs = abs(resk)
    do j = 1, 15
       jtw = j*2
       absc = hlgth*xgk(jtw)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtw) = fval1
       fv2(jtw) = fval2
       fsum = fval1+fval2
       resg = resg+wg(j)*fsum
       resk = resk+wgk(jtw)*fsum
       resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do
    do j = 1, 15
       jtwm1 = j*2-1
       absc = hlgth*xgk(jtwm1)
       fval1 = f(centr-absc)
       fval2 = f(centr+absc)
       fv1(jtwm1) = fval1
       fv2(jtwm1) = fval2
       fsum = fval1+fval2
       resk = resk+wgk(jtwm1)*fsum
       resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
     end do
     reskh = resk*0.5d+00
     resasc = wgk(31)*abs(fc-reskh)
     do j = 1, 30
        RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
     end do
     result = resk*hlgth
     resabs = resabs*dhlgth
     resasc = resasc*dhlgth
     abserr = abs((resk-resg)*hlgth)
     if (resasc.ne.0.0d+00.and.abserr.ne.0.0d+00) &
          abserr = resasc*min(0.1d+01,(0.2d+03*abserr/resasc)**1.5d+00)
     if (resabs.gt.uflow/(0.5d+02*epmach)) abserr = max &
          ((epmach*0.5d+02)*resabs,abserr)
   end subroutine dqk61

  module subroutine gauss_kronrod &
       (type, f, a, b, limit, result, abserr, epsabs, epsrel)
    integer, intent(in) :: type
    procedure(g_func) :: f
    real(default), intent(in) :: a, b
    real(default), intent(in) :: epsabs, epsrel
    real(default), intent(out) :: result, abserr
    integer, intent(in) :: limit
    real(default) :: area, errsum
    real(default) :: tolerance_g, round_off
    real(default) :: result0, abserr0, resabs0, resasc0
    type(int_workspace_t) :: work
    real(default) :: a1, b1, a2, b2, a_i, b_i, r_i, e_i
    real(default) :: area1, area2, area12, error1, error2, error12
    real(default) :: resasc1, resasc2, resabs1, resabs2, delta
    integer :: i, iteration
    integer :: roundoff_type1, roundoff_type2
    !!! epsilon(dbl)
    real(default), parameter :: dbl_prec = 2.2204460492503131e-16_default
    !!! tiny(dbl)
    real(default), parameter :: dbl_min = 2.2250738585072014e-308_default
    iteration = 0
    if (epsabs <= 0) then
       if (epsrel < 50 * dbl_prec .or. epsrel < 0.5e-28_default) then
          error stop ("Gauss_Kronrod: tolerance cannot be achieved with given " &
               // "epsabs and epsrel")
       end if
    end if
    call work%init (a, b, limit)
    ! Perform the first integration
    select case (type)
    case (GAUSS_KRONROD_61)
       call dqk61 (f, a, b, result0, abserr0, resabs0, resasc0)
    case default
       call dqk41 (f, a, b, result0, abserr0, resabs0, resasc0)
    end select
    call work%set_initial (result0, abserr0)
    ! Test on accuracy
    tolerance_g = max (epsabs, epsrel * abs(result0))
    round_off = 50 * dbl_prec * resabs0
    if (abserr0 <= round_off .and. abserr0 > tolerance_g) then
       result = result0
       abserr = abserr0
       error stop ("Gauss_Kronrod: cannot reach tolerance because of roundoff " &
            // "error on first attempt")
    else if ((abserr0 <= tolerance_g .and. abserr0 /= resasc0) .or. &
         abserr0 == 0.0_default) then
       result = result0
       abserr = abserr0
       return
    else if (limit == 1) then
       result = result0
       abserr = abserr0
       error stop ("Gauss_Kronrod: a maximum of one iteration was insufficient")
    end if
    area = result0
    errsum = abserr0
    roundoff_type1 = 0
    roundoff_type2 = 0
    ITERATION_LOOP: do iteration = 1, limit
       a_i = work%alist(work%i)
       b_i = work%blist(work%i)
       r_i = work%rlist(work%i)
       e_i = work%elist(work%i)
       a1 = a_i
       b1 = 0.5_default * (a_i + b_i)
       a2 = b1
       b2 = b_i
       select case (type)
       case (GAUSS_KRONROD_61)
          call dqk61 (f, a1, b1, area1, error1, resabs1, resasc1)
          call dqk61 (f, a2, b2, area2, error2, resabs2, resasc2)
       case default
          call dqk41 (f, a1, b1, area1, error1, resabs1, resasc1)
          call dqk41 (f, a2, b2, area2, error2, resabs2, resasc2)
       end select
       area12 = area1 + area2
       error12 = error1 + error2
       errsum = errsum + error12 - e_i
       area = area + area12 - r_i
       if (resasc1 /= error1 .and. resasc2 /= error2) then
          delta = r_i - area12
          if (abs(delta) <= 1.0e-5_default * abs (area12) .and. &
               error12 >= 0.99_default * e_i) then
             roundoff_type1 = roundoff_type1 + 1
          end if
          if (iteration >= 10 .and. error12 > e_i) then
             roundoff_type2 = roundoff_type2 + 1
          end if
       end if
       tolerance_g = max (epsabs, epsrel * abs(area))
       if (errsum > tolerance_g) then
          if (roundoff_type1 >= 6 .or. roundoff_type2 >= 20) &
               error stop ("Gauss_Kronrod: roundoff error prevents tolerance " &
               // "from being achieved")
          ! Set error flag in the case of bad integrand behavior at a
          ! point of the integration range
          if (subinterval_too_small (a1, a2, b2)) &
               error stop ("Gauss_Kronrod: bad integrand behavior found in " &
               // "integration interval")
       end if
       call work%update (a1, b1, area1, error1, a2, b2, area2, error2)
       if (errsum <= tolerance_g)  exit ITERATION_LOOP
    end do ITERATION_LOOP
    result = sum (work%rlist (1:work%size))
    abserr = errsum
    if (errsum > tolerance_g) &
       error stop ("Gauss_Kronrod: could not integrate function")
  contains
    function subinterval_too_small (a_1, a_2, b_2) result (flag)
      real(default), intent(in) :: a_1, a_2, b_2
      logical :: flag
      real(default) :: tmp
      tmp = (1._default + 100._default * dbl_prec) * (abs(a_2) + &
           1000._default * dbl_min)
      flag = abs(a_1) <= tmp .and. abs(b_2) <= tmp
    end function subinterval_too_small
  end subroutine gauss_kronrod

  elemental module subroutine pacify_real_default (x, tolerance)
    real(default), intent(inout) :: x
    real(default), intent(in) :: tolerance
    if (abs (x) < tolerance)  x = 0._default
  end subroutine pacify_real_default

  elemental module subroutine pacify_complex_default (x, tolerance)
    complex(default), intent(inout) :: x
    real(default), intent(in) :: tolerance
    if (abs (real (x)) < tolerance)   &
         x = cmplx (0._default, aimag (x), kind=default)
    if (abs (aimag (x)) < tolerance)  &
         x = cmplx (real (x), 0._default, kind=default)
  end subroutine pacify_complex_default


end submodule numeric_utils_s

