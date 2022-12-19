! divisions.f90 --
! Copyright (C) 1998 by Thorsten Ohl <ohl@hep.tu-darmstadt.de>
! 
! VAMP is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by 
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
! 
! VAMP is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This version of the source code of vamp has no comments and
! can be hard to understand, modify, and improve.  You should have
! received a copy of the literate `noweb' sources of vamp that
! contain the documentation in full detail.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module divisions
  use kinds
  use exceptions
  use vamp_stat
  use utils
  use iso_fortran_env
  implicit none
  private
  public :: create_division, create_empty_division
  public :: copy_division, delete_division
  public :: set_rigid_division, reshape_division
  public :: inject_division, inject_division_short
  public :: record_integral, record_variance, clear_integral_and_variance
  ! public :: record_efficiency
  public :: refine_division
  private :: rebinning_weights
  private :: rebin
  public :: probability
  public :: quadrupole_division
  public :: fork_division, join_division, sum_division
  private :: subdivide
  private :: distribute
  private :: collect
  public :: debug_division
  public :: dump_division
  public :: inside_division, stratified_division
  public :: volume_division, rigid_division, adaptive_division
  public :: copy_history, summarize_division
  private :: probabilities
  public :: print_history, write_history
  public :: write_division
  private :: write_division_unit, write_division_name
  public :: read_division
  private :: read_division_unit, read_division_name
  public :: write_division_raw
  private :: write_division_raw_unit, write_division_raw_name
  public :: read_division_raw
  private :: read_division_raw_unit, read_division_raw_name
  public :: marshal_division_size, marshal_division, unmarshal_division
  public :: marshal_div_history_size, marshal_div_history, unmarshal_div_history
  interface write_division
     module procedure write_division_unit, write_division_name
  end interface
  interface read_division
     module procedure read_division_unit, read_division_name
  end interface
  interface write_division_raw
     module procedure write_division_raw_unit, write_division_raw_name
  end interface
  interface read_division_raw
     module procedure read_division_raw_unit, read_division_raw_name
  end interface
  integer, private, parameter :: MIN_NUM_DIV = 3
  integer, private, parameter :: BUFFER_SIZE = 50
  character(len=*), parameter, private :: &
       descr_fmt =        "(1x,a)", &
       integer_fmt =      "(1x,a15,1x,i15)", &
       logical_fmt =      "(1x,a15,1x,l1)", &
       double_fmt =       "(1x,a15,1x,e30.22)", &
       double_array_fmt = "(1x,i15,1x,3(e30.22))"
  type, public :: division_t
  !   private
  !!! Avoiding a g95 bug 
     real(kind=default), dimension(:), pointer :: x => null ()
     real(kind=default), dimension(:), pointer :: integral => null ()
     real(kind=default), dimension(:), pointer &
                                        :: variance => null ()
  !                                      public :: variance => null ()
  !  real(kind=default), dimension(:), pointer :: efficiency => null ()
     real(kind=default) :: x_min, x_max
     real(kind=default) :: x_min_true, x_max_true
     real(kind=default) :: dx, dxg
     integer :: ng = 0
     logical :: stratified = .true.
  end type division_t
  type, public :: div_history
     private
     logical :: stratified
     integer :: ng, num_div
     real(kind=default) :: x_min, x_max, x_min_true, x_max_true
     real(kind=default) :: &
          spread_f_p, stddev_f_p, spread_p, stddev_p, spread_m, stddev_m
  end type div_history
  integer, parameter, private :: MAGIC_DIVISION = 11111111
  integer, parameter, private :: MAGIC_DIVISION_BEGIN = MAGIC_DIVISION + 1
  integer, parameter, private :: MAGIC_DIVISION_END = MAGIC_DIVISION + 2
contains
  elemental subroutine create_division &
       (d, x_min, x_max, x_min_true, x_max_true)
    type(division_t), intent(out) :: d
    real(kind=default), intent(in) :: x_min, x_max
    real(kind=default), intent(in), optional :: x_min_true, x_max_true
    allocate (d%x(0:1), d%integral(1), d%variance(1))
  ! allocate (d%efficiency(1))
    d%x(0) = 0.0
    d%x(1) = 1.0
    d%x_min = x_min
    d%x_max = x_max
    d%dx = d%x_max - d%x_min
    d%stratified = .false.
    d%ng = 1
    d%dxg = 1.0 / d%ng
    if (present (x_min_true)) then
       d%x_min_true = x_min_true
    else
       d%x_min_true = x_min
    end if
    if (present (x_max_true)) then
       d%x_max_true = x_max_true
    else
       d%x_max_true = x_max
    end if
  end subroutine create_division

  elemental subroutine create_empty_division (d)
    type(division_t), intent(out) :: d
    nullify (d%x, d%integral, d%variance)
  ! nullify (d%efficiency)
  end subroutine create_empty_division

  elemental subroutine set_rigid_division (d, ng)
    type(division_t), intent(inout) :: d
    integer, intent(in) :: ng
    d%stratified = ng > 1
    d%ng = ng
    d%dxg = real (ubound (d%x, dim=1), kind=default) / d%ng
  end subroutine set_rigid_division

  elemental subroutine reshape_division (d, max_num_div, ng, use_variance)
    type(division_t), intent(inout) :: d
    integer, intent(in) :: max_num_div
    integer, intent(in), optional :: ng
    logical, intent(in), optional :: use_variance
    real(kind=default), dimension(:), allocatable :: old_x, m
    integer :: num_div, equ_per_adap
    if (present (ng)) then
       if (max_num_div > 1) then
          d%stratified = ng > 1
       else
          d%stratified = .false.
       end if
    else
       d%stratified = .false.
    end if
    if (d%stratified) then
       d%ng = ng
       if (d%ng >= max_num_div / 2) then 
          d%stratified = .true.
          equ_per_adap = d%ng / max_num_div + 1
          num_div = d%ng / equ_per_adap
          if (num_div < 2) then
             d%stratified = .false.
             num_div = 2
             d%ng = 1
          else if (mod (num_div,2) == 1) then
             num_div = num_div - 1
             d%ng = equ_per_adap * num_div
          else
             d%ng = equ_per_adap * num_div
          end if
       else
          d%stratified = .false.
          num_div = max_num_div
          d%ng = 1
       end if
    else
       num_div = max_num_div
       d%ng = 1
    end if
    d%dxg = real (num_div, kind=default) / d%ng
    allocate (old_x(0:ubound(d%x,dim=1)), m(ubound(d%x,dim=1)))
    old_x = d%x
    if (present (use_variance)) then
       if (use_variance) then
          m = rebinning_weights (d%variance)
       else
          m = 1.0
       end if
    else
       m = 1.0
    end if
    if (ubound (d%x, dim=1) /= num_div) then
       deallocate (d%x, d%integral, d%variance)
    !  deallocate (d%efficiency)
       allocate (d%x(0:num_div), d%integral(num_div), d%variance(num_div))
    !  allocate (d%efficiency(num_div))
    end if
    d%x = rebin (m, old_x, num_div)
    deallocate (old_x, m)
  end subroutine reshape_division

  elemental subroutine inject_division (d, r, cell, x, x_mid, idx, wgt)
    type(division_t), intent(in) :: d
    real(kind=default), intent(in) :: r
    integer, intent(in) :: cell
    real(kind=default), intent(out) :: x, x_mid
    integer, intent(out) :: idx
    real(kind=default), intent(out) :: wgt
    real(kind=default) :: delta_x, xi
    integer :: i
    xi = (cell - r) * d%dxg + 1.0
    i = max (min (int (xi), ubound (d%x, dim=1)), 1)
    delta_x = d%x(i) - d%x(i-1)
    x = d%x_min + (d%x(i-1) + (xi - i) * delta_x) * d%dx
    wgt = delta_x * ubound (d%x, dim=1)
    idx = i
    x_mid = d%x_min + 0.5 * (d%x(i-1) + d%x(i)) * d%dx
  end subroutine inject_division

  elemental subroutine inject_division_short (d, r, x, idx, wgt)
    type(division_t), intent(in) :: d
    real(kind=default), intent(in) :: r
    integer, intent(out) :: idx
    real(kind=default), intent(out) :: x, wgt
    real(kind=default) :: delta_x, xi
    integer :: i
    xi = r * ubound (d%x, dim=1) + 1.0
    i = max (min (int (xi), ubound (d%x, dim=1)), 1)
    delta_x = d%x(i) - d%x(i-1)
    x = d%x_min + (d%x(i-1) + (xi - i) * delta_x) * d%dx
    wgt = delta_x * ubound (d%x, dim=1)
    idx = i
  end subroutine inject_division_short

  elemental subroutine record_integral (d, i, f)
    type(division_t), intent(inout) :: d
    integer, intent(in) :: i
    real(kind=default), intent(in) :: f
    d%integral(i) = d%integral(i) + f
    if (.not. d%stratified) then 
       d%variance(i) = d%variance(i) + f*f
    end if
  end subroutine record_integral

  elemental subroutine record_variance (d, i, var_f)
    type(division_t), intent(inout) :: d
    integer, intent(in) :: i
    real(kind=default), intent(in) :: var_f
    if (d%stratified) then 
       d%variance(i) = d%variance(i) + var_f
    end if
  end subroutine record_variance

  elemental subroutine clear_integral_and_variance (d)
    type(division_t), intent(inout) :: d
    d%integral = 0.0
    d%variance = 0.0
  ! d%efficiency = 0.0
  end subroutine clear_integral_and_variance
  elemental subroutine refine_division (d)
    type(division_t), intent(inout) :: d
    character(len=*), parameter :: FN = "refine_division"
    d%x = rebin (rebinning_weights (d%variance), d%x, size (d%variance))
  end subroutine refine_division
  pure function rebinning_weights (d) result (m)
    real(kind=default), dimension(:), intent(in) :: d
    real(kind=default), dimension(size(d)) :: m
    real(kind=default), dimension(size(d)) :: smooth_d
    real(kind=default), parameter :: ALPHA = 1.5
    integer :: nd
    if (any (d /= d)) then
       m = 1.0
       return
    end if
    nd = size (d)
    if (nd > 2) then
       smooth_d(1) = (d(1) + d(2)) / 2.0
       smooth_d(2:nd-1) = (d(1:nd-2) + d(2:nd-1) + d(3:nd)) / 3.0
       smooth_d(nd) = (d(nd-1) + d(nd)) / 2.0
    else
       smooth_d = d
    end if
    if (all (smooth_d < tiny (1.0_default))) then
       m = 1.0_default
    else
       smooth_d = smooth_d / sum (smooth_d)
       where (smooth_d < tiny (1.0_default))
          smooth_d = tiny (1.0_default)
       end where
       where (smooth_d /= 1._default)
          m = ((smooth_d - 1.0) / (log (smooth_d)))**ALPHA
       elsewhere
          m = 1.0_default
       endwhere
    end if
  end function rebinning_weights
  pure function rebin (m, x, num_div) result (x_new)
    real(kind=default), dimension(:), intent(in) :: m
    real(kind=default), dimension(0:), intent(in) :: x
    integer, intent(in) :: num_div
    real(kind=default), dimension(0:num_div) :: x_new
    integer :: i, k
    real(kind=default) :: step, delta
    step = sum (m) / num_div
    k = 0
    delta = 0.0
    x_new(0) = x(0)
    do i = 1, num_div - 1
       do
          if (step <= delta) then 
             exit
          end if
          k = k + 1
          delta = delta + m(k)
       end do
       delta = delta - step
       x_new(i) = x(k) - (x(k) - x(k-1)) * delta / m(k)
    end do
    x_new(num_div) = 1.0
  end function rebin
  elemental function probability (d, x) result (p)
    type(division_t), intent(in) :: d
    real(kind=default), intent(in) :: x
    real(kind=default) :: p
    real(kind=default) :: xi
    integer :: hi, mid, lo
    xi = (x - d%x_min) / d%dx
    if ((xi >= 0) .and. (xi <= 1)) then
       lo = lbound (d%x, dim=1)
       hi = ubound (d%x, dim=1)
       bracket: do
          if (lo >= hi - 1) then
             p = 1.0 / (ubound (d%x, dim=1) * d%dx * (d%x(hi) - d%x(hi-1)))
             return
          end if
          mid = (hi + lo) / 2
          if (xi > d%x(mid)) then
             lo = mid
          else
             hi = mid
          end if
       end do bracket
    else
       p = 0
    end if
  end function probability
  elemental function quadrupole_division (d) result (q)
    type(division_t), intent(in) :: d
    real(kind=default) :: q
    !!!   q = value_spread_percent (rebinning_weights (d%variance))
    q = standard_deviation_percent (rebinning_weights (d%variance))
  end function quadrupole_division
  pure subroutine fork_division (d, ds, sum_calls, num_calls, exc)
    type(division_t), intent(in) :: d
    type(division_t), dimension(:), intent(inout) :: ds
    integer, intent(in) :: sum_calls
    integer, dimension(:), intent(inout) :: num_calls
    type(exception), intent(inout), optional :: exc
    character(len=*), parameter :: FN = "fork_division"
    integer, dimension(size(ds)) :: n0, n1
    integer, dimension(0:size(ds)) :: n, ds_ng
    integer :: i, j, num_div, num_forks, nx
    real(kind=default), dimension(:), allocatable :: d_x, d_integral, d_variance
  ! real(kind=default), dimension(:), allocatable :: d_efficiency
    num_div = ubound (d%x, dim=1)
    num_forks = size (ds)
    if (d%ng == 1) then
       if (d%stratified) then
          call raise_exception (exc, EXC_FATAL, FN, &
                                "ng == 1 incompatiple w/ stratification")
       else
          call copy_division (ds, d)
          num_calls(2:) = ceiling (real (sum_calls) / num_forks)
          num_calls(1) = sum_calls - sum (num_calls(2:))
       end if
    else if (num_div >= num_forks) then
       if (modulo (d%ng, num_div) == 0) then
          n = (num_div * (/ (j, j=0,num_forks) /)) / num_forks
          n0(1:num_forks) = n(0:num_forks-1)
          n1(1:num_forks) = n(1:num_forks)
          do i = 1, num_forks
             call copy_array_pointer (ds(i)%x, d%x(n0(i):n1(i)), lb = 0)
             call copy_array_pointer (ds(i)%integral, d%integral(n0(i)+1:n1(i)))
             call copy_array_pointer (ds(i)%variance, d%variance(n0(i)+1:n1(i)))
          !  call copy_array_pointer (ds(i)%efficiency, d%efficiency(n0(i)+1:n1(i)))
             ds(i)%x = (ds(i)%x - ds(i)%x(0)) / (d%x(n1(i)) - d%x(n0(i)))
          end do
          ds%x_min = d%x_min + d%dx * d%x(n0)
          ds%x_max = d%x_min + d%dx * d%x(n1)
          ds%dx = ds%x_max - ds%x_min
          ds%x_min_true = d%x_min_true
          ds%x_max_true = d%x_max_true
          ds%stratified = d%stratified
          ds%ng = (d%ng * (n1 - n0)) / num_div
          num_calls = sum_calls !: this is a misnomer, it remains ``calls per cell'' here
          ds%dxg = real (n1 - n0, kind=default) / ds%ng
       else
          nx = lcm (d%ng / gcd (num_forks, d%ng), num_div)
          ds_ng = (d%ng * (/ (j, j=0,num_forks) /)) / num_forks
          n = (nx * ds_ng) / d%ng
          n0(1:num_forks) = n(0:num_forks-1)
          n1(1:num_forks) = n(1:num_forks)
          allocate (d_x(0:nx), d_integral(nx), d_variance(nx))
          ! allocate (d_efficiency(nx))
          call subdivide (d_x, d%x)
          call distribute (d_integral, d%integral)
          call distribute (d_variance, d%variance)
          ! call distribute (d_efficiency, d%efficiency)
          do i = 1, num_forks
             call copy_array_pointer (ds(i)%x, d_x(n0(i):n1(i)), lb = 0)
             call copy_array_pointer (ds(i)%integral, d_integral(n0(i)+1:n1(i)))
             call copy_array_pointer (ds(i)%variance, d_variance(n0(i)+1:n1(i)))
          !  call copy_array_pointer (ds(i)%efficiency, d_efficiency(n0(i)+1:n1(i)))
             ds(i)%x = (ds(i)%x - ds(i)%x(0)) / (d_x(n1(i)) - d_x(n0(i)))
          end do
          ds%x_min = d%x_min + d%dx * d_x(n0)
          ds%x_max = d%x_min + d%dx * d_x(n1)
          ds%dx = ds%x_max - ds%x_min
          ds%x_min_true = d%x_min_true
          ds%x_max_true = d%x_max_true
          ds%stratified = d%stratified
          ds%ng = ds_ng(1:num_forks) - ds_ng(0:num_forks-1)
          num_calls = sum_calls !: this is a misnomer, it remains ``calls per cell'' here
          ds%dxg = real (n1 - n0, kind=default) / ds%ng
          deallocate (d_x, d_integral, d_variance)
          ! deallocate (d_efficiency)
       end if
    else
       if (present (exc)) then
          call raise_exception (exc, EXC_FATAL, FN, "internal error")
       end if
       num_calls = 0
    end if
  end subroutine fork_division
  pure subroutine join_division (d, ds, exc)
    type(division_t), intent(inout) :: d
    type(division_t), dimension(:), intent(in) :: ds
    type(exception), intent(inout), optional :: exc
    character(len=*), parameter :: FN = "join_division"
    integer, dimension(size(ds)) :: n0, n1
    integer, dimension(0:size(ds)) :: n, ds_ng
    integer :: i, j, num_div, num_forks, nx
    real(kind=default), dimension(:), allocatable :: d_x, d_integral, d_variance
  ! real(kind=default), dimension(:), allocatable :: d_efficiency
    num_div = ubound (d%x, dim=1)
    num_forks = size (ds)
    if (d%ng == 1) then
       call sum_division (d, ds)
    else if (num_div >= num_forks) then
       if (modulo (d%ng, num_div) == 0) then
          n = (num_div * (/ (j, j=0,num_forks) /)) / num_forks
          n0(1:num_forks) = n(0:num_forks-1)
          n1(1:num_forks) = n(1:num_forks)
          do i = 1, num_forks
             d%integral(n0(i)+1:n1(i)) = ds(i)%integral
             d%variance(n0(i)+1:n1(i)) = ds(i)%variance
          !  d%efficiency(n0(i)+1:n1(i)) = ds(i)%efficiency
          end do
       else
          nx = lcm (d%ng / gcd (num_forks, d%ng), num_div)
          ds_ng = (d%ng * (/ (j, j=0,num_forks) /)) / num_forks
          n = (nx * ds_ng) / d%ng
          n0(1:num_forks) = n(0:num_forks-1)
          n1(1:num_forks) = n(1:num_forks)
          allocate (d_x(0:nx), d_integral(nx), d_variance(nx))
          ! allocate (d_efficiency(nx))
          do i = 1, num_forks
             d_integral(n0(i)+1:n1(i)) = ds(i)%integral
             d_variance(n0(i)+1:n1(i)) = ds(i)%variance
          !  d_efficiency(n0(i)+1:n1(i)) = ds(i)%efficiency
          end do
          call collect (d%integral, d_integral)
          call collect (d%variance, d_variance)
          ! call collect (d%efficiency, d_efficiency)
          deallocate (d_x, d_integral, d_variance)
          ! deallocate (d_efficiency)
       end if
    else
       if (present (exc)) then
          call raise_exception (exc, EXC_FATAL, FN, "internal error")
       end if
    end if
  end subroutine join_division
  pure subroutine subdivide (x, x0)
    real(kind=default), dimension(0:), intent(inout) :: x
    real(kind=default), dimension(0:), intent(in) :: x0
    integer :: i, n, n0
    n0 = ubound (x0, dim=1)
    n = ubound (x, dim=1) / n0
    x(0) = x0(0)
    do i = 1, n
       x(i::n) = x0(0:n0-1) * real (n - i) / n + x0(1:n0) * real (i) / n
    end do
  end subroutine subdivide
  pure subroutine distribute (x, x0)
    real(kind=default), dimension(:), intent(inout) :: x
    real(kind=default), dimension(:), intent(in) :: x0
    integer :: i, n
    n = ubound (x, dim=1) / ubound (x0, dim=1)
    do i = 1, n
       x(i::n) = x0 / n
    end do
  end subroutine distribute
  pure subroutine collect (x0, x)
    real(kind=default), dimension(:), intent(inout) :: x0
    real(kind=default), dimension(:), intent(in) :: x
    integer :: i, n, n0
    n0 = ubound (x0, dim=1)
    n = ubound (x, dim=1) / n0
    do i = 1, n0
       x0(i) = sum (x((i-1)*n+1:i*n))
    end do
  end subroutine collect
  pure subroutine sum_division (d, ds)
    type(division_t), intent(inout) :: d
    type(division_t), dimension(:), intent(in) :: ds
    integer :: i
    d%integral = 0.0
    d%variance = 0.0
  ! d%efficiency = 0.0
    do i = 1, size (ds)
       d%integral = d%integral + ds(i)%integral
       d%variance = d%variance + ds(i)%variance
  !    d%efficiency = d%efficiency + ds(i)%efficiency
    end do
  end subroutine sum_division
  subroutine debug_division (d, prefix)
    type(division_t), intent(in) :: d
    character(len=*), intent(in) :: prefix
    print "(1x,a,2(a,1x,i3,1x,f10.7))", prefix, ": d%x: ", &
         lbound(d%x,dim=1), d%x(lbound(d%x,dim=1)), &
         " ... ", &
         ubound(d%x,dim=1), d%x(ubound(d%x,dim=1))
    print "(1x,a,2(a,1x,i3,1x,f10.7))", prefix, ": d%i: ", &
         lbound(d%integral,dim=1), d%integral(lbound(d%integral,dim=1)), &
         " ... ", &
         ubound(d%integral,dim=1), d%integral(ubound(d%integral,dim=1))
    print "(1x,a,2(a,1x,i3,1x,f10.7))", prefix, ": d%v: ", &
         lbound(d%variance,dim=1), d%variance(lbound(d%variance,dim=1)), &
         " ... ", &
         ubound(d%variance,dim=1), d%variance(ubound(d%variance,dim=1))
  ! print "(1x,a,2(a,1x,i3,1x,f10.7))", prefix, ": d%e: ", &
  !      lbound(d%efficiency,dim=1), d%efficiency(lbound(d%efficiency,dim=1)), &
  !      " ... ", &
  !      ubound(d%efficiency,dim=1), d%efficiency(ubound(d%efficiency,dim=1))
  end subroutine debug_division
  subroutine dump_division (d, prefix)
    type(division_t), intent(in) :: d
    character(len=*), intent(in) :: prefix
  ! print "(2(1x,a),100(1x,f10.7))", prefix, ":x: ", d%x
    print "(2(1x,a),100(1x,f10.7))", prefix, ":x: ", d%x(1:)
    print "(2(1x,a),100(1x,e10.3))", prefix, ":i: ", d%integral
    print "(2(1x,a),100(1x,e10.3))", prefix, ":v: ", d%variance
  ! print "(2(1x,a),100(1x,e10.3))", prefix, ":e: ", d%efficiency
  end subroutine dump_division
  elemental function inside_division (d, x) result (theta)
    type(division_t), intent(in) :: d
    real(kind=default), intent(in) :: x
    logical :: theta
    theta = (x >= d%x_min_true) .and. (x <= d%x_max_true)
  end function inside_division
  elemental function stratified_division (d) result (yorn)
    type(division_t), intent(in) :: d
    logical :: yorn
    yorn = d%stratified
  end function stratified_division
  elemental function volume_division (d) result (vol)
    type(division_t), intent(in) :: d
    real(kind=default) :: vol
    vol = d%dx
  end function volume_division
  elemental function rigid_division (d) result (n)
    type(division_t), intent(in) :: d
    integer :: n
    n = d%ng
  end function rigid_division
  elemental function adaptive_division (d) result (n)
    type(division_t), intent(in) :: d
    integer :: n
    n = ubound (d%x, dim=1)
  end function adaptive_division
  elemental function summarize_division (d) result (s)
    type(division_t), intent(in) :: d
    type(div_history) :: s
    real(kind=default), dimension(:), allocatable :: p, m
    allocate (p(ubound(d%x,dim=1)), m(ubound(d%x,dim=1)))
    p = probabilities (d%x)
    m = rebinning_weights (d%variance)
    s%ng = d%ng
    s%num_div = ubound (d%x, dim=1)
    s%stratified = d%stratified
    s%x_min = d%x_min
    s%x_max = d%x_max
    s%x_min_true = d%x_min_true
    s%x_max_true = d%x_max_true
    s%spread_f_p = value_spread_percent (d%integral)
    s%stddev_f_p = standard_deviation_percent (d%integral)
    s%spread_p = value_spread_percent (p)
    s%stddev_p = standard_deviation_percent (p)
    s%spread_m = value_spread_percent (m)
    s%stddev_m = standard_deviation_percent (m)
    deallocate (p, m)
  end function summarize_division
  pure function probabilities (x) result (p)
    real(kind=default), dimension(0:), intent(in) :: x
    real(kind=default), dimension(ubound(x,dim=1)) :: p
    integer :: num_div
    num_div = ubound (x, dim=1)
    p = 1.0 / (x(1:num_div) - x(0:num_div-1))
    p = p / sum(p)
  end function probabilities
  subroutine print_history (h, tag)
    type(div_history), dimension(:), intent(in) :: h
    character(len=*), intent(in), optional :: tag
    call write_history (output_unit, h, tag)
    flush (output_unit)
  end subroutine print_history
  subroutine write_history (u, h, tag)
    integer, intent(in) :: u
    type(div_history), dimension(:), intent(in) :: h
    character(len=*), intent(in), optional :: tag
    character(len=BUFFER_SIZE) :: pfx
    character(len=1) :: s
    integer :: i
    if (present (tag)) then
       pfx = tag
    else
       pfx = "[vamp]"
    end if
    if ((minval (h%x_min) == maxval (h%x_min)) &
         .and. (minval (h%x_max) == maxval (h%x_max))) then
       write (u, "(1X,A11,1X,2X,1X,2(ES10.3,A4,ES10.3,A7))") pfx, &
            h(1)%x_min, " <= ", h(1)%x_min_true, &
            " < x < ", h(1)%x_max_true, " <= ", h(1)%x_max
    else
       do i = 1, size (h)
          write (u, "(1X,A11,1X,I2,1X,2(ES10.3,A4,ES10.3,A7))") pfx, &
               i, h(i)%x_min, " <= ", h(i)%x_min_true, &
               " < x < ", h(i)%x_max_true, " <= ", h(i)%x_max
       end do
    end if
    write (u, "(1X,A11,1X,A2,2(1X,A3),A1,6(1X,A8))") pfx, &
         "it", "nd", "ng", "", &
         "spr(f/p)", "dev(f/p)", "spr(m)", "dev(m)", "spr(p)", "dev(p)"
    iterations: do i = 1, size (h)
       if (h(i)%stratified) then
          s = "*"
       else
          s = ""
       end if
       write (u, "(1X,A11,1X,I2,2(1X,I3),A1,6(1X,F7.2,A1))") pfx, &
            i, h(i)%num_div, h(i)%ng, s, &
            h(i)%spread_f_p, "%", h(i)%stddev_f_p, "%", &
            h(i)%spread_m, "%", h(i)%stddev_m, "%", &
            h(i)%spread_p, "%", h(i)%stddev_p, "%"
    end do iterations
    flush (u)
  end subroutine write_history
  subroutine write_division_unit (d, unit, write_integrals)
    type(division_t), intent(in) :: d
    integer, intent(in) :: unit
    logical, intent(in), optional :: write_integrals
    logical :: write_integrals0
    integer :: i
    write_integrals0 = .false.
    if (present(write_integrals)) write_integrals0 = write_integrals
    write (unit = unit, fmt = descr_fmt) "begin type(division_t) :: d"
    write (unit = unit, fmt = integer_fmt) "ubound(d%x,1) = ", ubound (d%x, dim=1)
    write (unit = unit, fmt = integer_fmt) "d%ng = ", d%ng
    write (unit = unit, fmt = logical_fmt) "d%stratified = ", d%stratified
    write (unit = unit, fmt = double_fmt) "d%dx = ", d%dx
    write (unit = unit, fmt = double_fmt) "d%dxg = ", d%dxg
    write (unit = unit, fmt = double_fmt) "d%x_min = ", d%x_min
    write (unit = unit, fmt = double_fmt) "d%x_max = ", d%x_max
    write (unit = unit, fmt = double_fmt) "d%x_min_true = ", d%x_min_true
    write (unit = unit, fmt = double_fmt) "d%x_max_true = ", d%x_max_true
    write (unit = unit, fmt = descr_fmt) "begin d%x" 
    do i = 0, ubound (d%x, dim=1)
       if (write_integrals0 .and. i/=0) then
          write (unit = unit, fmt = double_array_fmt) &
               i, d%x(i), d%integral(i), d%variance(i)
       else
          write (unit = unit, fmt = double_array_fmt) i, d%x(i)
       end if
    end do
    write (unit = unit, fmt = descr_fmt) "end d%x"
    write (unit = unit, fmt = descr_fmt) "end type(division_t)"
  end subroutine write_division_unit
  subroutine read_division_unit (d, unit, read_integrals)
    type(division_t), intent(inout) :: d
    integer, intent(in) :: unit
    logical, intent(in), optional :: read_integrals
    logical :: read_integrals0
    integer :: i, idum, num_div
    character(len=80) :: chdum
    read_integrals0 = .false.
    if (present(read_integrals)) read_integrals0 = read_integrals
    read (unit = unit, fmt = descr_fmt) chdum
    read (unit = unit, fmt = integer_fmt) chdum, num_div
    if (associated (d%x)) then
       if (ubound (d%x, dim=1) /= num_div) then
          deallocate (d%x, d%integral, d%variance)
    !     deallocate (d%efficiency)
          allocate (d%x(0:num_div), d%integral(num_div), d%variance(num_div))
    !     allocate (d%efficiency(num_div))
       end if
    else
       allocate (d%x(0:num_div), d%integral(num_div), d%variance(num_div))
    !  allocate (d%efficiency(num_div))
    end if
    read (unit = unit, fmt = integer_fmt) chdum, d%ng
    read (unit = unit, fmt = logical_fmt) chdum, d%stratified
    read (unit = unit, fmt = double_fmt) chdum, d%dx
    read (unit = unit, fmt = double_fmt) chdum, d%dxg
    read (unit = unit, fmt = double_fmt) chdum, d%x_min
    read (unit = unit, fmt = double_fmt) chdum, d%x_max
    read (unit = unit, fmt = double_fmt) chdum, d%x_min_true
    read (unit = unit, fmt = double_fmt) chdum, d%x_max_true
    read (unit = unit, fmt = descr_fmt) chdum
    do i = 0, ubound (d%x, dim=1)
       if (read_integrals0 .and. i/=0) then
          read (unit = unit, fmt = double_array_fmt) &
               & idum, d%x(i), d%integral(i), d%variance(i)
       else
          read (unit = unit, fmt = double_array_fmt) idum, d%x(i)
       end if
    end do
    read (unit = unit, fmt = descr_fmt) chdum
    read (unit = unit, fmt = descr_fmt) chdum
    if (.not.read_integrals0) then
       d%integral = 0.0
       d%variance = 0.0
  !    d%efficiency = 0.0
    end if
  end subroutine read_division_unit
  subroutine write_division_name (d, name, write_integrals)
    type(division_t), intent(in) :: d
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: write_integrals
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "write", status = "replace", file = name)
    call write_division_unit (d, unit, write_integrals)
    close (unit = unit)
  end subroutine write_division_name
  subroutine read_division_name (d, name, read_integrals)
    type(division_t), intent(inout) :: d
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: read_integrals
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "read", status = "old", file = name)
    call read_division_unit (d, unit, read_integrals)
    close (unit = unit)
  end subroutine read_division_name
  subroutine write_division_raw_unit (d, unit, write_integrals)
    type(division_t), intent(in) :: d
    integer, intent(in) :: unit
    logical, intent(in), optional :: write_integrals
    logical :: write_integrals0
    integer :: i
    write_integrals0 = .false.
    if (present(write_integrals)) write_integrals0 = write_integrals
    write (unit = unit) MAGIC_DIVISION_BEGIN
    write (unit = unit) ubound (d%x, dim=1)
    write (unit = unit) d%ng
    write (unit = unit) d%stratified
    write (unit = unit) d%dx
    write (unit = unit) d%dxg
    write (unit = unit) d%x_min
    write (unit = unit) d%x_max
    write (unit = unit) d%x_min_true
    write (unit = unit) d%x_max_true
    do i = 0, ubound (d%x, dim=1)
       if (write_integrals0 .and. i/=0) then
          write (unit = unit) d%x(i), d%integral(i), d%variance(i)
       else
          write (unit = unit) d%x(i)
       end if
    end do
    write (unit = unit) MAGIC_DIVISION_END
  end subroutine write_division_raw_unit
  subroutine read_division_raw_unit (d, unit, read_integrals)
    type(division_t), intent(inout) :: d
    integer, intent(in) :: unit
    logical, intent(in), optional :: read_integrals
    logical :: read_integrals0
    integer :: i, num_div, magic
    character(len=*), parameter :: FN = "read_division_raw_unit"
    read_integrals0 = .false.
    if (present(read_integrals)) read_integrals0 = read_integrals
    read (unit = unit) magic
    if (magic /= MAGIC_DIVISION_BEGIN) then
       print *, FN, " fatal: expecting magic ", MAGIC_DIVISION_BEGIN, &
                    ", found ", magic
       stop
    end if
    read (unit = unit) num_div
    if (associated (d%x)) then
       if (ubound (d%x, dim=1) /= num_div) then
          deallocate (d%x, d%integral, d%variance)
    !     deallocate (d%efficiency)
          allocate (d%x(0:num_div), d%integral(num_div), d%variance(num_div))
    !     allocate (d%efficiency(num_div))
       end if
    else
       allocate (d%x(0:num_div), d%integral(num_div), d%variance(num_div))
    !  allocate (d%efficiency(num_div))
    end if
    read (unit = unit) d%ng
    read (unit = unit) d%stratified
    read (unit = unit) d%dx
    read (unit = unit) d%dxg
    read (unit = unit) d%x_min
    read (unit = unit) d%x_max
    read (unit = unit) d%x_min_true
    read (unit = unit) d%x_max_true
    do i = 0, ubound (d%x, dim=1)
       if (read_integrals0 .and. i/=0) then
          read (unit = unit) d%x(i), d%integral(i), d%variance(i)
       else
          read (unit = unit) d%x(i)
       end if
    end do
    if (.not.read_integrals0) then
       d%integral = 0.0
       d%variance = 0.0
  !    d%efficiency = 0.0
    end if
    read (unit = unit) magic
    if (magic /= MAGIC_DIVISION_END) then
       print *, FN, " fatal: expecting magic ", MAGIC_DIVISION_END, &
                    ", found ", magic
       stop
    end if
  end subroutine read_division_raw_unit
  subroutine write_division_raw_name (d, name, write_integrals)
    type(division_t), intent(in) :: d
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: write_integrals
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "write", status = "replace", &
          form = "unformatted", file = name)
    call write_division_unit (d, unit, write_integrals)
    close (unit = unit)
  end subroutine write_division_raw_name
  subroutine read_division_raw_name (d, name, read_integrals)
    type(division_t), intent(inout) :: d
    character(len=*), intent(in) :: name
    logical, intent(in), optional :: read_integrals
    integer :: unit
    call find_free_unit (unit)
    open (unit = unit, action = "read", status = "old", &
          form = "unformatted", file = name)
    call read_division_unit (d, unit, read_integrals)
    close (unit = unit)
  end subroutine read_division_raw_name
  pure subroutine marshal_division (d, ibuf, dbuf)
    type(division_t), intent(in) :: d
    integer, dimension(:), intent(inout) :: ibuf
    real(kind=default), dimension(:), intent(inout) :: dbuf
    integer :: num_div
    num_div = ubound (d%x, dim=1)
    ibuf(1) = d%ng
    ibuf(2) = num_div
    if (d%stratified) then
       ibuf(3) = 1
    else
       ibuf(3) = 0
    end if
    dbuf(1) = d%x_min
    dbuf(2) = d%x_max
    dbuf(3) = d%x_min_true
    dbuf(4) = d%x_max_true
    dbuf(5) = d%dx
    dbuf(6) = d%dxg
    dbuf(7:7+num_div) = d%x
    dbuf(8+  num_div:7+2*num_div) = d%integral
    dbuf(8+2*num_div:7+3*num_div) = d%variance
  ! dbuf(8+3*num_div:7+4*num_div) = d%efficiency
  end subroutine marshal_division
  pure subroutine marshal_division_size (d, iwords, dwords)
    type(division_t), intent(in) :: d
    integer, intent(out) :: iwords, dwords
    iwords = 3
    dwords = 7 + 3 * ubound (d%x, dim=1)
  ! dwords = 7 + 4 * ubound (d%x, dim=1)
  end subroutine marshal_division_size
  pure subroutine unmarshal_division (d, ibuf, dbuf)
    type(division_t), intent(inout) :: d
    integer, dimension(:), intent(in) :: ibuf
    real(kind=default), dimension(:), intent(in) :: dbuf
    integer :: num_div
    d%ng = ibuf(1)
    num_div = ibuf(2)
    d%stratified = ibuf(3) /= 0
    d%x_min = dbuf(1)
    d%x_max = dbuf(2)
    d%x_min_true = dbuf(3)
    d%x_max_true = dbuf(4)
    d%dx = dbuf(5)
    d%dxg = dbuf(6)
    if (associated (d%x)) then
       if (ubound (d%x, dim=1) /= num_div) then
          deallocate (d%x, d%integral, d%variance)
    !     deallocate (d%efficiency)
          allocate (d%x(0:num_div), d%integral(num_div), d%variance(num_div))
    !     allocate (d%efficiency(num_div))
       end if
    else
       allocate (d%x(0:num_div), d%integral(num_div), d%variance(num_div))
    !  allocate (d%efficiency(num_div))
    end if
    d%x = dbuf(7:7+num_div)
    d%integral = dbuf(8+  num_div:7+2*num_div)
    d%variance = dbuf(8+2*num_div:7+3*num_div)
  ! d%efficiency = dbuf(8+3*num_div:7+4*num_div)
  end subroutine unmarshal_division
  pure subroutine marshal_div_history (h, ibuf, dbuf)
    type(div_history), intent(in) :: h
    integer, dimension(:), intent(inout) :: ibuf
    real(kind=default), dimension(:), intent(inout) :: dbuf
    ibuf(1) = h%ng
    ibuf(2) = h%num_div
    if (h%stratified) then
       ibuf(3) = 1
    else
       ibuf(3) = 0
    end if
    dbuf(1) = h%x_min
    dbuf(2) = h%x_max
    dbuf(3) = h%x_min_true
    dbuf(4) = h%x_max_true
    dbuf(5) = h%spread_f_p
    dbuf(6) = h%stddev_f_p
    dbuf(7) = h%spread_p
    dbuf(8) = h%stddev_p
    dbuf(9) = h%spread_m
    dbuf(10) = h%stddev_m
  end subroutine marshal_div_history
  pure subroutine marshal_div_history_size (h, iwords, dwords)
    type(div_history), intent(in) :: h
    integer, intent(out) :: iwords, dwords
    iwords = 3
    dwords = 10
  end subroutine marshal_div_history_size
  pure subroutine unmarshal_div_history (h, ibuf, dbuf)
    type(div_history), intent(inout) :: h
    integer, dimension(:), intent(in) :: ibuf
    real(kind=default), dimension(:), intent(in) :: dbuf
    h%ng = ibuf(1)
    h%num_div = ibuf(2)
    h%stratified = ibuf(3) /= 0
    h%x_min = dbuf(1)
    h%x_max = dbuf(2)
    h%x_min_true = dbuf(3)
    h%x_max_true = dbuf(4)
    h%spread_f_p = dbuf(5)
    h%stddev_f_p = dbuf(6)
    h%spread_p = dbuf(7)
    h%stddev_p = dbuf(8)
    h%spread_m = dbuf(9)
    h%stddev_m = dbuf(10)
  end subroutine unmarshal_div_history
  elemental subroutine copy_division (lhs, rhs)
    type(division_t), intent(inout) :: lhs
    type(division_t), intent(in) :: rhs
    if (associated (rhs%x)) then
       call copy_array_pointer (lhs%x, rhs%x, lb = 0)
    else if (associated (lhs%x)) then
       deallocate (lhs%x)
    end if
    if (associated (rhs%integral)) then
       call copy_array_pointer (lhs%integral, rhs%integral)
    else if (associated (lhs%integral)) then
       deallocate (lhs%integral)
    end if
    if (associated (rhs%variance)) then
       call copy_array_pointer (lhs%variance, rhs%variance)
    else if (associated (lhs%variance)) then
       deallocate (lhs%variance)
    end if
  ! if (associated (rhs%efficiency)) then
  !    call copy_array_pointer (lhs%efficiency, rhs%efficiency)
  ! else if (associated (lhs%efficiency)) then
  !    deallocate (lhs%efficiency)
  ! end if
    lhs%dx = rhs%dx
    lhs%dxg = rhs%dxg
    lhs%x_min = rhs%x_min
    lhs%x_max = rhs%x_max
    lhs%x_min_true = rhs%x_min_true
    lhs%x_max_true = rhs%x_max_true
    lhs%ng = rhs%ng
    lhs%stratified = rhs%stratified
  end subroutine copy_division
  elemental subroutine delete_division (d)
    type(division_t), intent(inout) :: d
    if (associated (d%x)) then
       deallocate (d%x, d%integral, d%variance)
  !    deallocate (d%efficiency)
    end if
  end subroutine delete_division
  elemental subroutine copy_history (lhs, rhs)
    type(div_history), intent(out) :: lhs
    type(div_history), intent(in) :: rhs
    lhs%stratified = rhs%stratified
    lhs%ng = rhs%ng
    lhs%num_div = rhs%num_div
    lhs%x_min = rhs%x_min
    lhs%x_max = rhs%x_max
    lhs%x_min_true = rhs%x_min_true
    lhs%x_max_true = rhs%x_max_true
    lhs%spread_f_p = rhs%spread_f_p
    lhs%stddev_f_p = rhs%stddev_f_p
    lhs%spread_p = rhs%spread_p
    lhs%stddev_p = rhs%stddev_p
    lhs%spread_m = rhs%spread_m
    lhs%stddev_m = rhs%stddev_m
  end subroutine copy_history
end module divisions
