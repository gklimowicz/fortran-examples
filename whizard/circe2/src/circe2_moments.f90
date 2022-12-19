module sampling
  use kinds
  implicit none
  private
  type sample
    integer :: n = 0
    real(kind=default) :: w = 0
    real(kind=default) :: w2 = 0
  end type sample
  public :: sample
  public :: reset, record
  public :: mean, variance
contains
  elemental subroutine reset (s)
    type(sample), intent(inout) :: s
    s%n = 0
    s%w = 0
    s%w2 = 0
  end subroutine reset
  elemental subroutine record (s, w)
    type(sample), intent(inout) :: s
    real(kind=default), intent(in), optional :: w
    s%n = s%n + 1
    if (present (w)) then
      s%w = s%w + w
      s%w2 = s%w2 + w*w
    else
      s%w = s%w + 1
      s%w2 = s%w2 + 1
    endif
  end subroutine record
  elemental function mean (s)
    type(sample), intent(in) :: s
    real(kind=default) :: mean
    mean = s%w / s%n
  end function mean
  elemental function variance (s)
    type(sample), intent(in) :: s
    real(kind=default) :: variance
    variance = (s%w2 / s%n - mean(s)**2) / s%n
    variance = max (variance, epsilon (variance))
  end function variance
end module sampling
module circe2_moments_library
  use kinds
  use tao_random_objects !NODEP!
  use sampling !NODEP!
  use circe2
  implicit none
  private
    public :: generate_beta
  type moment
    integer, dimension(2) :: n, m
    type(sample) :: sample = sample (0, 0.0_default, 0.0_default)
  end type moment
  public :: moment
  public :: init_moments
  public :: reset_moment, record_moment
  public :: mean_moment, variance_moment
  public :: beta_moment
  type channel
    real(kind=default) :: w = 1
    real(kind=default), dimension(2) :: a = 1, b = 1
    logical, dimension(2) :: delta = .false.
  end type channel
  public :: channel
  public :: generate_beta_multi, beta_moments_multi
  public :: selftest
  public :: random2_seed
  public :: read_channels
  public :: report_results
  public :: results_ok
  public :: generate
  public :: compare
  public :: check
contains
    subroutine generate_beta (rng, x, xmin, xmax, a, b)
      class(rng_type), intent(inout) :: rng
      real(kind=default), intent(out) :: x
      real(kind=default), intent(in) :: xmin, xmax, a, b
      real(kind=default) :: t, p, u, umin, umax, w
            if (a >= 1 .or. b <= 1) then
               x = -1
               print *, 'ERROR: beta-distribution expects a<1<b'
               return
            end if
                  t = (1 - a) / (b + 1 - a)
            p = b*t / (b*t + a * (1 - t)**b)
            if (xmin <= 0) then
               umin = 0
            elseif (xmin < t) then
               umin = p * (xmin/t)**a
            elseif (xmin == t) then
               umin = p
            elseif (xmin < 1) then
               umin = 1 - (1 - p) * ((1 - xmin)/(1 - t))**b
            else
               umin = 1
            endif
            if (xmax >= 1) then
               umax = 1
            elseif (xmax > t) then
               umax = 1 - (1 - p) * ((1 - xmax)/(1 - t))**b
            elseif (xmax == t) then
               umax = p
            elseif (xmax > 0) then
               umax = p * (xmax/t)**a
            else
               umax = 0
            endif
            if (umax < umin) then
               x = -1
               return
            endif
      do 
               call rng%generate (u)
               u = umin + (umax - umin) * u
               if (u <= p) then
                  x = t * (u/p)**(1/a)
                  w = (1 - x)**(b-1)
               else
                  x = 1 - (1 - t) * ((1 - u)/(1 - p))**(1/b)
                  w = (x/t)**(a-1)
               end if
         call rng%generate (u)
         if (w > u) exit
      end do
    end subroutine generate_beta
  subroutine init_moments (moments)
    type(moment), dimension(0:,0:,0:,0:), intent(inout) :: moments
    integer :: nx, mx, ny, my
    do nx = lbound(moments,1), ubound(moments,1)
       do mx = lbound(moments,2), ubound(moments,2)
          do ny = lbound(moments,3), ubound(moments,3)
             do my = lbound(moments,4), ubound(moments,4)
                moments(nx,mx,ny,my) = moment([nx,ny],[mx,my])
             end do
          end do
       end do
    end do
    call reset_moment (moments)
  end subroutine init_moments
  elemental subroutine reset_moment (m)
    type(moment), intent(inout) :: m
    call reset (m%sample)
  end subroutine reset_moment
  elemental subroutine record_moment (m, x1, x2, w)
    type(moment), intent(inout) :: m
    real(kind=default), intent(in) :: x1, x2
    real(kind=default), intent(in), optional :: w
    real(kind=default) :: p
    p = pwr (x1, m%n(1)) * pwr (1-x1, m%m(1)) &
      * pwr (x2, m%n(2)) * pwr (1-x2, m%m(2))
    if (present (w)) p = p*w
    call record (m%sample, p)
  contains
    pure function pwr (x, n)
      real(kind=default), intent(in) :: x
      integer, intent(in) :: n
      real(kind=default) :: pwr
      if (n == 0) then
         pwr = 1
      else 
         pwr = x**n
      end if
    end function pwr
  end subroutine record_moment
  elemental function mean_moment (m)
    type(moment), intent(in) :: m
    real(kind=default) :: mean_moment
    mean_moment = mean (m%sample)
  end function mean_moment
  elemental function variance_moment (m)
    type(moment), intent(in) :: m
    real(kind=default) :: variance_moment
    variance_moment = variance (m%sample)
  end function variance_moment
  elemental function beta_moment (n, m, a, b)
    integer, intent(in) :: n, m
    real(kind=default), intent(in) :: a, b
    real(kind=default) :: beta_moment
    beta_moment = &
      gamma_ratio (a, n) * gamma_ratio (b, m) / gamma_ratio (a+b, n+m)
  end function beta_moment
  elemental function gamma_ratio (x, n)
    real(kind=default), intent(in) :: x
    integer, intent(in) :: n
    real(kind=default) :: gamma_ratio
    integer :: i
    gamma_ratio = 1
    do i = 0, n - 1
      gamma_ratio = gamma_ratio * (x + i)
    end do
  end function gamma_ratio
  subroutine generate_beta_multi (rng, x, channels)
    class(rng_type), intent(inout) :: rng
    real(kind=default), dimension(:), intent(out) :: x
    type(channel), dimension(:), intent(in) :: channels
    real(kind=default) :: u, accum
    integer :: i, n
    call rng%generate (u)
    u = u * sum (channels%w)
    accum = 0
    scan: do n = 1, size (channels) - 1
      accum = accum + channels(n)%w
      if (accum >= u) exit scan
    end do scan
    do i = 1, size (x)
      if (channels(n)%delta(i)) then
        x(i) = 1
      else
        if (channels(n)%a(i) == 1 .and. channels(n)%b(i) == 1) then
          call rng%generate (x(i))
        else if (channels(n)%b(i) < channels(n)%a(i)) then
          call generate_beta (rng, x(i), 0.0_default, 1.0_default, &
                              channels(n)%b(i), channels(n)%a(i))
          x(i) = 1 - x(i)
        else
          call generate_beta (rng, x(i), 0.0_default, 1.0_default, &
                              channels(n)%a(i), channels(n)%b(i))
        end if
      end if
    end do
  end subroutine generate_beta_multi
  pure function beta_moments_multi (n, m, channels)
    integer, intent(in), dimension(2) :: n, m
    type(channel), dimension(:), intent(in) :: channels
    real(kind=default) :: beta_moments_multi
    real(kind=default) :: w
    integer :: c, i
    beta_moments_multi = 0
    do c = 1, size (channels)
      w = channels(c)%w
      do i = 1, 2
        if (channels(c)%delta(i)) then
          if (m(i) > 0) w = 0
        else
          w = w * beta_moment (n(i), m(i), channels(c)%a(i), channels(c)%b(i))
        end if
      end do
      beta_moments_multi = beta_moments_multi + w
    end do
    beta_moments_multi = beta_moments_multi / sum (channels%w)
  end function beta_moments_multi
  subroutine selftest (rng, nevents)
    class(rng_type), intent(inout) :: rng
    integer, intent(in) :: nevents
    integer, parameter :: N = 1
    type(moment), dimension(0:N,0:N,0:N,0:N) :: moments
    integer :: i
    real(kind=default), dimension(2) :: x
    type(channel), dimension(:), allocatable :: channels
    call read_channels (channels)
    call init_moments (moments)
    do i = 1, nevents
       call generate_beta_multi (rng, x, channels)
       call record_moment (moments, x(1), x(2))
    end do
    call report_results (moments, channels)
  end subroutine selftest
  subroutine random2_seed (rng, seed)
    class(rng_tao), intent(inout) :: rng
    integer, intent(in), optional:: seed
    integer, dimension(8) :: date_time
    integer :: seed_value
    if (present (seed)) then
       seed_value = seed
    else
       call date_and_time (values = date_time)
       seed_value = product (date_time)
    endif
    call rng%init (seed_value)
  end subroutine random2_seed
  subroutine read_channels (channels)
    type(channel), dimension(:), allocatable, intent(out) :: channels
    type(channel), dimension(100) :: buffer
    real(kind=default) :: w
    real(kind=default), dimension(2) :: a, b
    logical, dimension(2) :: delta
    integer :: n, status
    do n = 1, size (buffer)
      read (*, *, iostat = status) w, a(1), b(1), a(2), b(2), delta
      if (status == 0) then
        buffer(n) = channel (w, a, b, delta)
      else
        exit
      end if
    end do
    allocate (channels(n-1))
    channels = buffer(1:n-1)
  end subroutine read_channels
  subroutine report_results (moments, channels)
    type(moment), dimension(0:,0:,0:,0:), intent(in) :: moments
    type(channel), dimension(:), intent(in) :: channels
    integer :: nx, mx, ny, my
    real(kind=default) :: truth, estimate, sigma, pull, eps
    do nx = lbound(moments,1), ubound(moments,1)
       do mx = lbound(moments,2), ubound(moments,2)
          do ny = lbound(moments,3), ubound(moments,3)
             do my = lbound(moments,4), ubound(moments,4)
                truth = beta_moments_multi ([nx, ny], [mx, my], channels)
                estimate = mean_moment (moments(nx,mx,ny,my))
                sigma = sqrt (variance_moment (moments(nx,mx,ny,my)))
                pull = estimate - truth
                eps = pull / max (epsilon (1.0_default), epsilon (1.0_double))
                if (sigma /= 0.0_default) pull = pull / sigma
                write (*, "(' x^', I1, ' (1-x)^', I1, &
                           &' y^', I1, ' (1-y)^', I1, &
                           &': ', F8.5, ': est = ', F8.5, &
                           &' +/- ', F8.5,&
                           &', pull = ', F8.2,&
                           &', eps = ', F8.2)") &
                     nx, mx, ny, my, truth, estimate, sigma, pull, eps
             end do
          end do
       end do
    end do
  end subroutine report_results
  function results_ok (moments, channels, threshold, fraction)
    ! use, intrinsic :: ieee_arithmetic
    type(moment), dimension(0:,0:,0:,0:), intent(in) :: moments
    type(channel), dimension(:), intent(in) :: channels
    real(kind=default), intent(in), optional :: threshold, fraction
    logical :: results_ok
    integer :: nx, mx, ny, my, failures
    real(kind=default) :: thr, frac, eps
    real(kind=default) :: truth, estimate, sigma
    ! we mut not expect to measure zero better than the
    ! double precision used in the ocaml code:
    eps = 200 * max (epsilon (1.0_default), epsilon (1.0_double))
    if (present(threshold)) then
       thr = threshold
    else
       thr = 5
    end if
    if (present(fraction)) then
       frac = fraction
    else
       frac = 0.01_default
    end if
    failures = 0
    do nx = lbound(moments,1), ubound(moments,1)
       do mx = lbound(moments,2), ubound(moments,2)
          do ny = lbound(moments,3), ubound(moments,3)
             do my = lbound(moments,4), ubound(moments,4)
                truth = beta_moments_multi ([nx, ny], [mx, my], channels)
                estimate = mean_moment (moments(nx,mx,ny,my))
                sigma = sqrt (variance_moment (moments(nx,mx,ny,my)))
                if (.not. (      ieee_is_normal (truth) &
                           .and. ieee_is_normal (estimate) &
                           .and. ieee_is_normal (sigma)) &
                    .or. abs (estimate - truth) > max (thr * sigma, eps)) then
                   failures = failures + 1
                end if
             end do
          end do
       end do
    end do
    if (failures >= frac * size (moments)) then
      results_ok = .false.
    else
      results_ok = .true.
    end if
  contains
    function ieee_is_normal (x)
      real(kind=default), intent(in) :: x
      logical :: ieee_is_normal
      ieee_is_normal = .not. (x /= x)
    end function ieee_is_normal
  end function results_ok
  subroutine generate (rng, nevents)
    class(rng_type), intent(inout) :: rng
    integer, intent(in) :: nevents
    type(channel), dimension(:), allocatable :: channels
    real(kind=default), dimension(2) :: x
    integer :: i
    call read_channels (channels)
    do i = 1, nevents
       call generate_beta_multi (rng, x, channels)
       write (*, "(3(5x,F19.17))") x, 1.0_default
    end do
  end subroutine generate
  subroutine compare (rng, nevents, file)
    class(rng_type), intent(inout) :: rng
    integer, intent(in) :: nevents
    character(len=*), intent(in) :: file
    type(channel), dimension(:), allocatable :: channels
    integer, parameter :: N = 1
    type(moment), dimension(0:N,0:N,0:N,0:N) :: moments
    real(kind=default), dimension(2) :: x
    character(len=128) :: design
    real(kind=default) :: roots
    integer :: ierror
    integer, dimension(2) :: p, h
    integer :: i
    type(circe2_state) :: c2s
    call read_channels (channels)
    call init_moments (moments)
    design = "CIRCE2/TEST"
    roots = 42
    p = [11, -11]
    h = 0
    call circe2_load (c2s, trim(file), trim(design), roots, ierror)
    do i = 1, nevents
       call circe2_generate (c2s, rng, x, p, h)
       call record_moment (moments, x(1), x(2))
    end do
    call report_results (moments, channels)
  end subroutine compare
  subroutine check (rng, nevents, file, distributions, fail)
    class(rng_type), intent(inout) :: rng
    integer, intent(in) :: nevents
    character(len=*), intent(in) :: file
    logical, intent(in), optional :: distributions, fail
    type(channel), dimension(:), allocatable :: channels
    type(channel), dimension(1) :: unit_channel
    integer, parameter :: N = 1
    type(moment), dimension(0:N,0:N,0:N,0:N) :: moments, unit_moments
    real(kind=default), dimension(2) :: x
    character(len=128) :: design
    real(kind=default) :: roots, weight
    integer :: ierror
    integer, dimension(2) :: p, h
    integer :: i
    logical :: generation_ok, distributions_ok
    logical :: check_distributions, expect_failure
    type(circe2_state) :: c2s
    if (present (distributions)) then
       check_distributions = distributions
    else
       check_distributions = .true.
    end if
    if (present (fail)) then
       expect_failure = fail
    else
       expect_failure = .false.
    end if
    call read_channels (channels)
    call init_moments (moments)
    if (check_distributions) call init_moments (unit_moments)
    design = "CIRCE2/TEST"
    roots = 42
    p = [11, -11]
    h =   0
    call circe2_load (c2s, trim(file), trim(design), roots, ierror)
    do i = 1, nevents
       call circe2_generate (c2s, rng, x, p, h)
       call record_moment (moments, x(1), x(2))
       if (check_distributions) then
          weight = circe2_distribution (c2s, p, h, x)
          call record_moment (unit_moments, x(1), x(2), w = 1 / weight)
       end if
    end do
    generation_ok = results_ok (moments, channels)
    if (check_distributions) then
       distributions_ok = results_ok (unit_moments, unit_channel)
    else
       distributions_ok = .not. expect_failure
    end if
    if (expect_failure) then
       if (generation_ok .and. distributions_ok) then
          print *, "FAIL: unexpected success"
       else
          if (.not. generation_ok) then
             print *, "OK: expected failure in generation"
          end if
          if (.not. distributions_ok) then
             print *, "OK: expected failure in distributions"
          end if
       end if
       call report_results (moments, channels)
    else
       if (generation_ok .and. distributions_ok) then
          print *, "OK"
       else
          if (.not. generation_ok) then
             print *, "FAIL: generation"
             call report_results (moments, channels)
          end if
          if (.not. distributions_ok) then
             print *, "FAIL: distributions"
             call report_results (unit_moments, unit_channel)
          end if
       end if
    end if
  end subroutine check
end module circe2_moments_library
program circe2_moments
  use circe2
  use circe2_moments_library !NODEP!
  use tao_random_objects !NODEP!
  implicit none
  type(rng_tao), save :: rng
  character(len=1024) :: mode, filename, buffer
  integer :: status, nevents, seed
  call get_command_argument (1, value = mode, status = status)
  if (status /= 0) mode = ""
  call get_command_argument (2, value = filename, status = status)
  if (status /= 0) filename = ""
  call get_command_argument (3, value = buffer, status = status)
  if (status == 0) then
     read (buffer, *, iostat = status) nevents
     if (status /= 0) nevents = 1000
  else
     nevents = 1000
  end if
  call get_command_argument (4, value = buffer, status = status)
  if (status == 0) then
     read (buffer, *, iostat = status) seed
     if (status == 0) then
        call random2_seed (rng, seed)
     else
        call random2_seed (rng)
     end if
  else
     call random2_seed (rng)
  end if
  select case (trim (mode))
  case ("check")
     call check (rng, nevents, trim (filename))
  case ("!check")
     call check (rng, nevents, trim (filename), fail = .true.)
  case ("check_generation")
     call check (rng, nevents, trim (filename), distributions = .false.)
  case ("!check_generation")
     call check (rng, nevents, trim (filename), fail = .true., &
                                                distributions = .false.)
  case ("compare")
     call compare (rng, nevents, trim (filename))
  case ("generate")
     call generate (rng, nevents)
  case ("selftest")
     call selftest (rng, nevents)
  case default
     print *, &
      "usage: circe2_moments " // &
      "[check|check_generation|generate|selftest] " // &
      "filename [events] [seed]"
  end select
end program circe2_moments
