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
module omega_testtools
  use kinds
  implicit none
  private
  real(kind=default), parameter, private :: ABS_THRESHOLD_DEFAULT = 1E-17
  real(kind=default), parameter, private :: THRESHOLD_DEFAULT = 0.6
  real(kind=default), parameter, private :: THRESHOLD_WARN = 0.8
  public :: agreement
  interface agreement
     module procedure agreement_real, agreement_complex, &
            agreement_real_complex, agreement_complex_real, &
            agreement_integer_complex, agreement_complex_integer, &
            agreement_integer_real, agreement_real_integer
  end interface
  private :: agreement_real, agreement_complex, &
       agreement_real_complex, agreement_complex_real, &
       agreement_integer_complex, agreement_complex_integer, &
       agreement_integer_real, agreement_real_integer
  public:: vanishes
  interface vanishes
     module procedure vanishes_real, vanishes_complex
  end interface
  private :: vanishes_real, vanishes_complex
  public :: expect
  interface expect
     module procedure expect_integer, expect_real, expect_complex, &
          expect_real_integer, expect_integer_real, &
          expect_complex_integer, expect_integer_complex, &
          expect_complex_real, expect_real_complex
  end interface
  private :: expect_integer, expect_real, expect_complex, &
       expect_real_integer, expect_integer_real, &
       expect_complex_integer, expect_integer_complex, &
       expect_complex_real, expect_real_complex
  public :: expect_zero
  interface expect_zero
     module procedure expect_zero_integer, expect_zero_real, expect_zero_complex
  end interface
  private :: expect_zero_integer, expect_zero_real, expect_zero_complex
  public :: print_matrix
contains
  elemental function agreement_real (x, y, base) result (a)
    real(kind=default) :: a
    real(kind=default), intent(in) :: x, y
    real(kind=default), intent(in), optional :: base
    real(kind=default) :: scale, dxy
    if (present (base)) then
       scale = max (abs (x), abs (y), abs (base))
    else
       scale = max (abs (x), abs (y))
    end if
    if (ieee_is_nan (x) .or. ieee_is_nan (y)) then
       a = 0
    else if (scale <= 0) then
       a = -1
    else
       dxy = abs (x - y) / scale
       if (dxy <= 0.0_default) then
          a = 1
       else
          a = log (dxy) / log (epsilon (scale))
          a = max (0.0_default, min (1.0_default, a))
          if (ieee_is_nan (a)) then
             a = 0
          end if
       end if
    end if
    if (ieee_is_nan (a)) then
       a = 0
    end if
  end function agreement_real
  elemental function ieee_is_nan (x) result (yorn)
    logical :: yorn
    real (kind=default), intent(in) :: x
    yorn = (x /= x)
  end function ieee_is_nan
  elemental function agreement_complex (x, y, base) result (a)
    real(kind=default) :: a
    complex(kind=default), intent(in) :: x, y
    real(kind=default), intent(in), optional :: base
    real(kind=default) :: scale, dxy
    if (present (base)) then
       scale = max (abs (x), abs (y), abs (base))
    else
       scale = max (abs (x), abs (y))
    end if
    if (      ieee_is_nan (real (x, kind=default)) .or. ieee_is_nan (aimag (x)) &
         .or. ieee_is_nan (real (y, kind=default)) .or. ieee_is_nan (aimag (y))) then
       a = 0
    else if (scale <= 0) then
       a = -1
    else
       dxy = abs (x - y) / scale
       if (dxy <= 0.0_default) then
          a = 1
       else
          a = log (dxy) / log (epsilon (scale))
          a = max (0.0_default, min (1.0_default, a))
          if (ieee_is_nan (a)) then
             a = 0
          end if
       end if
    end if
    if (ieee_is_nan (a)) then
       a = 0
    end if
  end function agreement_complex
  elemental function agreement_real_complex (x, y, base) result (a)
    real(kind=default) :: a
    real(kind=default), intent(in) :: x
    complex(kind=default), intent(in) :: y
    real(kind=default), intent(in), optional :: base
    a = agreement_complex (cmplx (x, kind=default), y, base)
  end function agreement_real_complex
  elemental function agreement_complex_real (x, y, base) result (a)
    real(kind=default) :: a
    complex(kind=default), intent(in) :: x
    real(kind=default), intent(in) :: y
    real(kind=default), intent(in), optional :: base
    a = agreement_complex (x, cmplx (y, kind=default), base)
  end function agreement_complex_real
  elemental function agreement_integer_complex (x, y, base) result (a)
    real(kind=default) :: a
    integer, intent(in) :: x
    complex(kind=default), intent(in) :: y
    real(kind=default), intent(in), optional :: base
    a = agreement_complex (cmplx (x, kind=default), y, base)
  end function agreement_integer_complex
  elemental function agreement_complex_integer (x, y, base) result (a)
    real(kind=default) :: a
    complex(kind=default), intent(in) :: x
    integer, intent(in) :: y
    real(kind=default), intent(in), optional :: base
    a = agreement_complex (x, cmplx (y, kind=default), base)
  end function agreement_complex_integer
  elemental function agreement_integer_real (x, y, base) result (a)
    real(kind=default) :: a
    integer, intent(in) :: x
    real(kind=default), intent(in) :: y
    real(kind=default), intent(in), optional :: base
    a = agreement_real (real(x, kind=default), y, base)
  end function agreement_integer_real
  elemental function agreement_real_integer (x, y, base) result (a)
    real(kind=default) :: a
    real(kind=default), intent(in) :: x
    integer, intent(in) :: y
    real(kind=default), intent(in), optional :: base
    a = agreement_real (x, real (y, kind=default), base)
  end function agreement_real_integer
  elemental function vanishes_real (x, scale) result (a)
    real(kind=default) :: a
    real(kind=default), intent(in) :: x
    real(kind=default), intent(in), optional :: scale
    real(kind=default) :: scaled_x
    if (x == 0.0_default) then
       a = 1
       return
    else if (ieee_is_nan (x)) then
       a = 0
       return
    end if
    scaled_x = x
    if (present (scale)) then
       if (scale /= 0) then
          scaled_x = x / abs (scale)
       else
          a = 0
          return
       end if
    else
    end if
    a = log (abs (scaled_x)) / log (epsilon (scaled_x))
    a = max (0.0_default, min (1.0_default, a))
    if (ieee_is_nan (a)) then
       a = 0
    end if
  end function vanishes_real
  elemental function vanishes_complex (x, scale) result (a)
    real(kind=default) :: a
    complex(kind=default), intent(in) :: x
    real(kind=default), intent(in), optional :: scale
    a = vanishes_real (abs (x), scale)
  end function vanishes_complex
  subroutine expect_integer (x, x0, msg, passed, quiet, buffer, unit)
    integer, intent(in) :: x, x0
    character(len=*), intent(in) :: msg
    logical, intent(inout), optional :: passed
    logical, intent(in), optional :: quiet
    character(len=*), intent(inout), optional :: buffer
    integer, intent(in), optional :: unit
    logical :: failed, verbose
    character(len=*), parameter :: fmt = "(1X,A,': ',A)"
    character(len=*), parameter :: &
         fmt_verbose = "(1X,A,': ',A,' [expected ',I6,', got ',I6,']')"
    failed = .false.
    verbose = .true.
    if (present (quiet)) then
       verbose = .not.quiet
    end if
    if (x == x0) then
       if (verbose) then
          if (.not. (present (buffer) .or. present (unit))) then
             write (unit = *, fmt = fmt) msg, "passed"
          end if
          if (present (unit)) then
             write (unit = unit, fmt = fmt) msg, "passed"
          end if
          if (present (buffer)) then
             write (unit = buffer, fmt = fmt) msg, "passed"
          end if
       end if
    else
       if (.not. (present (buffer) .or. present (unit))) then
          write (unit = *, fmt = fmt_verbose) msg, "failed", x0, x
       end if
       if (present (unit)) then
          write (unit = unit, fmt = fmt_verbose) msg, "failed", x0, x
       end if
       if (present (buffer)) then
          write (unit = buffer, fmt = fmt_verbose) msg, "failed", x0, x
       end if
       failed = .true.
    end if
    if (present (passed)) then
       passed = passed .and. .not.failed
    end if
  end subroutine expect_integer
  subroutine expect_real (x, x0, msg, passed, threshold, quiet, abs_threshold)
    real(kind=default), intent(in) :: x, x0
    character(len=*), intent(in) :: msg
    logical, intent(inout), optional :: passed
    real(kind=default), intent(in), optional :: threshold
    real(kind=default), intent(in), optional :: abs_threshold
    logical, intent(in), optional :: quiet
    logical :: failed, verbose
    real(kind=default) :: agreement_threshold, abs_agreement_threshold
    character(len=*), parameter :: fmt = "(1X,A,': ',A,' at ',I4,'%')"
    character(len=*), parameter :: fmt_verbose = "(1X,A,': ',A,' at ',I4,'%'," // &
         "' [expected ',E10.3,', got ',E10.3,']')"
    real(kind=default) :: a
    failed = .false.
    verbose = .true.
    if (present (quiet)) then
       verbose = .not.quiet
    end if
    if (x == x0) then
       if (verbose) then
          write (unit = *, fmt = fmt) msg, "passed", 100
       end if
    else
       if (x0 == 0) then
          a = vanishes (x)
       else
          a = agreement (x, x0)
       end if
       if (present (threshold)) then
          agreement_threshold = threshold
       else
          agreement_threshold = THRESHOLD_DEFAULT
       end if
       if (present (abs_threshold)) then
          abs_agreement_threshold = abs_threshold
       else
          abs_agreement_threshold = ABS_THRESHOLD_DEFAULT
       end if
       if (a >= agreement_threshold .or. &
           max(abs(x), abs(x0)) <= abs_agreement_threshold) then
          if (verbose) then
             if (a >= THRESHOLD_WARN) then
                write (unit = *, fmt = fmt) msg, "passed", int (a * 100)
             else
                write (unit = *, fmt = fmt_verbose) msg, "passed", int (a * 100), x0, x
             end if
          end if
       else
          failed = .true.
          write (unit = *, fmt = fmt_verbose) msg, "failed", int (a * 100), x0, x
       end if
    end if
    if (present (passed)) then
       passed = passed .and. .not. failed
    end if
  end subroutine expect_real
  subroutine expect_complex (x, x0, msg, passed, threshold, quiet, abs_threshold)
    complex(kind=default), intent(in) :: x, x0
    character(len=*), intent(in) :: msg
    logical, intent(inout), optional :: passed
    real(kind=default), intent(in), optional :: threshold
    real(kind=default), intent(in), optional :: abs_threshold
    logical, intent(in), optional :: quiet
    logical :: failed, verbose
    real(kind=default) :: agreement_threshold, abs_agreement_threshold
    character(len=*), parameter :: fmt = "(1X,A,': ',A,' at ',I4,'%')"
    character(len=*), parameter :: fmt_verbose = "(1X,A,': ',A,' at ',I4,'%'," // &
         "' [expected (',E10.3,',',E10.3,'), got (',E10.3,',',E10.3,')]')"
    character(len=*), parameter :: fmt_phase = "(1X,A,': ',A,' at ',I4,'%'," // &
         "' [modulus passed at ',I4,'%',', phases ',F5.3,' vs. ',F5.3,']')"
    real(kind=default) :: a, a_modulus
    failed = .false.
    verbose = .true.
    if (present (quiet)) then
       verbose = .not.quiet
    end if
    if (x == x0) then
       if (verbose) then
          write (unit = *, fmt = fmt) msg, "passed", 100
       end if
    else
       if (x0 == 0) then
          a = vanishes (x)
       else
          a = agreement (x, x0)
       end if
       if (present (threshold)) then
          agreement_threshold = threshold
       else
          agreement_threshold = THRESHOLD_DEFAULT
       end if
       if (present (abs_threshold)) then
          abs_agreement_threshold = abs_threshold
       else
          abs_agreement_threshold = ABS_THRESHOLD_DEFAULT
       end if
       if (a >= agreement_threshold .or. &
           max(abs(x), abs(x0)) <= abs_agreement_threshold) then
          if (verbose) then
             if (a >= THRESHOLD_WARN) then
                write (unit = *, fmt = fmt) msg, "passed", int (a * 100)
             else
                write (unit = *, fmt = fmt_verbose) msg, "passed", int (a * 100), x0, x
             end if
          end if
       else
          a_modulus = agreement (abs (x), abs (x0))
          if (a_modulus >= agreement_threshold) then
             write (unit = *, fmt = fmt_phase) msg, "failed", int (a * 100), &
                  int (a_modulus * 100), &
                  atan2 (real (x, kind=default), aimag (x)), &
                  atan2 (real (x0, kind=default), aimag (x0))
          else
             write (unit = *, fmt = fmt_verbose) msg, "failed", int (a * 100), x0, x
          end if
          failed = .true.
       end if
    end if
    if (present (passed)) then
       passed = passed .and. .not.failed
    end if
  end subroutine expect_complex
  subroutine expect_real_integer (x, x0, msg, passed, threshold, quiet)
    real(kind=default), intent(in) :: x
    integer, intent(in) :: x0
    character(len=*), intent(in) :: msg
    real(kind=default), intent(in), optional :: threshold
    logical, intent(inout), optional :: passed
    logical, intent(in), optional :: quiet
    call expect_real (x, real (x0, kind=default), msg, passed, threshold, quiet)
  end subroutine expect_real_integer
  subroutine expect_integer_real (x, x0, msg, passed, threshold, quiet)
    integer, intent(in) :: x
    real(kind=default), intent(in) :: x0
    character(len=*), intent(in) :: msg
    real(kind=default), intent(in), optional :: threshold
    logical, intent(inout), optional :: passed
    logical, intent(in), optional :: quiet
    call expect_real (real (x, kind=default), x0, msg, passed, threshold, quiet)
  end subroutine expect_integer_real
  subroutine expect_complex_integer (x, x0, msg, passed, threshold, quiet)
    complex(kind=default), intent(in) :: x
    integer, intent(in) :: x0
    character(len=*), intent(in) :: msg
    logical, intent(inout), optional :: passed
    real(kind=default), intent(in), optional :: threshold
    logical, intent(in), optional :: quiet
    call expect_complex (x, cmplx (x0, kind=default), msg, passed, threshold, quiet)
  end subroutine expect_complex_integer
  subroutine expect_integer_complex (x, x0, msg, passed, threshold, quiet)
    integer, intent(in) :: x
    complex(kind=default), intent(in) :: x0
    character(len=*), intent(in) :: msg
    logical, intent(inout), optional :: passed
    real(kind=default), intent(in), optional :: threshold
    logical, intent(in), optional :: quiet
    call expect_complex (cmplx (x, kind=default), x0, msg, passed, threshold, quiet)
  end subroutine expect_integer_complex
  subroutine expect_complex_real (x, x0, msg, passed, threshold, quiet)
    complex(kind=default), intent(in) :: x
    real(kind=default), intent(in) :: x0
    character(len=*), intent(in) :: msg
    logical, intent(inout), optional :: passed
    real(kind=default), intent(in), optional :: threshold
    logical, intent(in), optional :: quiet
    call expect_complex (x, cmplx (x0, kind=default), msg, passed, threshold, quiet)
  end subroutine expect_complex_real
  subroutine expect_real_complex (x, x0, msg, passed, threshold, quiet)
    real(kind=default), intent(in) :: x
    complex(kind=default), intent(in) :: x0
    character(len=*), intent(in) :: msg
    logical, intent(inout), optional :: passed
    real(kind=default), intent(in), optional :: threshold
    logical, intent(in), optional :: quiet
    call expect_complex (cmplx (x, kind=default), x0, msg, passed, threshold, quiet)
  end subroutine expect_real_complex
  subroutine expect_zero_integer (x, msg, passed)
    integer, intent(in) :: x
    character(len=*), intent(in) :: msg
    logical, intent(inout), optional :: passed
    call expect_integer (x, 0, msg, passed)
  end subroutine expect_zero_integer
  subroutine expect_zero_real (x, scale, msg, passed, threshold, quiet)
    real(kind=default), intent(in) :: x, scale
    character(len=*), intent(in) :: msg
    logical, intent(inout), optional :: passed
    real(kind=default), intent(in), optional :: threshold
    logical, intent(in), optional :: quiet
    logical :: failed, verbose
    real(kind=default) :: agreement_threshold
    character(len=*), parameter :: fmt = "(1X,A,': ',A,' at ',I4,'%')"
    character(len=*), parameter :: fmt_verbose = "(1X,A,': ',A,' at ',I4,'%'," // &
         "' [expected 0 (relative to ',E10.3,') got ',E10.3,']')"
    real(kind=default) :: a
    failed = .false.
    verbose = .true.
    if (present (quiet)) then
       verbose = .not.quiet
    end if
    if (x == 0) then
       if (verbose) then
          write (unit = *, fmt = fmt) msg, "passed", 100
       end if
    else
       a = vanishes (x, scale = scale)
       if (present (threshold)) then
          agreement_threshold = threshold
       else
          agreement_threshold = THRESHOLD_DEFAULT
       end if
       if (a >= agreement_threshold) then
          if (verbose) then
             if (a >= THRESHOLD_WARN) then
                write (unit = *, fmt = fmt) msg, "passed", int (a * 100)
             else
                write (unit = *, fmt = fmt_verbose) msg, "passed", int (a * 100), scale, x
             end if
          end if
       else
          failed = .true.
          write (unit = *, fmt = fmt_verbose) msg, "failed", int (a * 100), scale, x
       end if
    end if
    if (present (passed)) then
       passed = passed .and. .not.failed
    end if
  end subroutine expect_zero_real
  subroutine expect_zero_complex (x, scale, msg, passed, threshold, quiet)
    complex(kind=default), intent(in) :: x
    real(kind=default), intent(in) :: scale
    character(len=*), intent(in) :: msg
    logical, intent(inout), optional :: passed
    real(kind=default), intent(in), optional :: threshold
    logical, intent(in), optional :: quiet
    call expect_zero_real (abs (x), scale, msg, passed, threshold, quiet)
  end subroutine expect_zero_complex
  subroutine print_matrix (a)
    complex(kind=default), dimension(:,:), intent(in) :: a
    integer :: row
    do row = 1, size (a, dim=1)
       write (unit = *, fmt = "(10(tr2, f5.2, '+', f5.2, 'I'))") a(row,:)
    end do
  end subroutine print_matrix
end module omega_testtools
