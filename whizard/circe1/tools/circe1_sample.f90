! circe1_sample.f90 -- canonical beam spectra for linear collider physics
! 
! Copyright (C) 1999-2022 by 
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!     with contributions from
!     Christian Speckner <cnspeckn@googlemail.com>
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
! to the source 'circe1.nw'
module sample_routines
  use kinds
  use circe1 !NODEP!  

  implicit none
  private

  public :: sigma
  public :: d12
  public :: d1
  public :: d2
  public :: d12a
  public :: random
  public :: gauss1 
  public :: gaussx
  public :: gauss2

contains

  function sigma (s)
    real(kind=double) :: s, sigma
    sigma = 1d0 / s
  end function sigma

  function d12 (t1, t2)
    real(kind=double) :: d12, t1, t2, x1, x2
        real(kind=double), parameter :: EPS = 1d-6, PWR = 5d0
        real(kind=double), parameter :: KIREPS = 1D-6
    x1 = 1d0 - t1**PWR
    x2 = 1d0 - t2**PWR
    d12 = PWR*PWR * (t1*t2)**(PWR-1d0) &
           * sigma (x1*x2) * circee (x1, x2)
  end function d12

  function d1 (t1)
    real(kind=double) :: t1, x1, d1
        real(kind=double), parameter :: EPS = 1d-6, PWR = 5d0
        real(kind=double), parameter :: KIREPS = 1D-6
    x1 = 1d0 - t1**PWR
    d1 = PWR * t1**(PWR-1d0) * sigma (x1) * circee (x1, 1d0)
  end function d1

  function d2 (t2)
    real(kind=double) :: t2, x2, d2
        real(kind=double), parameter :: EPS = 1d-6, PWR = 5d0
        real(kind=double), parameter :: KIREPS = 1D-6
    x2 = 1d0 - t2**PWR
    d2 = PWR * t2**(PWR-1d0) * sigma (x2) * circee (1d0, x2)
  end function d2

  function d12a (x1, x2)
    real(kind=double) :: x1, x2, d12a
    d12a = sigma (x1*x2) * kirkee (x1, x2)
  end function d12a

  subroutine random (r)
    real(kind=double), intent(out) :: r
    integer :: m = 259200, a = 7141, c = 54773
    integer, save :: n = 0
    ! data n /0/
    n = mod(n*a+c,m)
    r = real (n, kind=double) / real (m, kind=double)
  end subroutine random

  function gauss1 (f, a, b, eps)
    real(kind=double) :: gauss1
      real(kind=double) :: f, a, b, eps
      external f
      real(kind=double), parameter :: Z1 = 1, HF = Z1/2, CST = 5*Z1/1000
      integer :: i
      real(kind=double) :: h, const, aa, bb, c1, c2, s8, s16, u
          real(kind=double), dimension(12), parameter :: &       
            x = (/ 9.6028985649753623d-1, &
                   7.9666647741362674d-1, &
                   5.2553240991632899d-1, &
                   1.8343464249564980d-1, &
                   9.8940093499164993d-1, &
                   9.4457502307323258d-1, &
                   8.6563120238783174d-1, &
                   7.5540440835500303d-1, &
                   6.1787624440264375d-1, &
                   4.5801677765722739d-1, &
                   2.8160355077925891d-1, &
                   9.5012509837637440d-2 /), &
            w = (/ 1.0122853629037626d-1, & 
                   2.2238103445337447d-1, &
                   3.1370664587788729d-1, &
                   3.6268378337836198d-1, &
                   2.7152459411754095d-2, &
                   6.2253523938647893d-2, &
                   9.5158511682492785d-2, &
                   1.2462897125553387d-1, &
                   1.4959598881657673d-1, &
                   1.6915651939500254d-1, &
                   1.8260341504492359d-1, &
                   1.8945061045506850d-1 /)
      h = 0
      if (b .eq. a) go to 99
      const = CST / dabs(b-a)
      bb = a
    1 continue
         aa = bb
         bb = b
    2 continue
         c1 = HF*(bb+aa)
         c2 = HF*(bb-aa)
         s8 = 0
         do i = 1, 4
            u = c2*x(i)
    s8 = s8 + w(i) * (f (c1+u) + f (c1-u))
         end do
         s16 = 0
         do i = 5, 12
            u = c2*x(i)
    s16 = s16 + w(i) * (f (c1+u) + f (c1-u))
             end do
             s16 = c2*s16
          if (dabs(s16-c2*s8) .le. eps*(1+dabs(s16))) then
             h = h + s16
             if (bb .ne. b) go to 1
          else
             bb = c1
             if (1 + const*dabs(c2) .ne. 1) go to 2
             h = 0
             print *, 'gauss: too high accuracy required'
             go to 99
          end if
       99 continue
    gauss1 = h
  end function gauss1

    function gaussx (f, y, a, b, eps)
      real(kind=double) :: y
      real(kind=double) :: gaussx
      real(kind=double) :: f, a, b, eps
      external f
      real(kind=double), parameter :: Z1 = 1, HF = Z1/2, CST = 5*Z1/1000
      integer :: i
      real(kind=double) :: h, const, aa, bb, c1, c2, s8, s16, u
          real(kind=double), dimension(12), parameter :: &       
            x = (/ 9.6028985649753623d-1, &
                   7.9666647741362674d-1, &
                   5.2553240991632899d-1, &
                   1.8343464249564980d-1, &
                   9.8940093499164993d-1, &
                   9.4457502307323258d-1, &
                   8.6563120238783174d-1, &
                   7.5540440835500303d-1, &
                   6.1787624440264375d-1, &
                   4.5801677765722739d-1, &
                   2.8160355077925891d-1, &
                   9.5012509837637440d-2 /), &
            w = (/ 1.0122853629037626d-1, & 
                   2.2238103445337447d-1, &
                   3.1370664587788729d-1, &
                   3.6268378337836198d-1, &
                   2.7152459411754095d-2, &
                   6.2253523938647893d-2, &
                   9.5158511682492785d-2, &
                   1.2462897125553387d-1, &
                   1.4959598881657673d-1, &
                   1.6915651939500254d-1, &
                   1.8260341504492359d-1, &
                   1.8945061045506850d-1 /)
      h = 0
      if (b .eq. a) go to 99
      const = CST / dabs(b-a)
      bb = a
    1 continue
         aa = bb
         bb = b
    2 continue
         c1 = HF*(bb+aa)
         c2 = HF*(bb-aa)
         s8 = 0
         do i = 1, 4
            u = c2*x(i)
      s8 = s8 + w(i) * (f (y, c1+u) + f (y, c1-u))
         end do
         s16 = 0
         do i = 5, 12
            u = c2*x(i)
      s16 = s16 + w(i) * (f (y, c1+u) + f (y, c1-u))
               end do
               s16 = c2*s16
            if (dabs(s16-c2*s8) .le. eps*(1+dabs(s16))) then
               h = h + s16
               if (bb .ne. b) go to 1
            else
               bb = c1
               if (1 + const*dabs(c2) .ne. 1) go to 2
               h = 0
               print *, 'gauss: too high accuracy required'
               go to 99
            end if
         99 continue
      gaussx = h
    end function gaussx

  function gauss2 (f, a, b, a1, b1, eps)   
    real(kind=double) :: a1, b1
    real(kind=double) :: gauss2
      real(kind=double) :: f, a, b, eps
      external f
      real(kind=double), parameter :: Z1 = 1, HF = Z1/2, CST = 5*Z1/1000
      integer :: i
      real(kind=double) :: h, const, aa, bb, c1, c2, s8, s16, u
          real(kind=double), dimension(12), parameter :: &       
            x = (/ 9.6028985649753623d-1, &
                   7.9666647741362674d-1, &
                   5.2553240991632899d-1, &
                   1.8343464249564980d-1, &
                   9.8940093499164993d-1, &
                   9.4457502307323258d-1, &
                   8.6563120238783174d-1, &
                   7.5540440835500303d-1, &
                   6.1787624440264375d-1, &
                   4.5801677765722739d-1, &
                   2.8160355077925891d-1, &
                   9.5012509837637440d-2 /), &
            w = (/ 1.0122853629037626d-1, & 
                   2.2238103445337447d-1, &
                   3.1370664587788729d-1, &
                   3.6268378337836198d-1, &
                   2.7152459411754095d-2, &
                   6.2253523938647893d-2, &
                   9.5158511682492785d-2, &
                   1.2462897125553387d-1, &
                   1.4959598881657673d-1, &
                   1.6915651939500254d-1, &
                   1.8260341504492359d-1, &
                   1.8945061045506850d-1 /)
      h = 0
      if (b .eq. a) go to 99
      const = CST / dabs(b-a)
      bb = a
    1 continue
         aa = bb
         bb = b
    2 continue
         c1 = HF*(bb+aa)
         c2 = HF*(bb-aa)
         s8 = 0
         do i = 1, 4
            u = c2*x(i)
    s8 = s8 + w(i) * (gaussx (f, c1+u, a1, b1, eps) &
                     + gaussx (f, c1-u, a1, b1, eps))
         end do
         s16 = 0
         do i = 5, 12
            u = c2*x(i)
    s16 = s16 + w(i) * (gaussx (f, c1+u, a1, b1, eps) &
                       + gaussx (f, c1-u, a1, b1, eps))
         end do
         s16 = c2*s16
      if (dabs(s16-c2*s8) .le. eps*(1+dabs(s16))) then
         h = h + s16
         if (bb .ne. b) go to 1
      else
         bb = c1
         if (1 + const*dabs(c2) .ne. 1) go to 2
         h = 0
         print *, 'gauss: too high accuracy required'
         go to 99
      end if
   99 continue
    gauss2 = h
  end function gauss2


end module sample_routines

  program circe1_sample
    use kinds
    use sample_routines
    use circe1

    implicit none

        integer, parameter :: SBAND  = 1
        integer, parameter :: TESLA  = 2
        integer, parameter :: XBAND  = 3
        integer, parameter :: JLCNLC = 3
        integer, parameter :: SBNDEE = 4
        integer, parameter :: TESLEE = 5
        integer, parameter :: XBNDEE = 6
        integer, parameter :: NLCH   = 7
        integer, parameter :: ILC    = 8
        integer, parameter :: CLIC   = 9
        integer, parameter :: NACC   = 9
          real(kind=double), parameter :: EPS = 1d-6, PWR = 5d0
          real(kind=double), parameter :: KIREPS = 1D-6
          real(kind=double) :: s       
          real(kind=double) :: w, s2, x1, x2
          integer, parameter :: NEVENT = 10000
          integer :: n
      integer :: acc, ver, i
      real(kind=double), dimension(9) :: roots(9) = &
         (/ 90D0,  170D0,  250D0,  350D0,  500D0, &
                 800D0, 1000D0, 1200D0, 1500D0 /)
      do acc = 1, NACC
      ! do acc = JLCNLC, NLCH, NLCH-JLCNLC
         do ver = 9, 9
            do i = 1, 9
               call circes (0d0, 0d0, roots(i), acc, ver, 20020328, 1)
               s = sigma (1d0) * circee (1d0, 1d0) &
                  + gauss1 (d1, 0d0, 1d0, EPS) & 
                  + gauss1 (d2, 0d0, 1d0, EPS) &
                  + gauss2 (d12, 0d0, 1d0, 0d0, 1d0, EPS)
               write (*, 1000) 'delta(sigma) (Gauss) =', (s-1d0)*100d0
               1000 format (1X, A22, 1X, F6.2, '%')
               s = gauss2 (d12a, 0d0, 1d0-KIREPS, 0d0, 1d0-KIREPS, EPS) &
                 + gauss2 (d12a, 0d0, 1d0-KIREPS, 1d0-KIREPS, 1d0, EPS) &
                 + gauss2 (d12a, 1d0-KIREPS, 1d0, 0d0, 1d0-KIREPS, EPS) &
                 + gauss2 (d12a, 1d0-KIREPS, 1d0, 1d0-KIREPS, 1d0, EPS)
               write (*, 1000) 'delta(sigma) (Gauss) =', (s-1d0)*100d0
               s = 0d0
               s2 = 0d0
               do n = 1, NEVENT
                 call gircee (x1, x2, random)
                 w = sigma (x1*x2)
                 s = s + w
                 s2 = s2 + w*w
               end do
               s = s / dble(NEVENT)
               s2 = s2 / dble(NEVENT)
               write (*, 1000) 'delta(sigma) (MC)    =', (s-1d0)*100d0
               write (*, 1000) '                   +/-', sqrt((s2-s*s)/dble(NEVENT))*100d0
            end do
         end do
      end do   
  end program circe1_sample

