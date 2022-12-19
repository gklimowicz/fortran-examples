! circe1_fit.f90 -- fitting for circe
! 
! Copyright (C) 1999-2022 by 
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
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
! to the source 'minuit.nw'

module fit_routines
  use kinds

  implicit none
  private

  public :: fct
  public :: gethst
  public :: fixerr
  public :: sumsqu
  public :: scale
  public :: phie
  public :: phig
  public :: pslice

contains
  subroutine fct (nx, df, f, a, mode, g)
    integer :: nx, mode
    real(kind=double) :: f, g
    real(kind=double), dimension(:) :: df, a     
        integer, parameter :: NDATA = 20
        real(kind=double) :: see, tee, dtee
        real(kind=double) :: seg, teg, dteg
        real(kind=double) :: sge, tge, dtge
        real(kind=double) :: sgg, tgg, dtgg
        real(kind=double), dimension(2,0:NDATA+1,0:NDATA+1) :: xee, xeg, &
             xge, xgg
        real(kind=double), dimension(0:NDATA+1,0:NDATA+1) :: fee, dfee, &
             feg, dfeg, fge, dfge, fgg, dfgg
        real(kind=double) :: pwr
        integer :: i, j
        real(kind=double) :: delta
            integer, parameter :: NPARAM = 6  
        real(kind=double), dimension(NPARAM) :: a1
        integer, parameter :: REV = 1
        integer :: ndof
        real(kind=double) :: ab
        real(kind=double) :: eplus, eminus, epara, corr
        integer :: n
    if (mode .eq. 1) then
                call gethst ('ee', NDATA, xee, fee, dfee, see, tee, pwr)
                call gethst ('eg', NDATA, xeg, feg, dfeg, seg, teg, pwr)
                call gethst ('ge', NDATA, xge, fge, dfge, sge, tge, pwr)
                call gethst ('gg', NDATA, xgg, fgg, dfgg, sgg, tgg, pwr)
                call fixerr (NDATA, dfee, 20d0, 30d0, 40d0)
                call fixerr (NDATA, dfeg, 15d0, 20d0,  0d0)
                call fixerr (NDATA, dfge, 15d0, 20d0,  0d0)
                call fixerr (NDATA, dfgg, 10d0,  0d0,  0d0)
                dtee = sumsqu (NDATA, dfee)
                dteg = sumsqu (NDATA, dfeg)
                dtge = sumsqu (NDATA, dfge)
                dtgg = sumsqu (NDATA, dfgg)
                call scale (NDATA, 1d0/tee, fee)
                call scale (NDATA, 1d0/tee, dfee)
                call scale (NDATA, 1d0/tee, feg)
                call scale (NDATA, 1d0/tee, dfeg)
                call scale (NDATA, 1d0/tee, fge)
                call scale (NDATA, 1d0/tee, dfge)
                call scale (NDATA, 1d0/tee, fgg)
                call scale (NDATA, 1d0/tee, dfgg)
    else if (mode .eq. 2) then
          print *, "ERROR: $\nabla f$ n.a."
          stop
    end if
        f = 0d0
        do i = 1, NDATA
           do j = 1, NDATA
              if (dfee(i,j) .gt. 0d0) then
                 f = f + ((phie(xee(1,i,j),a) * phie(xee(2,i,j),a) &
                            - fee(i,j)) / dfee(i,j))**2
              end if
              if (dfeg(i,j) .gt. 0d0) then
                 f = f + ((phie(xeg(1,i,j),a) * phig(xeg(2,i,j),a) &
                            - feg(i,j)) / dfeg(i,j))**2
              end if
              if (dfge(i,j) .gt. 0d0) then
                 f = f + ((phig(xge(1,i,j),a) * phie(xge(2,i,j),a) &
                            - fge(i,j)) / dfge(i,j))**2
              end if
              if (dfgg(i,j) .gt. 0d0) then
                 f = f + ((phig(xgg(1,i,j),a) * phig(xgg(2,i,j),a) &
                            - fgg(i,j)) / dfgg(i,j))**2
              end if
           end do
        end do
        if ((a(2) .le. -1d0) .or. (a(3) .le. -1d0/pwr)) then
           print *, "warning: discarding out-of-range a2/3: ", a(2), a(3)
               f = 1d100
        else
          delta = 1d0 - exp(a(1)) * beta(a(2)+1d0,a(3)+1d0/pwr) * dble(NDATA) / pwr
          if (delta .lt. 0d0) then
             print *, "warnimg: delta forced to 0 from ", delta
             delta = 0d0
          end if
        do i = 1, NDATA
           if (dfee(ndata+1,i) .gt. 0d0) then
              f = f + ((delta*phie(xee(2,ndata+1,i),a) &
                        - fee(ndata+1,i)) / dfee(ndata+1,i))**2
           end if
           if (dfeg(ndata+1,i) .gt. 0d0) then
              f = f + ((delta*phig(xeg(2,ndata+1,i),a) &
                        - feg(ndata+1,i)) / dfeg(ndata+1,i))**2
           end if
           if (dfee(i,ndata+1) .gt. 0d0) then
              f = f + ((delta*phie(xee(1,i,ndata+1),a) &
                        - fee(i,ndata+1)) / dfee(i,ndata+1))**2
           end if
           if (dfge(i,ndata+1) .gt. 0d0) then
              f = f + ((delta*phig(xge(1,i,ndata+1),a) &
                        - fge(i,ndata+1)) / dfge(i,ndata+1))**2
           end if
        end do
        if (dfee(ndata+1,ndata+1) .gt. 0d0) then
           f = f + ((delta*delta &
                     - fee(ndata+1,ndata+1)) / dfee(ndata+1,ndata+1))**2
        end if
    end if
    if (mode .eq. 3) then
           a1(1) = exp(a(1)) * dble(NDATA) / pwr
           a1(2) = a(2)
           a1(3) = a(3) - 1d0 + 1d0/pwr
           a1(4) = exp(a(4)) * dble(NDATA) / pwr
           a1(5) = a(5) - 1d0 + 1d0/pwr
           a1(6) = a(6)
           open (10, file = 'Parameters')
           write (10, 1000) REV, tee / 1D32
           write (10, 1001) REV, &
                 1d0 - a1(1) * beta(a1(2)+1d0,a1(3)+1d0), &
                 a1(1), a1(2), a1(3), a1(4), a1(5), a1(6), &
                 a1(4) * beta(a1(5)+1d0,a1(6)+1d0)
       1000 format ('      data xa5lum(@ENERGY@,@ACC@,', I2, ') / ', E12.5, ' /')
       1001 format ('      data (xa5(i,@ENERGY@,@ACC@,', I2 ,'),i=0,7) /', /, &
                    '     $  ', 4(E12.5,', '), /, &
                    '     $  ', 3(E12.5,', '), E12.5, ' /')
           close (10)
           delta = 1d0 - a1(1) * beta(a1(2)+1d0,a1(3)+1d0)
           print *, '< x_e > = ', delta + (1d0-delta)*(a1(2)+1d0)/(a1(2)+a1(3)+2d0)
           print *, '< x_g > = ', (a1(5)+1d0)/(a1(5)+a1(6)+2d0)
           ndof = 0
           do i = 0, ndata+1
              do j = 0, ndata+1
                 if (dfee(i,j) .gt. 0d0) ndof = ndof + 1
                 if (dfeg(i,j) .gt. 0d0) ndof = ndof + 1
                 if (dfge(i,j) .gt. 0d0) ndof = ndof + 1
                 if (dfgg(i,j) .gt. 0d0) ndof = ndof + 1
              end do
           end do
           print *, 'CHI2 = ', f / ndof
           open (10, file = 'Errors.tex')
           write (10, 1099) tee / 1d32, dtee / 1d32, dtee / 1d32
       1099 format ('$', F8.2, '_{-', F4.2, '}^{+', F4.2, '}$')
           call mnerrs (1,  eplus, eminus, epara, corr)
           ab = a1(1) * beta(a1(2)+1d0,a1(3)+1d0)
           write (10, 1100) ab, abs (ab*eminus), abs (ab*eplus)
       1100 format ('$', F8.4, '_{-', F6.4, '}^{+', F6.4, '}$')
           do i = 2, 3
              call mnerrs (i,  eplus, eminus, epara, corr)
              write (10, 1100) a1(i), abs (eminus), abs (eplus)
           end do
           call mnerrs (4,  eplus, eminus, epara, corr)
           ab = a1(4) * beta(a1(5)+1d0,a1(6)+1d0)
           write (10, 1100) ab, abs (ab*eminus), abs (ab*eplus)
           do i = 5, 6
              call mnerrs (i,  eplus, eminus, epara, corr)
              write (10, 1100) a1(i), abs (eminus), abs (eplus)
           end do
           close (10)
           do n = 1, 10
              call pslice ('ee','x',n,NDATA,xee,fee,dfee,phie,phie,a)
              call pslice ('eg','x',n,NDATA,xeg,feg,dfeg,phie,phig,a)
              call pslice ('ge','x',n,NDATA,xge,fge,dfge,phig,phie,a)
              call pslice ('gg','x',n,NDATA,xgg,fgg,dfgg,phig,phig,a)
              call pslice ('ee','y',n,NDATA,xee,fee,dfee,phie,phie,a)
              call pslice ('eg','y',n,NDATA,xeg,feg,dfeg,phie,phig,a)
              call pslice ('ge','y',n,NDATA,xge,fge,dfge,phig,phie,a)
              call pslice ('gg','y',n,NDATA,xgg,fgg,dfgg,phig,phig,a)
           end do
           call pslice ('ee','x',21,NDATA,xee,fee,dfee,phie,phie,a)
           call pslice ('eg','x',21,NDATA,xeg,feg,dfeg,phie,phig,a)
           call pslice ('ee','y',21,NDATA,xee,fee,dfee,phie,phie,a)
           call pslice ('ge','y',21,NDATA,xge,fge,dfge,phig,phie,a)
           open (10, file = 'Slices.mp4')
           write (10,*) "picture eslice[], gslice[];"
           do n = 1, NDATA
              write (10,*) 'eslice[', n, '] := ', &
                   'btex $x_{e^\\pm} = ', xee(1,n,1), '$ etex;'
              write (10,*) 'gslice[', n, '] := ', &
                   'btex $x_\\gamma = ', xgg(1,n,1), '$ etex;'
           end do
           close (10)
    end if
  end subroutine fct

  subroutine gethst (tag, ndata, x, f, df, s, t, pwr)
    character(len=2) :: tag
    integer :: ndata
    real(kind=double) :: s, t, pwr
    real(kind=double), dimension(2,0:ndata+1,0:ndata+1) :: x
    real(kind=double), dimension(0:ndata+1,0:ndata+1) :: f, df
    integer :: i, j
    open (10, file = 'lumidiff-'//tag//'.dat')
    read (10, *) pwr
    s = 0d0
        do i = 1, ndata
           do j = 1, ndata
              read (10, *) x(1,i,j), x(2,i,j), f(i,j), df(i,j)
              s = s + f(i,j)
           end do
        end do
    t = s
        do i = 1, ndata
           read (10, *) x(1,i,0), f(i,0), df(i,0), &
                                  f(i,ndata+1), df(i,ndata+1)
           x(1,i,ndata+1) = x(1,i,0)
           t = t + f(i,0) + f(i,ndata+1)
        end do
        do i = 1, ndata
           read (10, *) x(2,0,i), f(0,i), df(0,i), &
                                  f(ndata+1,i), df(ndata+1,i)
           x(2,ndata+1,i) = x(2,0,i)
           t = t + f(0,i) + f(ndata+1,i)
        end do
        read (10, *) f(0,0), df(0,0), f(0,ndata+1), df(0,ndata+1)
        t = t + f(0,0) + f(0,ndata+1)
        read (10, *) f(ndata+1,0), df(ndata+1,0), &
                     f(ndata+1,ndata+1), df(ndata+1,ndata+1)
        t = t + f(ndata+1,0) + f(ndata+1,ndata+1)
    close (10)
  end subroutine gethst

  subroutine fixerr (ndata, df, c, sd, dd)
    integer :: ndata
    real(kind=double) :: c, sd, dd
    real(kind=double), dimension(0:ndata+1,0:ndata+1) :: df
    integer :: i, j
    do i = 1, NDATA
       do j = 1, NDATA
          df(i,j) = c * df(i,j)
       end do
    end do
    do i = 1, NDATA
       df(0,i) = sd * df(0,i)
       df(i,0) = sd * df(i,0)
       df(ndata+1,i) = sd * df(ndata+1,i)
       df(i,ndata+1) = sd * df(i,ndata+1)
    end do
    df(0,0) = dd * df(0,0)
    df(ndata+1,0) = dd * df(ndata+1,0)
    df(0,ndata+1) = dd * df(0,ndata+1)
    df(ndata+1,ndata+1) = dd * df(ndata+1,ndata+1)
  end subroutine fixerr

  function sumsqu (ndata, f)
    integer :: ndata
    real(kind=double) :: sumsqu
    real(kind=double), dimension(0:ndata+1,0:ndata+1) :: f
    integer :: i, j
    real(kind=double) :: s2
    s2 = 0
    do i = 0, NDATA+1
       do j = 0, NDATA+1
          s2 = s2 + f(i,j)*f(i,j)
       end do
    end do
    sumsqu = sqrt (s2)
  end function sumsqu

  subroutine scale (ndata, s, f)
    integer :: ndata
    real(kind=double) :: s
    real(kind=double), dimension(0:ndata+1,0:ndata+1) :: f
    integer :: i, j
    do i = 0, NDATA+1
       do j = 0, NDATA+1
          f(i,j) = s * f(i,j)
       end do
    end do
  end subroutine scale

  function phie (x, a)
    real(kind=double) :: x, phie 
    real(kind=double), dimension(6) :: a
    phie = exp (a(1) + a(2)*log(x) + a(3)*log(1d0-x))
  end function phie

  function phig (x, a)
    real(kind=double) :: x, phig
    real(kind=double), dimension(6) :: a
    phig = exp (a(4) + a(5)*log(x) + a(6)*log(1d0-x))
  end function phig

  subroutine pslice (pp, xy, n, ndata, x, f, df, phi1, phi2, a)
    character(len=2) :: pp
    character(len=1) :: xy
    integer :: n, ndata
    real(kind=double), dimension(2,0:ndata+1,0:ndata+1) :: x
    real(kind=double), dimension(0:ndata+1,0:ndata+1) :: f, df
    real(kind=double), dimension(6) :: a
    real(kind=double) :: z
    real(kind=double) :: phi1, phi2, d, delta, pwr
    external phi1, phi2
    integer :: i
    character(len=2) digits
    write (digits, '(I2.2)') n
    open (10, file = 'lumidiff-'//pp//xy//digits//'.dat')
    open (11, file = 'lumidiff-'//pp//xy//digits//'.fit')
    open (12, file = 'lumidiff-'//pp//xy//digits//'.chi')
    if (n .eq. ndata+1) then
       pwr = 5d0
       delta = 1d0 - exp(a(1))*beta(a(2)+1d0,a(3)+1d0/pwr) &
                       * dble(NDATA) / pwr
    else
       delta = 0
    end if
    if (xy .eq. 'x') then
       do i = 1, ndata
          if (df(n,i) .gt. 0d0) then
             if (pp(2:2) .eq. 'g') then
                z = x(2,n,i)
             else
                z = 1d0 - x(2,n,i)
             endif
             if (n .eq. ndata+1) then
                d = delta*phi2(x(2,n,i),a)
             else
                d = phi1(x(1,n,i),a)*phi2(x(2,n,i),a)
             endif
             write (10,*) z, f(n,i), df(n,i)
             write (11,*) z, d
             write (12,*) z, (f(n,i) - d) / df(n,i)
          endif
       end do
    else if (xy .eq. 'y') then
       do i = 1, ndata
          if (df(i,n) .gt. 0d0) then
             if (pp(1:1) .eq. 'g') then
                z = x(1,i,n)
             else
                z = 1d0 - x(1,i,n)
             endif
             if (n .eq. ndata+1) then
                d = phi1(x(1,i,n),a)*delta
             else
                d = phi1(x(1,i,n),a)*phi2(x(2,i,n),a)
             endif
             write (10,*) z, f(i,n), df(i,n)
             write (11,*) z, d
             write (12,*) z, (f(i,n) - d) / df(i,n)
          endif
       end do
    endif
    close (10)
    close (11)
    close (12)
  end subroutine pslice

  function beta (a, b)
    real(kind=double) :: a, b, beta
    beta = exp (dlgama(a) + dlgama(b) - dlgama(a+b))
    contains
      function dlgama (x)
        real(kind=double) :: dlgama
        real(kind=double), dimension(7) :: p1, q1, p2, q2, p3, q3
        real(kind=double), dimension(5) :: c, xl
        real(kind=double) :: x, y, zero, one, two, half, ap, aq
        integer :: i
        data ZERO /0.0D0/, ONE /1.0D0/, TWO /2.0D0/, HALF /0.5D0/
        data XL /0.0D0,0.5D0,1.5D0,4.0D0,12.0D0/
        data p1 /+3.8428736567460D+0, +5.2706893753010D+1, &
                 +5.5584045723515D+1, -2.1513513573726D+2, &
                 -2.4587261722292D+2, -5.7500893603041D+1, &
             -2.3359098949513D+0/
        data q1 /+1.0000000000000D+0, +3.3733047907071D+1, &
             +1.9387784034377D+2, +3.0882954973424D+2, &
             +1.5006839064891D+2, +2.0106851344334D+1, &
             +4.5717420282503D-1/
        data p2 /+4.8740201396839D+0, +2.4884525168574D+2, &
             +2.1797366058896D+3, +3.7975124011525D+3, &
             -1.9778070769842D+3, -3.6929834005591D+3, &
             -5.6017773537804D+2/
        data q2 /+1.0000000000000D+0, +9.5099917418209D+1, &
             +1.5612045277929D+3, +7.2340087928948D+3, &
             +1.0459576594059D+4, +4.1699415153200D+3, &
             +2.7678583623804D+2/
        data p3 /-6.8806240094594D+3, -4.3069969819571D+5, &
                 -4.7504594653440D+6, -2.9423445930322D+6, &
             +3.6321804931543D+7, -3.3567782814546D+6, &
             -2.4804369488286D+7/
        data q3 /+1.0000000000000D+0, -1.4216829839651D+3, &
             -1.5552890280854D+5, -3.4152517108011D+6, &
             -2.0969623255804D+7, -3.4544175093344D+7, &
             -9.1605582863713D+6/
        data c / 1.1224921356561D-1,  7.9591692961204D-2, &
                -1.7087794611020D-3,  9.1893853320467D-1, &
                 1.3469905627879D+0/
        if (x .le. xl(1)) then
             print *, 'ERROR: DLGAMA non positive argument: ', X
                 dlgama = zero
        end if
        if (x .le. xl(2)) then
          y = x + one
          ap = p1(1)
          aq = q1(1)
          do i = 2, 7
            ap = p1(i) + y * ap
            aq = q1(i) + y * aq
          end do
          y = - log(x) + x * ap / aq
        else if (x .le. xl(3)) then
          ap = p1(1)
          aq = q1(1)
          do i = 2, 7
             ap = p1(i) + x * ap
             aq = q1(i) + x * aq
          end do
          y = (x - one) * ap / aq
        else if (x .le. xl(4)) then
          ap = p2(1)
          aq = q2(1)
          do i = 2, 7
             ap = p2(i) + x * ap
             aq = q2(i) + x * aq
          end do         
          y = (x-two) * ap / aq
        else if (x .le. xl(5)) then
          ap = p3(1)
          aq = q3(1)
          do i = 2, 7
             ap = p3(i) + x * ap
             aq = q3(i) + x * aq
          end do 
          y = ap / aq
        else
         y = one / x**2
         y = (x-half) * log(x) - x + c(4) + &
             (c(1) + y * (c(2) + y * c(3))) / ((c(5) + y) * x)
        end if
        dlgama = y
      end function dlgama
  end function beta

end module fit_routines

program fit
  use kinds       
  use fit_routines

  implicit none

  integer :: i, rcode
    integer, parameter :: NPARAM = 6  
  integer, dimension(NPARAM) :: pnum
  character(len=10), dimension(NPARAM) :: pname
  real(kind=double), dimension(NPARAM) :: pstart, pstep
  integer, parameter :: ARGC = 10
  real(kind=double), dimension(ARGC) :: argv

        data pnum   /     1,     2,       3,     4,     5,       6 /
        data pname  / '1_e', 'x_e', '1-x_e', '1_g', 'x_g', '1-x_g' /
        data pstart / -1.00, 20.00,    0.20, -1.00,  0.20,   20.00 /
        data pstep  /  0.01,  0.01,    0.01,  0.01,  0.01,    0.01 /
      call mninit (5, 6, 7)
        do i = 1, NPARAM
           call mnparm (pnum(i), pname(i), pstart(i), pstep (i), 0d0, 0d0, rcode)
           if (rcode .ne. 0) then
              print *, "fit: MINUIT won''t accept parameter ", pnum(i)
              stop
           endif
        end do
      call mnseti ('CIRCE: fit version 1     ')
      argv(1) = 1
      call mnexcm (fct, 'SET PRINTOUT        ', argv, 1, rcode, 0d0)
      argv(1) = 1
      call mnexcm (fct, 'CALL FCT            ', argv, 1, rcode, 0d0)
      call mnexcm (fct, 'MIGRAD              ', argv, 0, rcode, 0d0)
      call mnexcm (fct, 'MINOS               ', argv, 0, rcode, 0d0)
      argv(1) = 3
      call mnexcm (fct, 'CALL FCT            ', argv, 1, rcode, 0d0)
      call mnexcm (fct, 'STOP                ', argv, 0, rcode, 0d0)

end program fit

