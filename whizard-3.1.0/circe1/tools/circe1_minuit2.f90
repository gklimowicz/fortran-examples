! minuit2.f90 -- fitting for circe
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
module minuit2
  use kinds

  implicit none

  public :: fct
  public :: phi

contains

  subroutine fct (nx, df, f, a, mode, g)
    integer, intent(in) ::  nx, mode
    real(kind=double) :: f, g
    real(kind=double), dimension(:) :: df, a
        integer, parameter :: NDATA = 20
        real(kind=double) :: chi, chi2
        real(kind=double), dimension(2,NDATA,NDATA) :: xi
        real(kind=double), dimension(NDATA,NDATA) :: fi, dfi
        integer :: i, j, n
    if (mode .eq. 1) then
          open (10, file = 'minuit.data')
          do i = 1, NDATA
             do j = 1, NDATA
                read (10, *) xi(1,i,j), xi(2,i,j), fi(i,j), dfi(i,j)
                fi(i,j) = fi(i,j)/1d30
                dfi(i,j) = dfi(i,j)/1d30
             end do   
          end do
          close (10)
    else if (mode .eq. 2) then
          print *, "ERROR: $\nabla f$ n.a."
          stop
    end if
        f = 0d0
        do i = 1, NDATA
           do j = 1, NDATA
              if (dfi(i,j).gt.0d0) then
                 f = f + ((phi(xi(1,i,j),xi(2,i,j),a) &
                            - fi(i,j)) / dfi(i,j))**2
              end if
           end do
        end do
    if (mode .eq. 3) then
           chi2 = 0d0
           n = 0
           open (10, file = 'minuit.fit')
           do i = 1, NDATA
              do j = 1, NDATA
                 if (dfi(i,j).gt.0d0) then
                    chi = (phi(xi(1,i,j),xi(2,i,j),a)-fi(i,j))/dfi(i,j)
                    write (10,*) xi(1,i,j), xi(2,i,j), &
                                 1d30 * phi(xi(1,i,j),xi(2,i,j),a), &
                                 1d30 * fi(i,j), &
                                 chi
                    chi2 = chi2 + chi**2
                    n = n + 1
                 else
                    write (10,*) xi(1,i,j), xi(2,i,j), &
                                 1d30 * phi(xi(1,i,j),xi(2,i,j),a), &
                                 1d30 * fi(i,j)
                 end if
              end do
           end do
           close (10)
           print *, 'CHI2 = ', chi2/n
    end if
  end subroutine fct

  function phi (e1, e2, a)
    real(kind=double) :: e1, e2
    real(kind=double), dimension(33) :: a
    real(kind=double) :: phi
    real(kind=double) :: y1, y2
    y1 = e1 / 250d0
    y2 = e2 / 250d0
    phi = exp (                         &
         + a( 1) * 1d0                  &        
         + a( 2) * log(y1)              &
         + a( 3) * log(1d0-y1)          &
         + a( 4) * log(-log(y1))        &
         + a( 5) * log(-log(1d0-y1))    &
         + a( 6) * y1                   &
         + a( 7) * log(y1)**2           &
         + a( 8) * log(1d0-y1)**2       &
         + a( 9) * log(-log(y1))**2     &
         + a(10) * log(-log(1d0-y1))**2 &
         + a(11) * y1**2                &
         + a(12) / log(y1)              &
         + a(13) / log(1d0-y1)          &
         + a(14) / log(-log(y1))        &
         + a(15) / log(-log(1d0-y1))    &
         + a(16) / y1                   &
         + a(17) / (1d0-y1)             &
         + a(18) * log(y2)              &
         + a(19) * log(1d0-y2)          &
         + a(20) * log(-log(y2))        &
         + a(21) * log(-log(1d0-y2))    &
         + a(22) * y2                   &
         + a(23) * log(y2)**2           &
         + a(24) * log(1d0-y2)**2       &
         + a(25) * log(-log(y2))**2     &
         + a(26) * log(-log(1d0-y2))**2 &
         + a(27) * y2**2                &
         + a(28) / log(y2)              &
         + a(29) / log(1d0-y2)          &
         + a(30) / log(-log(y2))        &
         + a(31) / log(-log(1d0-y2))    &
         + a(32) / y2                   &
         + a(33) / (1d0-y2)             &
         )
    end function phi

end module minuit2

program fit
  use kinds
  use minuit2

  implicit none

  call minuit (fct, 0d0)
end program fit
