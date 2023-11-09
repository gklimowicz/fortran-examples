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
module omega_color
  use kinds
  implicit none
  private
  public :: omega_color_factor
  type omega_color_factor
     integer :: i1, i2
     real(kind=default) :: factor
  end type omega_color_factor
  public :: omega_color_sum
  public :: ovm_color_sum
  integer, parameter, public :: omega_color_2010_01_A = 0
contains
  
  function omega_color_sum (flv, hel, amp, cf) result (amp2)
    complex(kind=default) :: amp2
    integer, intent(in) :: flv, hel
    complex(kind=default), dimension(:,:,:), intent(in) :: amp
    type(omega_color_factor), dimension(:), intent(in) :: cf
    integer :: n
    amp2 = 0
    !$omp parallel do reduction(+:amp2)
    do n = 1, size (cf)
       amp2 = amp2 + cf(n)%factor * &
                     amp(flv,cf(n)%i1,hel) * conjg (amp(flv,cf(n)%i2,hel))
    end do
    !$omp end parallel do
  end function omega_color_sum
  
  function ovm_color_sum (flv, hel, amp, cf) result (amp2)
    real(kind=default) :: amp2
    integer, intent(in) :: flv, hel
    complex(kind=default), dimension(:,:,:), intent(in) :: amp
    type(omega_color_factor), dimension(:), intent(in) :: cf
    integer :: n
    amp2 = 0
    !$omp parallel do reduction(+:amp2)
    do n = 1, size (cf)
       if (cf(n)%i1 == cf(n)%i2) then
          amp2 = amp2 + cf(n)%factor * &
                 real(amp(flv,cf(n)%i1,hel) * conjg(amp(flv,cf(n)%i2,hel)))
       else
          amp2 = amp2 + cf(n)%factor * 2 * &
                 real(amp(flv,cf(n)%i1,hel) * conjg(amp(flv,cf(n)%i2,hel)))
       end if
    end do
    !$omp end parallel do
  end function ovm_color_sum
end module omega_color
