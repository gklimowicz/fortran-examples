! test_qed_eemm.f90 --
! driver.f90 -- O'Mega self test driver
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Copyright (C) 1999-2022 by 
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!      Christian Speckner <cnspeckn@googlemail.com>
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

program test2_qed_eemm

  use kinds
  use parameters_QED
  use amplitude_qed_eemm

  real(default), dimension(0:3,4) :: p
  complex(default) :: a
  integer :: h, n_flv, n_hel, n_col

  n_flv = number_flavor_states ()
  n_hel = number_spin_states ()
  n_col = number_color_flows ()

  if (n_flv /= 1) then
     print *, "unexpected # of flavor combinations"
     stop 1
  end if

  if (n_hel /= 16) then
     print *, "unexpected # of helicity combinations"
     stop 1
  end if

  if (n_col /= 1) then
     print *, "unexpected # of color flows"
     stop 1
  end if

  call init_parameters

  p(:,1) = (/ 100.0_default, 0.0_default,     0.0_default,   100.0_default /)
  p(:,2) = (/ 100.0_default, 0.0_default,     0.0_default, - 100.0_default /)
  p(:,3) = (/ 100.0_default, 0.0_default,   100.0_default,     0.0_default /)
  p(:,4) = (/ 100.0_default, 0.0_default, - 100.0_default,     0.0_default /)

  call new_event (p)

  do h = 1, n_hel
     a = get_amplitude (1, h, 1)
     ! print *, "HEL = ", h, ", AMP = ", a
     if      (h == 6 .OR. h == 10) then
        if (abs (a + (0.0_default, 0.09_default)) > epsilon (1.0_default)) then
           print *, "unexpected value"
           stop 1
        end if
     else if (h == 7 .OR. h == 11) then
        if (abs (a - (0.0_default, 0.09_default)) > epsilon (1.0_default)) then
           print *, "unexpected value"
           stop 1
        end if
     end if
  end do

  stop 0

end program test2_qed_eemm

