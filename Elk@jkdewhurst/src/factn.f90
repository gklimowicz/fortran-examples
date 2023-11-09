
! Copyright (C) 2022 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

elemental real(8) function factn(n)
implicit none
! arguments
integer, intent(in) :: n
! local variables
integer i
real(8), parameter :: f(24)=[ &
                       1.d0,                        2.d0,  &
                       6.d0,                       24.d0,  &
                     120.d0,                      720.d0,  &
                    5040.d0,                    40320.d0,  &
                  362880.d0,                  3628800.d0,  &
                39916800.d0,                479001600.d0,  &
              6227020800.d0,              87178291200.d0,  &
           1307674368000.d0,           20922789888000.d0,  &
         355687428096000.d0,         6402373705728000.d0,  &
      121645100408832000.d0,      2432902008176640000.d0,  &
    51090942171709440000.d0,   1124000727777607680000.d0,  &
 25852016738884976640000.d0, 620448401733239439360000.d0]
if (n.le.1) then
  factn=1.d0
else if (n.le.24) then
  factn=f(n)
else
  factn=f(24)
  do i=25,n
    factn=factn*dble(i)
  end do
end if
end function

