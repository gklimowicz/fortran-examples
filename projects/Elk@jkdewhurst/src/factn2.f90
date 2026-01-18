
! Copyright (C) 2022 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

elemental real(8) function factn2(n)
implicit none
! arguments
integer, intent(in) :: n
! local variables
integer i
real(8), parameter :: f(38)=[ &
                       1.d0,                        2.d0,  &
                       3.d0,                        8.d0,  &
                      15.d0,                       48.d0,  &
                     105.d0,                      384.d0,  &
                     945.d0,                     3840.d0,  &
                   10395.d0,                    46080.d0,  &
                  135135.d0,                   645120.d0,  &
                 2027025.d0,                 10321920.d0,  &
                34459425.d0,                185794560.d0,  &
               654729075.d0,               3715891200.d0,  &
             13749310575.d0,              81749606400.d0,  &
            316234143225.d0,            1961990553600.d0,  &
           7905853580625.d0,           51011754393600.d0,  &
         213458046676875.d0,         1428329123020800.d0,  &
        6190283353629375.d0,        42849873690624000.d0,  &
      191898783962510625.d0,      1371195958099968000.d0,  &
     6332659870762850625.d0,     46620662575398912000.d0,  &
   221643095476699771875.d0,   1678343852714360832000.d0,  &
  8200794532637891559375.d0,  63777066403145711616000.d0]
if (n.le.1) then
  factn2=1.d0
else if (n.le.38) then
  factn2=f(n)
else
  factn2=dble(n)
  do i=1,n/2-1
    factn2=factn2*dble(n-2*i)
  end do
end if
end function

