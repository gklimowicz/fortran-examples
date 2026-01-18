
! Copyright (C) 2011 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

!BOP
! !ROUTINE: mtdmin
! !INTERFACE:
pure subroutine mtdmin(is,js,dmin)
! !USES:
use modmain
! !INPUT/OUTPUT PARAMETERS:
!   is, js : species numbers (out,integer)
!   dmin   : minimum distance between muffin-tin surfaces (out,real)
! !DESCRIPTION:
!   Finds the atomic species pair for which the distance between the muffin-tin
!   surfaces is a minimum. This distance may be negative if the muffin-tins
!   overlap.
!
! !REVISION HISTORY:
!   Created October 2011 (JKD)
!EOP
!BOC
implicit none
! arguments
integer, intent(out) :: is,js
real(8), intent(out) :: dmin
! local variables
integer i1,i2,i3,ks,ka,ls,la
real(8) d2,d2m,r,r2,rm,dr,drm
real(8) v1,v2,v3,w1,w2,w3
is=1
js=1
drm=1.d8
do i1=-1,1; do i2=-1,1; do i3=-1,1
  v1=i1*avec(1,1)+i2*avec(1,2)+i3*avec(1,3)
  v2=i1*avec(2,1)+i2*avec(2,2)+i3*avec(2,3)
  v3=i1*avec(3,1)+i2*avec(3,2)+i3*avec(3,3)
  do ks=1,nspecies
    do ka=1,natoms(ks)
      w1=v1+atposc(1,ka,ks)
      w2=v2+atposc(2,ka,ks)
      w3=v3+atposc(3,ka,ks)
      do ls=1,nspecies
        r=rmt(ks)+rmt(ls)
        r2=r**2
        do la=1,natoms(ls)
          if ((i1.ne.0).or.(i2.ne.0).or.(i3.ne.0).or. &
           (ks.ne.ls).or.(ka.ne.la)) then
            d2=(w1-atposc(1,la,ls))**2 &
              +(w2-atposc(2,la,ls))**2 &
              +(w3-atposc(3,la,ls))**2
            dr=d2-r2
            if (dr.lt.drm-epslat) then
              is=ks
              js=ls
              rm=r
              d2m=d2
              drm=dr
            end if
          end if
        end do
      end do
    end do
  end do
end do; end do; end do
dmin=sqrt(d2m)-rm
end subroutine
!EOC

