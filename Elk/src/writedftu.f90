
! Copyright (C) 2007 F. Bultmark, F. Cricchio and L. Nordstrom.
! This file is distributed under the terms of the GNU Lesser General Public
! License. See the file COPYING for license details.

subroutine writedftu
use modmain
use moddftu
implicit none
! local variables
integer ispn,jspn,idu
integer is,ia,ias
integer k,l,ll,m1,m2,lm1,lm2
if (dftu.eq.0) return
! machine readable density matrix file
open(50,file='DMATMT'//trim(filext),form='FORMATTED',action='WRITE')
do idu=1,ndftu
  is=isldu(1,idu)
  l=isldu(2,idu)
  ll=l*(l+1)+1
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,*)
    write(50,*)
    write(50,'(3I4," : species, atom, l")') is,ia,l
    do ispn=1,nspinor
      do jspn=1,nspinor
        write(50,*)
        write(50,'(2I4," : ispn, jspn; m1, m2, dmatmt below")') ispn,jspn
        do m1=-l,l
          lm1=ll+m1
          do m2=-l,l
            lm2=ll+m2
            write(50,'(2I6," ",2G18.10)') m1,m2,dmatmt(lm1,ispn,lm2,jspn,ias)
          end do
        end do
      end do
    end do
  end do
end do
close(50)
! machine readable potential matrix file
open(50,file='VMATMT'//trim(filext),form='FORMATTED',action='WRITE')
do idu=1,ndftu
  is=isldu(1,idu)
  l=isldu(2,idu)
  ll=l*(l+1)+1
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    write(50,*)
    write(50,*)
    write(50,'(3I4," : species, atom, l")') is,ia,l
    do ispn=1,nspinor
      do jspn=1,nspinor
        write(50,*)
        write(50,'(2I4," : ispn, jspn; m1, m2, vmatmt below")') ispn,jspn
        do m1=-l,l
          lm1=ll+m1
          do m2=-l,l
            lm2=ll+m2
            write(50,'(2I6," ",2G18.10)') m1,m2,vmatmt(lm1,ispn,lm2,jspn,ias)
          end do
        end do
      end do
    end do
  end do
end do
close(50)
! Slater parameters
open(50,file='FDU'//trim(filext),form='FORMATTED',action='WRITE')
do idu=1,ndftu
  is=isldu(1,idu)
  l=isldu(2,idu)
  write(50,*)
  write(50,'(2I4," : species, l")') is,l
  do k=0,2*l,2
    write(50,'(G18.10," : F^(",I1,")")') fdu(k,idu),k
  end do
  write(50,'(G18.10," : U")') ujdu(1,idu)
  write(50,'(G18.10," : J")') ujdu(2,idu)
  if (inpdftu.ge.4) write(50,'(G18.10," : screening length λ")') lamdu(idu)
end do
close(50)
end subroutine

