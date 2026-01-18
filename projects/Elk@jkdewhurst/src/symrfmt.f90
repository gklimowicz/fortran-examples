
! Copyright (C) 2018 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine symrfmt(nr,nri,np,ld,rfmt)
use modmain
implicit none
! arguments
integer, intent(in) :: nr(nspecies),nri(nspecies),np(nspecies)
integer, intent(in) :: ld
real(8), intent(inout) :: rfmt(ld,natmtot)
! local variables
integer is,ia,ja,ias,jas
integer n,isym,lspl
real(8) t0
! automatic arrays
logical done(natmmax)
! allocatable arrays
real(8), allocatable :: rfmt1(:,:),rfmt2(:)
allocate(rfmt1(ld,natmmax),rfmt2(ld))
t0=1.d0/dble(nsymcrys)
do is=1,nspecies
  n=np(is)
! make a copy of the input function
  do ia=1,natoms(is)
    ias=idxas(ia,is)
    rfmt1(1:n,ia)=rfmt(1:n,ias)
  end do
  done(:)=.false.
! loop over atoms
  do ia=1,natoms(is)
    if (done(ia)) cycle
    ias=idxas(ia,is)
    rfmt(1:n,ias)=0.d0
! loop over crystal symmetries
    do isym=1,nsymcrys
! index to spatial rotation lattice symmetry
      lspl=lsplsymc(isym)
! equivalent atom index (symmetry rotates atom ja into atom ia)
      ja=ieqatom(ia,is,isym)
! apply the rotation to the muffin-tin function
      call rotrfmt(symlatc(:,:,lspl),nr(is),nri(is),rfmt1(:,ja),rfmt2)
! accumulate in original function array
      rfmt(1:n,ias)=rfmt(1:n,ias)+rfmt2(1:n)
    end do
! normalise
    rfmt(1:n,ias)=t0*rfmt(1:n,ias)
    done(ia)=.true.
! rotate into equivalent atoms
    do isym=1,nsymcrys
      ja=ieqatom(ia,is,isym)
      if (done(ja)) cycle
      jas=idxas(ja,is)
! inverse symmetry (which rotates atom ia into atom ja)
      lspl=isymlat(lsplsymc(isym))
! rotate symmetrised function into equivalent muffin-tin
      call rotrfmt(symlatc(:,:,lspl),nr(is),nri(is),rfmt(:,ias),rfmt(:,jas))
      done(ja)=.true.
    end do
  end do
end do
deallocate(rfmt1,rfmt2)
end subroutine

