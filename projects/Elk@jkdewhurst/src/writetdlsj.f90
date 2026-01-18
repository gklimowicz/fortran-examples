
! Copyright (C) 2019 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

subroutine writetdlsj
use modmain
use modtddft
use modmpi
implicit none
! local variables
integer is,ia,ias
real(8) xl(3),xs(3)
character(256) fext
! allocatable arrays
complex(8), allocatable :: dmat(:,:,:,:,:)
if (ntswrite(1).le.0) return
if (mod(itimes-1,ntswrite(1)).ne.0) return
if ((itimes.gt.1).and.(itimes.lt.ntswrite(2))) return
allocate(dmat(lmmaxo,nspinor,lmmaxo,nspinor,natmtot))
! generate the density matrix in each muffin-tin
call gendmat(.false.,.false.,0,lmaxo,lmmaxo,dmat)
if (mp_mpi) then
  write(fext,'("_TS",I8.8,".OUT")') itimes
  open(50,file='LSJ'//trim(fext),form='FORMATTED')
  do is=1,nspecies
    write(50,*)
    write(50,'("Species : ",I4," (",A,")")') is,trim(spsymb(is))
    do ia=1,natoms(is)
      ias=idxas(ia,is)
! calculate the expectation value of L and S
      call dmatls(dmat(:,:,:,:,ias),xl,xs)
      write(50,'(" atom : ",I4)') ia
      write(50,'("  L : ",3G18.10)') xl(:)
      write(50,'("  S : ",3G18.10)') xs(:)
      write(50,'("  J : ",3G18.10)') xl(:)+xs(:)
    end do
  end do
  close(50)
end if
end subroutine

