
! Copyright (C) 2020 J. K. Dewhurst and S. Sharma.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module moddelf

use modramdisk
use modmpi

contains

subroutine delfiles(evec,devec,eval,occ,pmat,epsi)
use modphonon
implicit none
! arguments
logical, optional, intent(in) :: evec,devec,eval,occ,pmat,epsi
! local variables
integer ios
character(256) fext,fname
if (present(evec)) then
! delete the first-variational eigenvector file
  fname=trim(scrpath)//'EVECFV'//trim(filext)
  if (mp_mpi) then
    open(202,file=fname,form='UNFORMATTED',iostat=ios)
    close(202,status='DELETE',iostat=ios)
  end if
  if (ramdisk) call delfrd(fname)
! delete the second-variational eigenvector file
  fname=trim(scrpath)//'EVECSV'//trim(filext)
  if (mp_mpi) then
    open(206,file=fname,form='UNFORMATTED',iostat=ios)
    close(206,status='DELETE',iostat=ios)
  end if
  if (ramdisk) call delfrd(fname)
end if
if (present(devec)) then
! construct the dynamical matrix file extension
  call dynfext(iqph,isph,iaph,ipph,fext)
! delete the eigenvector derivative files
  fname=trim(scrpath)//'DEVECFV'//trim(fext)
  if (mp_mpi) then
    open(222,file=fname,form='UNFORMATTED',iostat=ios)
    close(222,status='DELETE',iostat=ios)
  end if
  if (ramdisk) call delfrd(fname)
  fname=trim(scrpath)//'DEVECSV'//trim(fext)
  if (mp_mpi) then
    open(226,file=fname,form='UNFORMATTED',iostat=ios)
    close(226,status='DELETE',iostat=ios)
  end if
  if (ramdisk) call delfrd(fname)
end if
if (present(eval)) then
! delete the first-variational eigenvalue file
  fname='EVALFV'//trim(filext)
  if (mp_mpi) then
    open(200,file=fname,form='UNFORMATTED',iostat=ios)
    close(200,status='DELETE',iostat=ios)
  end if
  if (ramdisk) call delfrd(fname)
! delete the second-variational eigenvalue file
  if (mp_mpi) then
    open(204,file='EVALSV'//trim(filext),form='UNFORMATTED',iostat=ios)
    close(204,status='DELETE',iostat=ios)
  end if
end if
if (present(occ)) then
! delete the occupation number file
  if (mp_mpi) then
    open(208,file='OCCSV'//trim(filext),form='UNFORMATTED',iostat=ios)
    close(208,status='DELETE',iostat=ios)
  end if
end if
if (present(pmat)) then
! delete the momentum matrix elements file
  if (mp_mpi) then
    open(230,file='PMAT.OUT',form='UNFORMATTED',iostat=ios)
    close(230,status='DELETE',iostat=ios)
  end if
  if (ramdisk) call delfrd('PMAT.OUT')
end if
if (present(epsi)) then
! delete the inverse epsilon file
  if (mp_mpi) then
    open(245,file='EPSINV.OUT',form='UNFORMATTED',iostat=ios)
    close(245,status='DELETE',iostat=ios)
  end if
end if
end subroutine

end module

