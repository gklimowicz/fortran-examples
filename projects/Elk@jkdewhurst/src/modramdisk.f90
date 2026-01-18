
! Copyright (C) 2021 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modramdisk

! ramdisk is .true. if the RAM disk should be used
logical ramdisk
! maximum allowed number of files on the RAM disk
integer, parameter, private :: maxfiles=16
! maximum number of records per file
integer, parameter, private :: maxrec=32768
! wrtdsk is .true. if files should also be written to disk
logical wrtdsk,wrtdsk0

! record data stored as 4-byte words
type, private :: rec_t
  integer(4), allocatable :: dat(:)
end type

! RAM disk file consisting of the filename and an array of records
type, private :: file_t
  character(256) fname
  type(rec_t), allocatable :: rec(:)
end type

! arrays of files constituting the RAM disk
type(file_t), allocatable, private :: file(:)

! private subroutines
private findfile,openfile

contains

subroutine initrd
! this subroutine should not be called from a parallel region
implicit none
! local variables
integer i
if (allocated(file)) then
  write(*,*)
  write(*,'("Error(initrd): RAM disk already initialised")')
  write(*,*)
  stop
end if
allocate(file(maxfiles))
do i=1,maxfiles
  file(i)%fname=''
end do
end subroutine

subroutine findfile(fname,fnum)
implicit none
! arguments
character(*), intent(in) :: fname
integer, intent(out) :: fnum
! local variables
integer i
if (.not.allocated(file)) then
  write(*,*)
  write(*,'("Error(findfile): RAM disk not initialised")')
  write(*,*)
  stop
end if
fnum=0
do i=1,maxfiles
  if (file(i)%fname.eq.fname) then
    fnum=i
    return
  end if
end do
end subroutine

subroutine openfile(fname,fnum)
implicit none
! arguments
character(*), intent(in) :: fname
integer, intent(out) :: fnum
! local variables
integer i
! check to see if the file already exists
call findfile(fname,fnum)
if (fnum.ne.0) return
! use the first available file number
do i=1,maxfiles
  if (file(i)%fname.eq.'') then
! assign the filename
    file(i)%fname=fname
! allocate the record array
    allocate(file(i)%rec(maxrec))
    fnum=i
    return
  end if
end do
write(*,*)
write(*,'("Error(openfile): too many RAM disk files open : ",I8)') maxfiles
write(*,*)
stop
end subroutine

subroutine delfrd(fname)
! this subroutine should not be called from a parallel region
implicit none
! arguments
character(*), intent(in) :: fname
! local variables
integer fnum,i
! determine the file number
call findfile(fname,fnum)
! return if the file does not exist
if (fnum.eq.0) return
! erase filename
file(fnum)%fname=''
! deallocate associated arrays
do i=1,maxrec
  if (allocated(file(fnum)%rec(i)%dat)) then
    deallocate(file(fnum)%rec(i)%dat)
  end if
end do
deallocate(file(fnum)%rec)
end subroutine

subroutine putrd(fname,irec,n1,n2,n3,v1,v2,nrv,rva,nzv,zva)
! this subroutine should be called from an OpenMP critical section
implicit none
! arguments
character(*), intent(in) :: fname
integer, intent(in) :: irec
integer, optional, intent(in) :: n1,n2,n3
real(8), optional, intent(in) :: v1(3),v2(3)
integer, optional, intent(in) :: nrv
real(8), optional, intent(in) :: rva(*)
integer, optional, intent(in) :: nzv
complex(8), optional, intent(in) :: zva(*)
! local variables
integer fnum,n,i,j
! check that the record number does not exceed the maximum
if (irec.gt.maxrec) then
! fail safe by enabling disk writes
  wrtdsk=.true.
  return
end if
! find the record length in 4-byte words
n=0
if (present(n1)) n=n+1
if (present(n2)) n=n+1
if (present(n3)) n=n+1
if (present(v1)) n=n+6
if (present(v2)) n=n+6
if (present(rva)) then
  if (present(nrv)) then
    n=n+2*nrv
  else
    write(*,*)
    write(*,'("Error(putrd): missing argument nrv")')
    write(*,*)
    stop
  end if
end if
if (present(zva)) then
  if (present(nzv)) then
    n=n+4*nzv
  else
    write(*,*)
    write(*,'("Error(putrd): missing argument nzv")')
    write(*,*)
    stop
  end if
end if
! open the file
call openfile(fname,fnum)
! allocate the record data array if required
if (allocated(file(fnum)%rec(irec)%dat)) then
  if (size(file(fnum)%rec(irec)%dat).lt.n) then
    deallocate(file(fnum)%rec(irec)%dat)
  end if
end if
if (.not.allocated(file(fnum)%rec(irec)%dat)) then
  allocate(file(fnum)%rec(irec)%dat(n))
end if
i=1
if (present(n1)) then
  file(fnum)%rec(irec)%dat(i)=n1
  i=i+1
end if
if (present(n2)) then
  file(fnum)%rec(irec)%dat(i)=n2
  i=i+1
end if
if (present(n3)) then
  file(fnum)%rec(irec)%dat(i)=n3
  i=i+1
end if
if (present(v1)) then
  j=i+5
  file(fnum)%rec(irec)%dat(i:j)=transfer(v1(:),file(fnum)%rec(irec)%dat(i:j))
  i=i+6
end if
if (present(v2)) then
  j=i+5
  file(fnum)%rec(irec)%dat(i:j)=transfer(v2(:),file(fnum)%rec(irec)%dat(i:j))
  i=i+6
end if
if (present(rva)) then
  j=i+2*nrv-1
  file(fnum)%rec(irec)%dat(i:j)=transfer(rva(1:nrv), &
   file(fnum)%rec(irec)%dat(i:j))
  i=i+2*nrv
end if
if (present(zva)) then
  j=i+4*nzv-1
  file(fnum)%rec(irec)%dat(i:j)=transfer(zva(1:nzv), &
   file(fnum)%rec(irec)%dat(i:j))
end if
end subroutine

subroutine getrd(fname,irec,tgs,n1,n2,n3,v1,v2,nrv,rva,nzv,zva)
implicit none
! arguments
character(*), intent(in) :: fname
integer, intent(in) :: irec
logical, intent(out) :: tgs
integer, optional, intent(out) :: n1,n2,n3
real(8), optional, intent(out) :: v1(3),v2(3)
integer, optional, intent(in) :: nrv
real(8), optional, intent(out) :: rva(*)
integer, optional, intent(in) :: nzv
complex(8), optional, intent(out) :: zva(*)
! local variables
integer fnum,n,i
if (present(rva)) then
  if (.not.present(nrv)) then
    write(*,*)
    write(*,'("Error(getrd): missing argument nrv")')
    write(*,*)
    stop
  end if
end if
if (present(zva)) then
  if (.not.present(nzv)) then
    write(*,*)
    write(*,'("Error(getrd): missing argument nzv")')
    write(*,*)
    stop
  end if
end if
tgs=.false.
if (irec.gt.maxrec) return
! determine the file number
call findfile(fname,fnum)
! return unsuccessfully if file is not found or record is unavailable
if (fnum.eq.0) return
if (.not.allocated(file(fnum)%rec(irec)%dat)) return
n=size(file(fnum)%rec(irec)%dat)
i=1
if (present(n1)) then
  if (n.lt.1) return
  n1=file(fnum)%rec(irec)%dat(i)
  i=i+1
  n=n-1
end if
if (present(n2)) then
  if (n.lt.1) return
  n2=file(fnum)%rec(irec)%dat(i)
  i=i+1
  n=n-1
end if
if (present(n3)) then
  if (n.lt.1) return
  n3=file(fnum)%rec(irec)%dat(i)
  i=i+1
  n=n-1
end if
if (present(v1)) then
  if (n.lt.6) return
  v1(:)=transfer(file(fnum)%rec(irec)%dat(i:i+5),v1(:))
  i=i+6
  n=n-6
end if
if (present(v2)) then
  if (n.lt.6) return
  v2(:)=transfer(file(fnum)%rec(irec)%dat(i:i+5),v2(:))
  i=i+6
  n=n-6
end if
if (present(rva)) then
  if (n.lt.2*nrv) return
  rva(1:nrv)=transfer(file(fnum)%rec(irec)%dat(i:i+2*nrv-1),rva(1:nrv))
  i=i+2*nrv
  n=n-2*nrv
end if
if (present(zva)) then
  if (n.lt.4*nzv) return
  zva(1:nzv)=transfer(file(fnum)%rec(irec)%dat(i:i+4*nzv-1),zva(1:nzv))
end if
! flag the get operation as successful
tgs=.true.
end subroutine

subroutine rdstatus
! this subroutine should not be called from a parallel region
implicit none
! local variables
integer nf,nr,i,j
integer(8) m,n
write(*,*)
write(*,'("Info(rdstatus):")')
if (.not.allocated(file)) then
  write(*,'(" RAM disk not initialised")')
  return
end if
nf=0
n=0
do i=1,maxfiles
  if (file(i)%fname.ne.'') then
    write(*,*)
    write(*,'(" Filename : ",A)') trim(file(i)%fname)
    nf=nf+1
    nr=0
    m=0
    do j=1,maxrec
      if (allocated(file(i)%rec(j)%dat)) then
        nr=nr+1
        m=m+size(file(i)%rec(j)%dat)
      end if
    end do
    n=n+m
    write(*,'("  number of records : ",I8)') nr
    write(*,'("  total number of bytes : ",I14)') 4*m
  end if
end do
write(*,*)
write(*,'(" Number of files on RAM disk : ",I4)') nf
write(*,'(" Total number of bytes used by RAM disk : ",I14)') 4*n
end subroutine

end module

