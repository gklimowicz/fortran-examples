!@sum diffreport reports on differences between arrays in 2 netcdf files
!@auth M. Kelley
      program diffreport
      implicit none
      include 'netcdf.inc'
      integer :: status,vtype,varid1,varid2,fid1,fid2,yesno_vid,nargs
      character(len=256) :: file1,file2
      character(len=40) :: vname,pos_str,yesno,attname,yesno_list
      integer :: iargc,nvars,ndims,dsizes(7),iatt,natts,n,cnt0,cntn
      integer :: arrsize1,arrsize2
      real*8, dimension(:), allocatable :: arr1,arr2,arrdiff,fd,fd1
      integer :: nmax_abs(1),nmax_rel(1)
      real*8 :: max_abs,max_rel
      logical, dimension(:), allocatable :: dothisvar
      logical :: print_usage

c
c parse command line
c
      nargs = iargc()
      if(nargs.ne.2 .and. nargs.ne.3) then
        write(6,*) 'usage: diffreport file1 file2 [optional: YesNoList]'
        write(6,*)
        write(6,*)
     &       '   Difference reports for selected variables can be'
        write(6,*)
     &       '   suppressed via the attributes of variable YesNoList.'
        write(6,*)
     &       '   Example: if file1 contains a variable'
        write(6,*)
        write(6,*)
     &       '      double IsImportant ;'
        write(6,*)
     &       '             IsImportant:x = "no" ;'
        write(6,*)
     &       '             IsImportant:y = "no" ;'
        write(6,*)
     &       '             IsImportant:z = "yes" ;'
        write(6,*)
        write(6,*)
     &       '   diffreport file1 file2 IsImportant'
        write(6,*)
        write(6,*)
     &       '   will skip fields x and y, but not z.'
        stop
      endif

c
c open input files
c
      call getarg(1,file1)
      call getarg(2,file2)
      call handle_err( nf_open(file1,nf_nowrite,fid1),
     &     'nonexistent/non-netcdf input file '//trim(file1) )
      call handle_err( nf_open(file2,nf_nowrite,fid2),
     &     'nonexistent/non-netcdf input file '//trim(file2) )

      status = nf_inq_nvars(fid1,nvars)
      allocate(dothisvar(nvars))
      dothisvar = .true.

c
c check for a list of names to skip
c
      if(nargs.eq.3) then
        call getarg(3,yesno_list)
        status = nf_inq_varid(fid1,trim(yesno_list),yesno_vid)
        if(status.eq.nf_noerr) then
          dothisvar(yesno_vid) = .false.
          status = nf_inq_varnatts(fid1,yesno_vid,natts)
          do iatt=1,natts
            attname = ''
            status = nf_inq_attname(fid1,yesno_vid,iatt,attname)
            yesno = ''
            status = nf_get_att_text(fid1,yesno_vid,trim(attname),yesno)
            if(nf_inq_varid(fid1,trim(attname),varid1).eq.nf_noerr)
     &           dothisvar(varid1) = yesno.eq.'yes'
          enddo
        else
          write(6,*) 'warning: nonexistent yesno_list '//
     &         trim(yesno_list)
        endif
      endif

c
c loop over arrays
c

      do varid1=1,nvars
        if(.not.dothisvar(varid1)) cycle
c skip character arrays
        status = nf_inq_vartype(fid1,varid1,vtype)
        if(vtype.eq.nf_char) cycle
c get array name
        status = nf_inq_varname(fid1,varid1,vname)
c skip certain modelE arrays
        if(trim(vname).eq.'cputime') cycle
c
        status = nf_inq_varid(fid2,vname,varid2)
c skip arrays not present in both files
        if(status.ne.nf_noerr) then
          write(6,*) 'file1 array '//trim(vname)//
     &         ' is absent in file2: skipping'
          cycle
        endif
c get array sizes and allocate space
        call get_varsize(fid1,vname,arrsize1)
        call get_varsize(fid2,vname,arrsize2)
c if array sizes do not match, skip
        if(arrsize1.ne.arrsize2) then
          write(6,*) 'array '//trim(vname)//
     &         ' has different sizes in file1/file2: skipping'
          cycle
        endif
        allocate(arr1(arrsize1),arr2(arrsize2))
c read arrays
        status = nf_get_var_double(fid1,varid1,arr1)
        status = nf_get_var_double(fid2,varid2,arr2)
        if(any(arr1.ne.arr1)) then
          write(6,*) 'NaNs in file1 array '//trim(vname)
        endif
        if(any(arr2.ne.arr2)) then
          write(6,*) 'NaNs in file2 array '//trim(vname)
        endif
c check for differences
        if(any(arr1.ne.arr2.and.
     &      .not.arr1.ne.arr1.and.
     &      .not.arr2.ne.arr2)) then
          write(6,*) trim(vname)//' max diffs:'
          call get_vdimsizes(fid1,vname,ndims,dsizes)
          allocate(arrdiff(arrsize1))
c absolute
          where(arr1.eq.arr1) arrdiff = abs(arr1-arr2)
          max_abs = maxval(arrdiff)
          nmax_abs = maxloc(arrdiff)
          call get_pos_str(nmax_abs(1),ndims,dsizes,pos_str)
          write(6,*) '          abs: ',max_abs,trim(pos_str)
c relative
          where(arr1.eq.arr1) arrdiff = arrdiff/max(abs(arr1),abs(arr2))
          max_rel = maxval(arrdiff)
          nmax_rel = maxloc(arrdiff)
          call get_pos_str(nmax_rel(1),ndims,dsizes,pos_str)
          write(6,*) '          rel: ',max_rel,trim(pos_str)
          deallocate(arrdiff)
c fraction of differences containing more than 8 bits of information
c          allocate(fd(arrsize1),fd1(arrsize1))
c          fd = fraction(abs(arr1-arr2))
c          cnt0 = count(fd.ne.0.)
c          do n=1,8
c            fd1 = fd-.5**n
c            where(fd1.ge.0.) fd = fd1
c            cntn = count(fd.ne.0.)
c            if(cntn.eq.0) exit
c          enddo
c          write(6,*) '      % diffs > 8 bits wide: ',
c     &         nint(100.*real(cntn,kind=8)/real(cnt0,kind=8))
c          deallocate(fd,fd1)
        endif

c deallocate arrays
        deallocate(arr1,arr2)
      enddo

c
c report arrays present in file2 but not in file1
c
      status = nf_inq_nvars(fid2,nvars)
      do varid2=1,nvars
c skip character arrays
        status = nf_inq_vartype(fid2,varid2,vtype)
        if(vtype.eq.nf_char) cycle
c get array name
        status = nf_inq_varname(fid2,varid2,vname)
c skip certain modelE arrays
        if(trim(vname).eq.'cputime') cycle
c
        status = nf_inq_varid(fid1,vname,varid1)
        if(status.ne.nf_noerr) then
          write(6,*) 'file2 array '//trim(vname)//
     &         ' is absent in file1: skipping'
          cycle
        endif
      enddo

c
c close input files
c
      status = nf_close(fid1)
      status = nf_close(fid2)

      deallocate(dothisvar)

      end program diffreport

      subroutine get_pos_str(n,ndims,dsizes,pos_str)
      implicit none
      integer :: n,ndims
      integer :: dsizes(ndims)
      character(len=40) :: pos_str
      integer :: inds(7)
      integer :: idim,nn,denom
      if(ndims.le.0) then
        pos_str=''
      else
        nn = n-1
        denom = product(dsizes)
        do idim=ndims,1,-1
          denom = denom/dsizes(idim)
          inds(idim) = nn/denom
          nn = nn-denom*inds(idim)
        enddo
        write(pos_str,'(7i6)') 1+inds(1:ndims)
        pos_str='at pos. '//trim(pos_str)
      endif
      return
      end subroutine get_pos_str
