      program write_giss2d
!@sum Writes one or more variables in the input file as
!@+   [title*80,real*4 data] fortran records to
!@+   the output file.  If dimension names are specified on
!@+   the command line, only the variables possessing the two
!@+   dimensions are written; the two dimensions need not
!@+   be the first two of a given variable, nor consecutive.
!@+   Variables with more than two dimensions are
!@+   written one slab at a time.
!@+   If a variable name is specified on the command line, only
!@+   that variable is written.
!@+   The title is created from attributes long_name/units.
!@+   If long_name is absent, the netcdf variable name is used.
!@auth M. Kelley
      implicit none
      real*4, dimension(:), allocatable :: xout
      character(len=80) :: title,infile,outfile,lname
      character(len=40) :: vname,dimname1,dimname2,units
      integer :: dsiz1,dsiz2,nargs,lunit,k,n,nslab

      include 'netcdf.inc'
      integer :: fid,status,varid,varid2,nvars,ndims,did1,did2,idim,jdim
      integer :: vid1,vid2
      integer, dimension(7) :: dids,srt,cnt,dsizes,kmod,p1,p2
      character(len=30), dimension(7) :: dnames
      character(len=30) :: diminfo,varname
      character(len=1) :: str1
      character(len=6) :: str3
      character(len=6) :: ifmt(7)='(ix.x)'
      real*4, parameter :: undef=-1e30
      real*4 :: shnhgm(3)
      character(len=30) :: run_info
      character(len=132) :: xlabel
      logical :: is_modelE_output,single_pair_of_dims
      integer :: ind,linfo,lrem,l1

      nargs = iargc()
      if(nargs.lt.2 .or. nargs.gt.4) then
        write(6,*)
     &  'usage: write_giss2d infile outfile '//
     &       '[ varname OR dimname1 dimname2 ]'
        stop
      endif

      single_pair_of_dims = nargs.eq.4

      call getarg(1,infile)
      call getarg(2,outfile)

      call handle_err(nf_open(infile,nf_nowrite,fid),
     &     'opening '//trim(infile))

      if(single_pair_of_dims) then
        call getarg(3,dimname1)
        call getarg(4,dimname2)
c
c Get info on the two requested dimensions, allocate workspace
c
        call get_dimsize(fid,trim(dimname1),dsiz1)
        call get_dimsize(fid,trim(dimname2),dsiz2)
        status = nf_inq_dimid(fid,trim(dimname1),did1)
        status = nf_inq_dimid(fid,trim(dimname2),did2)
        allocate(xout(dsiz1*dsiz2))
      endif

c
c get the number of quantities in the file
c
      status = nf_inq_nvars(fid,nvars)

      if(nargs.eq.3) then ! optional request of only 1 variable
        call getarg(3,varname)
        call handle_err(nf_inq_varid(fid,trim(varname),vid1),
     &       'variable '//trim(varname)//' does not exist')
        vid2 = vid1
      else ! all variables
        vid1 = 1
        vid2 = nvars
      endif

c
c check whether this is a modelE output file and try to extract
c the time period from the filename
c
      xlabel = ''
      status = nf_get_att_text(fid,nf_global,'xlabel',xlabel)
      is_modelE_output = .false.
      if(status.eq.nf_noerr) then
        run_info = xlabel(1:index(xlabel,' ')-1)
        ind = index(infile,trim(run_info))
        ind = index(infile(1:ind-1),'.')
        if(ind.gt.0) then
          is_modelE_output = .true.
          run_info = infile(1:ind-1)//' '//run_info
          linfo = len_trim(run_info)
        endif
      endif

c
c open fortran output file
c
      lunit = 10
      open(lunit,file=trim(outfile),form='unformatted',
     &     convert='big_endian')

c
c Loop over quantities in the file
c
      do varid=vid1,vid2
        status = nf_inq_varndims(fid,varid,ndims)
        if(ndims.lt.2) cycle
        dids = -1
        status = nf_inq_vardimid(fid,varid,dids)
        if(single_pair_of_dims) then
          if(count(dids.eq.did1).ne.1) cycle
          if(count(dids.eq.did2).ne.1) cycle
        else
          did1 = dids(1)
          did2 = dids(2)
          status = nf_inq_dimlen(fid,dids(1),dsiz1)
          status = nf_inq_dimlen(fid,dids(2),dsiz2)
          allocate(xout(dsiz1*dsiz2))
        endif
        do n=1,ndims
          if(dids(n).eq.did1) idim=n
          if(dids(n).eq.did2) jdim=n
          status = nf_inq_dimlen(fid,dids(n),dsizes(n))
          status = nf_inq_dimname(fid,dids(n),dnames(n))
        enddo
        srt = 1
        cnt = 1
        cnt(idim) = dsizes(idim)
        cnt(jdim) = dsizes(jdim)
        nslab = product(dsizes(1:ndims))/(dsiz1*dsiz2)
        if(ndims.gt.2) then
          k = 1
          diminfo=''
          do n=1,ndims
            if(n.eq.idim .or. n.eq.jdim) cycle
            kmod(n) = k
            k = k*dsizes(n)
            write(str1,'(i1)') int(1.+log10(real(dsizes(n))))
            ifmt(n)(3:3) = str1
            ifmt(n)(5:5) = str1
            write(str3,ifmt(n)) srt(n)
            if(len_trim(diminfo).eq.0) then
              diminfo=trim(dnames(n))//'='//trim(str3)
            else
              diminfo=
     &             trim(diminfo)//' '//trim(dnames(n))//'='//trim(str3)
            endif
            p1(n) = len_trim(diminfo)-len_trim(str3)+1
            p2(n) = len_trim(diminfo)
          enddo
        endif
        vname = ''
        status = nf_inq_varname(fid,varid,vname)
        lname = ''
        status = nf_get_att_text(fid,varid,'long_name',lname)
        if(status.eq.nf_noerr) then
          do k=1,len_trim(lname) ! remove extra NULL characters
            if(iachar(lname(k:k)).eq.0) lname(k:k)=' '
          enddo
        else
          lname = vname
        endif
        units = ''
        status = nf_get_att_text(fid,varid,'units',units)
        if(status.eq.nf_noerr) then
          do k=1,len_trim(units) ! remove extra NULL characters
            if(iachar(units(k:k)).eq.0) units(k:k)=' '
          enddo
          units = '('//trim(units)//')'
        endif
        shnhgm = undef
        if(ndims.eq.2) then ! look for global means
          status = nf_inq_varid(fid,trim(vname)//'_hemis',varid2)
          if(status.eq.nf_noerr) then
            status = nf_get_var_real(fid,varid2,shnhgm)
          endif
        endif
        do k=1,nslab
          status = nf_get_vara_real(fid,varid,srt,cnt,xout)
          title = trim(lname)//' '//units
          if(ndims.gt.2) title=trim(title)//' '//trim(diminfo)
          lrem = 80-(len_trim(title)+1) ! space remaining for run info
          if(lrem>0 .and. is_modelE_output) then
            l1 = 80-min(lrem,linfo)+1
            title(l1:80)=run_info(1:min(lrem,linfo))
          endif
          if(shnhgm(3).eq.undef) then ! no global mean available
            write(lunit) title,xout
          else                        ! write with global mean
            write(lunit) title,xout
     &           ,(undef,n=1,dsiz2)   ! have to write means at each lat
     &           ,shnhgm(3)           ! before the global mean
          endif
          if(ndims.gt.2) then
            do n=1,ndims        ! increment the start vector
              if(n.eq.idim .or. n.eq.jdim) cycle
              if(mod(k,kmod(n)).eq.0) then
                srt(n) = srt(n) + 1
                if(srt(n).gt.dsizes(n)) srt(n)=1
                write(str3,ifmt(n)) srt(n)
                diminfo(p1(n):p2(n))=trim(str3)
              endif
            enddo
          endif
        enddo
        if(.not.single_pair_of_dims) deallocate(xout)
      enddo

c
c deallocate workspace
c
      if(single_pair_of_dims) deallocate(xout)

c
c close fortran output file
c
      close(lunit)

      end program write_giss2d
