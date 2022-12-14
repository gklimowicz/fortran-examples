      program write_2d_as_giss4d
!@sum Writes one or more 2D variables in the input file as
!@+   "giss4d" fortran records to the output file.  If dimension
!@+   names are specified on the command line, only the variables
!@+   possessing the two dimensions are written.
!@+   If a variable name is specified on the command line, only
!@+   that variable is written.
!@+   The title is created from attributes long_name/units.
!@+   If long_name is absent, the netcdf variable name is used.
!@auth M. Kelley
      implicit none
      real*4, dimension(:), allocatable :: coord1,coord2,vmean
      real*4, dimension(:,:), allocatable :: xout,shnhgm
      character(len=80) :: title,infile,outfile,lname
      character(len=40) :: vname,dimname1,dimname2,units,varname
      integer :: dsiz1,dsiz2,nargs,lunit,k,n
      character*16 :: cx,cy
      character*16, parameter :: cblank = '                '
      real*4, parameter :: undef=-1e30,one=1.
      include 'netcdf.inc'
      integer :: fid,status,varid,varid2,nvars,ndims,did1,did2,idim,jdim
      integer :: vid1,vid2
      integer, dimension(2) :: dids
      logical :: single_pair_of_dims

      nargs = iargc()
      if(nargs.lt.2 .or. nargs.gt.4) then
        write(6,*)
     &       'usage: write_2d_as_giss4d infile outfile '//
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
        allocate(coord1(dsiz1),coord2(dsiz2))
        allocate(xout(dsiz1,dsiz2),vmean(dsiz1+3),shnhgm(3,dsiz2))
        call get_var_real(fid,trim(dimname1),coord1)
        call get_var_real(fid,trim(dimname2),coord2)
        status = nf_inq_varid(fid,trim(dimname1),varid)
        cx=''
        cy=''
        if(nf_get_att_text(fid,varid,'giss_name',cx).ne.nf_noerr)
     &       cx = dimname1
        status = nf_inq_varid(fid,trim(dimname2),varid)
        if(nf_get_att_text(fid,varid,'giss_name',cy).ne.nf_noerr)
     &       cy = dimname2
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
        if(ndims.ne.2) cycle
        dids = -1
        status = nf_inq_vardimid(fid,varid,dids)
        vname =''
        status = nf_inq_varname(fid,varid,vname)
        if(index(vname,'_hemis').gt.0) cycle
        if(single_pair_of_dims) then
          if(count(dids.eq.did1).ne.1) cycle
          if(count(dids.eq.did2).ne.1) cycle
          if(dids(1).ne.did1) then
            write(6,*) 'dimensions for variable ',trim(vname),
     &           ' match but are reversed: skipping'
            cycle
          endif
        else
          did1 = dids(1)
          did2 = dids(2)
          status = nf_inq_dimlen(fid,dids(1),dsiz1)
          status = nf_inq_dimlen(fid,dids(2),dsiz2)
          allocate(coord1(dsiz1),coord2(dsiz2))
          allocate(xout(dsiz1,dsiz2),vmean(dsiz1+3),shnhgm(3,dsiz2))
          dimname1 = ''
          status = nf_inq_dimname(fid,did1,dimname1)
          dimname2 = ''
          status = nf_inq_dimname(fid,did2,dimname2)
          call get_var_real(fid,trim(dimname1),coord1)
          call get_var_real(fid,trim(dimname2),coord2)
          status = nf_inq_varid(fid,trim(dimname1),varid2)
          cx=''
          cy=''
          if(nf_get_att_text(fid,varid2,'giss_name',cx).ne.nf_noerr)
     &         cx = dimname1
          status = nf_inq_varid(fid,trim(dimname2),varid2)
          if(nf_get_att_text(fid,varid2,'giss_name',cy).ne.nf_noerr)
     &         cy = dimname2
        endif

        lname = vname
        status = nf_get_att_text(fid,varid,'long_name',lname)
        do k=1,len_trim(lname)  ! remove extra NULL characters
          if(iachar(lname(k:k)).eq.0) lname(k:k)=' '
        enddo
        units = ''
        status = nf_get_att_text(fid,varid,'units',units)
        if(status.eq.nf_noerr) then
          do k=1,len_trim(units) ! remove extra NULL characters
            if(iachar(units(k:k)).eq.0) units(k:k)=' '
          enddo
          units = ' ('//trim(units)//') '
        endif
        title = trim(lname)//units
        status = nf_get_var_real(fid,varid,xout)
c look for horizontal and vertical means
        shnhgm = undef; vmean = undef
        if(nf_inq_varid(fid,trim(vname)//'_hemis',varid2).eq.nf_noerr)
     &       status = nf_get_var_real(fid,varid2,shnhgm)
        if(nf_inq_varid(fid,trim(vname)//'_vmean',varid2).eq.nf_noerr)
     &       status = nf_get_var_real(fid,varid2,vmean)
c write the record
        write(lunit) title,
     &       dsiz1,dsiz2,1,1,       ! dimension sizes
     &       xout,                  ! the field
     &       coord1,coord2,one,one, ! coordinate axes
     &       cx,cy,cblank,cblank,   ! names of coordinate axes
     &       'NASAGISS',            !
     &       vmean,                 ! vertical mean
     &       shnhgm                 ! horizontal means
        if(.not.single_pair_of_dims) then
          deallocate(xout,coord1,coord2,vmean,shnhgm)
        endif
      enddo

c
c deallocate workspace
c
      if(single_pair_of_dims) then
        deallocate(xout,coord1,coord2,vmean,shnhgm)
      endif

c
c close fortran output file
c
      close(lunit)

      end program write_2d_as_giss4d
