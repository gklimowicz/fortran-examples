      module cdl_mod
      implicit none
      private

      integer, parameter :: cdl_strlen=100
      integer, parameter :: ndims_max=100
      character(len=4) :: indent='    '
      character(len=8) :: indent2='        '

      public :: cdl_type
      type cdl_type
        private
        character(len=30) :: name ! name of this cdl for error reporting
        integer :: nlines     ! current number of lines
        integer :: ndims,ncoordlines,nvarlines,ndatalines
        character(len=cdl_strlen), dimension(ndims_max) :: dims
        character(len=cdl_strlen), dimension(3*ndims_max) :: coords
        character(len=cdl_strlen), dimension(:), allocatable ::
     &       vars,datavalues,text
      end type cdl_type

      public :: cdl_strlen,init_cdl_type,print_cdl,add_dim,add_coord,
     &     add_var,add_varline,add_dataline,add_vardata,add_unlimdim,
     &     copy_dims,copy_coord_vars,
     &     copy_data,merge_cdl,assemble_cdl,defvar_cdl,write_cdl
      

      interface add_vardata
        module procedure add_vardata_r8_1d
        module procedure add_vardata_int_1d
        module procedure add_vardata_1d_array_of_strings
      end interface add_vardata

      contains

      subroutine init_cdl_type(cdl_name,cdl)
      character(len=*), intent(in) :: cdl_name
      type(cdl_type), intent(inout) :: cdl
      cdl%name = cdl_name
      cdl%ndims = 0
      cdl%ncoordlines = 0
      cdl%nvarlines = 0
      cdl%ndatalines = 0
      return
      end subroutine init_cdl_type

      subroutine add_dim(cdl,dimname,dimsize)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: dimname
      integer, intent(in) :: dimsize
      character(len=10) :: dimstr
      integer :: k
      call checkname(trim(dimname))
      write(dimstr,'(i10)') dimsize
      k = cdl%ndims + 1
      cdl%dims(k) = indent//trim(dimname)//' = '//
     &     trim(adjustl(dimstr))//' ;'
      cdl%ndims  = k
      return
      end subroutine add_dim

      subroutine add_unlimdim(cdl,dimname)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: dimname
      integer :: k
      call checkname(trim(dimname))
      k = cdl%ndims + 1
      cdl%dims(k) = indent//trim(dimname)//
     &     ' = UNLIMITED ; // (0 currently)'
      cdl%ndims  = k
      return
      end subroutine add_unlimdim

      subroutine copy_dims(cdl_in,cdl_out)
      type(cdl_type), intent(in) :: cdl_in
      type(cdl_type), intent(inout) :: cdl_out
      integer :: k,kout
      kout = cdl_out%ndims
      do k=1,cdl_in%ndims
        kout = kout + 1
        cdl_out%dims(kout) = cdl_in%dims(k)
      enddo
      cdl_out%ndims = kout
      return
      end subroutine copy_dims

      subroutine add_coord(cdl,coordname,coordsize,
     &     units,long_name,coordvalues)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: coordname
      character(len=*), intent(in), optional :: units,long_name
      integer, intent(in) :: coordsize
      real*8, dimension(coordsize), intent(in), optional :: coordvalues
      integer :: k
      call add_dim(cdl,coordname,coordsize)
      k = cdl%ncoordlines + 1
      cdl%coords(k) = indent//'float '//trim(coordname)//'('//
     &     trim(coordname)//') ;'
      if(present(units)) then
        k = k + 1
        cdl%coords(k) = indent2//trim(coordname)//':units = "'//
     &     trim(units)//'" ;'
      endif
      if(present(long_name)) then
        k = k + 1
        cdl%coords(k) = indent2//trim(coordname)//':long_name = "'//
     &     trim(long_name)//'" ;'
      endif
      cdl%ncoordlines  = k
      if(present(coordvalues)) then
        call add_vardata(cdl,trim(coordname),coordvalues)
      endif
      return
      end subroutine add_coord

      subroutine copy_coord_vars(cdl_in,cdl_out)
      type(cdl_type), intent(in) :: cdl_in
      type(cdl_type), intent(inout) :: cdl_out
      integer :: k,kout
      kout = cdl_out%ncoordlines
      do k=1,cdl_in%ncoordlines
        kout = kout + 1
        cdl_out%coords(kout) = cdl_in%coords(k)
      enddo
      cdl_out%ncoordlines = kout
      return
      end subroutine copy_coord_vars

      subroutine copy_data(cdl_in,cdl_out)
      type(cdl_type), intent(in) :: cdl_in
      type(cdl_type), intent(inout) :: cdl_out
      integer :: k,kout
      call alloc_datavalues(cdl_out,
     &     cdl_out%ndatalines+cdl_in%ndatalines)
      kout = cdl_out%ndatalines
      do k=1,cdl_in%ndatalines
        kout = kout + 1
        cdl_out%datavalues(kout) = cdl_in%datavalues(k)
      enddo
      cdl_out%ndatalines = kout
      return
      end subroutine copy_data

      subroutine add_var(cdl,varstr,units,long_name,auxvar_string,
     &     set_miss,make_timeaxis)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: varstr
      character(len=*), intent(in), optional ::
     &     units,long_name,auxvar_string
      logical, intent(in), optional :: set_miss,make_timeaxis
      character(len=cdl_strlen) :: varname,tmpstr
      character(len=1024) :: thisline
      integer :: k,n1,n2,l,lmax
      logical :: make_timeaxis_
      make_timeaxis_ = .false.
      if(present(make_timeaxis)) make_timeaxis_ = make_timeaxis
      tmpstr = adjustl(varstr)
      call get_varname(tmpstr,varname)
      if(make_timeaxis_) call formvar_time(tmpstr)
      k = cdl%nvarlines + 1
      call alloc_vars(cdl,k)
      cdl%vars(k) = indent//tmpstr
      if(present(units)) then
        if(len_trim(units).gt.0) then
          k = k + 1
          cdl%vars(k) = indent2//trim(varname)//':units = "'//
     &         trim(units)//'" ;'
        endif
      endif
      if(present(long_name)) then
        thisline = indent2//trim(varname)//':long_name = "'//
     &     trim(long_name)//'" ;'
        ! Break this line into parts if necessary.  No
        ! attempt yet made to break at whitespaces or punctuation.
        ! Should put the line-breaking into a routine.
        lmax = len(cdl%vars(1))
        l = len_trim(thisline)
        do while(l.gt.0)
          k = k + 1
          cdl%vars(k) = thisline(1:min(l,lmax))
          if(l.gt.lmax) thisline = thisline(lmax+1:l)
          l = l - lmax
        enddo
      endif
      if(present(set_miss)) then
        if(set_miss) then
          k = k + 1
          cdl%vars(k) = indent2//trim(varname)//
     &         ':missing_value = -1.e30f ;'
        endif
      endif
      if(present(auxvar_string)) then
        tmpstr = adjustl(auxvar_string)
        call get_varname(tmpstr,varname)
        if(make_timeaxis_) call formvar_time(tmpstr)
        k = k + 1
        cdl%vars(k) = indent//tmpstr
      endif
      cdl%nvarlines  = k
      return
      end subroutine add_var

      subroutine add_varline(cdl,varstr)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: varstr
      character(len=cdl_strlen) :: varname,attname,tmpstr
      integer :: k,nc,ne,lt
      tmpstr = adjustl(varstr)
      lt = len_trim(tmpstr)
      nc = index(tmpstr,':')
      ne = index(tmpstr,'=')
      if(nc.eq.0 .or. nc.eq.ne-1 .or. ne.eq.0 .or. ne.eq.lt) then
        write(6,*) 'add_varline: syntax error: ',trim(varstr)
        call stop_model('add_varline: syntax error',255)
      else
        if(nc.gt.1) then
          varname=tmpstr(1:nc-1)
          call checkname(varname)
        endif
        attname=adjustl(tmpstr(nc+1:ne-1))
        call checkname(attname)
      endif
      k = cdl%nvarlines + 1
      call alloc_vars(cdl,k)
      cdl%vars(k) = indent2//trim(varstr)
      cdl%nvarlines  = k
      return
      end subroutine add_varline

      subroutine get_varname(varstr,varname)
      character(len=cdl_strlen) :: varstr,varname
      integer :: n1,n2
      n1 = index(varstr,' ')+1
      n2 = index(varstr,'(')-1
      if(n2.gt.0) then
        varname = adjustl(varstr(n1:n2))
      else
        n2 = index(varstr,';')-1
        varname = adjustl(varstr(n1:n2))
      endif
      call checkname(trim(varname))
      call checktype(varstr(1:n1-1),trim(varname))
      end subroutine get_varname

      subroutine formvar_time(varstr)
      character(len=cdl_strlen) :: varstr
      integer :: n1,n2
      n2 = index(varstr,';')
      n1 = index(varstr,'(')
      if(n1.gt.0) then
        varstr = varstr(1:n1-1)//'(time,'//varstr(n1+1:n2)
      else
        varstr = trim(varstr(1:n2-1))//'(time);'
      endif
      return
      end subroutine formvar_time

      subroutine checkname(namestr)
c Two rules:
c All characters must be A-Z|a-z|0-9|_|+|-|.|@
c First character must be A-Z|a-z
      character(len=*) :: namestr
      integer :: i,ascii,r
      logical :: is_bad
      integer, parameter :: num_ranges=7
      integer, dimension(2,num_ranges) :: bds=reshape((/
c          A-Z    a-z     0-9    +      -   .   @      _
     &     65,90, 97,122, 48,57, 43,43, 45,46, 64,64, 95,95
     &     /),(/2,num_ranges/))
      do i=1,len_trim(namestr)
        ascii = iachar(namestr(i:i))
        is_bad = .true.
        do r=1,num_ranges
          if(ascii.ge.bds(1,r) .and. ascii.le.bds(2,r)) then
            is_bad = .false.
            exit
          endif
          if(i.eq.1 .and. r.eq.2) exit ! 1st char A-Z|a-z
        enddo
        if(is_bad) exit
      enddo
      if(is_bad) then
        write(6,*) 'bad character at position ',i,
     &       'in name ',trim(namestr)
        call stop_model('invalid netcdf name',255)
      endif
      end subroutine checkname

      subroutine checktype(typestr,varname)
      character(len=*) :: typestr,varname
      select case(trim(typestr))
      case ('float','double','int','char','short')
      case default
        write(6,*) 'bad type '//trim(typestr)//
     &       ' for variable '//trim(varname)
        call stop_model('invalid netcdf type',255)
      end select
      end subroutine checktype

      subroutine add_dataline(cdl,varstr)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: varstr
      integer :: k
      k = cdl%ndatalines + 1
      call alloc_datavalues(cdl,k)
      cdl%datavalues(k) = indent//trim(varstr)
      cdl%ndatalines  = k
      return
      end subroutine add_dataline

      subroutine add_vardata_r8_1d(cdl,varname,values,fmtstr)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: varname
      real*8, dimension(:), intent(in) :: values
      character(len=*), intent(in), optional :: fmtstr
      integer :: k,line,i1,i2,pos,varsize,nl
      character(len=32) :: fmtstr_
      integer, parameter :: npl=6 ! 6 values per line
      if(present(fmtstr)) then
        fmtstr_ = '(6('//trim(fmtstr)//',","))'
      else
        fmtstr_ = '(6(1pe13.5,","))'
      endif
      varsize = size(values)
      k = cdl%ndatalines + 1
      nl = (varsize+npl-1)/npl
      call alloc_datavalues(cdl,k+nl)
      cdl%datavalues(k) = indent//trim(varname)//' = '
      do line=1,nl
        i1 = 1 + npl*(line-1)
        i2 = min(varsize,i1+npl-1)
        k = k + 1
        cdl%datavalues(k) = ''
        write(cdl%datavalues(k),fmtstr_) values(i1:i2)
        if(i2.eq.varsize) then
          pos = len_trim(cdl%datavalues(k))
          cdl%datavalues(k)(pos:pos) = ';'
        endif
      enddo
      cdl%ndatalines  = k
      return
      end subroutine add_vardata_r8_1d

      subroutine add_vardata_int_1d(cdl,varname,values)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: varname
      integer, dimension(:), intent(in) :: values
      integer :: k,line,i1,i2,pos,varsize,nl
      integer, parameter :: npl=8 ! 8 values per line
      varsize = size(values)
      k = cdl%ndatalines + 1
      nl = (varsize+npl-1)/npl
      call alloc_datavalues(cdl,k+nl)
      cdl%datavalues(k) = indent//trim(varname)//' = '
      do line=1,nl
        i1 = 1 + npl*(line-1)
        i2 = min(varsize,i1+npl-1)
        k = k + 1
        cdl%datavalues(k) = ''
        write(cdl%datavalues(k),'(8(i10,","))') values(i1:i2)
        if(i2.eq.varsize) then
          pos = len_trim(cdl%datavalues(k))
          cdl%datavalues(k)(pos:pos) = ';'
        endif
      enddo
      cdl%ndatalines  = k
      return
      end subroutine add_vardata_int_1d

      subroutine add_vardata_1d_array_of_strings(cdl,varname,strings)
      type(cdl_type), intent(inout) :: cdl
      character(len=*), intent(in) :: varname
      character(len=*), dimension(:), intent(in) :: strings(:)
      integer :: k,line,varsize
      character*1 punct
      varsize = size(strings)
      k = cdl%ndatalines + 1
      call alloc_datavalues(cdl,k+varsize)
      cdl%datavalues(k) = indent//trim(varname)//' = '
      punct = ','
      do line=1,varsize
        k = k + 1
        if(line.eq.varsize) punct=';'
        cdl%datavalues(k) = ''
        write(cdl%datavalues(k),'(a)')
     &       '"'//strings(line)//'"'//punct
      enddo
      cdl%ndatalines  = k
      return
      end subroutine add_vardata_1d_array_of_strings

      subroutine assemble_cdl(cdl)
      type(cdl_type), intent(inout) :: cdl
      integer :: k,n
      logical :: needs_time
      needs_time = .false.
      do n=1,cdl%nvarlines
        if(index(cdl%vars(n),'(time').gt.0) then
          needs_time = .true.
          exit
        endif
      enddo
      if(allocated(cdl%text)) deallocate(cdl%text)
      n = cdl%ndims + cdl%ncoordlines
     &  + cdl%nvarlines + cdl%ndatalines + 1000
      allocate(cdl%text(n))
      cdl%text = ''
      cdl%text(1) = 'netcdf xxx { '
      cdl%text(2) = 'dimensions:  '
      k = 2
      do n=1,cdl%ndims
        k = k + 1
        cdl%text(k) = cdl%dims(n)
        if(index(cdl%dims(n),'time =').gt.0) then
          needs_time = .false.
        endif
      enddo
      if(needs_time) then
        k = k + 1
        cdl%text(k) = indent//'time = UNLIMITED;'
      endif
      k = k + 1
      cdl%text(k) = 'variables:  '
      if(needs_time) then
        k = k + 1
        cdl%text(k) = indent//'double time(time);'
        k = k + 1
        cdl%text(k) = indent2//'time:units = "'//
     &       'years since 0000-01-01 00:00 UTC" ;'
      endif
      do n=1,cdl%ncoordlines
        k = k + 1
        cdl%text(k) = cdl%coords(n)
      enddo
      do n=1,cdl%nvarlines
        k = k + 1
        cdl%text(k) = cdl%vars(n)
      enddo
      if(needs_time .or. cdl%ndatalines.gt.0) then
        k = k + 1
        cdl%text(k) = 'data:  '
      endif
      if(needs_time) then
        k = k + 1
        cdl%text(k) = 'time = 0;'
      endif
      if(cdl%ndatalines.gt.0) then
        do n=1,cdl%ndatalines
          k = k + 1
          cdl%text(k) = cdl%datavalues(n)
        enddo
      endif
      k = k + 1
      cdl%text(k) = '}'
      cdl%nlines = k
      return
      end subroutine assemble_cdl

      subroutine merge_cdl(cdl1,cdl2,cdl)
      type(cdl_type), intent(in) :: cdl1,cdl2
      type(cdl_type), intent(inout) :: cdl
      integer :: k,n
c dims
      k = 0
      do n=1,cdl1%ndims
        k = k + 1
        cdl%dims(k) = cdl1%dims(n)
      enddo
      do n=1,cdl2%ndims
        k = k + 1
        cdl%dims(k) = cdl2%dims(n)
      enddo
      cdl%ndims = k
c coords
      k = 0
      do n=1,cdl1%ncoordlines
        k = k + 1
        cdl%coords(k) = cdl1%coords(n)
      enddo
      do n=1,cdl2%ncoordlines
        k = k + 1
        cdl%coords(k) = cdl2%coords(n)
      enddo
      cdl%ncoordlines = k
c vars
      call alloc_vars(cdl,cdl1%nvarlines+cdl2%nvarlines)
      k = 0
      do n=1,cdl1%nvarlines
        k = k + 1
        cdl%vars(k) = cdl1%vars(n)
      enddo
      do n=1,cdl2%nvarlines
        k = k + 1
        cdl%vars(k) = cdl2%vars(n)
      enddo
      cdl%nvarlines = k
c data
      call alloc_datavalues(cdl,cdl1%ndatalines+cdl2%ndatalines)
      k = 0
      do n=1,cdl1%ndatalines
        k = k + 1
        cdl%datavalues(k) = cdl1%datavalues(n)
      enddo
      do n=1,cdl2%ndatalines
        k = k + 1
        cdl%datavalues(k) = cdl2%datavalues(n)
      enddo
      cdl%ndatalines = k
      return
      end subroutine merge_cdl

      subroutine print_cdl(cdl)
      type(cdl_type), intent(in) :: cdl
      integer :: k
      do k=1,cdl%nlines
        write(6,*) cdl%text(k)
      enddo
      end subroutine print_cdl

      subroutine defvar_cdl(grid,fid,cdl,varstr)
      use dd2d_utils, only : dist_grid
      use pario, only : defvar
      type(dist_grid) :: grid
      integer :: fid
      type(cdl_type) :: cdl
      character(len=*) :: varstr
      call assemble_cdl(cdl)
      call defvar(grid,fid,cdl%text(1:cdl%nlines),trim(varstr))
      return
      end subroutine defvar_cdl

      subroutine write_cdl(grid,fid,varstr,cdl)
      use dd2d_utils, only : dist_grid
      use pario, only : write_data
      type(dist_grid) :: grid
      integer :: fid
      type(cdl_type) :: cdl
      character(len=*) :: varstr
      call assemble_cdl(cdl)
      call write_data(grid,fid,trim(varstr),cdl%text(1:cdl%nlines))
      return
      end subroutine write_cdl

      subroutine alloc_vars(cdl,k)
      type(cdl_type), intent(inout) :: cdl
      integer, intent(in) :: k
      integer :: n
      character(len=cdl_strlen), dimension(:), allocatable :: tmp_
      if(.not.allocated(cdl%vars)) then
        allocate(cdl%vars(max(k,1000)))
      else
        n = size(cdl%vars)
        if(k+100.gt.n) then
          allocate(tmp_(n))
          tmp_(1:n) = cdl%vars(1:n)
          deallocate(cdl%vars)
          allocate(cdl%vars(max(k+100,2*n)))
          cdl%vars(1:n) = tmp_(1:n)
          deallocate(tmp_)
        endif
      endif
      end subroutine alloc_vars

      subroutine alloc_datavalues(cdl,k)
      type(cdl_type), intent(inout) :: cdl
      integer, intent(in) :: k
      integer :: n
      character(len=cdl_strlen), dimension(:), allocatable :: tmp_
      if(.not.allocated(cdl%datavalues)) then
        allocate(cdl%datavalues(max(k,1000)))
      else
        n = size(cdl%datavalues)
        if(k+100.gt.n) then
          allocate(tmp_(n))
          tmp_(1:n) = cdl%datavalues(1:n)
          deallocate(cdl%datavalues)
          allocate(cdl%datavalues(max(k+100,2*n)))
          cdl%datavalues(1:n) = tmp_(1:n)
          deallocate(tmp_)
        endif
      endif
      end subroutine alloc_datavalues

      end module cdl_mod
