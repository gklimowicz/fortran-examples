! @sum sumfiles is a generic summation program for modelE acc-files.
! @+ See conventions.txt for documentation.
! @auth M. Kelley

      program sumfiles
      implicit none
      include 'netcdf.inc'
      integer :: status,ofid,ivarid,varid
      integer, dimension(:), allocatable :: fids,reduc_varids
      character(len=4096) :: ifile,ofile
      character(len=4096), dimension(:), allocatable :: ifiles
      integer :: n,nfiles,nfirst,nlast,iargc,nvars,nvars_reduc
      integer :: ndims,ivar
      integer :: itbeg,itend,itnow,n1Dif,n1To
      integer*8, dimension(:), allocatable :: reduc_varsizes
      integer*8 :: j,jj,jnext,accsize
      integer :: dsizes(7)
      character(len=6) :: ayear0,ayear
      character(len=3) :: amon0,amon
      real*4 :: dif_total,dif
      character(len=16) :: acc_period
      character(len=32) :: runid,vname
      character(len=100) :: fromto,fromto0
      character(len=132) :: xlabel
      integer, dimension(12) :: monacc,monacc1
      real*8, dimension(:), allocatable :: acc,acc_part
      real*8 :: bynfiles
      character(len=30) :: reduction
      character(len=30), dimension(:), allocatable :: reduc_list
      logical :: show_dif = .true.

c
c get the number of input files
c
      nfiles = iargc()

      if(nfiles.le.1) then
        write(6,*)
     &     'usage: sumfiles files_to_be_summed'
        write(6,*)
     &     '(works on modelE acc files, not on scaleacc_outputs)'
        stop
      endif
      allocate(fids(nfiles),ifiles(nfiles))
      bynfiles = 1d0/real(nfiles,kind=8)

c
c open each input file and find the min/max itime
c
      itbeg=+huge(itbeg)
      itend=-huge(itend)
      monacc(:) = 0  ;   dif_total = 0 ; dif = 0
      do n=1,nfiles
        call getarg(n,ifile)
        status = nf_open(trim(ifile),nf_nowrite,fids(n))
        ifiles(n) = ifile
        if(status.ne.nf_noerr) then
          write(6,*) 'nonexistent/non-netcdf input file ',trim(ifile)
          stop
        endif
        call get_var_int(fids(n),'itime0',itnow)
        if(itnow.lt.itbeg) then
          itbeg = itnow
          nfirst = n
        endif
        call get_var_int(fids(n),'itime',itnow)
        if(itnow.gt.itend) then
          itend = itnow
          nlast = n
        endif
        call get_var_int(fids(n),'monacc',monacc1)
        monacc(:) = monacc(:) + monacc1(:)
        status = nf_get_att_text(fids(n),nf_global,'fromto',fromto)
        n1Dif = index(fromto,'Dif:')+4
        if (fromto(n1Dif:n1Dif+6) == '       ') show_dif = .false.
        if (show_dif) read(fromto(n1Dif:n1Dif+7),*) dif
        dif_total = dif_total + dif
      enddo

c
c determine the name for the averaging period and the fromto-string
c
      status = nf_get_att_text(fids(nfirst),nf_global,'fromto',fromto0)
      status = nf_get_att_text(fids(nlast ),nf_global,'fromto',fromto )
      read(fromto0(6:16),'(a6,2x,a3)') ayear0,amon0
      n1To = index(fromto,'To:')+3
      read(fromto (n1To:n1To+10),'(a6,2x,a3)') ayear,amon
      call aperiod(monacc,ayear0,ayear,acc_period,amon0,amon)

      fromto(1:index(fromto0,'To:')) = fromto0(1:index(fromto0,'To:'))
      n1Dif = index(fromto,'Dif:')+4
      fromto(n1Dif:n1Dif+6) = '       ' ! if "Dif" is too large
      if (show_dif) then
        if (dif_total < 10000.) then
          write(fromto(n1Dif:n1Dif+6),'(f7.2)') dif_total
        else if (dif_total < 1000000.) then
          write(fromto(n1Dif:n1Dif+6),'(f7.0)') dif_total
!!      else
!!        show_dif = .false.
        end if
      end if

c
c copy the structure of the latest input file to the output file
c
      xlabel=''
      status = nf_get_att_text(fids(nlast),nf_global,'xlabel',xlabel)
      runid = xlabel(1:index(xlabel,' ')-1)
      ofile = trim(acc_period)//'.acc'//trim(runid)//'.nc'
      do n=1,nfiles
        call getarg(n,ifile)
        if(trim(ifile) .eq. trim(ofile)) then
          write(6,*) 'error: output file name ',trim(ofile),
     &         ' is an input file'
          stop
        endif
      enddo
      write(6,*) 'creating ',trim(ofile)
      status = nf_create(trim(ofile),nf_64bit_offset,ofid)
      call copy_file_structure(fids(nlast),ofid)

c
c make a list of the fields to be reduced
c
      status = nf_inq_nvars(ofid,nvars)
      allocate(reduc_varids(nvars),reduc_varsizes(nvars))
      allocate(reduc_list(nvars))
      nvars_reduc = 0
      do varid=1,nvars
        reduction=''
        status=nf_get_att_text(ofid,varid,'reduction',reduction)
        if(status.ne.nf_noerr) cycle
        nvars_reduc = nvars_reduc + 1
        reduc_varids(nvars_reduc) = varid
        status = nf_inq_varname(ofid,varid,vname)
        call get_varsize8(ofid,vname,accsize)
        reduc_varsizes(nvars_reduc) = accsize
        reduc_list(nvars_reduc) = reduction
      enddo

c
c collect the reduced data into the array acc
c

c Initialize acc
      accsize = sum(reduc_varsizes(1:nvars_reduc))
      allocate(acc(accsize))
      acc(:) = 0. ! default initialization
      j = 1
      do ivar=1,nvars_reduc
        jnext = j + reduc_varsizes(ivar)
        select case (trim(reduc_list(ivar)))
        case ('min')
          do jj=j,jnext-1
            acc(jj) = +1d30
          enddo
        case ('max')
          do jj=j,jnext-1
            acc(jj) = -1d30
          enddo
        end select
        j = jnext
      enddo
      allocate(acc_part(accsize))

c loop over files
      do n=1,nfiles

c loop over the fields to be reduced
        jnext = 1
        do ivar=1,nvars_reduc
          j = jnext
          jnext = j + reduc_varsizes(ivar)
          varid = reduc_varids(ivar)
          vname = ''
          status = nf_inq_varname(ofid,varid,vname)
          status = nf_inq_varid(fids(n),trim(vname),varid)
          if(status.ne.nf_noerr) then
            write(6,*) 'warning: skipping missing variable '//
     &           trim(vname)//' in '//trim(ifiles(n))
            cycle
          endif
          call get_varsize8(fids(n),vname,accsize)
          if(reduc_varsizes(ivar) .ne. accsize) then
            write(6,*) 'warning: skipping size-mismatched variable '//
     &           trim(vname)//' in '//trim(ifiles(n))
            cycle
          endif
          status = nf_get_var_double(fids(n),varid,acc_part(j))
          select case (trim(reduc_list(ivar)))
          case ('min')
            do jj=j,jnext-1
              acc(jj) = min(acc(jj),acc_part(jj))
            enddo
          case ('max')
            do jj=j,jnext-1
              acc(jj) = max(acc(jj),acc_part(jj))
            enddo
          case default ! 'sum' and 'avg' (postpone division)
            do jj=j,jnext-1
              acc(jj) = acc(jj) + acc_part(jj)
            enddo
          end select
        enddo
      enddo

      j = 1
      do ivar=1,nvars_reduc
        jnext = j + reduc_varsizes(ivar)
        select case (trim(reduc_list(ivar)))
        case ('avg')
          do jj=j,jnext-1
            acc(jj) = acc(jj)*bynfiles
          enddo
        end select
        j = jnext
      enddo

c
c fill in the output file data from acc or the last input file
c
      call put_shared_vars(fids(nlast),ofid,
     &     acc,nvars_reduc,reduc_varids,reduc_varsizes )

c
c correct accumulation period info: start, months included, "fromto"
c
      call put_var_int(ofid,'itime0',itbeg)
      call put_var_int(ofid,'monacc',monacc)
      status = nf_put_att_text(ofid,nf_global,'fromto'
     &     ,int(len_trim(fromto),kind=8),fromto)

c
c close input and output files
c
      status = nf_close(ofid)
      do n=1,nfiles
        status = nf_close(fids(n))
      enddo

      deallocate(acc,acc_part)

      deallocate(fids)

      end program sumfiles

      subroutine put_shared_vars(ifid,ofid,
     &     acc,nvars_reduc,reduc_varids,reduc_varsizes )
      implicit none
      include 'netcdf.inc'
      integer :: ifid,ofid
      real*8 :: acc(1)
      integer :: nvars_reduc
      integer, dimension(nvars_reduc) :: reduc_varids
      integer*8, dimension(nvars_reduc) :: reduc_varsizes
c
      integer :: status,ivarid,ovarid,nvars,vtype
      integer :: vsize,vsize_input
      character(len=32) :: vname
      real*4, dimension(:), allocatable :: v4
      real*8, dimension(:), allocatable :: v8
      integer*8 :: j
      integer :: k

      j = 1

      status = nf_inq_nvars(ofid,nvars)
      ovarid_loop: do ovarid=1,nvars
         do k=1,nvars_reduc
           if(reduc_varids(k).eq.ovarid) then
             status = nf_put_var_double(ofid,ovarid,acc(j))
             j = j + reduc_varsizes(k)
             cycle ovarid_loop
           endif
         enddo
         status = nf_inq_varname(ofid,ovarid,vname)
         status = nf_inq_varid(ifid,vname,ivarid)
         if(status.ne.nf_noerr) cycle
         status = nf_inq_vartype(ofid,ovarid,vtype)
         call get_varsize(ofid,vname,vsize)
         call get_varsize(ifid,vname,vsize_input)
         if(vsize.ne.vsize_input) cycle
         if(vtype.eq.nf_double) then ! real*8 data
           allocate(v8(vsize))
           status = nf_get_var_double(ifid,ivarid,v8)
           status = nf_put_var_double(ofid,ovarid,v8)
           deallocate(v8)
         elseif(vtype.eq.nf_char) then ! character data
           allocate(v4(vsize)) ! larger than necessary
           status = nf_get_var_text(ifid,ivarid,v4)
           status = nf_put_var_text(ofid,ovarid,v4)
           deallocate(v4)
         else                   ! floats and integers
           allocate(v4(vsize))
           status = nf_get_var_real(ifid,ivarid,v4)
           status = nf_put_var_real(ofid,ovarid,v4)
           deallocate(v4)
         endif
      enddo ovarid_loop
      return
      end subroutine put_shared_vars

      subroutine aperiod(monacc,ayr0,ayr1,acc_period,amon0,amon)
! Find appropriate name for the accumulation period e.g. MonYear1-Year2
! For NonEarth calendars we assume the months are called A..,B..,C..,..
!      Notes about the output file names:
! restrictions: months_per_year is currently restricted to 12,
!               the months are called JAN FEB ... or Jan Feb ...
! ambivalence:  different periods may have the same name X-Z...
!               e.g. A-M may be Aug-Mar or Aug-May
!               the name of any 12-month period is ANN
      implicit none
      integer :: monacc(12),mon0,yr0,yr1
      character(len=3) :: amon0,amon
      character(len=6) :: ayr0,ayr1
      character(len=16) :: acc_period
      character(len=26), parameter :: alphabet =
     *  'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
      character(len=36) :: emonth='JANFEBMARAPRMAYJUNJULAUGSEPOCTNOVDEC'
      character(len=36) :: emnthx='JanFebMarAprMayJunJulAugSepOctNovDec'
      character(len=12) :: mostr
      integer :: m,mm,nmo,ninc,ndec
      logical :: incyr0, earth = .true.

      if(minval(monacc,mask=monacc>0).ne.maxval(monacc,mask=monacc>0))
     &     stop 'unequal numbers of months'

      mon0 = index(emonth,amon0)                 ! all caps
      if(mon0 < 0) mon0 = index(emnthx,amon0)    ! mixed mode
      if(mod(mon0-1,3) .ne. 0) mon0 = -1         ! not likely
      if(mon0 > 0) mon0 = 1 + (mon0-1)/3
      if(mon0 < 0) then
        mon0 = index(alphabet,amon0(1:1))
        earth = .false.
      end if
      if(mon0 .le. 0) then
        write(*,'(a)') 'Unknown month: '//amon0
        stop
      end if

      mostr=''
      nmo = 0 ; incyr0 = .false.
      do mm=mon0,mon0+11             ! mon0+months_in_year - 1
        m = mm ; if(m > 12) m = m-12
        if(monacc(m).eq.0) cycle
        if(mm>12) incyr0 = .true.
        nmo = nmo + 1
        if(nmo.eq.1) then
          mostr(1:3)=amon0
        else
          if (earth) then
            mostr(nmo:nmo)=emonth(1+3*(m-1):1+3*(m-1))
          else
            mostr(nmo:nmo)=alphabet(m:m)
          end if
        end if
      enddo

      if(minval(monacc)>0) then
        mostr='ANN'
      elseif(nmo.eq.2) then
        mostr(3:3) = mostr(2:2)
        mostr(2:2) = '+'
      elseif(nmo.gt.3) then
        mostr(3:3) = mostr(nmo:nmo)
        mostr(2:2) = '-'
        if(earth) write(*,*) 'labeling may be ambiguous !!'
      end if

      read(ayr0,*) yr0 ; ayr0='      '
      if(incyr0) yr0 = yr0 + 1
      if(yr0.lt.10000) then
        write(ayr0(1:4),'(i4.4)') yr0
      else
        write(ayr0,'(i6)') yr0
      end if

      read(ayr1,*) yr1 ; ayr1='      '
      if(amon=='JAN') yr1 = yr1 - 1
      if(.not.earth.and.amon(1:1)=='A') yr1 = yr1 - 1
      if(yr1.lt.10000) then
        write(ayr1(1:4),'(i4.4)') yr1
      else
        write(ayr1,'(i6)') yr1
      end if

      acc_period=''
      acc_period = mostr(1:3)//trim(adjustl(ayr0))
      if(yr1.gt.yr0)  acc_period = mostr(1:3)//trim(adjustl(ayr0))//
     *   '-'//trim(adjustl(ayr1))

c check for gaps
      ninc = count(monacc(2:12).gt.monacc(1:11))
      ndec = count(monacc(2:12).lt.monacc(1:11))
      if(ninc+ndec.gt.2) then
        write(6,*) 'gap'
      end if

      return
      end subroutine aperiod
