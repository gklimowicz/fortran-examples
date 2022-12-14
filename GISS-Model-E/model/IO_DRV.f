#include "rundeck_opts.h"

      SUBROUTINE io_rsf(fname,it,iaction_in,ioerr)
!@sum  io_rsf manages the reading/writing of restart and acc files
!@auth M. Kelley
!@ver  beta.

      use iso_c_binding, only: C_CHAR, C_NULL_CHAR
C**** For all iaction < 0  ==> WRITE, For all iaction > 0  ==> READ
      USE DOMAIN_DECOMP_ATM, only : grid,am_i_root
      USE MODEL_COM, only : ioread_acc,iowrite_single,irerun,
     *                      ioread,iowrite,iowrite_mon
     &     ,itimei,rsf_file_name,kcopy
      use pario, only : par_open,par_close,par_enddef
#ifdef TRACERS_OceanBiology
      use obio_diag, only: new_io_obio_diag
#endif

      IMPLICIT NONE

      interface
        ! see: man 2 rename
        subroutine rename_c(oldpath,newpath) bind(C, name="rename")
          use iso_c_binding, only: c_char
          character(kind=c_char) :: oldpath(*)
          character(kind=c_char) :: newpath(*)
        end subroutine rename_c
      end interface



!@var fname name of file to be read or written
      character(len=*) :: fname
!@var iaction_in flag for reading or writing rsf file
      INTEGER, INTENT(IN) :: iaction_in
!@var it hour of model run
      INTEGER, INTENT(INOUT) :: it
!@var IOERR (1,0,-1) if there (is, is maybe, is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
      integer :: fid,iorw,iaction
      logical :: do_io_prog,do_io_acc,do_io_longacc,r4
      character(len=200) :: tmpname

      iaction = iaction_in

      call set_ioptrs_acc_default

      if(iaction.eq.iowrite_single) then
        call calc_derived_acc
        call set_ioptrs_acc_extended
      endif

      do_io_prog = .true.
      if(iaction.eq.iowrite_single) do_io_prog = .false.
      if(iaction.eq.ioread_acc)     do_io_prog = .false.

      do_io_acc = .false.
c this logic would be much simpler if there were fewer
c iaction possibilities related to reruns.  reading
c arrays by name will eliminate a number of the
c rerun cases.
#ifndef DEFER_ACC_READ
      if(iaction.eq.ioread)         do_io_acc = .true.
#endif
      if(iaction.eq.iowrite)        do_io_acc = .true.
      if(iaction.eq.iowrite_single) do_io_acc = .true.
      if(iaction.eq.ioread_acc)     do_io_acc = .true.

      do_io_longacc = do_io_acc
      if(iaction.eq.iowrite_mon)    do_io_longacc = .true.
      if(iaction.eq.irerun)         do_io_longacc = .true.

c for routines designed to accept only iowrite or ioread:
      if(iaction.le.iowrite) then
        iorw = iowrite
      else
        iorw = ioread
      endif
      if(iaction.eq.ioread_acc) iaction = ioread

      ioerr=-1

      tmpname = trim(fname)//'.nc'
      if(iorw.eq.ioread) then   ! open input file
        if(trim(fname) == 'AIC') tmpname=fname ! symbolic link w/o .nc suffix
        fid = par_open(grid,trim(tmpname),'read')
      else                      ! define contents of output file
        if(iaction.eq.iowrite) tmpname = 'checkpoint.nc'
        fid = par_open(grid,trim(tmpname),'create')
        call def_rsf_label(fid)
        if(do_io_prog) call def_rsf_prog(fid)
        r4 = .false.
        if(iaction.eq.iowrite_single) r4=.true. ! real*4 disk storage
        if(do_io_acc) call def_acc_all(fid,r4) 
        if(do_io_longacc) call def_rsf_longacc(fid,r4)
        if(iaction.eq.iowrite_single) call def_acc_meta(fid)
        call par_enddef(grid,fid)
      endif

      call new_io_label  (fid,iaction,it)

c
c prognostic arrays
c
      if(do_io_prog) then
        call new_io_atmvars(fid,iorw)
        call new_io_ocean  (fid,iorw)
        call new_io_seaice (fid,iorw)
      end if

c
c acc arrays
c
      if(do_io_acc) then
        if(iaction.eq.iowrite_single) call write_acc_meta(fid)
        call new_io_acc(fid,iaction)
#ifndef STANDALONE_OCEAN
        call new_io_glaacc(fid,iaction)
#endif
        call new_io_ocdiag(fid,iaction)
#ifdef TRACERS_OceanBiology
        call new_io_obio_diag(fid, iaction)
#endif
#ifndef STANDALONE_HYCOM
        call new_io_icdiag(fid,iorw)
#endif
      endif
C**** KCOPY > 2 : ALSO SAVE THE OCEAN DATA TO INITIALIZE DEEP OCEAN RUNS
      if(kcopy.gt.2 .and. iaction.eq.iowrite_mon) then
        call new_io_acc(fid,iaction)
      endif

c
c long-period acc arrays
c
      if(do_io_longacc) call new_io_longacc(fid,iorw)

c
c close input/output file
c
      call par_close(grid,fid)

      if(iaction.eq.iowrite .and. am_i_root()) then
        call rename_c(
     &      C_CHAR_'checkpoint.nc'//C_NULL_CHAR,
     &      C_CHAR_''//trim(fname)//'.nc'//C_NULL_CHAR)
      endif

      RETURN
      END SUBROUTINE io_rsf

      subroutine find_later_rsf(kdisk)
!@sum set kdisk such that Itime in rsf_file_name(kdisk) is
!@+   the larger of the Itimes in rsf_file_name(1:2).
!@auth  M. Kelley
      use model_com, only : rsf_file_name
      use domain_decomp_atm, only : grid
      use pario, only : read_data,par_open,par_close
      implicit none
      integer, intent(out) :: kdisk
      integer :: itimes(2),fid,k
      itimes(:) = -1
      do k=1,2
        fid = par_open(grid,trim(rsf_file_name(k))//'.nc','read')
        call read_data(grid,fid,'itime',itimes(k),bcast_all=.true.)
        call par_close(grid,fid)
      enddo
      if(maxval(itimes).eq.-1) call stop_model(
     &     'FIND_LATER_RSF: ERRORS ON BOTH RESTART FILES',255)
      kdisk = 1
      if(itimes(2).gt.itimes(1)) kdisk=2
      return
      end subroutine find_later_rsf

      subroutine def_rsf_prog(fid)
!@sum  def_rsf_prog defines prognostic array structure in restart files
!@auth M. Kelley
!@ver  beta
      implicit none
      integer :: fid
      call def_rsf_atmvars(fid)
      call def_rsf_ocean  (fid)
      call def_rsf_seaice (fid)
      return
      end subroutine def_rsf_prog

      subroutine def_acc_meta(fid)
!@sum  def_acc_meta defines metadata in acc files
!@auth M. Kelley
!@ver  beta
      use model_com, only: modelEclock
      USE MODEL_COM, only : jhour0,jdate0,amon,amon0,
     &     jyear0,itime,itime0,nday
      use domain_decomp_atm, only : grid
      use pario, only : write_attr
#ifdef TRACERS_OceanBiology
      use obio_diag, only: def_meta_obio_diag
#endif
      implicit none
      integer :: fid         !@var fid file id
      character(len=100) :: fromto
      real*8 :: days
      integer :: year, hour, date

      call modelEclock%get(year=year, hour=hour, date=date)

      days=(itime-itime0)/float(nday)
      write(fromto,902) jyear0,amon0,jdate0,jhour0,
     &     year,amon,date,hour,itime,days
      call write_attr(grid,fid,'global','fromto',fromto)

c idacc(5) is not additive
      call write_attr(grid,fid,'idacc','reduction','sum')

      call def_meta_atmacc(fid)
#ifndef STANDALONE_OCEAN
      call def_meta_glaacc(fid)
#endif
      call def_meta_ocdiag(fid)
#ifdef TRACERS_OceanBiology
      call def_meta_obio_diag(fid)
#endif
#ifndef STANDALONE_HYCOM
      call def_meta_icdiag(fid)
#endif

      return
  902 FORMAT ('From:',I6,A6,I2,',  Hr',I3,
     *  6X,'To:',I6,A6,I2,', Hr',I3,'  Model-Time:',I9,5X,
     *  'Dif:',F7.2,' Days')
      end subroutine def_acc_meta

      subroutine write_acc_meta(fid)
!@sum  def_acc_meta writes metadata to acc files
!@auth M. Kelley
!@ver  beta
#ifdef TRACERS_OceanBiology
      use obio_diag, only: write_meta_obio_diag
#endif
      implicit none
      integer :: fid         !@var fid file id
      call write_meta_atmacc(fid)
#ifndef STANDALONE_OCEAN
      call write_meta_glaacc(fid)
#endif
      call write_meta_ocdiag(fid)
#ifdef TRACERS_OceanBiology
      call write_meta_obio_diag(fid)
#endif
#ifndef STANDALONE_HYCOM
      call write_meta_icdiag(fid)
#endif
      return
      end subroutine write_acc_meta

      subroutine def_acc_all(fid,r4_on_disk)
!@sum  def_acc_all defines acc array structure in restart files
!@auth M. Kelley
!@ver  beta
#ifdef TRACERS_OceanBiology
      use obio_diag, only: def_rsf_obio_diag
#endif
      implicit none
      integer :: fid         !@var fid file id
      logical :: r4_on_disk  !@var r4_on_disk if true, real*8 stored as real*4
      call def_rsf_acc(fid,r4_on_disk)
#ifndef STANDALONE_OCEAN
      call def_rsf_glaacc(fid,r4_on_disk)
#endif
      call def_rsf_ocdiag(fid,r4_on_disk)
#ifndef STANDALONE_HYCOM
      call def_rsf_icdiag(fid,r4_on_disk)
#endif
#ifdef TRACERS_OceanBiology
      call def_rsf_obio_diag(fid, r4_on_disk)
#endif
      return
      end subroutine def_acc_all

      subroutine def_rsf_label(fid)
!@sum  def_rsf_label defines control info in restart/acc files
!@auth M. Kelley
!@ver  beta
      use model_com, only : xlabel,itime,itimee,itime0,itimei,
     &     iyear1,nday,iowrite
      use timings, only : ntimemax,ntimeacc,timestr,timing
      use domain_decomp_atm, only : grid
      use pario, only : defvar,write_attr
      implicit none
      integer fid   !@var fid file id
      integer :: intdum
      call write_attr(grid,fid,'global','xlabel',xlabel)
      call defvar(grid,fid,itime,'itime')
      call write_caldate(fid)
      call defvar(grid,fid,itimee,'itimee')
      call defvar(grid,fid,itime0,'itime0')
      call defvar(grid,fid,itimei,'itimei')
      call defvar(grid,fid,iyear1,'iyear1')
      call defvar(grid,fid,nday,'nday')
      call defvar(grid,fid,timing(1:ntimemax),'cputime(ntimemax)')
      call defvar(grid,fid,ntimeacc,'ntimeacc')
      call io_cputime(fid,iowrite)
c rparam, iparam, cparam are dummy vars which hold the
c parameter database in their attributes.
      call defvar(grid,fid,intdum,'rparam')
      call defvar(grid,fid,intdum,'iparam')
      call defvar(grid,fid,intdum,'cparam')
      call new_io_param(fid,iowrite,.false.)
      return
      end subroutine def_rsf_label

      subroutine new_io_label(fid,iaction,ihrX)
!@sum  new_io_label read/write control info from/to restart,acc files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com, only : xlabel,itime,itimee,itime0,itimei,
     &     iyear1,nday,iowrite,ioread,irsfic,irsficnt,irsficno
      use domain_decomp_atm, only: grid
      use pario, only : write_data,read_data,write_attr,read_attr
      use timings, only : ntimemax,ntimeacc,timestr,timing
      use Dictionary_mod
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      integer :: ihrX !@var ihrX a dummy itime read from IC files
      integer :: idum,nday_dummy
      logical :: is_ic
!@dbparam keep_params if =1, params are read in from a *.rsf* file
!@+                   if =0, defaults and rundeck values are used
!@+       this is used only in istart=8 starts
      integer :: keep_params=0 ! .false.

      select case (iaction)
      case (:iowrite) ! output to restart or acc file
        call write_data(grid, fid, 'itime', itime)
        call write_data(grid, fid, 'itimee', itimee)
        call write_data(grid, fid, 'itime0', itime0)
        call write_data(grid, fid, 'itimei', itimei)
        call write_data(grid, fid, 'iyear1', iyear1)
        call write_data(grid, fid, 'nday', nday)
        call write_data(grid, fid, 'ntimeacc', ntimeacc)
        call write_data(grid, fid, 'cputime', timing(1:ntimemax))
      case (ioread:) ! input from restart or acc file
        call read_attr(grid,fid,'global','xlabel',idum,xlabel)
        call read_data(grid,fid,'nday',nday_dummy,bcast_all=.true.)
        is_ic = .false.
        if(iaction.eq.irsfic)   is_ic = .true.
        if(iaction.eq.irsficnt) is_ic = .true.
        if(iaction.eq.irsficno) is_ic = .true.
        if(.not.is_ic) then
          nday = nday_dummy
          call read_data(grid,fid,'itime', itime, bcast_all=.true.)
          call read_data(grid,fid,'itimee',itimee,bcast_all=.true.)
          call read_data(grid,fid,'itime0',itime0,bcast_all=.true.)
          call read_data(grid,fid,'itimei',itimei,bcast_all=.true.)
          if(iyear1.lt.0) ! this check not present for ioread_acc!!!
     &    call read_data(grid,fid,'iyear1',iyear1,bcast_all=.true.)
          call read_data(grid,fid,'ntimeacc',ntimeacc,
     &         bcast_all=.true.)
          call read_data(grid,fid, 'cputime', timing(1:ntimemax),
     &         bcast_all=.true.)
          call io_cputime(fid,ioread)
          call new_io_param(fid,ioread,.false.)
        else
          call read_data(grid,fid,'itime', IhrX, bcast_all=.true.)
!!        IhrX = IhrX*24/nday_dummy
          IhrX=nint(IhrX*(24.d0/nday_dummy)) ! to prevent overflow
          call sync_param('keep_params',keep_params)
          if(keep_params==1) call new_io_param(fid,ioread,.false.)
        endif
      end select
      end subroutine new_io_label

      subroutine write_caldate(fid)
c write a text version of the date to a restart/acc file
      use model_com
      use pario, only : write_attr
      use domain_decomp_atm, only: grid
      implicit none
      integer :: fid
      character(len=19) :: caldate
      character(len=2) :: cmo,cday
      character(len=4) :: cyr
      character(len=5) :: chr
      integer :: month

      write(cmo,'(i2.2)') modelEclock%getMonth()
      write(cday,'(i2.2)') modelEclock%getDate()
      write(cyr,'(i4.4)') modelEclock%getYear()
!      write(chr,'(f4.1)') real(jhour)
      write(chr,'(f5.2)') (real(mod(itime,nday))/real(nday))*24.
      caldate=cmo//'/'//cday//'/'//cyr//' hr '//chr
      call write_attr(grid,fid,'itime','caldate',caldate)
      return
      end subroutine write_caldate

      subroutine io_cputime(fid,iaction)
c manage the reading/writing of timing information. could be done better
      use model_com, only : iowrite,ioread,nday,itime,itime0
      use timings, only : ntimemax,ntimeacc,timestr,timing
      use domain_decomp_atm, only: grid,broadcast
      use pario, only : write_attr,read_attr
      implicit none
      integer :: fid,iaction
      character*2 :: c2
      character(len=5) :: cfrac
      character(len=6) :: fracnames(ntimemax)
      character(len=17) :: timestr_pct
      integer :: n,idum
      real*8 :: tsum,min_per_day
      do n=1,ntimemax
        write(c2,'(i2.2)') n
        fracnames(n) = 'frac'//c2
      enddo
      select case (iaction)
      case (:iowrite)
        call broadcast(grid,timing)
        tsum = sum(timing(1:ntimeacc))+1d-20
        if(itime.gt.itime0) then
          min_per_day = (tsum/60.)*nday/dble(itime-itime0)
        else
          min_per_day = 0.
        endif
        min_per_day = int(min_per_day*100d0)*.01d0
        call write_attr(grid,fid,'cputime','minutes_per_day',
     &       min_per_day)
        do n=1,ntimeacc
          write(cfrac,'(f5.1)') 100.*timing(n)/tsum
          timestr_pct = timestr(n)//cfrac
          call write_attr(grid,fid,'cputime',fracnames(n),
     &         timestr_pct)
        enddo
      case (ioread:)
        do n=1,ntimeacc
          call read_attr(grid,fid,'cputime',fracnames(n),idum,
     &         timestr_pct)
          timestr(n) = timestr_pct(1:12)
        enddo
      end select
      return
      end subroutine io_cputime

      subroutine new_io_param(fid,iaction,ovrwrt)
!@sum  new_io_param read/write parameter database from/to restart files
!@auth M. Kelley
!@ver  beta new_ prefix avoids name clash with the default version
      use model_com
      use domain_decomp_atm, only : grid
      use pario, only : write_attr, read_attr, get_natts
      use Dictionary_mod
      implicit none
      integer fid   !@var fid file id
      integer iaction !@var iaction flag for reading or writing to file
      logical :: ovrwrt
      integer :: l,m,n,plen,strlen,niparam,nrparam,ncparam
      character*1 :: ptype
      integer :: pvali(1000)
      real*8 :: pvalr(1000)
      character(len=128) :: pname,pvalc(1000)
      integer, parameter :: lcstr=1500
      character(len=lcstr) :: cstr
      select case (iaction)
      case (iowrite) ! output
        niparam = 0; nrparam = 0; ncparam = 0;
        n = 0
        do
          n = n + 1
          call query_param( n, pname, plen, ptype )
          if(pname(1:5).eq.'EMPTY') exit
          select case(ptype)
          case ('i')
            call get_param(pname,pvali,plen,update_access_flag=.false.)
            call write_attr(grid,fid,'iparam',trim(pname),
     &           pvali(1:plen))
            niparam = niparam + 1
          case ('r')
            call get_param(pname,pvalr,plen,update_access_flag=.false.)
            call write_attr(grid,fid,'rparam',trim(pname),
     &           pvalr(1:plen))
            nrparam = nrparam + 1
          case ('c')
            call get_param(pname,pvalc,plen,update_access_flag=.false.)
            cstr = trim(pvalc(1))
c Arrays of strings are not supported here.  Separate the
c strings with the character "|".
            do m=2,plen
              cstr=trim(cstr)//'|'//trim(pvalc(m))
            enddo
            if(len_trim(cstr).gt.0) then
              call write_attr(grid,fid,'cparam',trim(pname),cstr)
              ncparam = ncparam + 1
            endif
          end select
        enddo
      case (ioread) ! input
        call get_natts(grid,fid,'iparam',niparam)
        call get_natts(grid,fid,'rparam',nrparam)
        call get_natts(grid,fid,'cparam',ncparam)
        do n=1,niparam
          call read_attr(grid,fid,'iparam',pname,plen,pvali,
     &         attnum=n)
          if ( (.not. is_set_param(pname)) .or. ovrwrt )
     &         call set_param( trim(pname), pvali, plen, 'o' )
        enddo
        do n=1,nrparam
          call read_attr(grid,fid,'rparam',pname,plen,pvalr,
     &         attnum=n)
          if ( (.not. is_set_param(pname)) .or. ovrwrt )
     &         call set_param( trim(pname), pvalr, plen, 'o' )
        enddo
        do n=1,ncparam
          cstr(lcstr:lcstr)='x' ! workaround for length checking
          call read_attr(grid,fid,'cparam',pname,strlen,cstr,
     &         attnum=n)
c Arrays of strings are not supported here.  Strings containing
c the character "|" are separated into arrays of strings.
          l = 0
          plen = 1
          pvalc = ''
          do m=1,strlen
            l = l + 1
            if(cstr(m:m).eq.'|') then
              l = 0
              plen = plen + 1
            else
              pvalc(plen)(l:l) = cstr(m:m)
            endif
          enddo
          if ( (.not. is_set_param(pname)) .or. ovrwrt )
     &         call set_param( trim(pname), pvalc, plen, 'o' )
        enddo
      end select
      return
      end subroutine new_io_param

      subroutine set_ioptrs_acc_default
c point i/o pointers for accumlated quantities to the
c instances of the arrays used during normal operation.
      implicit none
      call set_ioptrs_atmacc_default
      call set_ioptrs_ocnacc_default
      call set_ioptrs_iceacc_default
      return
      end subroutine set_ioptrs_acc_default

      subroutine set_ioptrs_acc_extended
c point i/o pointers for accumlated quantities to the
c instances of the arrays holding derived quantities as well
      implicit none
      call set_ioptrs_atmacc_extended
      call set_ioptrs_ocnacc_extended
c      call set_ioptrs_iceacc_extended
      return
      end subroutine set_ioptrs_acc_extended

      subroutine calc_derived_acc
      implicit none
      call calc_derived_acc_atm
      call diag_ocean_prep
      return
      end subroutine calc_derived_acc
