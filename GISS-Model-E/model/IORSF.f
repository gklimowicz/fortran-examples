#include "rundeck_opts.h"

      SUBROUTINE io_rsf(filenm,it,iaction,ioerr)
!@sum   io_rsf controls the reading and writing of the restart files
!@auth  Gavin Schmidt
!@calls io_model,io_ocean,io_lakes,io_seaice,io_earth,io_soils,io_snow
!@+     io_landice,io_bldat,io_pbl,io_clouds,io_somtq,io_rad,io_diags
!@+     io_ocdiag,io_icedyn,io_icdiag
      USE FILEMANAGER, only : openunit,closeunit
      USE DOMAIN_DECOMP_1D, only : am_i_root !REWIND_PARALLEL
      USE ATM_COM, only: Kradia
      USE MODEL_COM, only : ioread_single,iowrite_single
     *                     ,ioread,ioread_nodiag,iowrite

      IMPLICIT NONE
!@var filenm name of file to be read or written
      character(len=*) :: filenm
!@var iaction flag for reading or writing rsf file
      INTEGER, INTENT(IN) :: iaction
!@var it hour of model run
      INTEGER, INTENT(INOUT) :: it
!@var IOERR (1,0,-1) if there (is, is maybe, is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var IT1 hour for correct reading check
!@var ITM maximum hour for post-processing
!@var kunit Fortran unit number of file i/o
      INTEGER IT1,itm,iact,kunit
      logical skip_diag

      iact=iaction   ; skip_diag=.false.
      if(iaction==ioread_nodiag) then
         iact=ioread ; skip_diag=.true.
      end if

      ioerr=-1

      if(iaction.le.iowrite) then
c open output files with status = 'UNKNOWN'
        if(am_i_root()) ! only root PE can write
     &       call openunit(trim(filenm),kunit,.true.,.false.)
      elseif(iaction.ge.ioread) then
c open input files with status = 'OLD'
c all PEs need to read the label records, so no am_i_root check
        call openunit(trim(filenm),kunit,.true.,.true.)
      endif
c      if (kunit.gt.0) CALL REWIND_PARALLEL( kunit )

C**** For all iaction < 0  ==> WRITE, For all iaction > 0  ==> READ
C**** Particular values may produce variations in indiv. i/o routines

C**** Calls to individual i/o routines
      call io_label  (kunit,it,itm,iact,ioerr)
      it1=it
      if (Kradia>9) go to 10 ! used for testing daily stuff
      if (Kradia.gt.0) then
        if (iaction.ne.ioread_single .and.
     *   iaction.ne.iowrite_single) call io_rad (kunit,iact,ioerr)
        call io_diags  (kunit,it,iact,ioerr)
        if(it1.ne.it .or. ioerr.eq.1)
     &       call stop_model('restart problem',255)
        go to 10
      end if
      if(iaction.ne.ioread_single.and.iaction.ne.iowrite_single) then
        call io_model  (kunit,iact,ioerr)
        call io_strat  (kunit,iact,ioerr)
        call io_ocean  (kunit,iact,ioerr)
        call io_lakes  (kunit,iact,ioerr)
        call io_seaice (kunit,iact,ioerr)
        call io_earth  (kunit,iact,ioerr)
        call io_soils  (kunit,iact,ioerr)
        call io_vegetation  (kunit,iact,ioerr)
        !!! actually not sure if this call is needed
        !!! (seems like it is duplicated in io_vegetation...)
        call io_veg_related  (kunit,iact,ioerr)
        !call io_ent    (kunit,iaction,ioerr)
        call io_snow   (kunit,iact,ioerr)
        call io_landice(kunit,iact,ioerr)
        call io_bldat  (kunit,iact,ioerr)
        call io_pbl    (kunit,iact,ioerr)
        call io_clouds (kunit,iact,ioerr)
#ifdef AUTOTUNE_LIGHTNING
        call io_lightning(kunit,iact,ioerr)
#endif
        call io_somtq  (kunit,iact,ioerr)
        call io_rad    (kunit,iact,ioerr)
        call io_icedyn (kunit,iact,ioerr)
#ifdef CALCULATE_FLAMMABILITY
        call io_flammability(kunit,iact,ioerr)
#endif
#ifdef TRACERS_ON
        call io_tracer (kunit,iact,ioerr)
#endif
      end if
      if(skip_diag) go to 10
      call io_diags  (kunit,it,iaction,ioerr)
      call io_ocdiag (kunit,it,iaction,ioerr)
      call io_icdiag (kunit,it,iaction,ioerr)
#ifdef TRACERS_ON
      call io_trdiag (kunit,it,iaction,ioerr)
#endif
      if (it1.ne.it) THEN
        WRITE(6,*) "TIMES DO NOT MATCH READING IN RSF FILE",it,it1
        ioerr=1
      END IF
      if (ioerr.eq.1) WRITE(6,*) "I/O ERROR IN RESTART FILE: KUNIT="
     *     ,kunit

C**** return maximum time
   10 it=itm

      if(am_i_root() .or. iaction.ge.ioread) call closeunit(kunit)

      RETURN
      END SUBROUTINE io_rsf

      subroutine read_ground_ic
!@sum   read_ground_ic read initial conditions file for
!@+     sea ice, land surface, land ice.  Transplanted from INPUT.
!@auth  M. Kelley
!@calls io_seaice,io_earth,io_soils,io_landice
      use model_com, only : ioreadnt
      use filemanager, only : openunit,closeunit
      use domain_decomp_1d, only : am_i_root
      implicit none
      integer :: iu_GIC,ioerr
      call openunit("GIC",iu_GIC,.true.,.true.)
      ioerr=-1
      read(iu_GIC) ! ignore first line (ocean ic done in init_OCEAN)
      call io_seaice (iu_GIC,ioreadnt,ioerr)
      call io_earth  (iu_GIC,ioreadnt,ioerr)
      call io_soils  (iu_GIC,ioreadnt,ioerr)
      call io_landice(iu_GIC,ioreadnt,ioerr)
      if (ioerr.eq.1) then
        if (AM_I_ROOT())
     *       WRITE(6,*) "I/O ERROR IN GIC FILE: KUNIT=",iu_GIC
        call stop_model("INPUT: GIC READ IN ERROR",255)
      end if
      call closeunit (iu_GIC)
      return
      end subroutine read_ground_ic

      subroutine find_later_rsf(kdisk)
!@sum set kdisk such that Itime in rsf_file_name(kdisk) is
!@+   the larger of the Itimes in rsf_file_name(1:2).
!@+   Transplanted from INPUT.
      USE FILEMANAGER, only : openunit,closeunit
      use model_com, only : rsf_file_name
      implicit none
      integer, intent(out) :: kdisk
      integer :: Itime1,Itime2,kunit
      Itime1=-1
      call openunit(rsf_file_name(1),kunit,.true.,.true.)
      READ (kunit,ERR=410) Itime1
 410  continue
      call closeunit(kunit)
      Itime2=-1
      call openunit(rsf_file_name(2),kunit,.true.,.true.)
      READ (kunit,ERR=420) Itime2
 420  continue
      call closeunit(kunit)
      if (Itime1+Itime2.LE.-2) then
        call stop_model(
     &       'FIND_LATER_RSF: ERRORS ON BOTH RESTART FILES',255)
      endif
      KDISK=1
      IF (Itime2.GT.Itime1) KDISK=2
      return
      end subroutine find_later_rsf

      SUBROUTINE io_label(kunit,it,itm,iaction,ioerr)
!@sum  io_label reads and writes label/parameters to file
!@auth Gavin Schmidt
      USE ATM_COM, only: Kradia
      USE MODEL_COM
      USE RESOLUTION, only : IM,JM,LM
      USE RESOLUTION, only : LS1
      USE DOMAIN_DECOMP_1D, only : AM_I_ROOT
      USE TIMINGS, only : ntimemax,ntimeacc,timestr,timing
      USE Dictionary_mod
      IMPLICIT NONE

      INTEGER kunit   !@var kunit unit number of read/write
      INTEGER iaction !@var iaction flag for reading or writing to file
!@var IOERR 1 (or -1) if there is (or is not) an error in i/o
      INTEGER, INTENT(INOUT) :: IOERR
!@var it input/ouput value of hour
!@var itm maximum hour returned (different from it if post-processing)
      INTEGER, INTENT(INOUT) :: it,itm
!@var LABEL2 content of record 2
      CHARACTER*80 :: LABEL2
!@var NTIM1,TSTR1,TIM1 timing related dummy arrays
      INTEGER NTIM1,ITIM1(NTIMEMAX)
      REAL*8 ::      TIM1(NTIMEMAX)
      CHARACTER*12  TSTR1(NTIMEMAX)
!@var ITmin,ITmax minimal/maximal time in acc periods to be combined
      INTEGER, SAVE :: ITmax=-1, ITmin=-1 ! to protect against long runs
      INTEGER nd1,iy1,iti1,ite1,it01,im0,jm0,lm0,ls10

      SELECT CASE (IACTION)
      CASE (:IOWRITE)           ! output to end-of-month restart file
        IF (AM_I_ROOT()) THEN
          WRITE (kunit,err=10) it,XLABEL,nday,iyear1,itimei,itimee,
     *         itime0,NTIMEACC,TIMING(1:NTIMEACC),TIMESTR(1:NTIMEACC)
C**** doc line: basic model parameters
          write(label2,'(a13,4i4,a)') 'IM,JM,LM,LS1=',im,jm,lm,ls1,' '
          label2(74:80) = 'LABEL01'
          WRITE (kunit,err=10) LABEL2
C**** write parameters database here
          call write_param(kunit)
        END IF
      CASE (IOREAD:)          ! label always from input file
!****   Determine whether timing numbers were saved as integers or reals
        read(kunit) ; read(kunit) label2 ; rewind kunit
        if (label2(80:80)==' ') then         ! integers, convert to seconds
          READ (kunit,err=10) it,XLABEL,nd1,iy1,iti1,ite1,it01,
     *        NTIM1,ITIM1(1:NTIM1),TSTR1(1:NTIM1)
          tim1(1:NTIM1) = 1d-2*itim1(1:NTIM1)
        else                                 ! real*8
          READ (kunit,err=10) it,XLABEL,nd1,iy1,iti1,ite1,it01,
     *        NTIM1,TIM1(1:NTIM1),TSTR1(1:NTIM1)
        end if
C**** use doc-record to check the basic model parameters
        READ (kunit,err=10) LABEL2
        READ (label2,'(13x,4i4)',err=10) im0,jm0,lm0,ls10
        if (im.ne.im0.or.jm.ne.jm0.or.lm.ne.lm0.or.ls10.ne.ls1) then
          ioerr = 0   ! warning
        end if
        SELECT CASE (IACTION)   ! set model common according to iaction
        CASE (ioread)           ! parameters from rundeck & restart file
          call read_param(kunit,.false.)
          nday=nd1 ; itimei=iti1      ! changeable only at new starts
          itimee=ite1 ; itime0=it01   ! are changed later if appropriate
          if (iyear1.lt.0) iyear1=iy1 ! rarely changes on restarts
          NTIMEACC=NTIM1
          TIMESTR(1:NTIM1)=TSTR1(1:NTIM1)
          TIMING(1:NTIM1)=TIM1(1:NTIM1)
        CASE (IRSFIC,irsficnt,IRSFICNO) ! rundeck/defaults except label
          read(kunit,err=10)          ! skip parameters, dates
          it=it*24/nd1                ! switch itime to ihour
        CASE (IRERUN)           ! parameters from rundeck & restart file
          call read_param(kunit,.false.)
          nday=nd1 ; itimei=iti1      ! changeable only at new starts
          itimee=ite1 ; itime0=it01   ! is changed later if appropriate
          if (iyear1.lt.0) iyear1=iy1 ! rarely changes on restart/reruns
        CASE (IOREAD_SINGLE)    ! parameters/label from 1-many acc files
          call read_param(kunit,.false.)  ! use rundeck
          call sync_param( "kradia",kradia)
          nday=nd1 ; iyear1=iy1 ; itime0=it01
          NTIMEACC=NTIM1                 ! use timing from current file
          TIMESTR(1:NTIM1)=TSTR1(1:NTIM1)
          TIMING(1:NTIM1)=TIMING(1:NTIM1)+TIM1(1:NTIM1)

C**** keep track of min/max time over the combined diagnostic period
          if (it.gt.ITmax)                   ITmax = it
          if (ITmin.lt.0 .or. it01.lt.ITmin) ITmin = it01
          itime0 = ITmin

        END SELECT
      END SELECT
      itm = max(it,ITmax)
      RETURN
 10   IOERR=1
      RETURN
      END SUBROUTINE io_label

