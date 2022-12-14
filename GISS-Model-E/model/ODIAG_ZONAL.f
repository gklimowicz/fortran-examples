#include "rundeck_opts.h"

      module odiag_zonal
!@sum odiag_zonal (per-basin) zonal-mean diagnostics of transport and state
      USE OCEAN, only : im,jm,lmo
      USE MDIAG_COM, only : sname_strlen,units_strlen,lname_strlen
      use cdl_mod

!@var do_odiag_zonal whether to calculate zonal-mean diags,
!@+   currently determined by presence of KBASIN file
      logical :: do_odiag_zonal=.false.

!@var kbasin integer index of which basin a particular ocean point is in
      INTEGER, DIMENSION(:,:), allocatable :: KBASIN,KBASIN_glob

!@var NBAS number of ocean basins 
!@var nqty number of output qtys zonally averaged over basins
      INTEGER, PARAMETER :: NBAS=4,nqty=3,NCIRC=3
!@var BASIN names of ocean basins for diag output
      CHARACTER*16, DIMENSION(NBAS) :: BASIN=
     *     (/"Atlantic","Pacific ","Indian  ","Global  "/)
      character(len=4), dimension(nqty), parameter :: qtyname=(/
     &     'Mass','Heat','Salt' /)
      character(len=9), dimension(nqty), parameter :: qtyflxunit=(/
     &     '10^9 kg/s','10^15 W  ','10^6 kg/s' /)
      character(len=3), dimension(ncirc), parameter :: circstr=(/
     &     'moc','gmf','gyr' /)
      character(len=7), dimension(ncirc), parameter :: circname=(/
     &     'Overtrn', 'GM flx ', 'Hor gyr'/)

!@var KOJL,OJL (number of qtys having) zonal sums/means over basins
      INTEGER, PARAMETER :: KOJL=5*NBAS
      REAL*8, DIMENSION(JM,LMO,NBAS,KOJL/NBAS) :: OJL
      REAL*8, DIMENSION(JM,LMO,KOJL) :: OJL_out

!@var JL_xxx indices for qtys in OJL
      INTEGER :: JL_M,JL_PT,JL_S,JL_SF,JL_SFB

!@var lname_ojl Long names for OJL diagnostics
      CHARACTER(len=lname_strlen), DIMENSION(KOJL) :: LNAME_OJL
!@var sname_ojl Short names for OJL diagnostics
      CHARACTER(len=sname_strlen), DIMENSION(KOJL) :: SNAME_OJL
!@var units_ojl Units for OJL diagnostics
      CHARACTER(len=units_strlen), DIMENSION(KOJL) :: UNITS_OJL
!@var ia_ojl IDACC numbers for OJL diagnostics
      INTEGER, DIMENSION(KOJL) :: IA_OJL
!@var denom_ojl denominators for OJL diagnostics
      INTEGER, DIMENSION(KOJL) :: DENOM_OJL
!@var scale_ojl scales for OJL diagnostics
      REAL*8, DIMENSION(KOJL) :: SCALE_OJL
!@var [jl]grid_ojl Grid descriptors for OJL diagnostics
      INTEGER, DIMENSION(KOJL) :: JGRID_OJL,LGRID_OJL

!@var NSEC number of lat/lon sections for diags
      INTEGER, PARAMETER :: NSEC=3
!@var SEC_LAT, SEC_LON lat/lon for sectional tracer profiles
      REAL*8, PARAMETER :: SEC_LAT(NSEC) = (/-64.,0.,48./),
     *     SEC_LON(NSEC) = (/-165.0,-30.,65./)
C**** OTJ is integrated flux array OTJ(LATITUDE,BASIN,KQ)
C****   KQ 1   Mass (kg)
C****      2   Heat (J)
C****      3   Salt (kg)
C****  OTJCOMP OTJ(LATITUDE,BASIN,COMP,KQ) (KQ = 2,3)
C**** COMP 1   advected by overturning
C****      2   flux from GM
C****      3   advected by horizontal gyres (residual)
C****
      REAL*8 OTJ(0:JM,4,3),OTJCOMP(0:JM,4,3,3)
      integer, parameter :: kotj=4*4*3

!@var OTJ_out reshaped combination of OTJ and OTJCOMP
      REAL*8 OTJ_out(JM,kotj)

!@var ia_otj IDACC numbers for OTJ diagnostics
      integer, dimension(kotj) :: ia_otj
!@var scale_otj scales for OTJ diagnostics
      real*8, dimension(kotj) :: scale_otj
!@var sname_otj short names for OTJ diagnostics
      character(len=sname_strlen), dimension(kotj) :: sname_otj
!@var lname_otj Long names for OTJ diagnostics
      character(len=lname_strlen), dimension(kotj) :: lname_otj
!@var units_otj units for OTJ diagnostics
      character(len=units_strlen), dimension(kotj) :: units_otj

!@var SFM meridional overturning stream function for each basin
      REAL*8, DIMENSION(JM,0:LMO,4) :: SFM!,SFS SFS is for salt

      type(cdl_type) :: cdl_ojl,cdl_otj

      end module odiag_zonal

      subroutine init_ODIAG_zonal
!@sum  init_ODIAG initialises ocean basin/zonal diagnostics
      USE MDIAG_COM, only : ia_src=>ia_cpl
      USE MDIAG_COM, only : make_timeaxis
      use odiag, only : cdl_olats,cdl_odepths
      use odiag_zonal
      use filemanager, only : file_exists
      IMPLICIT NONE
c
      INTEGER k,kb,kq,kc,kk,n
      character(len=10) :: ystr,zstr
      logical :: set_miss


      do_odiag_zonal = file_exists('KBASIN')
      if(.not.do_odiag_zonal) return

C**** Initialise ocean basins
      CALL OBASIN

C**** Metadata for northward transports
      scale_otj(:) = 1.
      ia_otj(:) = ia_src
      kk = 0
      do kq=1,3
      do kb=1,4
        kk = kk + 1
        sname_otj(kk)='nt_'//trim(qtyname(kq))//'_'//basin(kb)(1:3)
        lname_otj(kk)='North. Trans. of '//trim(qtyname(kq))//
     &       ' in the '//trim(basin(kb))//' basin'
        units_otj(kk)=qtyflxunit(kq)
        do kc=1,3
          kk = kk + 1
          sname_otj(kk)='nt_'//trim(qtyname(kq))//'_'//basin(kb)(1:3)//
     &       '_'//trim(circstr(kc))
          lname_otj(kk)='North. Trans. of '//trim(qtyname(kq))//
     &       ' in the '//trim(basin(kb))//' basin by '//
     &       trim(circname(kc))
          units_otj(kk)=qtyflxunit(kq)
        enddo
      enddo
      enddo

C**** set properties for OJL diagnostics
      do k=1,kojl
        sname_ojl(k) = 'unused'
        ia_ojl(k) = ia_src
        denom_ojl(k) = 0
        scale_ojl(k) = 1.
        lname_ojl(k) = 'no output'
        units_ojl(k) = 'no output'
        jgrid_ojl(k) = 1
        lgrid_ojl(k) = 1
      enddo
c
      k=0
      kk=0
c
      k=k+1
      JL_M = k
      do n=1,nbas
        kk = kk + 1
        sname_ojl(kk) = 'mo_'//basin(n)(1:3)
      enddo
c
      k=k+1
      JL_PT = k
      do n=1,nbas
        kk = kk + 1
        denom_ojl(kk) = JL_M +n-1
        sname_ojl(kk) = 'temp_'//basin(n)(1:3)
        units_ojl(kk) = 'C'
        lname_ojl(kk) = 'Temperature, '//trim(basin(n))//' Basin'
      enddo
c
      k=k+1
      JL_S = k
      do n=1,nbas
        kk = kk + 1
        denom_ojl(kk) = JL_M +n-1
        sname_ojl(kk) = 'salt_'//basin(n)(1:3)
        units_ojl(kk) = 'psu'
        lname_ojl(kk) = 'Salinity, '//trim(basin(n))//' Basin'
      enddo
c
      k=k+1
      JL_SF = k
      do n=1,nbas
        kk = kk + 1
        sname_ojl(kk) = 'sf_'//basin(n)(1:3)
        units_ojl(kk) = 'Sv'
        lname_ojl(kk) = 'Total Stream Function, '//
     &       trim(basin(n))//' Basin'
        jgrid_ojl(kk) = 2
        lgrid_ojl(kk) = 2
      enddo
c
      k=k+1
      JL_SFB = k
      do n=1,nbas
        kk = kk + 1
        sname_ojl(kk) = 'sfb_'//basin(n)(1:3)
        units_ojl(kk) = 'Sv'
        lname_ojl(kk) = 'Bolus SF, '//trim(basin(n))//' Basin'
        jgrid_ojl(kk) = 2
        lgrid_ojl(kk) = 2
      enddo

      call merge_cdl(cdl_olats,cdl_odepths,cdl_ojl)
      do k=1,kojl
        if(trim(sname_ojl(k)).eq.'unused') cycle
        if(trim(lname_ojl(k)).eq.'no output') cycle
        ystr='lato) ;'
        if(jgrid_ojl(k).eq.2) ystr='lato2) ;'
        zstr='(zoc,'
        if(lgrid_ojl(k).eq.2) zstr='(zoce,'
        set_miss = denom_ojl(k).ne.0
        call add_var(cdl_ojl,
     &       'float '//trim(sname_ojl(k))//trim(zstr)//trim(ystr),
     &       long_name=trim(lname_ojl(k)),
     &       units=trim(units_ojl(k)),
     &       set_miss=set_miss,
     &       make_timeaxis=make_timeaxis)
      enddo

      cdl_otj = cdl_olats
      do k=1,kotj
        if(trim(sname_otj(k)).eq.'unused') cycle
        call add_var(cdl_otj,
     &       'float '//trim(sname_otj(k))//'(lato2) ;',
     &       long_name=trim(lname_otj(k)),
     &       units=trim(units_otj(k)),
     &       make_timeaxis=make_timeaxis)
      enddo

      end subroutine init_ODIAG_zonal

      subroutine def_meta_ocdiag_zonal(fid)
!@sum  def_meta_ocdiag defines basin/zonal diag metadata in ocean acc files
      use odiag_zonal
      use pario, only : defvar,write_attr
      USE OCEANR_DIM, only : grid=>ogrid
      use cdl_mod, only : defvar_cdl
      implicit none
      integer :: fid         !@var fid file id

      if(.not.do_odiag_zonal) return

      call defvar(grid,fid,ojl_out,'ojl(jmo,lmo,kojl)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'ojl','reduction','sum')
      call write_attr(grid,fid,'ojl','split_dim',3)
      call defvar(grid,fid,ia_ojl,'ia_ojl(kojl)')
      call defvar(grid,fid,denom_ojl,'denom_ojl(kojl)')
      call defvar(grid,fid,scale_ojl,'scale_ojl(kojl)')
      call defvar(grid,fid,sname_ojl,'sname_ojl(sname_strlen,kojl)')
      call defvar_cdl(grid,fid,cdl_ojl,
     &     'cdl_ojl(cdl_strlen,kcdl_ojl)')

      call defvar(grid,fid,otj_out,'otj(jmo,kotj)',
     &     r4_on_disk=.true.)
      call write_attr(grid,fid,'otj','reduction','sum')
      call write_attr(grid,fid,'otj','split_dim',2)
      call defvar(grid,fid,ia_otj,'ia_otj(kotj)')
      call defvar(grid,fid,scale_otj,'scale_otj(kotj)')
      call defvar(grid,fid,sname_otj,'sname_otj(sname_strlen,kotj)')
      call defvar_cdl(grid,fid,cdl_otj,
     &     'cdl_otj(cdl_strlen,kcdl_otj)')

      end subroutine def_meta_ocdiag_zonal

      subroutine write_meta_ocdiag_zonal(fid)
!@sum  write_meta_ocdiag write ocean basin/zonal accumulations to file
      use odiag_zonal
      use pario, only : write_dist_data,write_data
      USE OCEANR_DIM, only : grid=>ogrid
      use cdl_mod, only : write_cdl
      implicit none
      integer :: fid         !@var fid file id

      if(.not.do_odiag_zonal) return

      call write_data(grid,fid,'ojl',ojl_out)
      call write_data(grid,fid,'ia_ojl',ia_ojl)
      call write_data(grid,fid,'denom_ojl',denom_ojl)
      call write_data(grid,fid,'scale_ojl',scale_ojl)
      call write_data(grid,fid,'sname_ojl',sname_ojl)
      call write_cdl(grid,fid,'cdl_ojl',cdl_ojl)

      call write_data(grid,fid,'otj',otj_out)
      call write_data(grid,fid,'ia_otj',ia_otj)
      call write_data(grid,fid,'scale_otj',scale_otj)
      call write_data(grid,fid,'sname_otj',sname_otj)
      call write_cdl(grid,fid,'cdl_otj',cdl_otj)

      end subroutine write_meta_ocdiag_zonal

      SUBROUTINE STRMJL (OIJL,FAC,OLNST,FACST,SF)
!@sum  STRMJL calculates the latitude by layer stream function
!@auth G. Russell/G. Schmidt
C****
C**** Input:
!@var OIJL  = west-east and south-north tracer fluxes (kg/s)
!@var FAC   = global scaling factor
!@var OLNST = strait mass flux (kg/s)
!@var FACST = global scaling factor for straits
C**** Output:
!@var    SF = stream function (kg/s)
C****
      USE CONSTANT, only : undef
      USE OCEAN, only : im,jm,lmo
      USE STRAITS, only : nmst
      USE ODIAG_zonal, only : kbasin,kbasin_glob
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : pack_dataj,am_i_root
      IMPLICIT NONE
      REAL*8, INTENT(IN) :: FAC,FACST
      REAL*8, INTENT(IN) ::
     &     OIJL(IM,grid%j_strt_halo:grid%j_stop_halo,LMO,2)
      REAL*8, INTENT(IN), DIMENSION(LMO,NMST) :: OLNST
      REAL*8, INTENT(OUT), DIMENSION(JM,0:LMO,4) :: SF
      REAL*8, DIMENSION(4) :: SUMB
      INTEGER :: KB,I,J,L,K
      REAL*8 SFx(grid%j_strt_halo:grid%j_stop_halo,0:LMO,4)
      REAL*8 :: SF_str(JM,0:LMO,4)
C****
C**** Zero out the Stream Function at the ocean bottom
C****
      SFX(:,:,:) = 0.
C****
C**** Integrate the Stream Function upwards from the ocean bottom
C****
      DO L=LMO-1,0,-1
      DO J=grid%j_strt,min(grid%j_stop,jm-1)
        DO KB=1,4
          SFX(J,L,KB) = SFX(J,L+1,KB)
        END DO
        DO I=1,IM
          KB = KBASIN(I,J+1)
          IF(KB.gt.0) SFX(J,L,KB) = SFX(J,L,KB) + OIJL(I,J,L+1,2)*FAC
        END DO
      END DO
      END DO

      call pack_dataj(grid,sfx,sf)

      if(am_i_root()) then
C****
C**** Add strait flow to the Stream Function
C****
      sf_str = 0.
      DO L=LMO-1,0,-1
        CALL STRMJL_STRAITS(L,SF_str,OLNST,FACST)
        sf_str(:,l,:) = sf_str(:,l,:) + sf_str(:,l+1,:)
      END DO
      sf = sf + sf_str

C****
C**** Calculate global Stream Function by summing it over 3 oceans
C****
      SF(:,:,4) = SF(:,:,1) + SF(:,:,2) + SF(:,:,3) + SF(:,:,4)
C****
C**** Mask streamfunction so that topography is put to skip value
C****
      DO J=1,JM-1
        SUMB = 0
        DO I=1,IM
          IF (KBASIN_glob(I,J).gt.0) SUMB(KBASIN_glob(I,J))=1.
        END DO
        SUMB(4)=SUMB(1)+SUMB(2)+SUMB(3)+SUMB(4)
        DO K=1,4
          IF (SUMB(K).eq.0) THEN
            SF(J,0:LMO,K) = UNDEF   ! UNDEF areas of no integral
          ELSE
            DO L=LMO,1,-1
              IF(SF(J,L,K).eq.0d0 .and. SF(J,L-1,K).eq.0d0)
     &             SF(J,L,K) = UNDEF
              IF(SF(J,L-1,K).ne.0d0) EXIT
            ENDDO
          END IF
        END DO
      END DO
      endif ! am_i_root
C****
      RETURN
      END

      SUBROUTINE OBASIN
!@sum  OBASIN Read in KBASIN: 0=continent,1=Atlantic,2=Pacific,3=Indian
      USE Constant, only : UNDEF_VAL
      USE OCEAN, only : IM,JM,focean
      USE ODIAG_zonal, only : kbasin,kbasin_glob
      use pario, only : par_open,par_close,read_dist_data
      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : pack_data,halo_update,am_i_root
      IMPLICIT NONE
      INTEGER i,j,k,fid, j_0h,j_1h
      character(len=3), parameter :: basins(3)=(/'atl','pac','ind'/)
      real*8, dimension(im,grid%j_strt_halo:grid%j_stop_halo) ::
     &     zeroone
C****
C**** read in basin data
C****

      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo


      allocate(kbasin(im,j_0h:j_1h))
      if(am_i_root()) then
        allocate(kbasin_glob(im,jm))
      endif

      fid = par_open(grid,'KBASIN','read')
      kbasin = 0
      zeroone = UNDEF_VAL
      do k=1,3
        call read_dist_data(grid,fid,'mask_'//basins(k),zeroone)
        call halo_update(grid,zeroone)
        where(zeroone==1d0) kbasin = k
      enddo
      call par_close(grid,fid)

      do j=max(1,j_0h),min(jm,j_1h) !j=j_0h,j_1h
      do i=1,im
        if(focean(i,j).gt.0. .and. kbasin(i,j).eq.0) then
          kbasin(i,j) = 4
        endif
      enddo
      enddo

      call pack_data(grid,kbasin,kbasin_glob)

      RETURN
      END SUBROUTINE OBASIN

      subroutine basin_prep
c
c Calculate zonal sums for ocean basins
c
      use constant, only : undef
      use ocean, only : im,jm,lmo,imaxj,focean,dxypo
      USE OCEAN, only : dts
      USE STRAITS, only : nmst

      use odiag, only : oijl=>oijl_loc
     &     ,ijl_mo,ijl_g0m,ijl_s0m,ijl_mfu,ijl_mfv,ijl_mfvb
     &     ,ijl_gflx,ijl_sflx
     &     ,ijl_ggmfl,ijl_sgmfl
     &     ,olnst,ln_mflx,ln_gflx,ln_sflx

      use odiag_zonal, only : do_odiag_zonal,kbasin
     &     ,ojl,ojl_out,nbas,nqty,jl_m,jl_pt,jl_s,jl_sf,jl_sfb
     &     ,otj,otjcomp,otj_out
     &     ,sfm

      use oceanr_dim, only : grid=>ogrid
      use domain_decomp_1d, only : pack_dataj,am_i_root
      implicit none
      integer i,j,l,kb
      real*8 gos,sos,temgs
      integer :: j_0,j_1, j_0h,j_1h
      real*8 :: ojlx(grid%j_strt_halo:grid%j_stop_halo,lmo,nbas,4)
      INTEGER K,KQ,NOIJL,NOIJLGM,NOLNST,NS,KK,KC
      REAL*8 SOIJL,SOIJLGM,FAC,FACST
      REAL*8 SCALEM(3),SCALES(3),SOLNST(NMST),MV(4),MT(4)
      INTEGER IS(4)
      REAL*8 ::
     &     X(grid%j_strt_halo:grid%j_stop_halo,4,3)
     &    ,XCOMP(grid%j_strt_halo:grid%j_stop_halo,4,3,3)
      real*8 dumarr(1,1)

      if(.not.do_odiag_zonal) return

      j_0 = grid%j_strt
      j_1 = grid%j_stop
      j_0h = grid%j_strt_halo
      j_1h = grid%j_stop_halo

      ojlx = 0.
      do l=1,lmo
      do j=j_0,j_1
      do i=1,imaxj(j)
        if(focean(i,j).lt..5) cycle
        kb=kbasin(i,j)
        ojlx(j,l,kb,jl_m) = ojlx(j,l,kb,jl_m) +oijl(i,j,l,ijl_mo)
        ojlx(j,l,kb,jl_pt) = ojlx(j,l,kb,jl_pt) +oijl(i,j,l,ijl_g0m)
        ojlx(j,l,kb,jl_s) = ojlx(j,l,kb,jl_s) +oijl(i,j,l,ijl_s0m)
      enddo
      enddo
      enddo

      ojlx(:,:,4,:) = sum(ojlx(:,:,1:4,:),dim=3)

      do kb=1,4
      do l=1,lmo
      do j=j_0,j_1
        ojlx(j,l,kb,jl_m) = ojlx(j,l,kb,jl_m)*dxypo(j)
        if(ojlx(j,l,kb,jl_m).eq.0.) cycle
        gos = ojlx(j,l,kb,jl_pt)/ojlx(j,l,kb,jl_m)
        sos = ojlx(j,l,kb,jl_s )/ojlx(j,l,kb,jl_m)
        ojlx(j,l,kb,jl_pt) = temgs(gos,sos)*ojlx(j,l,kb,jl_m)
        ojlx(j,l,kb,jl_s) = ojlx(j,l,kb,jl_s)*1d3
      enddo
      enddo
      enddo

      call pack_dataj(grid,ojlx,ojl)

C****
C**** Calculate Mass Stream Function
C****
      FAC   = -1d-9/DTS
      FACST = -1d-9/DTS
      if(nmst.gt.0) then
        CALL STRMJL (OIJL(1,j_0h,1,IJL_MFU),FAC,
     &       OLNST(1,1,LN_MFLX),FACST,SFM)
      else
        CALL STRMJL (OIJL(1,j_0h,1,IJL_MFU),FAC,
     &       DUMARR,FACST,SFM)
      endif
      if(am_i_root()) then
        ojl(1,:,:,jl_sf) = 0.
        ojl(2:jm,:,:,jl_sf) = sfm(1:jm-1,1:lmo,:)
      endif

C****
C**** Zonally sum the bolus streamfunction and add it to
C**** the resolved streamfunction
C****
      FAC   = -1d-9/DTS
      FACST = 0.
      if(nmst.gt.0) then
        CALL STRMJL (OIJL(1,j_0h,1,IJL_MFVB-1),FAC,
     &       OLNST(1,1,LN_MFLX),FACST,SFM)
      else
        CALL STRMJL (OIJL(1,j_0h,1,IJL_MFVB-1),FAC,
     &       DUMARR,FACST,SFM)
      endif
      if(am_i_root()) then
        ojl(1,:,:,jl_sfb) = 0.
        ojl(2:jm,:,:,jl_sfb) = sfm(1:jm-1,1:lmo,:)
        where( ojl(2:jm,:,:,jl_sf) .ne. undef )
     &  ojl(2:jm,:,:,jl_sf ) = sfm(1:jm-1,1:lmo,:)
     & +ojl(2:jm,:,:,jl_sf )
      endif

c copy to j,l,n array
      if(am_i_root()) ojl_out = reshape(ojl,shape(ojl_out))


c
c Transports
c

C****
C****  N          Contents of OIJL(I,J,L,N)       Scaling factor
C****  -          -------------------------       --------------
C**** IJL_MV      South-north mass flux (kg/s)    2/IDACC(1)*NDYNO
C**** IJL_GFLX+1  South-north heat flux (J/s)     1/IDACC(1)*DTS
C**** IJL_SFLX+1  South-north salt flux (kg/s)    1/IDACC(1)*DTS
C****
C****  N       Contents of OLNST(L,N)             Scaling factor
C****  -       ----------------------             --------------
C**** LN_MFLX  Strait mass flux (kg/s)            1/IDACC(1)*DTS
C**** LN_GFLX  Strait heat flux (J/s)            .5/IDACC(1)*DTS
C**** LN_SFLX  Strait salt flux (kg/s)           .5/IDACC(1)*DTS
C****
      SCALEM(1) = 1.D- 9/(DTS)
      SCALEM(2) = 1.D-15/(DTS)
      SCALEM(3) = 1.D- 6/(DTS)
      SCALES(1) = 1.D- 9/(DTS)
      SCALES(2) = .5D-15/(DTS)
      SCALES(3) = .5D- 6/(DTS)
C****
      X = 0. ; XCOMP =0.
C****
C**** Loop over quantites
C****
      DO KQ=1,3
C****
C**** Calculate the transport entering each domain from the
C**** southern boundary
C****
      SELECT CASE (KQ)
      CASE (1)
        NOIJL=IJL_MFV
        NOIJLGM=0
      CASE (2)
        NOIJL=IJL_GFLX+1
        NOIJLGM=IJL_GGMFL+1
      CASE (3)
        NOIJL=IJL_SFLX+1
        NOIJLGM=IJL_SGMFL+1
      END SELECT

      DO J=J_0,min(J_1,jm-1)
        DO I=1,IM
          IF(OIJL(I,J,1,IJL_MFV).eq.0) CYCLE
          KB = KBASIN(I,J+1)
          SOIJL = 0.
          SOIJLGM = 0.
          DO L=1,LMO
            SOIJL = SOIJL + OIJL(I,J,L,NOIJL)
            IF (NOIJLGM.gt.0) SOIJLGM = SOIJLGM + OIJL(I,J,L,NOIJLGM)
          END DO
          X(J,KB,KQ) = X(J,KB,KQ) + (SOIJL+SOIJLGM)*SCALEM(KQ)
          IF (KQ.ne.1) THEN
c        X(J,KB,5) = X(J,KB,5) + SOIJLGM*SCALEM(KQ)
c        X(J, 4,5) = X(J, 4,5) + SOIJLGM*SCALEM(KQ)
            XCOMP(J,KB,2,KQ) = XCOMP(J,KB,2,KQ) + SOIJLGM*SCALEM(KQ)
          END IF
        END DO
        X(J,4,KQ) = SUM(X(J,1:4,KQ))
        XCOMP(J,4,2,KQ) = SUM(XCOMP(J,1:4,2,KQ))
      END DO
C****
C**** Accumulate northward transports by overturning and by
C**** horizontal gyre
C**** 
      DO J=J_0,min(J_1,jm-1)
        DO L=1,LMO
          IS=0 ; MV=0 ; MT=0
          DO I=1,IM
            IF(OIJL(I,J,L,IJL_MFV).eq.0)  CYCLE
            KB = KBASIN(I,J+1)
            IS(KB) = IS(KB) + 1
            MV(KB) = MV(KB) + OIJL(I,J,L,IJL_MFV)
            IF (NOIJL.gt.0) MT(KB) = MT(KB) + OIJL(I,J,L,NOIJL)/OIJL(I,J
     *           ,L,IJL_MFV)
          END DO
          IS(4) = SUM(IS(1:4))
          MV(4) = SUM(MV(1:4))
          IF (NOIJL.gt.0) MT(4) = SUM(MT(1:4))
          DO KB=1,4
            IF(IS(KB).eq.0) CYCLE
c     X(J,KB,4) = X(J,KB,4) + SCALEM(2)*MV(KB)*MT(KB)/IS(KB)
            IF (KQ.eq.1) THEN
              XCOMP(J,KB,1,KQ) = XCOMP(J,KB,1,KQ) +SCALEM(KQ)*MV(KB)
            ELSE
              XCOMP(J,KB,1,KQ) = XCOMP(J,KB,1,KQ) +SCALEM(KQ)*MV(KB)
     *             *MT(KB)/IS(KB)
            END IF
          END DO
        END DO
        DO KB=1,4
c     X(J,KB,6) = X(J,KB,2) - X(J,KB,4) - X(J,KB,5)
          XCOMP(J,KB,3,KQ)=
     &         X(J,KB,KQ)-XCOMP(J,KB,1,KQ)-XCOMP(J,KB,2,KQ)
        END DO
      END DO

      END DO ! end loop over kq

      call pack_dataj(grid,x,otj(1:jm,:,:))
      call pack_dataj(grid,xcomp,otjcomp(1:jm,:,:,:))

      if(am_i_root()) then
      DO KQ=1,3
C****
C**** Calculate transport through straits from latitude to another
C**** within the same basin
C****
      SELECT CASE (KQ)
      CASE (1)
        NOLNST=LN_MFLX
      CASE (2)
        NOLNST=LN_GFLX
      CASE (3)
        NOLNST=LN_SFLX
      END SELECT

      DO NS=1,NMST
        SOLNST(NS) = 0.
        DO L=1,LMO
          SOLNST(NS) = SOLNST(NS) + OLNST(L,NS,NOLNST)
        END DO
      END DO

      CALL OTJ_STRAITS(OTJ,SOLNST,SCALES(KQ),KQ)

C**** Fill in south polar value
      OTJ(0,1:4,KQ)  = OTJ(1,1:4,KQ)

      END DO ! end loop over kq

C**** Replace salt by salt - .035 times mass
      DO KB=1,4
        DO J=0,JM
          OTJ(J,KB,3) = OTJ(J,KB,3) - .035d3*OTJ(J,KB,1)
          OTJCOMP(J,KB,:,3) = OTJCOMP(J,KB,:,3)-.035d3*OTJCOMP(J,KB,:,1)
        END DO
      END DO

C**** combine OTJ and OTJCOMP into a single output array
      otj_out(1,:) = 0.
      kk = 0
      do kq=1,3
      do kb=1,4
        kk = kk + 1
        otj_out(2:jm,kk) = otj(1:jm-1,kb,kq)
        do kc=1,3
          kk = kk + 1
          otj_out(2:jm,kk) = otjcomp(1:jm-1,kb,kc,kq)
        enddo
      enddo
      enddo

      endif ! am_i_root

      return
      end subroutine basin_prep

      Subroutine STRMJL_STRAITS (L,SF,OLNST,FACST)
C****
C**** Add strait flow to JxL Strean Function
C****
C**** Input:   L = ocean model layer index
C****      OLNST = ocean strait mass flux (kg/s)
C****      FACST = global scaling factor for ocean straits
C**** Output: SF = JxL stream function (kg/s)
C****
      Use OCEAN,   Only: JM,LMO
      Use STRAITS, Only: NMST, IST,JST
      Use ODIAG_zonal,   Only: KBASIN=>kbasin_glob
      Implicit None

      Integer*4,Intent(In) :: L
      Real*8,Intent(InOut) :: SF(JM,0:LMO,4)
      Real*8,Intent(In)    :: OLNST(LMO,NMST),FACST

C**** Local variables
      Integer*4 N, I1,J1,K1, I2,J2,K2

      Do 30 N=1,NMST
      I1 = IST(N,1)       ;  I2 = IST(N,2)
      J1 = JST(N,1)       ;  J2 = JST(N,2)
      K1 = KBASIN(I1,J1)  ;  K2 = KBASIN(I1,J1)
      If (J2 - J1) 10,30,20

C**** JST(N,2) < JST(N,1)
   10 If (K1 == K2)  Then
         SF(J2:J1-1,L,K1) = SF(J2:J1-1,L,K1) - OLNST(L+1,N)*FACST
      Else
         SF(J2:J1-1,L,K1) = SF(J2:J1-1,L,K1) - OLNST(L+1,N)*FACST*.5
         SF(J2:J1-1,L,K2) = SF(J2:J1-1,L,K2) - OLNST(L+1,N)*FACST*.5
      EndIf
      GoTo 30

C**** JST(N,2) > JST(N,1)
   20 If (K1 == K2)  Then
         SF(J1:J2-1,L,K1) = SF(J1:J2-1,L,K1) + OLNST(L+1,N)*FACST
      Else
         SF(J1:J2-1,L,K1) = SF(J1:J2-1,L,K1) + OLNST(L+1,N)*FACST*.5
         SF(J1:J2-1,L,K2) = SF(J1:J2-1,L,K2) + OLNST(L+1,N)*FACST*.5
      EndIf
   30 Continue
      Return
      EndSubroutine STRMJL_STRAITS

      Subroutine OTJ_STRAITS (X,SOLNST,SCALE,KQ)
C****
C**** Calculate transport through straits from one latitude to another
C****
C**** Input: SOLNST = vertically integrated flux through each strait
C****         SCALE = strait scaling factor
C****            KQ = 1: mass flux; 2: heat flux; 3: salt flux
C**** Output:     X = northward flux as a function of J and basin
C****
      Use OCEAN,   Only: JM
      Use STRAITS, Only: NMST, IST,JST
      Use ODIAG_zonal,   Only: KBASIN=>kbasin_glob
      Implicit None

      Real*8,Intent(InOut) :: X(0:JM,4,3)
      Real*8,Intent(In)    :: SOLNST(NMST),SCALE
      Integer*4,Intent(In) :: KQ

C**** Local variables
      Integer*4 N, I1,J1,K1, I2,J2,K2

      Do 30 N=1,NMST
      I1 = IST(N,1)       ;  I2 = IST(N,2)
      J1 = JST(N,1)       ;  J2 = JST(N,2)
      K1 = KBASIN(I1,J1)  ;  K2 = KBASIN(I2,J2)
      If (J2 - J1) 10,30,20

C**** JST(N,2) < JST(N,1)
   10 If (K1 == K2)  Then
        X(J2:J1-1,K1,KQ) = X(J2:J1-1,K1,KQ) - SOLNST(N)*SCALE
      Else
        X(J2:J1-1,K1,KQ) = X(J2:J1-1,K1,KQ) - SOLNST(N)*SCALE*.5
        X(J2:J1-1,K2,KQ) = X(J2:J1-1,K2,KQ) - SOLNST(N)*SCALE*.5
      Endif
        X(J2:J1-1, 4,KQ) = X(J2:J1-1, 4,KQ) - SOLNST(N)*SCALE
      GoTo 30

C**** JST(N,2) > JST(N,1)
   20 If (K1 == K2)  Then
        X(J1:J2-1,K1,KQ) = X(J1:J2-1,K1,KQ) + SOLNST(N)*SCALE
      Else
        X(J1:J2-1,K1,KQ) = X(J1:J2-1,K1,KQ) + SOLNST(N)*SCALE*.5
        X(J1:J2-1,K2,KQ) = X(J1:J2-1,K2,KQ) + SOLNST(N)*SCALE*.5
      Endif
        X(J1:J2-1, 4,KQ) = X(J1:J2-1, 4,KQ) + SOLNST(N)*SCALE
   30 Continue
      Return
      EndSubroutine OTJ_STRAITS
