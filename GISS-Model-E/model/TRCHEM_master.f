#include "rundeck_opts.h"
      SUBROUTINE masterchem
!@sum masterchem main chemistry routine
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!@calls fastj2_drv,Crates,Oxfam,HOxfam,NOxfam,chemstep
C
C**** GLOBAL parameters and variables:
c
!!    use precision_mod, only : reduce_precision 
      USE Dictionary_mod, only : get_param, is_set_param
      USE SOMTQ_COM, only   : qmom
      USE DOMAIN_DECOMP_1D, only : PACK_DATA ! for DU_O3
      USE DOMAIN_DECOMP_ATM,only: GRID,getDomainBounds,AM_I_ROOT,
     &                        GLOBALSUM,GLOBALMAX,
     &                        write_parallel,writet8_column,
     &                        writet_parallel
      USE RESOLUTION, only  : ls1=>ls1_nominal,plbot
      USE RESOLUTION, only  : IM,JM
      USE ATM_COM, only     : T,Q
      use model_com, only: modelEclock
      use model_com, only: itime, itimeI, itime0
      use TimeConstants_mod, only: HOURS_PER_DAY
      USE TRACER_COM, only  : ntm
      USE TRACER_COM, only  : COUPLED_CHEM
      USE RAD_COM, only     : o2x
      USE CONSTANT, only    : radian,gasc,mair,mb2kg,pi,avog,rgas,pO2,
     &                        bygrav,lhe,undef,teeny,byavog
      USE ATM_COM, only     : pedn,PMIDL00,LTROPO
      USE FILEMANAGER, only : openunit,closeunit,nameunit
      USE RAD_COM, only     : COSZ1,alb,rcloudfj=>rcld,
     &                        rad_to_chem,chem_tracer_save,H2ObyCH4,
     &                        SRDN,rad_to_file,ghg_yr,clim_interact_chem
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
     &                        ,stratO3_tracer_save
#endif
      USE GEOM, only        : BYAXYP, AXYP, LAT2D_DG, IMAXJ, LAT2D
      USE FLUXES, only      : tr3Dsource
      use OldTracer_mod, only: tr_wd_type, nWater
      USE TRACER_COM, only  : ntm_chem_beg, ntm_chem_end
      USE TRACER_COM, only  : n_Ox,n_NOx,n_N2O5,n_HNO3,n_H2O2,
     &                      n_HCHO,n_HO2NO2,n_CO,n_CH4,
     &                      n_Isoprene,n_AlkylNit,n_Alkenes,n_stratOx,
     &                      n_Terpenes,n_SO4,n_H2O2_s,oh_live,no3_live,
     &                      ntm_chem,n_DMS,n_MSA,n_SO2,
     &                      trm,trmom,nChemistry,nOverwrite,
     &                      n_isopp1g,n_isopp1a,n_isopp2g,n_isopp2a,
     &                      n_apinp1g,n_apinp1a,n_apinp2g,n_apinp2a,
     &                      rsulf1,rsulf2,rsulf3,rsulf4,
     &                      n_HBr,n_HOCl,n_HCl,n_ClONO2,n_ClOx,
     &                      n_BrOx,n_BrONO2,n_CFC,n_N2O,n_HOBR
#ifdef TRACERS_dCO
     &                     ,n_d13Calke
     &                     ,n_dHCH17O,n_dHCH18O,n_dH13CHO
     &                     ,n_dC17O,n_dC18O,n_d13CO
      use tracers_dCO, only: dacetone_fact, dalke_IC_fact
#endif  /* TRACERS_dCO */
#ifdef TRACERS_AMP
      USE TRACER_COM, only  : n_M_AKK_SU,n_M_ACC_SU,n_M_DD1_SU,
     &                        n_M_DS1_SU,n_M_DD2_SU,n_M_DS2_SU,
     &                        n_M_SSA_SU,n_M_OCC_SU,n_M_BC1_SU,
     &                        n_M_BC2_SU,n_M_BC3_SU,n_M_DBC_SU,
     &                        n_M_BOC_SU,n_M_BCS_SU,n_M_MXX_SU
#endif
#ifdef TRACERS_TOMAS
      USE TRACER_COM, only  : n_ASO4,nbins
#endif  
      use OldTracer_mod, only: tr_mm, mass2vol, vol2mass, trname
#ifdef TRACERS_HETCHEM
      USE TRACER_COM, only  : 
     &                        krate,n_N_d1,n_N_d2,n_N_d3
#endif
#ifdef INTERACTIVE_WETLANDS_CH4
      USE TRACER_SOURCES, only: avg_model,n__sw
#endif
      USE TRDIAG_COM, only    : taijs=>taijs_loc,taijls=>taijls_loc
     &     ,ijlt_NO3,jls_COp,jls_COd,jls_Oxp,jls_N2O5sulf,jls_O3vmr
     &     ,jls_Oxd,jls_OxpT,jls_OxdT,ijs_NO2_1030,ijs_NO2_1030c
     &     ,ijlt_COp,ijlt_COd,ijlt_Oxd,ijlt_Oxp,ijlt_pO1D,ijs_O3mass
     &     ,ijlt_pOH,ijlt_OxpHO2,ijlt_OxpCH3O2,ijlt_OxlHO2,ijlt_OxlALK
     &     ,ijlt_OxlOH,ijs_NO2_1330,ijs_NO2_1330c,ijlt_NO2vmr,ijlt_NOvmr
     &     ,ijlt_JO1D,ijlt_JNO2,ijlt_JH2O2,ijlt_O3ppbv,ijlt_O3cmatm
     &     ,jls_ClOcon,jls_H2Ocon
      USE TRCHEM_Shindell_COM
#ifdef TRACERS_AEROSOLS_SOA
      USE TRACERS_SOA, only: soa_aerosolphase,voc2nox,soa_apart,
     &                       whichsoa,apartmolar,LM_soa
#endif  /* TRACERS_AEROSOLS_SOA */
      use zonalmean_mod, only : zonalmean_ij2ij

      use TRACER_COM, only: nn_CH4,  nn_N2O, nn_Ox,   nn_NOx, 
     &      nn_N2O5,   nn_HNO3,  nn_H2O2,  nn_HCHO,
     &      nn_HO2NO2, nn_H2O17,             
     &      nn_Isoprene, nn_AlkylNit, nn_Alkenes,
     &      nn_stratOx, nn_Terpenes,
     &      nn_isopp1g,nn_isopp1a,nn_isopp2g,nn_isopp2a,         
     &      nn_apinp1g,nn_apinp1a,nn_apinp2g,nn_apinp2a,         
     &      nn_ClOx,   nn_BrOx,  nn_HCl,   nn_HOCl,   nn_ClONO2,  
     &      nn_HBr,    nn_HOBr,  nn_BrONO2,nn_CFC,    nn_GLT
#ifdef TRACERS_dCO
     &     ,nn_d13Calke
     &     ,nn_dHCH17O,nn_dHCH18O,nn_dH13CHO
#endif  /* TRACERS_dCO */
#ifdef CACHED_SUBDD
      use subdd_mod, only : subdd_groups,subdd_type,subdd_ngroups
     &     ,inc_subdd,find_groups, LmaxSUBDD
#endif
      use photolysis, only: fastj2_drv,o3_fastj,rj
     &                     ,sza,szamax,zj,jpnl,sf3_fact,sf2_fact

      IMPLICIT NONE

#ifdef CACHED_SUBDD
      integer :: igrp,ngroups,grpids(subdd_ngroups)
      type(subdd_type), pointer :: subdd
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,
     &                  LM) :: mrno,mrno2,mro3,OH_conc,HO2_conc,
     &                         JO1D_rate,JNO2_rate
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: O3col
#endif
C**** Local parameters and variables and arguments:
!@param by35 1/35 used for spherical geometry constant
      REAL*8, PARAMETER  :: by35=1.d0/35.d0
      REAL*8, PARAMETER  :: bymair = 1.d0/mair
!@var FASTJ_PFACT temp factor for vertical pressure-weighting
!@var FACT1,2,3 temp variable for strat overwrite
!@var bydtsrc reciprocal of the timestep dtsrc
!@var local logical for error checking 
!@var PRES2 local nominal pressure for verticle interpolations
!@var thick thickness of each layer in various units
!@var ClTOT total chlorine in all forms (reactive and reservoir)
!@var BrTOT total bromine in all forms (reactive and reservoir)
!@var colmO2, colmO3 are overhead oxygen and ozone columns
!@var CH4FACT, r179 for setting CH4 ICs and strat distribution
!@var changeClONO2,changeClOx,changeHOCl,changeHCl nighttime changes
!@+ also changeBrOx,changeBrONO2,changeBrOx2,changeHBr
!@var changehetClONO2 nighttime het change in ClONO2 (on sulfate)
!@var chgHT3,chgHT4,chgHT5 reaction rates for het rxns on pscs
!@var rmrClOx,rmrBrOx dummy vars with mixing ratios of halogens
!@var rmv dummy variable for halogne removal in trop vs height
!@var PIfact strat-overwrite scaling
!@var pfactor to convert units on species chemical changes
!@var bypfactor to convert units on species chemical changes
!@var dNO3,gwprodHNO3,gwprodN2O5,changeAldehyde,
!@+   changeAlkenes,changeIsoprene,changeHCHO,changeAlkylNit,
!@+   changeTerpenes,changeisopp1g,changeisopp2g,changeapinp1g,
!@+   changeapinp2g,changeHNO3,changeNOx,changeN2O5,wprodHCHO
!@+   working variables to calculate nighttime chemistry changes
!@var rlossN,rprodN,ratioN variables for nitrogen conservation
!@var I,J,L,N,igas,inss,LL,L2,n2 dummy loop variables
!@var avgTT_CH4 Itime avg CH4 # density at LTROPO between 20N and 20S
!@var avgTT_H2O Itime avg H2O # density at LTROPO between 20N and 20S
!@var countTT # of points between 20N and 20S on LTROPO plane
!@var maxT chosen tropopause 0=LTROPO(I,J), 1=LS1-1. Or the top layer
!@+ of chemistry in the unlikely event that is lower. Note in that case,
!@+ loops like L=maxT+1,topLevelOfChemistry will do nothing.
!@var sumOx for summing regional Ox tracers
!@var bysumOx reciprocal of sum of regional Ox tracers
      REAL*8, DIMENSION(NTM) :: PIfact
      REAL*8, DIMENSION(LM) :: PRES2 ! keep LM; based on PMIDL00(:)
      REAL*8 :: FACT1,FACT2,FACT3,FACT4,FACT5,FACT6,FACT7,fact_so4,
     &  FASTJ_PFACT,bydtsrc,CH4FACT,r179,rlossN,
     &  rprodN,ratioN,pfactor,bypfactor,gwprodHNO3,
#ifdef TRACERS_dCO
     &  changed17Oald,changed18Oald,changed13Cald,
     &  rd17OaldplusNO3,rd18OaldplusNO3,rd13CaldplusNO3,
     &  gwprodHNO3dHCH17O,gwprodHNO3dHCH18O,gwprodHNO3dH13CHO,
#endif  /* TRACERS_dCO */
     &  gwprodN2O5,wprod_sulf,wprodCO,dNO3,wprodHCHO,prod_sulf,
     &  RVELN2O5,changeAldehyde,changeAlkenes,changeAlkylNit,
     &  changeIsoprene,changeHCHO,changeHNO3,changeNOx,changeN2O5,
     &  rNO3plusNO2,rN2O5decomp,rHCHOplusNO3,rAldplusNO3,rIsopplusNO3,
     &  rClOplusNO2,rDMSplusNO3,wlossN2O5,wlossNOx,
     &  wprodNOx,rlossNO3,rprodNO3,wprodN2O5,pNO3temp,rBrOplusNO2,
     &  wlossClONO2,wprodClONO2,changeN,sphericalCorrectionReg1,
     &  sphericalCorrectionReg2,sphericalCorrectionReg3,
     &  sphericalCorrectionReg4,
     &  changeTerpenes,rTerpplusNO3,changeisopp1g,changeisopp2g,
     &  changeapinp1g,changeapinp2g,changeOx,fraQ,
     &  thick,changeCO,changeN_d1,changeN_d2,changeN_d3,changeNO3p,
#ifdef TRACERS_dCO
     &  rdHCH17OplusNO3,rdHCH18OplusNO3,rdH13CHOplusNO3,
     &  changed13Calke,
     &  changedHCH17O,changedHCH18O,changedH13CHO,
     &  changedC17O,changedC18O,changed13CO,
#endif  /* TRACERS_dCO */
     &  temp_SW,BRTOT,CLTOT,colmO2,colmO3,changeClONO2,changeClOx,
     &  changeHOCl,changeHCl,changehetClONO2,chgHT3,albedoToUse,
     &  chgHT4,chgHT5,rmrClOx,rmrBrOx,rmv,rmrOx,avgTT_H2O,avgTT_CH4,
     &  countTT,bHNO3,mHNO3,HNO3_thresh,Ttemp,changeBrOx,changeBrONO2,
     &  changeBrOx2,changeHBr,tempChangeNOx,ss27x2,ss27x2_c,OHpptv,
     &  HO2pptv,ObyO3,NO2byNO,ClbyClO,voc2nox_denom,tempChangeOx,pNOloc
      integer :: igas,LL,I,J,L,N,inss,L2,n2,Jqq,Iqq,
     & maxT,iu,itemp_iter,ih1330e,ih1030e,ih1030,ih1330,m,istep,index1,
     & index2,nb
      LOGICAL                   :: error, jay, daylight
      CHARACTER*4               :: ghg_name
      CHARACTER*80              :: ghg_file,title
      character(len=300)        :: out_line

      real*8 :: ghg_out(LM,GRID%I_STRT_HALO:GRID%I_STOP_HALO, ! remains LM
     &                     GRID%J_STRT_HALO:GRID%J_STOP_HALO)

      real*8, dimension(JM)     :: DU_O3_glob

      real*8, dimension(GRID%I_STRT_HALO:GRID%I_STOP_HALO,
     &                  GRID%J_STRT_HALO:GRID%J_STOP_HALO) ::
     &     avgTT_CH4_part,avgTT_H2O_part,countTT_part,
     &     surfIsop,zonalIsop

      INTEGER :: J_0, J_1, J_0S, J_1S, J_0H, J_1H, I_0, I_1
      integer :: initial_GHG_setup
      LOGICAL :: HAVE_SOUTH_POLE, HAVE_NORTH_POLE     
      real*8 :: qsat ! this is a function in UTILDBL.f
#if (defined TRACERS_TOMAS) || (defined CACHED_SUBDD)
      integer :: k
#endif
      integer :: hour, idx

      call modelEclock%get(hour=hour)

      call getDomainBounds(grid, 
     &               J_STRT    =J_0,  J_STOP    =J_1,
     &               I_STRT    =I_0,  I_STOP    =I_1,
     &               J_STRT_SKP=J_0S, J_STOP_SKP=J_1S,
     &               J_STRT_HALO=J_0H, J_STOP_HALO=J_1H,
     &               HAVE_SOUTH_POLE = have_south_pole,
     &               HAVE_NORTH_POLE = have_north_pole)
      
! calculate what longitudes to accumulate for 10:30am/1:30pm NO2 diags:
! Um... Does use of Jhour here assume starting the model at midnight?
      istep = NINT(real(IM)/HOURS_PER_DAY) ! number of boxes per hour
      ! ih1030/1330 are westmost I index that hour (careful: int arith.)
      ih1030 = istep*(10-hour)+IM/2+NINT(real(istep)/2.)-(istep-1)/2
      ih1330 = istep*(13-hour)+IM/2+NINT(real(istep)/2.)-(istep-1)/2
      if(ih1030 < 0) ih1030 = IM+ih1030  
      if(ih1330 < 0) ih1330 = IM+ih1330  
      if(ih1030 > IM) ih1030 = ih1030-IM
      if(ih1330 > IM) ih1330 = ih1330-IM
      ih1330e=ih1330+istep-1
      ih1030e=ih1030+istep-1
      if(ih1030e < 0) ih1030e = IM+ih1030e  
      if(ih1330e < 0) ih1330e = IM+ih1330e  
      if(ih1030e > IM) ih1030e = ih1030e-IM
      if(ih1330e > IM) ih1330e = ih1330e-IM

! meanwhile, initialize the instantaneous SUBDD of NO2 column:
! I think it will be overwritten for all i,j, so this can be
! a temporary check:
      save_NO2column(I_0:I_1,J_0:J_1)=undef
 
      if (is_set_param('initial_GHG_setup')) then
        call get_param('initial_GHG_setup', initial_GHG_setup)
        if (initial_GHG_setup == 1 .and. itime == itimeI) then
C--------special section for ghg runs ---------
          write(out_line,*)'Warning: INITIAL_GHG_SETUP is on!'
          call write_parallel(trim(out_line))
          if(use_rad_ch4>0 .or. use_rad_n2o>0 .or. use_rad_cfc>0)then
            rad_to_file(1,:,I_0:I_1,J_0:J_1)=
     &           rad_to_chem(1,:,I_0:I_1,J_0:J_1)
            rad_to_file(2,:,I_0:I_1,J_0:J_1)=
     &           rad_to_chem(2,:,I_0:I_1,J_0:J_1)
            do j=J_0,J_1
              do i=I_0,I_1
                rad_to_file(3,:,i,j)=rad_to_chem(3,:,i,j)*2.69e20*
     &               byavog*axyp(i,j)*tr_mm(n_N2O) ! i.e. in trm units now!
                rad_to_file(4,:,i,j)=rad_to_chem(4,:,i,j)*2.69e20*
     &               byavog*axyp(i,j)*tr_mm(n_CH4) ! i.e. in trm units now!
                rad_to_file(5,:,i,j)=rad_to_chem(5,:,i,j)*2.69e20*
     &               byavog*axyp(i,j)*tr_mm(n_CFC)*fact_CFC ! i.e. in trm units now!
              enddo
            enddo 
            if(ghg_yr/=0)then; write(ghg_name,'(I4.4)')ghg_yr
            else; write(ghg_name,'(I4.4)')modelEclock%getYear(); endif
            ghg_file='GHG_IC_'//ghg_name
            call openunit(ghg_file,iu,.true.,.false.)
            do m=1,5
             ghg_out(:,I_0:I_1,J_0:J_1)=rad_to_file(m,:,I_0:I_1,J_0:J_1)
             CALL WRITET8_COLUMN(grid,iu,NAMEUNIT(iu),GHG_OUT,ghg_file)
            enddo
            call closeunit(iu)          
            if(AM_I_ROOT( ))then
              write(6,*)'Kludge in masterchem to output inital'
              write(6,*)'conditions for ghgs.: ',trim(ghg_file)
              write(6,*)'First time step has used default values.'
              write(6,*)'If you wish to produce a correct'
              write(6,*)'first time step, then redo setup for this'
              write(6,*)' rundeck with initial_GHG_setup set to 0.'
              write(6,*)'Address questions to G. Faluvegi or T. Clune.'
              write(6,*)'Thanks.'
            end if
          end if
        end if
      end if

C Some INITIALIZATIONS :
      bydtsrc = 1.d0/dtsrc
      BYFJM   = 1.d0/real(JM)
      PRES2(1:LM) = PMIDL00(1:LM)

      if(H2ObyCH4 /= 0. .and. clim_interact_chem > 0)                
     &call stop_model('H2ObyCH4.ne.0 .and. clim_interact_chem > 0',13)

#ifdef INTERACTIVE_WETLANDS_CH4
C Not really related to chemistry, but convenient place to update
C running-averages for interactive wetlands CH4:
      do J=J_0,J_1; do I=I_0,IMAXJ(J)
        temp_SW=ALB(I,J,1)*(SRDN(I,J)+1.d-20)*COSZ1(I,J)
        call running_average(temp_SW,I,J,1.d0,n__sw)
      end do      ; end do
#endif

C Calculation of gas phase reaction rates for sulfur chemistry:
C Now called from tracer_3Dsource
c      CALL GET_SULF_GAS_RATES
      
#ifdef TRACERS_HETCHEM
c Calculation of removal rates on dust surfaces:
      CALL HETCDUST
#endif

      ! Note to self: move all Itime==ItimeI things to TRCHEM_init.f
      if(Itime==ItimeI)then
        if(use_rad_n2o > 0)then
          write(out_line,*) 'Warning:use_rad_n2o overrides PIfact_N2O'
          call write_parallel(trim(out_line))
        endif
        if(use_rad_cfc > 0)then
          write(out_line,*) 'Warning:use_rad_cfc overrides PIfact_CFC'
          call write_parallel(trim(out_line))
        endif
      endif

c Set "chemical time step". Really this is a method of applying only
c a fraction of the chemistry change to the tracer mass for the first
c 30 hours.  That fraction is: dt2/dtscr.  E.g. in the first hour it
c is (dtsrc/24)/dtsrc = 1/24th of the chemistry change is applied.
c This is to work around initial instabilities.

      if(allowSomeChemReinit == 1)then
        if(Itime-ItimeI <= 3)then
          dt2=dtsrc/24.d0          ! e.g. 150.
        elseif(Itime-ItimeI > 3 .and. Itime-ItimeI <= 6)then
          dt2=dtsrc/12.d0          ! e.g. 300.
        elseif(Itime-ItimeI > 6 .and. Itime-ItimeI <= 11)then
          dt2=dtsrc/6.d0           ! e.g. 600.
        elseif(Itime-ItimeI > 11 .and. Itime-ItimeI <= 30)then
          dt2=dtsrc/2.4d0          ! e.g. 1500.
        elseif(Itime-ItimeI > 30)then
          dt2=dtsrc                ! e.g. 3600
        endif
      else 
        dt2=dtsrc
      endif

c Calculate new photolysis rates every n_phot main timesteps:

C CALCULATE TX, THE REAL TEMPERATURE:
C (note this section is already done in DIAG.f)
      IF(HAVE_SOUTH_POLE) THEN
        DO L=1,LM
          TX(1,1,L)=T(1,1,L)*PK(L,1,1)
          TX(I_0:I_1,1,L)=TX(1,1,L)
        END DO
      ENDIF  
      IF(HAVE_NORTH_POLE) THEN
        DO L=1,LM
          TX(1,JM,L)=T(1,JM,L)*PK(L,1,JM)
          TX(I_0:I_1,JM,L)=TX(1,JM,L)
        END DO
      ENDIF
      DO L=1,LM
        DO J=J_0,J_1
          TX(I_0:I_1,J,L)=T(I_0:I_1,J,L)*PK(L,I_0:I_1,J)
        END DO
      END DO

C info to set strat H2O based on tropical tropopause H2O and CH4:
      if(Itime == ItimeI .and. allowSomeChemReinit == 1 )then
        avgTT_H2O_part(I_0:I_1,J_0:J_1)=0.d0
        avgTT_CH4_part(I_0:I_1,J_0:J_1)=0.d0
        countTT_part(I_0:I_1,J_0:J_1)=0.d0
        do J=J_0,J_1
          do I=I_0,IMAXJ(J)
            if(LAT2D_DG(I,J) >= -20. .and. LAT2D_DG(I,J) <= 20.)then
              avgTT_H2O_part(I,J) = Q(I,J,LTROPO(I,J))*MWabyMWw
              if(use_rad_ch4 > 0) then
                avgTT_CH4_part(I,J) =
     &          rad_to_chem(4,LTROPO(I,J),I,J)
     &          *2.69d20*byavog*mair*byMA(LTROPO(I,J),I,J)
              else
                avgTT_CH4_part(I,J) =
     &          trm(I,J,LTROPO(I,J),n_CH4)
     &          *mass2vol(n_CH4)*BYAXYP(I,J)*byMA(LTROPO(I,J),I,J)
              endif
              countTT_part(I,J) = 1.d0
            end if
          end do
        end do
        CALL GLOBALSUM(grid, avgTT_CH4_part, avgTT_CH4, all=.true.)
        CALL GLOBALSUM(grid, avgTT_H2O_part, avgTT_H2O, all=.true.)
        CALL GLOBALSUM(grid, countTT_part,   countTT,   all=.true.)
        if(countTT <= 0.)call stop_model('countTT.le.0',255)
      end if

! calculate Isoprene zonal mean, to be used for acetone
      do j=J_0,J_1
        do i=I_0,IMAXJ(j)
          surfIsop(i,j)=trm(i,j,1,n_Isoprene)*mass2vol(n_Isoprene)*
     &         byaxyp(i,j)*byMA(1,i,j)
        enddo
      enddo
      call zonalmean_ij2ij(surfIsop,zonalIsop)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      j_loop: DO J=J_0,J_1          ! ===> MAIN J LOOP BEGINS <===

      i_loop: DO I=I_0,IMAXJ(J)     ! ===> MAIN I LOOP BEGINS <===
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      y = 0.d0

      select case(which_trop)
      case(0); maxT=min(ltropo(I,J),topLevelOfChemistry)
      case(1); maxT=min(ls1-1,topLevelOfChemistry)
      case default; call stop_model('which_trop problem 4',255)
      end select

! Define acetone in terms of Isoprene:
!kt Terpenes should also be included here in the future
      do L=1,maxT
        acetone(L)=max(0.d0, ! in molec/cm3
     &  (1.25d0*(
     &    zonalIsop(i,j)-trm(i,j,L,n_Isoprene)*mass2vol(n_Isoprene)*
     &    byaxyp(i,j)*byMA(L,i,j)))*PMID(L,i,j)/(TX(i,j,L)*cboltz))
      enddo
      do L=maxT+1,topLevelOfChemistry
        acetone(L)=0.d0
      enddo
#ifdef TRACERS_dCO
      do L=1,topLevelOfChemistry
        d17Oacetone(L)=acetone(L)*dacetone_fact
        d18Oacetone(L)=acetone(L)*dacetone_fact
        d13Cacetone(L)=acetone(L)*dacetone_fact
      enddo

#endif  /* TRACERS_dCO */
      DO L=1,topLevelOfChemistry
c Initialize the 2D change variable:
       changeL(L,:)=0.d0  
c Save presure, temperature, thickness, rel. hum. in local arrays:
       ta(L)=TX(I,J,L)
       rh(L)=Q(i,j,l)/min(1.d0,QSAT(ta(L),lhe,pmid(L,i,j)))
       bythick(L)=1.d0/
     & (rgas*bygrav*TX(i,j,L)*LOG(pedn(L,i,j)/pedn(L+1,i,j)))
c Calculate M and set fixed ratios for O2 & H2:
       y(nM,L)=pmid(L,i,j)/(ta(L)*cboltz)
       y(nO2,L)=y(nM,L)*pO2*o2x
       if(pres2(l) > 20.d0)then
         y(nH2,L)=y(nM,L)*pfix_H2
       else
         ! Was: y(nH2,L)=y(nM,L)*pfix_H2*7.d1/(7.d1+L-maxT+1)
         ! Now: a drop of 0.3 molec/cm3 per 12 hPa decrease:
         y(nH2,L)=y(nM,L)*pfix_H2 + 2.5d-2*(pres2(L)-20.d0)
       endif
       CLTOT=0.d0

      ! check trm 1 vs 2
c Tracers (converted from mass to number density):
       do igas=1,ntm_chem
         idx=igas+ntm_chem_beg-1
         y(igas,L)=trm(I,J,L,idx)*y(nM,L)*mass2vol(idx)*
     &   BYAXYP(I,J)*byMA(L,I,J)
       enddo

! If we are fixing methane for chemistry purposes set it's y here:
! 0.55866d0 below is 1/1.79 (HALOE observations)
       if(fix_CH4_chemistry == 1) then
         if(lat2d_dg(i,j) < 0.) then         ! Southern Hemisphere
           if(L<LS1)y(nn_CH4,L)=y(nM,L)*ch4_init_sh*1.d-6 !troposphere
           if(abs(lat2d_dg(i,j)) > 30.) then ! extratropics
             if(L>=LS1)then
               y(nn_CH4,L)=                   ! stratosphere
     &         y(nM,L)*ch4_init_sh*0.55866d0*1.d-6*CH4altX(L)
             end if
           else                              ! tropics
             if(L>=LS1)then
               y(nn_CH4,L)=                   ! stratosphere
     &         y(nM,L)*ch4_init_sh*0.55866d0*1.d-6*CH4altT(L)
             end if
           end if
         else                                ! Northern Hemisphere
           if(L<LS1)y(nn_CH4,L)=y(nM,L)*ch4_init_sh*1.d-6 !troposphere
           if(abs(lat2d_dg(i,j)) > 30.) then ! extratropics
             if(L>=LS1)then
               y(nn_CH4,L)=                   ! stratosphere
     &         y(nM,L)*ch4_init_nh*0.55866d0*1.d-6*CH4altX(L)
             end if
           else                              ! tropics
             if(L>=LS1)then
               y(nn_CH4,L)=                   ! stratosphere
     &         y(nM,L)*ch4_init_nh*0.55866d0*1.d-6*CH4altT(L)
             end if
           end if
         end if
       end if

#if defined(TRACERS_AEROSOLS_Koch) || defined(TRACERS_AMP) || \
    defined(TRACERS_TOMAS)
C Concentrations of DMS and SO2 for sulfur chemistry:
       if (coupled_chem == 1) then
         ydms(i,j,L)=trm(i,j,L,n_dms)*y(nM,L)*(28.0D0/62.0D0)*
     &   BYAXYP(I,J)*byMA(L,I,J)
         yso2(i,j,L)=trm(i,j,L,n_so2)*y(nM,L)*(28.0D0/64.0D0)*
     &   BYAXYP(I,J)*byMA(L,I,J)
       else
         ! Convert from pptv to molecule cm-3:
         ydms(i,j,L)=dms_offline(i,j,L)*1.0d-12*y(nM,L)
         yso2(i,j,L)=so2_offline(i,j,L)*1.0d-12*y(nM,L)
       end if
#endif /* TRACERS_{AEROSOLS_Koch,AMP,TOMAS} */

c Save initial ClOx amount for use in ClOxfam:
       ClOx_old(L)=trm(I,J,L,n_ClOx)*y(nM,L)*mass2vol(n_ClOx)*
     & BYAXYP(I,J)*byMA(L,I,J)

c Limit N2O5 number density:
       if(y(nn_N2O5,L) < 1.) y(nn_N2O5,L)=1.d0
c Set H2O, based on Q:
       y(nH2O,L)=Q(I,J,L)*MWabyMWw*y(nM,L)
c Initialize stratospheric y(H2O) & GCM Q variable (!),
c based on tropical tropopause H2O and CH4:
       if(allowSomeChemReinit == 1) then
         if(Itime == ItimeI .and. L > LTROPO(I,J)) then
           y(nH2O,L) =  y(nM,L)*(avgTT_H2O/countTT +
     &     2.d0*(avgTT_CH4/countTT-y(nn_CH4,L)/y(nM,L)))
           if(clim_interact_chem > 0)then 
             fraQ=(y(nH2O,L)/(y(nM,L)*MWabyMWw))/Q(I,J,L)
             Q(I,J,L)=y(nH2O,L)/(y(nM,L)*MWabyMWw)
             if(fraQ < 1.)qmom(:,i,j,L)=qmom(:,i,j,L)*fraQ
#ifdef TRACERS_WATER
C**** Add water to relevant tracers as well
             do n=1,ntm
               select case (tr_wd_type(n))
               case (nWater)       ! water: initialise tracers
                 trm(i,j,L,n) = trm(i,j,L,n)*fraQ
                 if(fraQ < 1.)trmom(:,i,j,L,n) = trmom(:,i,j,L,n)*fraQ
               end select
             end do
#endif
           end if
         end if
       end if

c Initialize various other species:
c - set [NO]=0 (where?) for first HOx calc, NO2 = NOx:
c - set reactive species for use in family chemistry & nighttime NO2:

       if(Itime == ItimeI .and. allowSomeChemReinit == 1)then 
         y(nAldehyde,L)=y(nM,L)*pfix_Aldehyde
#ifdef TRACERS_dCO
         y(nd17Oald,L)=y(nM,L)*pfix_Aldehyde*dalke_IC_fact
         y(nd18Oald,L)=y(nM,L)*pfix_Aldehyde*dalke_IC_fact
         y(nd13Cald,L)=y(nM,L)*pfix_Aldehyde*dalke_IC_fact
#endif  /* TRACERS_dCO */
       else
         y(nAldehyde,L)=yAldehyde(I,J,L)
#ifdef TRACERS_dCO
         y(nd17Oald,L)=yd17Oald(I,J,L)
         y(nd18Oald,L)=yd18Oald(I,J,L)
         y(nd13Cald,L)=yd13Cald(I,J,L)
#endif  /* TRACERS_dCO */
       endif
       yNO3(I,J,L)   =pNO3(I,J,L)*y(nn_NOx,L)
       y(nNO2,L)     =y(nn_NOx,L)*pNOx(I,J,L)
       y(nNO,L)      =y(nn_NOx,L)-(y(nNO2,L)+yNO3(I,J,L))
       if(y(nNO,L) < 1.d0)y(nNO,L)=1.d0
       y(nO3,L)      =pOx(I,J,L)*y(nn_Ox,L)
       y(nCH3O2,L)   =yCH3O2(I,J,L)
#ifdef TRACERS_dCO
       y(ndCH317O2,L)=ydCH317O2(I,J,L)
       y(ndCH318O2,L)=ydCH318O2(I,J,L)
       y(nd13CH3O2,L)=yd13CH3O2(I,J,L)
#endif  /* TRACERS_dCO */
       y(nC2O3,L)    =yC2O3(I,J,L)
#ifdef TRACERS_dCO
       y(ndC217O3,L) =ydC217O3(I,J,L)
       y(ndC218O3,L) =ydC218O3(I,J,L)
       y(nd13C2O3,L) =yd13C2O3(I,J,L)
#endif  /* TRACERS_dCO */
       y(nXO2,L)     =yXO2(I,J,L)
       y(nXO2N,L)    =yXO2N(I,J,L)
       y(nRXPAR,L)   =yRXPAR(I,J,L)
#ifdef TRACERS_dCO
       y(nd13CXPAR,L)   =yd13CXPAR(I,J,L)
#endif  /* TRACERS_dCO */
       y(nROR,L)     =yROR(I,J,L)
#ifdef TRACERS_dCO
       y(nd17OROR,L) =yd17OROR(I,J,L)
       y(nd18OROR,L) =yd18OROR(I,J,L)
       y(nd13CROR,L) =yd13CROR(I,J,L)
#endif  /* TRACERS_dCO */
       y(nCl2,L)     =yCl2(I,J,L)
       y(nCl2O2,L)   =yCl2O2(I,J,L)
       y(nOClO,L)    =y(nn_ClOx,L)*pOClOx(I,J,L)
       y(nClO,L)     =y(nn_ClOx,L)*pClOx(I,J,L)
       y(nCl,L)      =y(nn_ClOx,L)*pClx(I,J,L)
       y(nBr,L)      =y(nn_BrOx,L)*(1.d0-pBrOx(I,J,L))
       y(nBrO,L)     =y(nn_BrOx,L)*pBrOx(I,J,L)
      END DO ! L

C For solar zenith angle, we use the arccosine of the COSZ1
C from the radiation code, which is the cosine of the solar zenith 
C angle averaged over the physics time step.
C If the solar zenith angle (sza) from the radiation code is > 90 deg,
C (and hence COSZ1 is set to 0), recalculate it with get_sza routine:
      IF(COSZ1(I,J) == 0.d0) THEN
        call get_sza(I,J,sza)
      ELSE
        sza = acos(COSZ1(I,J))*byradian
      END IF

C SUNLIGHT criteria:
      albedoToUse=ALB(I,J,1)
      ! previously: daylight=((ALB(I,J,1)/=0.d0).and.(sza<szamax))
      daylight=(sza<szamax)
      if(daylight)then
        if(albedoToUse/=0.d0)then
          mostRecentNonZeroAlbedo(I,J)=albedoToUse
        else
          albedoToUse=mostRecentNonZeroAlbedo(I,J)
        end if
      end if

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                 BEGIN PHOTOLYSIS                               C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      if(daylight)then

c Pass O3 array (in ppmv; here seems to be ppv) to fastj. Above these
C levels fastj2 uses Nagatani climatological O3, read in by chem_init: 
        DO L=1,topLevelOfChemistry
          O3_FASTJ(L)=y(nO3,L)/y(nM,L)
        END DO

! calculate photolysis rates
        call fastj2_drv(I, J, ta, rh, albedoToUse)
        call photo_acetone(I,J,sza*radian) ! simpler calculation for acetone

C Define and alter resulting photolysis coefficients (zj --> ss):

        ! Set above-chemistry-top O2 and O3 columns. Initial hardcoded numbers here
        ! were for the 0.1 model top. Scaling this linearly in pressure now. Note
        ! that, in the fastj2_init routine, the model will stop if the pressure at
        ! min(JPNL,topLevelOfChemistry) level top encroaches on the top of the ozone
        ! layer (to remind user that a rethink of below formula - and many other
        ! things - would be needed). In next two lines, 5.6d21 is really 5.6d20/0.1
        ! and 5.0d17 is 5.0d16/0.1:
        colmO2=5.6d21*plbot(min(JPNL,topLevelOfChemistry)+1)
        colmO3=5.0d17*plbot(min(JPNL,topLevelOfChemistry)+1)

        ! Using MAX() here because 
        ! letting this spherical corrections get too small (0?) causes NaNs
        ! e.g. in the SOA code. Also sza can be > 90, so COS could go negative:
        sphericalCorrectionReg1=
     &    (MAX(1.d-2,DCOS(sza*radian)))**reg1Power_SpherO2andN2Ocorr
        sphericalCorrectionReg2=
     &    (MAX(1.d-2,DCOS(sza*radian)))**reg2Power_SpherO2andN2Ocorr
        sphericalCorrectionReg3=
     &    (MAX(1.d-2,DCOS(sza*radian)))**reg3Power_SpherO2andN2Ocorr
        sphericalCorrectionReg4=
     &    (MAX(1.d-2,DCOS(sza*radian)))**reg4Power_SpherO2andN2Ocorr

        DO L=min(JPNL,topLevelOfChemistry),1,-1
          do inss=1,n_rj
            ss(inss,L,I,J)=zj(L,inss)
#ifndef SHINDELL_SKIP_WINDOW_TUNE /* note NOT defined */
            !reduce rates for gases that photolyze in window region (~200nm):
            if(inss == rj%O2__O_O .or.
     &         inss == rj%N2O__M_O1D) then ! for O2 and N2O reactions:
              ! Apply spherical corrections:
              if(pres2(L)>reg1TopPres_SpherO2andN2Ocorr)then
                ss(inss,L,I,J)=ss(inss,L,I,J)*sphericalCorrectionReg1
              else if(pres2(L)>reg2TopPres_SpherO2andN2Ocorr .and.
     &                pres2(L).le.reg1TopPres_SpherO2andN2Ocorr) then
                ss(inss,L,I,J)=ss(inss,L,I,J)*sphericalCorrectionReg2
              else if(pres2(L)>reg3TopPres_SpherO2andN2Ocorr .and.
     &                pres2(L).le.reg2TopPres_SpherO2andN2Ocorr) then
                ss(inss,L,I,J)=ss(inss,L,I,J)*sphericalCorrectionReg3
              else
                ss(inss,L,I,J)=ss(inss,L,I,J)*sphericalCorrectionReg4
              end if
              ! Then apply linear corrections for same reactions:
              if(inss == rj%O2__O_O) then
                ss(inss,L,I,J)=ss(inss,L,I,J)*windowO2corr
              else if(inss == rj%N2O__M_O1D) then
                ss(inss,L,I,J)=ss(inss,L,I,J)*windowN2Ocorr
              end if
#ifdef TRACERS_dCO
#ifndef TRACERS_dCO_bin_reprod
            else if(inss == rj%d17Oald__dHCH17O_CO
     &         .or. inss == rj%d17Oald__HCHO_dC17O
     &         .or. inss == rj%d18Oald__dHCH18O_CO
     &         .or. inss == rj%d18Oald__HCHO_dC18O
     &         .or. inss == rj%d13Cald__dH13CHO_CO
     &         .or. inss == rj%d13Cald__HCHO_d13CO
     &             ) then
              ! the yield is one third, since one isotopically labeled atom
              ! is assumed to exist in each aldehyde, not two
              ss(inss,L,I,J)=ss(inss,L,I,J)/2.d0
#endif  /* not TRACERS_dCO_bin_reprod */
#endif  /* TRACERS_dCO */
            end if
#endif /* not defined to skip */
          enddo
          taijls(i,j,L,ijlt_JO1D)=taijls(i,j,L,ijlt_JO1D)
     &      +ss(rj%O3__O1D_O2,L,i,j)
          taijls(i,j,L,ijlt_JNO2)=taijls(i,j,L,ijlt_JNO2)
     &      +ss(rj%NO2__NO_O,L,i,j)
          taijls(i,j,L,ijlt_JH2O2)=taijls(i,j,L,ijlt_JH2O2)
     &      +ss(rj%H2O2__OH_OH,L,i,j)
          thick=
     &    1.d-3*rgas*bygrav*TX(I,J,L)*LOG(PEDN(L,i,j)/PEDN(L+1,i,j))
          colmO2=colmO2+y(nO2,L)*thick*1.d5
          colmO3=colmO3+y(nO3,L)*thick*1.d5
! SF3 is photolysis of water in Schumann-Runge bands based on:
! Nicolet, Pl. Space Sci., p 871, 1984.
! SF3_fact is, if x[ ] = bin4_flux[ ]:
! {(x[present] - x[1988]) / (x[1991] - x[1988])} * 0.1E-6
! This gets ADDED to the 1.3E-6 factor in the SF3 calculation. Here,
! bin4_flux is a proxy for the flux from all 175-200nm bins. 
! (Drew says the ratio would be the same.)
          if(SF2_fact == 0.)call stop_model('SF2_fact=0 in master',255)
          if(pres2(L) <= 10.)then
            if((SF3_FACT+1.3d-6) < 0.)call stop_model
     &      ('(SF3_FACT+1.3d-6) < 0 in master',255)
            SF3(I,J,L)=(SF3_FACT+1.3d-6)*EXP(-1.d-7*colmO2**.35)
     &      *by35*SQRT(1.224d3*COSZ1(I,J)**2.+1.d0)
            ! SF3(I,J,L)=SF3(I,J,L)*5.d-2
          else
            SF3(I,J,L)=0.d0
          endif
! SF2 is photlysis of NO in bands (0-0) and (1-0) based on Nicolet,
! Pl. Space Sci., p 111, 1980. SF2_fact is a ratio representative of
! bands (0-0) and (1-0); =  bin5_flux[present] / bin5_flux[1988] :
          if(colmO2 > 2.d19)then
            SF2(I,J,L)=4.5d-6*EXP(-(1.d-8*colmO2**.38+5.d-19*colmO3))
          else
            SF2(I,J,L)=4.75d-6*EXP(-1.5d-20*colmO2)
          endif
          SF2(I,J,L)=SF2(I,J,L)*SF2_fact*
     &    by35*SQRT(1.224d3*COSZ1(I,J)**2.+1.d0)
        END DO

      endif ! (sunlight)
      
CCCCCCCCCCCCCCCCC END PHOTOLYSIS SECTION CCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                 BEGIN CHEMISTRY                                C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ! Check for PSC's existance. 
      pscX(:)=.false.
      do L=1,topLevelOfChemistry
        if(pres2(L) <= 250.d0 .and. pres2(L) >= 3.d0)then    ! pres crit for now
          if(lat2D_dg(I,J)<=PSClatS.or.lat2D_dg(I,J)>=PSClatN)then! lat crit for now
            if(lat2d_dg(I,J)<=PSClatS)then
              Ttemp=ta(L)+Tpsc_offset_S
            else if(lat2d_dg(I,J)>=PSClatN)then
              Ttemp=ta(L)+Tpsc_offset_N
            endif
            if(Ttemp <= T_thresh)then                         ! cold enough
              bHNO3=38.9855d0-11397.d0/Ttemp+0.009179d0*Ttemp ! H2O and HNO3
              mHNO3= -2.7836d0 - 8.8d-4*Ttemp                 ! criteria from
              HNO3_thresh=2.69d19/760.d0*10.d0**              ! Hanson+Mauersberger
     &        (mHNO3*log10(y(nH2O,L)*760.d0/2.69d19)+bHNO3)   ! 1988
              if(y(nn_HNO3,L) >= HNO3_thresh) pscX(L)=.true.! <-- yes PSC
            endif ! temperature
          endif   ! lat
        endif     ! pressure
      enddo       ! L
     
c Calculate the chemical reaction rates:
      call Crates(I,J)

#ifdef TRACERS_AEROSOLS_SOA
      voc2nox(:)=0.d0
#endif  /* TRACERS_AEROSOLS_SOA */

      if(daylight)then
CCCCCCCCCCCCCCCCCCCC   SUNLIGHT   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCC FAMILY PARTITIONING CCCCCCCCCCCCCCCCCCCCCCCCCC

       call Oxfam(topLevelOfChemistry,I,J)
       call HOxfam(topLevelOfChemistry,I,J)
       call NOxfam(topLevelOfChemistry,I,J)
       call BrOxfam(topLevelOfChemistry,I,J)
       ! ClOxfam done a bit below here

CCCCCCCCCCCCCCCCC NON-FAMILY CHEMISTRY CCCCCCCCCCCCCCCCCCCCCCCC
#ifdef TRACERS_AEROSOLS_SOA
! calculate voc2nox for SOA precursor chemistry
      do L=1,topLevelOfChemistry
        voc2nox_denom=(rr(rrbi%XO2_NO__NO2_M,L)*y(nNO,L)+
     &                 rr(rrbi%XO2_HO2__CH3OOH_O2,L)*y(nHO2,L)+
     &                 rr(rrbi%XO2_XO2__M_M,L)*yXO2(I,J,L))
        if (voc2nox_denom==0.d0) then
          voc2nox(L)=0.d0
        else
          voc2nox(L)=rr(rrbi%XO2_NO__NO2_M,L)*y(nNO,L)/
     &               voc2nox_denom
        end if
      end do
#endif  /* TRACERS_AEROSOLS_SOA */

      call chemstep(topLevelOfChemistry,I,J)

C Save 3D radical arrays to pass to aerosol code:
      if(coupled_chem == 1) then
        do L=1,topLevelOfChemistry
          oh_live(i,j,L)=y(nOH,L)
          no3_live(i,j,L)=yNO3(i,j,L)
        end do
        do L=topLevelOfChemistry+1,LM
          oh_live(i,j,L)=0.d0
          no3_live(i,j,L)=0.d0
        end do
      end if

      call ClOxfam(topLevelOfChemistry,I,J) ! needed something from chemstep.

      call printDaytimeChemistryDiags()
      
      else

CCCCCCCCCCCCCCCCCCCC END SUNLIGHT CCCCCCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCC   DARKNESS  CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C*****************************************************************
C               ABOUT: N2O5 sink on sulfate aerosol 
C                REACTION PROBABLITY FORMULATION:
C
C To evaluate 'k' for N2O5 + H2O -> 2HNO3, assume k = GAMMA*SA*v/4
C (Dentener and Crutzen, 1993). GAMMA = RGAMMASULF defined below,
c SA = given, v = molecular velocity (cm/s) where v = SQRT (8.Kb.T / PI*M);
C Kb = 1.38062E-23; T = Temp (K); M = mass N2O5 (Kg)
C
C Off-line sulfate fields to run in uncoupled mode give SO4 in
C cm2/cm3 'surface area density'.
C
C On-line sulfate in coupled mode must be converted to aerosol
C surface (cm2 aerosol/cm3 air) via aerosol volume fraction
C (cm3 aerosol/cm3 air). Assume a monodispersed aerosol with diameter
C of 0.078 microns. Specific aerosol density = 1.1g/cm3 (1.7g/cm3 ?)
C See Dentener and Crutzen, 1993 for details. Mr[sulfate] = 96.0g;
C Mr[(NH4)HSO4] = 115.0gC
C*****************************************************************

CCCCCCCCCCCCCCCC NIGHTTIME CCCCCCCCCCCCCCCCCCCCCC
#ifdef TRACERS_AEROSOLS_SOA
      call soa_apart ! calculate current apartmolar factors
#ifdef SOA_DIAGS
     &              (I,J)
#endif  /* SOA_DIAGS */
#endif  /* TRACERS_AEROSOLS_SOA */

      DO L=1,topLevelOfChemistry

! ----------------- RGAMMASULF ---------------------------------------
! Until more sophisticated method arrives, or when aerosol tracers are
! off, use method recommended by Faye, based on Kane et al., JPC, 2001
! and Hallquist et al., PCCP, 2003:
! For RH>50%, gamma=0.015. For RH <=50%:
! T<=290K, gamma=0.052 - 2.79d-4*RH [RH in percent]
! T>290K, gamma=above - log10(T-290)*0.05 [minimum gamma=0.001]
!
        if(rh(L)>0.5)then
          rgammasulf = 1.5d-2
        else
          rgammasulf = 5.2d-2 - 2.79d-4*100.d0*rh(L)
          if(ta(L)>290.) rgammasulf=
     &    max(1.d-3,rgammasulf-log10(ta(L)-290.d0)*5.d-2)
        end if
! --------------------------------------------------------------------
        if (coupled_chem == 1) then
          ! Convert SO4 from mass (kg) to aerosol surface per grid box:
          ! Here there is a factor of 1d-3  that converts kg/m3 to g/cm3
          ! and 1.76d5 is cm2/g from Dentener and Crutzen, 1993.
          ! So 1.d-3*1.76d5=1.76d2, and that value is for a relative
          ! humidity of 0.75 (1/0.75 = 1.33333d0 below). Reciprocal
          ! layer thickness below is in 1/m units:
           sulfate(i,j,l)=0.0
#ifdef TRACERS_AMP
          sulfate(i,j,L)=
     &    trm(i,j,L,n_M_AKK_SU)+trm(i,j,L,n_M_ACC_SU)+
     &    trm(i,j,L,n_M_DD1_SU)+trm(i,j,L,n_M_DS1_SU)+
     &    trm(i,j,L,n_M_DD2_SU)+trm(i,j,L,n_M_DS2_SU)+
     &    trm(i,j,L,n_M_SSA_SU)+trm(i,j,L,n_M_OCC_SU)+
     &    trm(i,j,L,n_M_BC1_SU)+trm(i,j,L,n_M_BC2_SU)+
     &    trm(i,j,L,n_M_BC3_SU)+trm(i,j,L,n_M_DBC_SU)+
     &    trm(i,j,L,n_M_BOC_SU)+trm(i,j,L,n_M_BCS_SU)+
     &    trm(i,j,L,n_M_MXX_SU)
#elif (defined TRACERS_AEROSOLS_Koch)
          sulfate(i,j,L)=trm(i,j,L,n_SO4)
#elif (defined TRACERS_TOMAS)
          do k=1,nbins
            sulfate(i,j,L)=sulfate(i,j,l)+trm(i,j,L,n_ASO4(k))
          end do
#endif
          sulfate(i,j,l)=sulfate(i,j,l)
     &      *1.76d2*byaxyp(i,j)*bythick(L)
     &      *max(0.1d0,rh(L)*1.33333d0)
          ! just in case loop changes (b/c sulfate is defined to LM):
          if(L>topLevelOfChemistry)sulfate(i,j,L)=0.d0
        end if

        pfactor=axyp(I,J)*MA(L,I,J)/y(nM,L)
        bypfactor=1.D0/pfactor
        RVELN2O5=SQRT(TX(I,J,L)*RKBYPIM)*100.d0
C       Calculate sulfate sink
c       in troposphere loss is rxn on sulfate, in strat rxn w PSC or sulfate
        wprod_sulf=
     &  dt2*sulfate(I,J,L)*y(nn_N2O5,L)*RGAMMASULF*RVELN2O5*0.25d0
        if(pres2(L)>5.d0)then
c       if there is reaction on strat particulate (in Crates), use that
          if(rr(rrhet%N2O5_H2O__HNO3_HNO3,L)>1.0d-25)
     &      wprod_sulf=dt2*y(nn_N2O5,L)*rr(rrhet%N2O5_H2O__HNO3_HNO3,L)
        else
          wprod_sulf=0.d0
        end if
        ! we used to limit (w)prod_sulf here w.r.t. N2O5

C*****************************************************************
c g signifies gas phase
c while prod_sulf and wprod_sulf are sulfate rxn
c wprods are in molecules/cm3/s
c prods/mass2vol are in mass units to add to tracers
c
c NO3 amounts are a function of reactions:
c    NO2+O3-->NO3+O2
c    NO2+NO3-->NO+NO2
c    NO3+NO3-->NO2+NO2
c leaving out, as outside of NOx family:
c    NO3+HCHO-->HNO3+CO
c    90% of Isoprene+NO3-->HO2+Alkenes,
c    N2O5+M-->NO3+NO2
c    NO3+NO2-->N2O5+M
c For NO2, similarly leave out:
c    PAN+M-->C2O3+NO2
c    HO2NO2+M-->HO2+NO2
c    N2O5+M-->NO3+NO2
c Keep NOx unchanged as this is only intrafamily.
C*****************************************************************

c       calculate NO3 vs NO2 (assume no NO at night)
        do itemp_iter=1,5
          rprodNO3=rr(rrbi%NO2_O3__NO3_O2,L)*y(nNO2,L)*y(nn_Ox,L)
     &      *pOx(I,J,L)
          rlossNO3=(2.d0*rr(rrbi%NO3_NO3__NO2_NO2,L)*yNO3(I,J,L))
     &        *yNO3(I,J,L)
     &      -(rr(rrbi%Alkenes_NO3__HCHO_NO2,L)*y(nn_Alkenes,L))
     &        *yNO3(I,J,L)
          if((rlossNO3+rprodNO3)==0.)
     &    call stop_model('(rlossNO3+rprodNO3)=0',255)
          pNO3temp=rprodNO3/(rlossNO3+rprodNO3)
          if(pNO3temp > 0.99d0)pNO3temp=0.95d0
          if(pNO3temp < 0.01d0)pNO3temp=0.05d0
          yNO3(I,J,L)=pNO3temp*y(nn_NOx,L)
          y(nNO2,L)=y(nn_NOx,L)-yNO3(I,J,L)
        end do
        pNOx(I,J,L)=y(nNO2,L)/y(nn_NOx,L)
        pNO3(I,J,L)=yNO3(I,J,L)/y(nn_NOx,L)

c       set reaction rates, then limit any uniformly across all
c       paths if lead to negative conc:
        rNO3plusNO2=
     &    rr(rrtri%NO3_NO2__N2O5_M,L)*y(nNO2,L)*yNO3(I,J,L)*dt2
        rN2O5decomp=rr(rrmono%N2O5_M__NO3_NO2,L)*y(nn_N2O5,L)*dt2
        chgHT5=rr(rrhet%N2O5_HCl__Cl_HNO3,L)*y(nn_N2O5,L)*dt2
        rHCHOplusNO3=
     &    y(nn_HCHO,L)*rr(rrbi%NO3_HCHO__HNO3_CO,L)*yNO3(I,J,L)*dt2
#ifdef TRACERS_dCO
        rdHCH17OplusNO3=
     &    y(nn_dHCH17O,L)*rr(rrbi%NO3_dHCH17O__HNO3_dC17O,L)
     &      *yNO3(I,J,L)*dt2
        rdHCH18OplusNO3=
     &    y(nn_dHCH18O,L)*rr(rrbi%NO3_dHCH18O__HNO3_dC18O,L)
     &      *yNO3(I,J,L)*dt2
        rdH13CHOplusNO3=
     &    y(nn_dH13CHO,L)*rr(rrbi%NO3_dH13CHO__HNO3_d13CO,L)
     &      *yNO3(I,J,L)*dt2
#endif  /* TRACERS_dCO */
        rAldplusNO3=2.5d-15*yAldehyde(I,J,L)*yNO3(I,J,L)*dt2
#ifdef TRACERS_dCO
        rd17OaldplusNO3=2.5d-15*yd17Oald(I,J,L)*yNO3(I,J,L)*dt2
        rd18OaldplusNO3=2.5d-15*yd18Oald(I,J,L)*yNO3(I,J,L)*dt2
        rd13CaldplusNO3=2.5d-15*yd13Cald(I,J,L)*yNO3(I,J,L)*dt2
#endif  /* TRACERS_dCO */
        rIsopplusNO3=rr(rrbi%Isoprene_NO3__HO2_Alkenes,L)
     &    *y(nn_Isoprene,L)*yNO3(I,J,L)*dt2
#ifdef TRACERS_TERP
        rTerpplusNO3=
     &    rr(rrbi%Terpenes_NO3__HO2_Alkenes,L)*y(nn_Terpenes,L)
     &      *yNO3(I,J,L)*dt2
#endif  /* TRACERS_TERP */
        rClOplusNO2=
     &    y(nClO,L)*rr(rrtri%ClO_NO2__ClONO2_M,L)*y(nNO2,L)*dt2
        rDMSplusNO3=ydms(i,j,L)*rsulf3(i,j,L)*yNO3(I,J,L)*dt2
        rBrOplusNO2=rr(rrtri%BrO_NO2__BrONO2_M,L)*y(nNO2,L) 
     &      *y(nn_BrOx,L)*pBrOx(I,J,L)*dt2
        chgHT3=rr(rrhet%ClONO2_HCl__Cl_HNO3,L)*y(nn_ClONO2,L)*dt2
        changehetClONO2=
     &    -1.d0*(rr(rrhet%ClONO2_H2O__HOCl_HNO3,L)*y(nn_ClONO2,L))*dt2

c       Examine prod and loss of N2O5
        wlossN2O5=rN2O5decomp+wprod_sulf+chgHT5
        wprodN2O5=rNO3plusNO2
        if(wlossN2O5>0.99d0*(y(nn_N2O5,L)+wprodN2O5))then
          if(wlossN2O5==0.)call stop_model('wprodN2O5=0',255)
          ratioN=0.99d0*(y(nn_N2O5,L)+wprodN2O5)/wlossN2O5
          rN2O5decomp=rN2O5decomp*ratioN
          wprod_sulf=wprod_sulf*ratioN
          chgHT5=chgHT5*ratioN
        end if

c       Examine prod and loss of NOx
        wlossNOx=rHCHOplusNO3+rAldplusNO3+2.d0*rNO3plusNO2
     &   +0.9d0*rIsopplusNO3+rClOplusNO2+rDMSplusNO3+rBrOplusNO2
#ifdef TRACERS_TERP
     &   +0.9d0*rTerpplusNO3
#endif  /* TRACERS_TERP */
        wprodNOx=2.d0*rN2O5decomp
        if(wlossNOx>0.99d0*(y(nn_NOx,L)+wprodNOx))then
          if(wlossNOx==0.)call stop_model('wlossNOx=0',255)
          ratioN=0.99d0*(y(nn_NOx,L)+wprodNOx)/wlossNOx
          rHCHOplusNO3=rHCHOplusNO3*ratioN
          rAldplusNO3=rAldplusNO3*ratioN
#ifdef TRACERS_dCO
          rdHCH17OplusNO3=rdHCH17OplusNO3*ratioN
          rdHCH18OplusNO3=rdHCH18OplusNO3*ratioN
          rdH13CHOplusNO3=rdH13CHOplusNO3*ratioN
          rd17OaldplusNO3=rd17OaldplusNO3*ratioN
          rd18OaldplusNO3=rd18OaldplusNO3*ratioN
          rd13CaldplusNO3=rd13CaldplusNO3*ratioN
#endif  /* TRACERS_dCO */
          rNO3plusNO2=rNO3plusNO2*ratioN
          rIsopplusNO3=rIsopplusNO3*ratioN
#ifdef TRACERS_TERP
          rTerpplusNO3=rTerpplusNO3*ratioN
#endif  /* TRACERS_TERP */
          rClOplusNO2=rClOplusNO2*ratioN
          rDMSplusNO3=rDMSplusNO3*ratioN
          rBrOplusNO2=rBrOplusNO2*ratioN
        end if

        if(rClOplusNO2 >= y(nClO,L))rClOplusNO2=0.8d0*y(nClO,L)
        changeClONO2=rClOplusNO2
        changeClOx=-1.d0*changeClONO2

c       Examine prod and loss of ClONO2
        wlossClONO2=changehetClONO2+chgHT3
        wprodClONO2=rClOplusNO2
        if(wlossClONO2>0.99d0*(y(nn_ClONO2,L)+wprodClONO2))then
          if(wlossClONO2==0.)call stop_model('wprodClONO2=0',255)
          ratioN=0.99d0*(y(nn_ClONO2,L)+wprodClONO2)/wlossClONO2
          changehetClONO2=changehetClONO2*ratioN
          chgHT3=chgHT3*ratioN
        end if

        prod_sulf=wprod_sulf*pfactor

C       LOWER LIMIT ON N2O5:
        if(y(nn_N2O5,L) <= 1.d0) y(nn_N2O5,L)=1.d0

C Calculate and limit gaseous changes to HNO3, HCHO, N2O5, Aldehyde,
C Alkenes, Isoprene, Terpenes (if used) and AlkylNit:

        gwprodHNO3=rHCHOplusNO3+rAldplusNO3
        if(gwprodHNO3 > y(nn_HCHO,L))gwprodHNO3=y(nn_HCHO,L)
#ifdef TRACERS_dCO
        gwprodHNO3dHCH17O=rdHCH17OplusNO3+rd17OaldplusNO3
        if(gwprodHNO3dHCH17O > y(nn_dHCH17O,L))
     &    gwprodHNO3dHCH17O=y(nn_dHCH17O,L)
        gwprodHNO3dHCH18O=rdHCH18OplusNO3+rd18OaldplusNO3
        if(gwprodHNO3dHCH18O > y(nn_dHCH18O,L))
     &    gwprodHNO3dHCH18O=y(nn_dHCH18O,L)
        gwprodHNO3dH13CHO=rdH13CHOplusNO3+rd13CaldplusNO3
        if(gwprodHNO3dH13CHO > y(nn_dH13CHO,L))
     &    gwprodHNO3dH13CHO=y(nn_dH13CHO,L)
#endif  /* TRACERS_dCO */

        changeAldehyde=(rr(rrbi%Alkenes_NO3__HCHO_NO2,L)*y(nn_Alkenes,L)
     &      +rr(rrbi%Isoprene_NO3__HO2_Alkenes,L)*y(nn_Isoprene,L)
     &        *0.12d0
#ifdef TRACERS_TERP
     &      +rr(rrbi%Terpenes_NO3__HO2_Alkenes,L)*y(nn_Terpenes,L)
     &        *0.12d0
#endif  /* TRACERS_TERP */
     &      -2.5d-15*yAldehyde(I,J,L)
     &    )*yNO3(I,J,L)*dt2
        if(-changeAldehyde > 0.75d0*yAldehyde(I,J,L))changeAldehyde=
     &  -0.75d0*yAldehyde(I,J,L)

#ifdef TRACERS_dCO
        changed17Oald=(rr(rrbi%Alkenes_NO3__HCHO_NO2,L)*y(nn_Alkenes,L)
     &      +rr(rrbi%Isoprene_NO3__HO2_Alkenes,L)*y(nn_Isoprene,L)
     &        *0.12d0
#ifdef TRACERS_TERP
     &      +rr(rrbi%Terpenes_NO3__HO2_Alkenes,L)*y(nn_Terpenes,L)
     &        *0.12d0
#endif  /* TRACERS_TERP */
     &      -2.5d-15*yd17Oald(I,J,L)
     &    )*yNO3(I,J,L)*dt2
        if(-changed17Oald > 0.75d0*yd17Oald(I,J,L))changed17Oald=
     &  -0.75d0*yd17Oald(I,J,L)

        changed18Oald=(rr(rrbi%Alkenes_NO3__HCHO_NO2,L)*y(nn_Alkenes,L)
     &      +rr(rrbi%Isoprene_NO3__HO2_Alkenes,L)*y(nn_Isoprene,L)
     &        *0.12d0
#ifdef TRACERS_TERP
     &      +rr(rrbi%Terpenes_NO3__HO2_Alkenes,L)*y(nn_Terpenes,L)
     &        *0.12d0
#endif  /* TRACERS_TERP */
     &      -2.5d-15*yd18Oald(I,J,L)
     &    )*yNO3(I,J,L)*dt2
        if(-changed18Oald > 0.75d0*yd18Oald(I,J,L))changed18Oald=
     &  -0.75d0*yd18Oald(I,J,L)

        changed13Cald=(rr(rrbi%d13Calke_NO3__dH13CHO_NO2,L)
     &        *y(nn_d13Calke,L)
     &      +rr(rrbi%Isoprene_NO3__HO2_d13Calke,L)*y(nn_Isoprene,L)
     &        *0.12d0
#ifdef TRACERS_TERP
     &      +rr(rrbi%Terpenes_NO3__HO2_d13Calke,L)*y(nn_Terpenes,L)
     &        *0.12d0
#endif  /* TRACERS_TERP */
     &      -2.5d-15*yd13Cald(I,J,L)
     &    )*yNO3(I,J,L)*dt2
        if(-changed13Cald > 0.75d0*yd13Cald(I,J,L))changed13Cald=
     &  -0.75d0*yd13Cald(I,J,L)
#endif  /* TRACERS_dCO */

        changeAlkenes=(rr(rrbi%Isoprene_NO3__HO2_Alkenes,L)
     &        *y(nn_Isoprene,L)*0.45d0
#ifdef TRACERS_TERP
     &      +rr(rrbi%Terpenes_NO3__HO2_Alkenes,L)*y(nn_Terpenes,L)
     &        *0.45d0
#endif  /* TRACERS_TERP */
     &      -rr(rrbi%Alkenes_NO3__HCHO_NO2,L)*y(nn_Alkenes,L)
     &    )*yNO3(I,J,L)*dt2
     &    +(rr(rrbi%Isoprene_O3__HCHO_Alkenes,L)*y(nn_Isoprene,L)*0.55d0
#ifdef TRACERS_TERP
     &    +rr(rrbi%Terpenes_O3__HCHO_Alkenes,L)*y(nn_Terpenes,L)*0.55d0
#endif  /* TRACERS_TERP */
     &    -rr(rrbi%Alkenes_O3__HCHO_CO,L)*y(nn_Alkenes,L)
     &    )*y(nO3,L)*dt2
        if(-changeAlkenes > 0.75d0*y(nn_Alkenes,L))changeAlkenes=
     &  -0.75d0*y(nn_Alkenes,L)

#ifdef TRACERS_dCO
        changed13Calke=(rr(rrbi%Isoprene_NO3__HO2_d13Calke,L)
     &        *y(nn_Isoprene,L)*0.45d0
#ifdef TRACERS_TERP
     &      +rr(rrbi%Terpenes_NO3__HO2_d13Calke,L)*y(nn_Terpenes,L)
     &        *0.45d0
#endif  /* TRACERS_TERP */
     &      -rr(rrbi%d13Calke_NO3__dH13CHO_NO2,L)*y(nn_d13Calke,L)
     &    )*yNO3(I,J,L)*dt2
     &    +(rr(rrbi%Isoprene_O3__dH13CHO_d13Calke,L)*y(nn_Isoprene,L)
     &      *0.55d0
#ifdef TRACERS_TERP
     &    +rr(rrbi%Terpenes_O3__dH13CHO_d13Calke,L)*y(nn_Terpenes,L)
     &      *0.55d0
#endif  /* TRACERS_TERP */
     &    -rr(rrbi%d13Calke_O3__dH13CHO_d13CO,L)*y(nn_d13Calke,L)
     &    )*y(nO3,L)*dt2
        if(-changed13Calke > 0.75d0*y(nn_d13Calke,L))changed13Calke=
     &  -0.75d0*y(nn_d13Calke,L)
#endif  /* TRACERS_dCO */

#ifdef TRACERS_AEROSOLS_SOA
! WARNING!!!
! No nighttime production from OH reactions, since no nighttime OH exist.
! This should be improved in the future.
        if (L<=LM_soa) then
          changeisopp1g=apartmolar(L,whichsoa(n_isopp1a))
     &      *(rr(rrbi%Isoprene_O3__HCHO_Alkenes,L)*y(nO3,L))
     &      *y(nn_Isoprene,L)*dt2
          if(-changeisopp1g > 0.75d0*y(nn_isopp1g,L))changeisopp1g=
     &    -0.75d0*y(nn_isopp1g,L)
          changeisopp2g=apartmolar(L,whichsoa(n_isopp2a))
     &      *(rr(rrbi%Isoprene_O3__HCHO_Alkenes,L)*y(nO3,L))
     &      *y(nn_Isoprene,L)*dt2
          if(-changeisopp2g > 0.75d0*y(nn_isopp2g,L))changeisopp2g=
     &    -0.75d0*y(nn_isopp2g,L)
#ifdef TRACERS_TERP
          changeapinp1g=apartmolar(L,whichsoa(n_apinp1a))
     &      *(rr(rrbi%Terpenes_O3__HCHO_Alkenes,L)*y(nO3,L))
     &      *y(nn_Terpenes,L)*dt2
          if(-changeapinp1g > 0.75d0*y(nn_apinp1g,L))changeapinp1g=
     &    -0.75d0*y(nn_apinp1g,L)
          changeapinp2g=apartmolar(L,whichsoa(n_apinp2a))
     &      *(rr(rrbi%Terpenes_O3__HCHO_Alkenes,L)*y(nO3,L))
     &      *y(nn_Terpenes,L)*dt2
          if(-changeapinp2g > 0.75d0*y(nn_apinp2g,L))changeapinp2g=
     &    -0.75d0*y(nn_apinp2g,L)
#endif  /* TRACERS_TERP */
        else
          changeisopp1g=0.d0
          changeisopp2g=0.d0
#ifdef TRACERS_TERP
          changeapinp1g=0.d0
          changeapinp2g=0.d0
#endif  /* TRACERS_TERP */
        end if
#endif  /* TRACERS_AEROSOLS_SOA */

        changeIsoprene=
     &    -(rr(rrbi%Isoprene_NO3__HO2_Alkenes,L)*yNO3(I,J,L)
     &    +rr(rrbi%Isoprene_O3__HCHO_Alkenes,L)*y(nO3,L))
     &      *y(nn_Isoprene,L)*dt2
        if(-changeIsoprene > 0.75d0*y(nn_Isoprene,L))changeIsoprene=
     &  -0.75d0*y(nn_Isoprene,L)

#ifdef TRACERS_TERP
        changeTerpenes=-(
     &      rr(rrbi%Terpenes_NO3__HO2_Alkenes,L)*yNO3(I,J,L)
     &      +rr(rrbi%Terpenes_O3__HCHO_Alkenes,L)*y(nO3,L)
     &    )*y(nn_Terpenes,L)*dt2
        if(-changeTerpenes > 0.75d0*y(nn_Terpenes,L))changeTerpenes=
     &  -0.75d0*y(nn_Terpenes,L)
#endif  /* TRACERS_TERP */

        changeHCHO=(rr(rrbi%Alkenes_NO3__HCHO_NO2,L)*y(nn_Alkenes,L)
     &      +rr(rrbi%Isoprene_NO3__HO2_Alkenes,L)*y(nn_Isoprene,L)
     &        *0.03d0
#ifdef TRACERS_TERP
     &      +rr(rrbi%Terpenes_NO3__HO2_Alkenes,L)*y(nn_Terpenes,L)
     &        *0.03d0
#endif  /* TRACERS_TERP */
     &    )*yNO3(I,J,L)*dt2
     &    -gwprodHNO3
     &    +(rr(rrbi%Isoprene_O3__HCHO_Alkenes,L)*y(nn_Isoprene,L)*0.9d0
#ifdef TRACERS_TERP
     &    +rr(rrbi%Terpenes_O3__HCHO_Alkenes,L)*y(nn_Terpenes,L)*0.9d0
#endif  /* TRACERS_TERP */
     &    +rr(rrbi%Alkenes_O3__HCHO_CO,L)*y(nn_Alkenes,L)
     &    )*y(nO3,L)*0.64d0*dt2

#ifdef TRACERS_dCO
        changedHCH17O=(
     *      rr(rrbi%Alkenes_NO3__dHCH17O_NO2,L)*y(nn_Alkenes,L)
     &      +rr(rrbi%Isoprene_NO3__HO2_Alkenes,L)*y(nn_Isoprene,L)
     &        *0.03d0
#ifdef TRACERS_TERP
     &      +rr(rrbi%Terpenes_NO3__HO2_Alkenes,L)*y(nn_Terpenes,L)
     &        *0.03d0
#endif  /* TRACERS_TERP */
     &    )*yNO3(I,J,L)*dt2
     &    -gwprodHNO3dHCH17O
     &    +(rr(rrbi%Isoprene_O3__dHCH17O_Alkenes,L)*y(nn_Isoprene,L)
     &      *0.9d0
#ifdef TRACERS_TERP
     &    +rr(rrbi%Terpenes_O3__dHCH17O_Alkenes,L)*y(nn_Terpenes,L)
     &      *0.9d0
#endif  /* TRACERS_TERP */
     &    +rr(rrbi%Alkenes_O3__dHCH17O_CO,L)*y(nn_Alkenes,L))*y(nO3,L)
     &      *0.64d0*dt2

        changedHCH18O=(
     &      rr(rrbi%Alkenes_NO3__dHCH18O_NO2,L)*y(nn_Alkenes,L)
     &      +rr(rrbi%Isoprene_NO3__HO2_Alkenes,L)*y(nn_Isoprene,L)
     &        *0.03d0
#ifdef TRACERS_TERP
     &      +rr(rrbi%Terpenes_NO3__HO2_Alkenes,L)*y(nn_Terpenes,L)
     &        *0.03d0
#endif  /* TRACERS_TERP */
     &    )*yNO3(I,J,L)*dt2
     &    -gwprodHNO3dHCH18O
     &    +(rr(rrbi%Isoprene_O3__dHCH18O_Alkenes,L)*y(nn_Isoprene,L)
     &      *0.9d0
#ifdef TRACERS_TERP
     &    +rr(rrbi%Terpenes_O3__dHCH18O_Alkenes,L)*y(nn_Terpenes,L)
     &      *0.9d0
#endif  /* TRACERS_TERP */
     &    +rr(rrbi%Alkenes_O3__dHCH18O_CO,L)*y(nn_Alkenes,L))*y(nO3,L)
     &      *0.64d0*dt2

        changedH13CHO=(
     *      rr(rrbi%d13Calke_NO3__dH13CHO_NO2,L)*y(nn_d13Calke,L)
     &      +rr(rrbi%Isoprene_NO3__HO2_d13Calke,L)*y(nn_Isoprene,L)
     &        *0.03d0
#ifdef TRACERS_TERP
     &      +rr(rrbi%Terpenes_NO3__HO2_d13Calke,L)*y(nn_Terpenes,L)
     &        *0.03d0
#endif  /* TRACERS_TERP */
     &    )*yNO3(I,J,L)*dt2
     &    -gwprodHNO3dH13CHO
     &    +(rr(rrbi%Isoprene_O3__dH13CHO_d13Calke,L)*y(nn_Isoprene,L)
     *      *0.9d0
#ifdef TRACERS_TERP
     &    +rr(rrbi%Terpenes_O3__dH13CHO_d13Calke,L)*y(nn_Terpenes,L)
     *      *0.9d0
#endif  /* TRACERS_TERP */
     &    +rr(rrbi%d13Calke_O3__dH13CHO_d13CO,L)*y(nn_d13Calke,L)
     &    )*y(nO3,L)*0.64d0*dt2
#endif  /* TRACERS_dCO */

        changeAlkylNit=rIsopplusNO3*0.9d0
#ifdef TRACERS_TERP
     &                +rTerpplusNO3*0.9d0
#endif  /* TRACERS_TERP */

c Convert some changes to molecules/cm3/s:
        changeHNO3=gwprodHNO3+2.d0*wprod_sulf  !always positive

        wlossNOx=rHCHOplusNO3+rAldplusNO3+2.d0*rNO3plusNO2
     &   +0.9d0*rIsopplusNO3+rClOplusNO2+rDMSplusNO3
#ifdef TRACERS_TERP
     &   +0.9d0*rTerpplusNO3
#endif  /* TRACERS_TERP */
        wprodNOx=2.d0*rN2O5decomp

        changeNOx=wprodNOx-wlossNOx

        wlossN2O5=rN2O5decomp+wprod_sulf
        gwprodN2O5=rNO3plusNO2
        changeN2O5=gwprodN2O5-wlossN2O5

c       Nighttime changes in Bromine-containing species
        changeBrOx=-1.d0*rBrOplusNO2
        if(rBrOplusNO2>0.5d0*y(nn_BrOx,L))
     &   changeBrOx=-0.5d0*y(nn_BrOx,L)
        changeBrONO2=-changeBrOx
        changeNOx=changeNOx+changeBrOx

c       Br+H2O2 converts to HBr+HO2. HO2 assumed to revert to H2O2
        changeBrOx2=-rr(rrbi%Br_H2O2__HBr_HO2,L)*y(nn_H2O2,L)
     &    *y(nn_BrOx,L)*(1.d0-pBrOx(I,J,L))*dt2
        if(-1.d0*changeBrOx2>0.2d0*y(nn_BrOx,L))
     &    changeBrOx2=-0.2d0*y(nn_BrOx,L)
        changeBrOx=changeBrOx+changeBrOx2
        changeHBr=-1.d0*changeBrOx2

c Heterogeneous reaction ClONO2+H2O on sulfate (and PSCs if present):
        if(rr(rrhet%ClONO2_H2O__HOCl_HNO3,L) > 2.d-35)then
          changeClONO2=changeClONO2+changehetClONO2
          changeHOCl=-changehetClONO2
          changeHNO3=changeHNO3-changehetClONO2
        else
          changeHOCl=0.d0
        end if
        
        changeHCl=0.d0

c Polar Stratospheric Clouds (PSC) Chemistry:
c 106 N2O5    +H2O     -->HNO3    +HNO3  (calculated above)
c 107 ClONO2  +H2O     -->HOCl    +HNO3  (calculated above)
c 108 ClONO2  +HCl     -->Cl      +HNO3  !really makes Cl2
c 109 HOCl    +HCl     -->Cl      +H2O   !really makes Cl2
c 110 N2O5    +HCl     -->Cl      +HNO3  !really makes ClNO2 (calc above)
        if(pscX(L)) then  ! PSCs exist
          if(chgHT3 >= 0.2d0*y(nn_HCl,L))chgHT3=0.2d0*y(nn_HCl,L)
          chgHT4=rr(rrhet%HOCl_HCl__Cl_H2O,L)*y(nn_HOCl,L)*dt2
          if(chgHT4 >= 0.2d0*y(nn_HOCl,L))chgHT4=0.2d0*y(nn_HOCl,L)
          if(chgHT4 >= 0.2d0*y(nn_HCl,L))chgHT4=0.2d0*y(nn_HCl,L)
          if(chgHT5 >= 0.5d0*y(nn_HCl,L))chgHT5=0.5d0*y(nn_HCl,L)
          changeClONO2=changeClONO2-chgHT3
          changeHOCl=changeHOCl-chgHT4
          changeN2O5=changeN2O5-chgHT5
          changeNOx=changeNOx+chgHT5
c         Note that really the following 3 produce Cl2, not ClOx, and Cl2
C         at night is stable and doesn't go back into ClONO2, so
C         should eventually keep track of Cl2/ClOx partitioning!
          changeHCl=changeHCl-chgHT3-chgHT4-chgHT5
          changeHNO3=changeHNO3+chgHT3+chgHT5
          changeClOx=changeClOx+2.d0*(chgHT3+chgHT4+chgHT5)
          ! Here we USED TO remove some of the HNO3 formed heterogeneously,
          ! as it doesn't come back to the gas phase. No longer.
        end if

        if(-1.d0*changeNOx>y(nn_NOx,L))changeNOx=-1.d0*y(nn_NOx,L)

#ifdef TRACERS_HETCHEM
C       Include reactions on dust for HNO3:
        changeHNO3 = changeHNO3 - krate(i,j,L,1,1)*y(nn_HNO3,L)*dt2
        changeN_d1 = krate(i,j,L,2,1) * y(nn_HNO3,L) *dt2
        changeN_d2 = krate(i,j,L,3,1) * y(nn_HNO3,L) *dt2
        changeN_d3 = krate(i,j,L,4,1) * y(nn_HNO3,L) *dt2
#endif

C Apply Alkenes, AlkyNit, and Aldehyde changes here:
        y(nn_Alkenes,L)  =y(nn_Alkenes,L)  +changeAlkenes
#ifdef TRACERS_dCO
        y(nn_d13Calke,L)  =y(nn_d13Calke,L)  +changed13Calke
#endif  /* TRACERS_dCO */
        y(nn_AlkylNit,L) =y(nn_AlkylNit,L) +changeAlkylNit
        yAldehyde(I,J,L)=yAldehyde(I,J,L)+changeAldehyde
#ifdef TRACERS_dCO
        yd17Oald(I,J,L)=yd17Oald(I,J,L)+changed17Oald
        yd18Oald(I,J,L)=yd18Oald(I,J,L)+changed18Oald
        yd13Cald(I,J,L)=yd13Cald(I,J,L)+changed13Cald
#endif  /* TRACERS_dCO */

#ifdef TRACERS_AEROSOLS_SOA
        y(nn_isopp1g,L)  =y(nn_isopp1g,L)  +changeisopp1g
        y(nn_isopp2g,L)  =y(nn_isopp2g,L)  +changeisopp2g
#ifdef TRACERS_TERP
        y(nn_apinp1g,L)  =y(nn_apinp1g,L)  +changeapinp1g
        y(nn_apinp2g,L)  =y(nn_apinp2g,L)  +changeapinp2g
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */

C Note: the lower limit of minKG placed on the resulting tracer mass
C from the following changes is to prevent negative tracer mass:

C -- HCHO --
c       Gas phase NO3 + HCHO -> HNO3 + CO yield of HCHO & CO:
        changeL(L,n_HCHO)=changeHCHO*pfactor*vol2mass(n_HCHO)
        if(-changeL(L,n_HCHO) > trm(I,J,L,n_HCHO))then
          changeL(L,n_HCHO)=-.95d0*trm(I,J,L,n_HCHO)
          changeHCHO=changeL(L,n_HCHO)*mass2vol(n_HCHO)*bypfactor
        endif
        IF((trm(i,j,l,n_HCHO)+changeL(l,n_HCHO)) < minKG) THEN
          changeL(l,n_HCHO) = minKG - trm(i,j,l,n_HCHO)
          changeHCHO=changeL(L,n_HCHO)*mass2vol(n_HCHO)*bypfactor
        ENDIF
        wprodHCHO=changeHCHO
#ifdef TRACERS_dCO
        changeL(L,n_dHCH17O)=changedHCH17O*pfactor*vol2mass(n_dHCH17O)
        if(-changeL(L,n_dHCH17O) > trm(I,J,L,n_dHCH17O))then
          changeL(L,n_dHCH17O)=-.95d0*trm(I,J,L,n_dHCH17O)
          changedHCH17O=changeL(L,n_dHCH17O)*mass2vol(n_dHCH17O)*
     &                  bypfactor
        endif
        IF((trm(i,j,l,n_dHCH17O)+changeL(l,n_dHCH17O)) < minKG) THEN
          changeL(l,n_dHCH17O) = minKG - trm(i,j,l,n_dHCH17O)
          changedHCH17O=changeL(L,n_dHCH17O)*mass2vol(n_dHCH17O)*
     &                  bypfactor
        ENDIF
        changeL(L,n_dHCH18O)=changedHCH18O*pfactor*vol2mass(n_dHCH18O)
        if(-changeL(L,n_dHCH18O) > trm(I,J,L,n_dHCH18O))then
          changeL(L,n_dHCH18O)=-.95d0*trm(I,J,L,n_dHCH18O)
          changedHCH18O=changeL(L,n_dHCH18O)*mass2vol(n_dHCH18O)*
     &                  bypfactor
        endif
        IF((trm(i,j,l,n_dHCH18O)+changeL(l,n_dHCH18O)) < minKG) THEN
          changeL(l,n_dHCH18O) = minKG - trm(i,j,l,n_dHCH18O)
          changedHCH18O=changeL(L,n_dHCH18O)*mass2vol(n_dHCH18O)*
     &                  bypfactor
        ENDIF
        changeL(L,n_dH13CHO)=changedH13CHO*pfactor*vol2mass(n_dH13CHO)
        if(-changeL(L,n_dH13CHO) > trm(I,J,L,n_dH13CHO))then
          changeL(L,n_dH13CHO)=-.95d0*trm(I,J,L,n_dH13CHO)
          changedH13CHO=changeL(L,n_dH13CHO)*mass2vol(n_dH13CHO)*
     &                  bypfactor
        endif
        IF((trm(i,j,l,n_dH13CHO)+changeL(l,n_dH13CHO)) < minKG) THEN
          changeL(l,n_dH13CHO) = minKG - trm(i,j,l,n_dH13CHO)
          changedH13CHO=changeL(L,n_dH13CHO)*mass2vol(n_dH13CHO)*
     &                  bypfactor
        ENDIF
#endif  /* TRACERS_dCO */
C -- CO --
        changeL(L,n_CO)=rHCHOplusNO3*pfactor*vol2mass(n_CO)
        changeCO=changeL(L,n_CO)*mass2vol(n_CO)*bypfactor
        if((trm(i,j,l,n_CO)+changeL(l,n_CO)) < minKG)then
          changeL(l,n_CO) = minKG - trm(i,j,l,n_CO)
          changeCO=changeL(L,n_CO)*mass2vol(n_CO)*bypfactor
        endif
        wprodCO=rHCHOplusNO3   ! <-- note
        if(changeL(L,n_CO) >= 0.) then  
          CALL INC_TAJLS(I,J,L,jls_COp,changeL(L,n_CO))
#ifdef ACCMIP_LIKE_DIAGS
          taijls(i,j,L,ijlt_COp)=taijls(i,j,L,ijlt_COp)+changeCO*cpd
#endif
        else
          CALL INC_TAJLS(I,J,L,jls_COd,changeL(L,n_CO))
#ifdef ACCMIP_LIKE_DIAGS
          taijls(i,j,L,ijlt_COd)=taijls(i,j,L,ijlt_COd)+changeCO*cpd
#endif
        end if       
#ifdef TRACERS_dCO
C -- dC17O --
        changeL(L,n_dC17O)=rdHCH17OplusNO3*pfactor*vol2mass(n_dC17O)
        changedC17O=changeL(L,n_dC17O)*mass2vol(n_dC17O)*bypfactor
        if((trm(i,j,l,n_dC17O)+changeL(l,n_dC17O)) < minKG)then
          changeL(l,n_dC17O) = minKG - trm(i,j,l,n_dC17O)
          changedC17O=changeL(L,n_dC17O)*mass2vol(n_dC17O)*bypfactor
        endif
C -- dC18O --
        changeL(L,n_dC18O)=rdHCH18OplusNO3*pfactor*vol2mass(n_dC18O)
        changedC18O=changeL(L,n_dC18O)*mass2vol(n_dC18O)*bypfactor
        if((trm(i,j,l,n_dC18O)+changeL(l,n_dC18O)) < minKG)then
          changeL(l,n_dC18O) = minKG - trm(i,j,l,n_dC18O)
          changedC18O=changeL(L,n_dC18O)*mass2vol(n_dC18O)*bypfactor
        endif
C -- d13CO --
        changeL(L,n_d13CO)=rdH13CHOplusNO3*pfactor*vol2mass(n_d13CO)
        changed13CO=changeL(L,n_d13CO)*mass2vol(n_d13CO)*bypfactor
        if((trm(i,j,l,n_d13CO)+changeL(l,n_d13CO)) < minKG)then
          changeL(l,n_d13CO) = minKG - trm(i,j,l,n_d13CO)
          changed13CO=changeL(L,n_d13CO)*mass2vol(n_d13CO)*bypfactor
        endif
#endif  /* TRACERS_dCO */
C -- HNO3 --  (HNO3 from gas and het phase rxns )
        changeL(L,n_HNO3)=changeHNO3*pfactor*vol2mass(n_HNO3)
        IF((trm(i,j,L,n_HNO3)+changeL(L,n_HNO3)) < minKG) THEN
          changeL(L,n_HNO3) = minKG - trm(i,j,L,n_HNO3)
          changeHNO3=changeL(L,n_HNO3)*mass2vol(n_HNO3)*bypfactor
        END IF
#ifdef TRACERS_HETCHEM
#ifdef TRACERS_NITRATE
        changeL(L,n_N_d1)=changeN_d1*pfactor*vol2mass(n_N_d1)
!       if(i==36.and.j==28.and.l==1) then
!         write(out_line,*)'Mchange L 2 ', changeL(L,n_N_d1),changeN_d1
!         call write_parallel(trim(out_line),crit=.true.)
!       endif
        IF((trm(i,j,l,n_N_d1)+changeL(l,n_N_d1)) < minKG) THEN
          changeL(l,n_N_d1) = minKG - trm(i,j,l,n_N_d1)
          changeN_d1=changeL(L,n_N_d1)*mass2vol(n_N_d1)*bypfactor
        END IF
        changeL(L,n_N_d2)=changeN_d2*pfactor*vol2mass(n_N_d2)
        IF((trm(i,j,l,n_N_d2)+changeL(l,n_N_d2)) < minKG) THEN
          changeL(l,n_N_d2) = minKG - trm(i,j,l,n_N_d2)
          changeN_d2=changeL(L,n_N_d2)*mass2vol(n_N_d2)*bypfactor
        END IF
        changeL(L,n_N_d3)=changeN_d3*pfactor*vol2mass(n_N_d3)
        IF((trm(i,j,l,n_N_d3)+changeL(l,n_N_d3)) < minKG) THEN
          changeL(l,n_N_d3) = minKG - trm(i,j,l,n_N_d3)
          changeN_d3=changeL(L,n_N_d3)*mass2vol(n_N_d3)*bypfactor
        END IF
#endif  /* TRACERS_NITRATE */
#endif  /* TRACERS_HETCHEM */
C -- N2O5 --  (N2O5 from gas and het phase rxns)
        changeL(L,n_N2O5)=changeN2O5*pfactor*vol2mass(n_N2O5)
        IF((trm(i,j,l,n_N2O5)+changeL(l,n_N2O5)) < minKG) THEN
          changeL(l,n_N2O5) = minKG - trm(i,j,l,n_N2O5)
          changeN2O5=changeL(L,n_N2O5)*mass2vol(n_N2O5)*bypfactor
        END IF
c -- NOx --   (NOx from gas phase rxns)
        changeL(L,n_NOx)=changeNOx*pfactor*vol2mass(n_NOx)
        IF((trm(i,j,l,n_NOx)+changeL(l,n_NOx)) < minKG) THEN
          changeL(l,n_NOx) = minKG - trm(i,j,l,n_NOx)
          changeNOx=changeL(L,n_NOx)*mass2vol(n_NOx)*bypfactor
        END IF
C -- Alkenes --  (Alkenes from gas phase rxns)
        changeL(L,n_Alkenes)=
     &  changeAlkenes*pfactor*vol2mass(n_Alkenes)
        IF((trm(i,j,l,n_Alkenes)+changeL(l,n_Alkenes)) < minKG)THEN
          changeL(l,n_Alkenes) = minKG - trm(i,j,l,n_Alkenes)
          changeAlkenes=changeL(L,n_Alkenes)*mass2vol(n_Alkenes)
     &    *bypfactor
        END IF
#ifdef TRACERS_dCO
C -- d13Calke --  (d13Calke from gas phase rxns)
        changeL(L,n_d13Calke)=
     &  changed13Calke*pfactor*vol2mass(n_d13Calke)
        IF((trm(i,j,l,n_d13Calke)+changeL(l,n_d13Calke)) < minKG)THEN
          changeL(l,n_d13Calke) = minKG - trm(i,j,l,n_d13Calke)
          changed13Calke=changeL(L,n_d13Calke)*mass2vol(n_d13Calke)
     &    *bypfactor
        END IF
#endif  /* TRACERS_dCO */
#ifdef TRACERS_AEROSOLS_SOA
C -- isopp1g --  (isopp1g from gas phase rxns)
        changeL(L,n_isopp1g)=
     &  changeisopp1g*pfactor*vol2mass(n_isopp1g)
        IF((trm(i,j,l,n_isopp1g)+changeL(l,n_isopp1g)) < minKG)THEN
          changeL(l,n_isopp1g) = minKG - trm(i,j,l,n_isopp1g)
          changeisopp1g=changeL(L,n_isopp1g)*mass2vol(n_isopp1g)
     &    *bypfactor
        END IF
C -- isopp2g --  (isopp2g from gas phase rxns)
        changeL(L,n_isopp2g)=
     &  changeisopp2g*pfactor*vol2mass(n_isopp2g)
        IF((trm(i,j,l,n_isopp2g)+changeL(l,n_isopp2g)) < minKG)THEN
          changeL(l,n_isopp2g) = minKG - trm(i,j,l,n_isopp2g)
          changeisopp2g=changeL(L,n_isopp2g)*mass2vol(n_isopp2g)
     &    *bypfactor
        END IF
#ifdef TRACERS_TERP
C -- apinp1g --  (apinp1g from gas phase rxns)
        changeL(L,n_apinp1g)=
     &  changeapinp1g*pfactor*vol2mass(n_apinp1g)
        IF((trm(i,j,l,n_apinp1g)+changeL(l,n_apinp1g)) < minKG)THEN
          changeL(l,n_apinp1g) = minKG - trm(i,j,l,n_apinp1g)
          changeapinp1g=changeL(L,n_apinp1g)*mass2vol(n_apinp1g)
     &    *bypfactor
        END IF
C -- apinp2g --  (apinp2g from gas phase rxns)
        changeL(L,n_apinp2g)=
     &  changeapinp2g*pfactor*vol2mass(n_apinp2g)
        IF((trm(i,j,l,n_apinp2g)+changeL(l,n_apinp2g)) < minKG)THEN
          changeL(l,n_apinp2g) = minKG - trm(i,j,l,n_apinp2g)
          changeapinp2g=changeL(L,n_apinp2g)*mass2vol(n_apinp2g)
     &    *bypfactor
        END IF
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
c -- Isoprene -- (Isoprene from gas phase rxns)
        changeL(L,n_Isoprene)=
     &  changeIsoprene*pfactor*vol2mass(n_Isoprene)
        IF((trm(i,j,l,n_Isoprene)+changeL(l,n_Isoprene)) < minKG)
     &  THEN
          changeL(l,n_Isoprene) = minKG - trm(i,j,l,n_Isoprene)
          changeIsoprene=changeL(L,n_Isoprene)*mass2vol(n_Isoprene)
     &    *bypfactor
        END IF
#ifdef TRACERS_TERP
c -- Terpenes -- (Terpenes from gas phase rxns)
        changeL(L,n_Terpenes)=
     &  changeTerpenes*pfactor*vol2mass(n_Terpenes)
        IF((trm(i,j,l,n_Terpenes)+changeL(l,n_Terpenes)) < minKG)
     &  THEN
          changeL(l,n_Terpenes) = minKG - trm(i,j,l,n_Terpenes)
          changeTerpenes=changeL(L,n_Terpenes)*mass2vol(n_Terpenes)
     &    *bypfactor
        END IF
#endif  /* TRACERS_TERP */
c -- AlkylNit -- (AlkylNit from gas phase rxns)
        changeL(L,n_AlkylNit)=
     &  changeAlkylNit*pfactor*vol2mass(n_AlkylNit)
        IF((trm(i,j,l,n_AlkylNit)+changeL(l,n_AlkylNit)) < minKG)
     &  THEN
          changeL(l,n_AlkylNit) = minKG - trm(i,j,l,n_AlkylNit)
          changeAlkylNit=changeL(L,n_AlkylNit)*mass2vol(n_AlkylNit)
     &    *bypfactor
        END IF

C Save 3D radical arrays to pass to aerosol code:
C Make sure we get the nightime values; Set OH to zero for now:
        if(coupled_chem == 1) then
          oh_live(i,j,L)=0.d0
          no3_live(i,j,L)=yNO3(i,j,L)
        end if

c --  Ox --   ( Ox from gas phase rxns)
        changeOx=-1.d0*rr(rrbi%NO2_O3__NO3_O2,L)*y(nNO2,L)*y(nn_Ox,L)
     &    *pOx(I,J,L)*dt2
        changeL(L,n_Ox)=changeOx*pfactor*vol2mass(n_Ox)
        IF((trm(i,j,L,n_Ox)+changeL(L,n_Ox)) < minKG) THEN
          changeL(L,n_Ox) = minKG - trm(i,j,L,n_Ox)
          changeOx=changeL(L,n_Ox)*mass2vol(n_Ox)*bypfactor
        END IF
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
        if(trm(i,j,L,n_Ox)==0.)call stop_model('zero Ox denom',255)
        changeL(L,n_stratOx)=changeL(L,n_Ox)*
     &  trm(i,j,L,n_stratOx)/trm(i,j,L,n_Ox)
        if((trm(i,j,L,n_stratOx)+changeL(L,n_stratOx)) < minKG)
     &  changeL(L,n_stratOx) = minKG - trm(i,j,L,n_stratOx)
#endif
        ! then come diags:
        if(changeL(L,n_Ox) >= 0.) then
          CALL INC_TAJLS(I,J,L,jls_Oxp,changeL(L,n_Ox))
          if(L<=maxT)CALL INC_TAJLS(I,J,L,jls_OxpT,changeL(L,n_Ox))
#ifdef ACCMIP_LIKE_DIAGS
          taijls(i,j,L,ijlt_Oxp)=taijls(i,j,L,ijlt_Oxp)+changeOx*cpd
#endif
        else
          CALL INC_TAJLS(I,J,L,jls_Oxd,changeL(L,n_Ox))
          if(L<=maxT)CALL INC_TAJLS(I,J,L,jls_OxdT,changeL(L,n_Ox))
#ifdef ACCMIP_LIKE_DIAGS
          taijls(i,j,L,ijlt_Oxd)=taijls(i,j,L,ijlt_Oxd)+changeOx*cpd
#endif
        end if
c -- ClONO2 --   (ClONO2 from gas and het phase rxns)
        changeL(L,n_ClONO2)=changeClONO2*pfactor*
     &  vol2mass(n_ClONO2)
        IF((trm(i,j,l,n_ClONO2)+changeL(l,n_ClONO2)) < minKG) THEN
          changeL(l,n_ClONO2) = minKG - trm(i,j,l,n_ClONO2)
          changeClONO2=changeL(L,n_ClONO2)*mass2vol(n_ClONO2)*
     &    bypfactor
        END IF
c -- ClOx --   (ClOx from gas and het phase rxns)
        changeL(L,n_ClOx)=changeClOx*pfactor*vol2mass(n_ClOx)
        IF((trm(i,j,l,n_ClOx)+changeL(l,n_ClOx)) < minKG) THEN
          changeL(l,n_ClOx) = minKG - trm(i,j,l,n_ClOx)
          changeClOx=changeL(L,n_ClOx)*mass2vol(n_ClOx)*bypfactor
        END IF
        if(pscX(L))then
c -- HOCl --   (HOCl from het phase rxns)
          changeL(L,n_HOCl)=changeHOCl*pfactor*vol2mass(n_HOCl)
          IF((trm(i,j,l,n_HOCl)+changeL(l,n_HOCl)) < minKG) THEN
            changeL(l,n_HOCl) = minKG - trm(i,j,l,n_HOCl)
            changeHOCl=changeL(L,n_HOCl)*mass2vol(n_HOCl)*bypfactor
          END IF
c -- HCl --   (HCl from het phase rxns)
          changeL(L,n_HCl)=changeHCl*pfactor*vol2mass(n_HCl)
          IF((trm(i,j,l,n_HCl)+changeL(l,n_HCl)) < minKG) THEN
            changeL(l,n_HCl) = minKG - trm(i,j,l,n_HCl)
            changeHCl=changeL(L,n_HCl)*mass2vol(n_HCl)*bypfactor
          END IF
c -- HBr --   (HBr from gas phase rxns)
          changeL(L,n_HBr)=changeHBr*pfactor*vol2mass(n_HBr)
          IF((trm(i,j,l,n_HBr)+changeL(l,n_HBr)) < minKG) THEN
            changeL(l,n_HBr) = minKG - trm(i,j,l,n_HBr)
            changeHBr=changeL(L,n_HBr)*mass2vol(n_HBr)*bypfactor
          END IF
c -- BrOx --   (BrOx from gas phase rxns)
          changeL(L,n_BrOx)=changeBrOx*pfactor*vol2mass(n_BrOx)
          IF((trm(i,j,l,n_BrOx)+changeL(l,n_BrOx)) < minKG) THEN
            changeL(l,n_BrOx) = minKG - trm(i,j,l,n_BrOx)
            changeBrOx=changeL(L,n_BrOx)*mass2vol(n_BrOx)*bypfactor
          END IF
c -- BrONO2 --   (BrONO2 from gas phase rxns)
          changeL(L,n_BrONO2)=changeBrONO2*pfactor*vol2mass(n_BrONO2)
          IF((trm(i,j,l,n_BrONO2)+changeL(l,n_BrONO2)) < minKG) THEN
            changeL(l,n_BrONO2) = minKG - trm(i,j,l,n_BrONO2)
           changeBrONO2=changeL(L,n_BrONO2)*mass2vol(n_BrONO2)*bypfactor
          END IF
        endif  ! PSCs exist

        call printNightChemistryDiags()
        call checkNighttimeTolerances()

C       ACCUMULATE 3D NO3 diagnostic: 
        if (yNO3(I,J,L) > 0.d0 .and. yNO3(I,J,L) < 1.d20)
     &  taijls(i,j,L,ijlt_NO3)=taijls(i,j,L,ijlt_NO3)+yNO3(i,j,L)

        if (y(nClO,L) > 0.d0 .and. y(nClO,L) < 1.d20)
     &  CALL INC_TAJLS2(I,J,L,jls_ClOcon,y(nClO,L)/y(nM,L))
        if (y(nH2O,L) > 0.d0 .and. y(nH2O,L) < 1.d20)
     &  CALL INC_TAJLS2(I,J,L,jls_H2Ocon,y(nH2O,L)/y(nM,L))
     
       end do  ! L loop <===========


       ! aerosol code uses radicals up to LM, so fill in above
       ! chemistry (nighttime case):
       if(coupled_chem == 1) then
         do L=topLevelOfChemistry+1,LM
           oh_live(i,j,L)=0.d0
           no3_live(i,j,L)=0.d0
         end do
       end if

CCCCCCCCCCCCCCCC END NIGHTTIME CCCCCCCCCCCCCCCCCCCC

      end if
CCCCCCCCCCCCCCCCCCCC END DARKNESS CCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      save_NO2column(i,j)=0.d0 ! initialize sum outside L loop.

      DO L=1,topLevelOfChemistry

C Lower limit on HO2NO2 : 
        if(trm(i,j,l,n_HO2NO2)+changeL(l,n_HO2NO2) < minKG)
     &  changeL(l,n_HO2NO2) = minKG - trm(i,j,l,n_HO2NO2)

c Tropospheric halogen sink Br & Cl :
        if(y(nH2O,L)/y(nM,L) > 10.d-6)then ! sink by wet removal in trop
          rmv=0.5d0
          changeL(L,n_ClOx)  =-trm(I,J,L,n_ClOx)  *rmv
          changeL(L,n_HCl)   =-trm(I,J,L,n_HCl)   *rmv
          changeL(L,n_HOCl)  =-trm(I,J,L,n_HOCl)  *rmv
          changeL(L,n_ClONO2)=-trm(I,J,L,n_ClONO2)*rmv
          changeL(L,n_BrOx)  =-trm(I,J,L,n_BrOx)  *rmv
          changeL(L,n_HBr)   =-trm(I,J,L,n_HBr)   *rmv
          changeL(L,n_HOBr)  =-trm(I,J,L,n_HOBr)  *rmv
          changeL(L,n_BrONO2)=-trm(I,J,L,n_BrONO2)*rmv
        else
c Set CLTOT based on CFCs (2.4 ppbv yield from complete oxidation of
c 1.8 ppbv CFC plus 0.8 ppbv background which is tied to methane) :
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !WARNING: RESETTING SOME Y's HERE; SO DON'T USE THEM BELOW!     
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          y(nn_ClOx,L)=(trm(I,J,L,n_ClOx)+changeL(L,n_ClOx))*y(nM,L)*
     &    mass2vol(n_ClOx)*BYAXYP(I,J)*byMA(L,I,J)
          y(nn_HCl,L)= (trm(I,J,L,n_HCl)+changeL(L,n_HCl))*y(nM,L)*
     &    mass2vol(n_HCl)*BYAXYP(I,J)*byMA(L,I,J)
          y(nn_HOCl,L)=(trm(I,J,L,n_HOCl)+changeL(L,n_HOCl))*y(nM,L)*
     &    mass2vol(n_HOCl)*BYAXYP(I,J)*byMA(L,I,J)
          y(nn_ClONO2,L)=(trm(I,J,L,n_ClONO2)+changeL(L,n_ClONO2))*
     &    y(nM,L)*mass2vol(n_ClONO2)*BYAXYP(I,J)*byMA(L,I,J)
          CLTOT=((y(nn_CFC,1)/y(nM,1) -
     &         y(nn_CFC,L)/y(nM,L))*(3.0d0/1.8d0)*
     &    y(nn_CFC,1)/(1.8d-9*y(nM,1)))
          CLTOT=CLTOT+0.8d-9*(y(nn_CH4,1)/y(nM,1)-y(nn_CH4,L)/y(nM,L))/
     &    (y(nn_CH4,1)/y(nM,1))
          CLTOT=CLTOT*y(nM,L)/
     &    (y(nn_ClOx,L)+y(nn_HCl,L)+y(nn_HOCl,L)+y(nn_ClONO2,L))
          if(prnchg.and.J == jprn.and.I == iprn.and.L == lprn)then  
            write(out_line,'("CLTOT = ",F20.5)') CLTOT
            call write_parallel(trim(out_line),crit=jay)
          end if
          IF(CLTOT <= 0.999d0 .OR. CLTOT >= 1.001d0) THEN
            changeL(L,n_ClOx)=changeL(L,n_ClOx)*CLTOT+
     &      trm(I,J,L,n_ClOx)*(CLTOT-1.D0)
            changeL(L,n_HCl)=changeL(L,n_HCl)*CLTOT+
     &      trm(I,J,L,n_HCl)*(CLTOT-1.D0)
            changeL(L,n_HOCl)=changeL(L,n_HOCl)*CLTOT+
     &      trm(I,J,L,n_HOCl)*(CLTOT-1.D0)
c           Conserve N wrt ClONO2 once inital Cl changes past:
            if(Itime-ItimeI >= 6 .OR. allowSomeChemReinit .NE. 1)then ! note logic
              changeL(L,n_NOx)=changeL(L,n_NOx)-
     &        (trm(I,J,L,n_ClONO2)+changeL(L,n_ClONO2))*
     &        (CLTOT-1.D0)*tr_mm(n_NOx)/tr_mm(n_ClONO2)
              if(-changeL(L,n_NOx) > trm(I,J,L,n_NOx))changeL(L,n_NOx)=
     &        -0.8d0*trm(I,J,L,n_NOx)
            endif
            changeL(L,n_ClONO2)=changeL(L,n_ClONO2)*CLTOT+
     &      trm(I,J,L,n_ClONO2)*(CLTOT-1.D0)
          END IF

c Set Total Bromine based on CFCs (4.5 pptv yield
C from complete oxidation of 1.8 ppbv CFC plus 0.5 pptv background) :
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !WARNING: RESETTING SOME Y's HERE; SO DON'T USE THEM BELOW!     
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          y(nn_BrOx,L)=(trm(I,J,L,n_BrOx)+changeL(L,n_BrOx))*y(nM,L)*
     &    mass2vol(n_BrOx)*BYAXYP(I,J)*byMA(L,I,J)
          y(nn_HBr,L)= (trm(I,J,L,n_HBr)+changeL(L,n_HBr))*y(nM,L)*
     &    mass2vol(n_HBr)*BYAXYP(I,J)*byMA(L,I,J)
          y(nn_HOBr,L)=(trm(I,J,L,n_HOBr)+changeL(L,n_HOBr))*y(nM,L)*
     &    mass2vol(n_HOBr)*BYAXYP(I,J)*byMA(L,I,J)
          y(nn_BrONO2,L)=(trm(I,J,L,n_BrONO2)+changeL(L,n_BrONO2))*
     &    y(nM,L)*mass2vol(n_BrONO2)*BYAXYP(I,J)*byMA(L,I,J)
     
          BRTOT=((y(nn_CFC,1)/y(nM,1) - 
     &         y(nn_CFC,L)/y(nM,L))*(4.5d-3/1.8d0)
     &    *y(nn_CFC,1)/(1.8d-9*y(nM,1)))
          BRTOT=BRTOT+0.5d-12*(y(nn_CH4,1)/y(nM,1)-y(nn_CH4,L)/y(nM,L))/
     &    (y(nn_CH4,1)/y(nM,1))
          BRTOT=BRTOT*y(nM,L)/
     &    (y(nn_BrOx,L)+y(nn_HBr,L)+y(nn_HOBr,L)+y(nn_BrONO2,L))
          if(prnchg.and.J == jprn.and.I == iprn.and.L == lprn)then  
            write(out_line,'("BrTOT = ",F20.5)') BRTOT
            call write_parallel(trim(out_line),crit=jay)
          end if
          IF(BRTOT <= 0.999d0 .OR. BRTOT >= 1.001d0) THEN
            changeL(L,n_BrOx)=changeL(L,n_BrOx)*BRTOT+
     &      trm(I,J,L,n_BrOx)*(BRTOT-1.D0)
            changeL(L,n_HBr)=changeL(L,n_HBr)*BRTOT+
     &      trm(I,J,L,n_HBr)*(BRTOT-1.D0)
            changeL(L,n_HOBr)=changeL(L,n_HOBr)*BRTOT+
     &      trm(I,J,L,n_HOBr)*(BRTOT-1.D0)
c           Conserve N wrt BrONO2 once inital Br changes past:
            if(Itime-ItimeI >= 6 .OR. allowSomeChemReinit .NE. 1)then ! note logic
              changeL(L,n_NOx)=changeL(L,n_NOx)-
     &        (trm(I,J,L,n_BrONO2)+changeL(L,n_BrONO2))*
     &        (BRTOT-1.D0)*tr_mm(n_NOx)/tr_mm(n_BrONO2)
              if(-changeL(L,n_NOx) > trm(I,J,L,n_NOx))changeL(L,n_NOx)=
     &        -0.8d0*trm(I,J,L,n_NOx)
            endif
            changeL(L,n_BrONO2)=changeL(L,n_BrONO2)*BRTOT+
     &      trm(I,J,L,n_BrONO2)*(BRTOT-1.D0)
          END IF
        end if ! i.e. y(nH2O,L)/y(nM,L) <= 10.d-6 

#ifdef TRACERS_AEROSOLS_SOA
        pfactor=axyp(I,J)*MA(L,I,J)/y(nM,L)
        bypfactor=1.D0/pfactor
        call soa_aerosolphase(I,J,L,changeL,bypfactor)
#endif  /* TRACERS_AEROSOLS_SOA */

        tempChangeNOx= ! this needed for several diags below:
     &  changeL(L,n_NOx)*mass2vol(n_NOx)*y(nM,L)/(axyp(I,J)*MA(L,I,J))

        tempChangeOx=
     &  changeL(L,n_Ox)*mass2vol(n_Ox)*y(nM,L)/(axyp(I,J)*MA(L,I,J))

! Accumulate NO2 10:30am/1:30pm tropo column diags:
! -- moved from sunlight/darkness sections because needed changeNOx
! -- saved here in molecules/cm2
        if(L<=min(maxT,LTROPO(I,J)))then
          if(daylight)then

            index1=0 ; index2=0

            if(ih1030 < ih1030e)then ! normal case
              if(i>=ih1030.and.i<=ih1030e)then
                index1=ijs_NO2_1030; index2=ijs_NO2_1030c
              end if 
            else                     ! crossing date line
              if(i<=ih1030.or.i>=ih1030e)then
                index1=ijs_NO2_1030; index2=ijs_NO2_1030c
              end if 
            end if
            if(ih1330 < ih1330e)then ! normal case
              if(i>=ih1330.and.i<=ih1330e)then
                index1=ijs_NO2_1330; index2=ijs_NO2_1330c
              end if
            else                     ! crossing date line
              if(i<=ih1330.or.i>=ih1330e)then
                index1=ijs_NO2_1330; index2=ijs_NO2_1330c
              end if
            end if

            if(index1/=0 .and. index2/=0)then
              thick= ! layer thickness in cm
     &        1.d2*rgas*bygrav*TX(I,J,L)*LOG(PEDN(L,i,j)/PEDN(L+1,i,j))
              taijs(i,j,index1)=taijs(i,j,index1)+thick*
     &        pNOx(i,j,L)*(y(nn_NOx,L)+tempChangeNOx)
              if(L==1)taijs(i,j,index2)=taijs(i,j,index2)+1.d0
            end if

          end if ! sunlight criteria
        end if ! troposphere criterion

! Also save instantaneous NO2 tropospheric column for SUBDDiag:
! Conversion is only from molecules/cm3 to molecules/cm2:
! save_NO2column is initialized to 0 outside this L loop.
! [note: we should consolodate all these "thick/byThick" guys.]
        if(L<=min(maxT,LTROPO(I,J)))then
          thick= ! layer thickness in cm
     &    1.d2*rgas*bygrav*TX(I,J,L)*LOG(PEDN(L,i,j)/PEDN(L+1,i,j))
          save_NO2column(i,j) = save_NO2column(i,j)+
     &    thick*pNOx(i,j,L)*(y(nn_NOx,L)+tempChangeNOx)
        end if

#ifdef ACCMIP_LIKE_DIAGS
! accumulate some 3D diagnostics in moles/m3/s units:
        ! chemical_production_of_O1D_from_ozone:
        taijls(i,j,l,ijlt_pO1D)=taijls(i,j,l,ijlt_pO1D)+
     &  ss(rj%O3__O1D_O2,l,i,j)*y(nO3,l)*cpd

        ! chemical_production_of_OH_from_O1D_plus_H2O:
        taijls(i,j,l,ijlt_pOH)=taijls(i,j,l,ijlt_pOH)+
     &  2.d0*rr(rrbi%O1D_H2O__OH_OH,l)*y(nH2O,l)*y(nO1D,l)*cpd

        ! chemical_production_rate_of_ozone_by_HO2_plus_NO:
        taijls(i,j,l,ijlt_OxpHO2)=taijls(i,j,l,ijlt_OxpHO2)+
     &  rr(rrbi%HO2_NO__OH_NO2,l)*y(nHO2,l)*y(nNO,l)*cpd
   
        ! chemical_production_rate_of_ozone_by_CH3O2_plus_NO:
        taijls(i,j,l,ijlt_OxpCH3O2)=taijls(i,j,l,ijlt_OxpCH3O2)+
     &  rr(rrbi%CH3O2_NO__HCHO_NO2,l)*y(nCH3O2,l)*y(nNO,l)*cpd
    
        ! chemical_destruction_rate_of_ozone_by_OH:
        taijls(i,j,l,ijlt_OxlOH)=taijls(i,j,l,ijlt_OxlOH)+
     &  rr(rrbi%OH_O3__HO2_O2,l)*y(nOH,l)*y(nO3,l)*cpd ! (positive)

        !chemical_destruction_rate_of_ozone_by_HO2:
        taijls(i,j,l,ijlt_OxlHO2)=taijls(i,j,l,ijlt_OxlHO2)+
     &  rr(rrbi%HO2_O3__OH_O2,l)*y(nOH,l)*y(nO3,l)*cpd ! (positive)

        !chemical_destruction_rate_of_ozone_by_Alkenes:
        taijls(i,j,l,ijlt_OxlALK)=taijls(i,j,l,ijlt_OxlALK)+
     &  rr(rrbi%Alkenes_O3__HCHO_CO,l)*y(nn_Alkenes,l)*y(nO3,l)*cpd ! (positive)

        !Save 3D NO2 separately from NOx (pppv here):
        ! need to add NOx change to match the NOx tracer diag:
        taijls(i,j,l,ijlt_NO2vmr)=taijls(i,j,l,ijlt_NO2vmr)+
     &  pNOx(i,j,l)*(y(nn_NOx,l)+tempChangeNOx)/y(nM,l)

        !Save 3D NO separately from NOx (pppv here):
        ! need to add NOx change to match the NOx tracer diag:
        pNOloc=1.d0-pNOx(i,j,L)-pNO3(i,j,L)
        if(pNOloc > 0.d0)then
          taijls(i,j,l,ijlt_NOvmr)=taijls(i,j,l,ijlt_NOvmr)+
     &    pNOloc*(y(nn_NOx,l)+tempChangeNOx)/y(nM,l)
        else
          ! avoid accumulating small negatives due to round-off
          ! of local pNO (e.g. at night)
          continue
        end if
#endif

        ! Below there is a 3D O3 diagnostic in cm-atm units for more
        ! direct NINT input. Here try to save it in vmr(ppbv) for humans,
        ! (both JL and IJL). Above top of chemistry accumulate Ox (which is
        ! likely, depending on your deck settings actually O3 from NINT
        ! input anyway) further on in code below:
        if(L <= topLevelOfChemistry) then
          taijls(i,j,L,ijlt_O3ppbv)=taijls(i,j,L,ijlt_O3ppbv)+
     &    1.e9*pOx(i,j,L)*(y(nn_Ox,L)+tempChangeOx)/y(nM,L)
          CALL INC_TAJLS2  ! (V/V air)
     &    (I,J,L,jls_O3vmr,pOx(i,j,L)*(y(nn_Ox,L)+tempChangeOx)/y(nM,L))
        end if

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Save chemistry changes for applying in apply_tracer_3Dsource.  C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

        do n=ntm_chem_beg,ntm_chem_end
          tr3Dsource(i,j,L,nChemistry,n) = changeL(L,n) * bydtsrc
        end do

#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
        tr3Dsource(i,j,L,nChemistry,n_stratOx)=
     &  changeL(L,n_stratOx)*bydtsrc
#endif

        ! save NO2 volume mixing ratio for sub-daily diagnosic:
        ! Note that for a long time the NON-Cached version of this used to
        ! neglect +tempChangeNOx. That term is needed to match the
        ! NOx tracer values:
        mNO2(i,j,L)=pNOx(i,j,L)*(y(nn_NOx,L)+tempChangeNOx)/y(nM,L)
#ifdef CACHED_SUBDD
        mrno2(i,j,L)=pNOx(i,j,L)*(y(nn_NOx,L)+tempChangeNOx)/y(nM,L)
        pNOloc=1.d0-pNOx(i,j,L)-pNO3(i,j,L)
        if(pNOloc > 0.d0)then
          mrno(i,j,L)=pNOloc*(y(nn_NOx,L)+tempChangeNOx)/y(nM,L)
        else
          ! avoid accumulating small negatives due to round-off
          ! of local pNO (e.g. at night)
          mrno(i,j,L)=0.d0
        end if
        mro3(i,j,L)=pOx(i,j,L)*(y(nn_Ox,L)+tempChangeOx)/y(nM,L)
        OH_conc(i,j,l)=y(nOH,L)
        HO2_conc(i,j,l)=y(nHO2,L)
        JO1D_rate(i,j,l)=zj(l,rj%O3__O1D_O2)
        JNO2_rate(i,j,l)=zj(l,rj%NO2__NO_O)
#endif
     
#ifdef TRACERS_HETCHEM
#ifdef TRACERS_NITRATE
        tr3Dsource(i,j,l,nChemistry,n_N_d1) = changeL(l,n_N_d1) *bydtsrc
        tr3Dsource(i,j,l,nChemistry,n_N_d2) = changeL(l,n_N_d2) *bydtsrc
        tr3Dsource(i,j,l,nChemistry,n_N_d3) = changeL(l,n_N_d3) *bydtsrc
#endif  /* TRACERS_NITRATE */
#endif  /* TRACERS_HETCHEM */

      END DO ! end current altitude loop

      call printSS27x2Etc()

#ifdef CACHED_SUBDD
      ! I guess set these to 0 above the chemistry... though that is
      ! not quite right...
      do L=topLevelOfChemistry+1,LM
        mrno2(i,j,L)=0.d0
        mrno(i,j,L)=0.d0
        mro3(i,j,L)=0.d0
        OH_conc(i,j,L)=0.d0
        HO2_conc(i,j,L)=0.d0
        JO1D_rate(i,j,L)=0.d0
        JNO2_rate(i,j,L)=0.d0
      end do
#endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      END DO i_loop ! ===> MAIN I LOOP ENDS <===

      END DO j_loop ! ===> MAIN J LOOP ENDS <===
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

#ifdef CACHED_SUBDD
      call find_groups('taijlh',grpids,ngroups)
      do igrp=1,ngroups
        subdd => subdd_groups(grpids(igrp))
        do k=1,subdd%ndiags
          select case (subdd%name(k))
          case ('MRNO2')
            call inc_subdd(subdd,k,mrno2)
          case ('MRNO')
            call inc_subdd(subdd,k,mrno)
          case ('MRO3')
            call inc_subdd(subdd,k,mro3)
          case ('OH_conc')
            call inc_subdd(subdd,k,OH_conc)
          case ('HO2_conc')
            call inc_subdd(subdd,k,HO2_conc)
          case ('JO1D')
            call inc_subdd(subdd,k,JO1D_rate)
          case ('JNO2')
            call inc_subdd(subdd,k,JNO2_rate)
          end select
        enddo ! k
      enddo ! igroup

      call find_groups('taijph',grpids,ngroups)
      do igrp=1,ngroups
        subdd => subdd_groups(grpids(igrp))
        do k=1,subdd%ndiags
          select case (subdd%name(k))
          case ('MRNO2cp')
            call inc_subdd(subdd,k,mrno2)
          case ('MRNOcp')
            call inc_subdd(subdd,k,mrno)
          case ('MRO3cp')
            call inc_subdd(subdd,k,mro3)
          case ('OH_conccp')
            call inc_subdd(subdd,k,OH_conc)
          case ('HO2_conccp')
            call inc_subdd(subdd,k,HO2_conc)
          end select
        enddo ! k
      enddo ! igroup

      call find_groups('taijh',grpids,ngroups)
      do igrp=1,ngroups
        subdd => subdd_groups(grpids(igrp))
        do k=1,subdd%ndiags
          select case (subdd%name(k))
          case ('MRNO2l1')
            call inc_subdd(subdd,k,mrno2(:,:,1))
          case ('MRNOl1')
            call inc_subdd(subdd,k,mrno(:,:,1))
          case ('MRO3l1max')
            call inc_subdd(subdd,k,mro3(:,:,1))
          end select
        enddo ! k
      enddo ! igroup
#endif

CCCCCCCCCCCCCCCCCC END CHEMISTRY SECTION CCCCCCCCCCCCCCCCCCCCCCCCC

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                 BEGIN OVERWRITING                              C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

C If fix_CH4_chemistry is turned on, reset the CH4 tracer everywhere
C to initial conditions and set the chemistry change to zero...
      if(fix_CH4_chemistry == 1)then
        tr3Dsource(:,:,:,nChemistry,n_CH4) = 0.d0 
        call get_CH4_IC(1)
      end if 

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Special cases of overwriting, when doing stratospheric chemistry C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC     
C N2O, CFC, and optional CH4 L=1 overwriting: with all these "fact"s
C this looks complicated, but basically, you are either converting
C from mixing ratio to KG (normal case) or from cm-atm to KG  
C (interactive radiation case - for more on that conversion, see
C the notes on O3MULT in the TRCHEM_Shindell_COM program):
      PIfact(:)=1.d0     
      if(PI_run == 1)then
        PIfact(n_NOx)=PIratio_N
        if(use_rad_n2o == 0) PIfact(n_N2O)=PIratio_N2O
        if(use_rad_cfc == 0) PIfact(n_CFC)=PIratio_CFC
      endif
      fact2=n2o_pppv  ! default N2O mixing ratio overwrite
      fact3=cfc_pppv  ! default CFC mixing ratio overwrite
      fact7=fact_cfc
      if(use_rad_cfc == 0)fact7=1.d0
      do j=J_0,J_1
        do i=I_0,IMAXJ(j)
          fact6=2.69d20*axyp(i,j)*byavog
          fact1=bymair*MA(1,i,j)*axyp(i,j)
          fact5=fact6 
          fact4=fact6
          if(use_rad_n2o == 0)fact4=fact1 
          if(use_rad_cfc == 0)fact5=fact1
          if(use_rad_n2o > 0)fact2=rad_to_chem(3,1,i,j)
          if(use_rad_cfc > 0)fact3=rad_to_chem(5,1,i,j)
          tr3Dsource(i,j,1,nOverwrite,n_N2O)=(fact2*fact4*
     &    tr_mm(n_N2O)*PIfact(n_N2O) - (trm(i,j,1,n_N2O)+ 
     &    tr3Dsource(i,j,1,nChemistry,n_N2O)*dtsrc))*bydtsrc
          tr3Dsource(i,j,1,nOverwrite,n_CFC)=(fact3*fact5*fact7*
     &    tr_mm(n_CFC)*PIfact(n_CFC) - (trm(i,j,1,n_CFC)+
     &    tr3Dsource(i,j,1,nChemistry,n_CFC)*dtsrc))*bydtsrc
          if(use_rad_ch4 > 0)then
            tr3Dsource(i,j,1,nOverwrite,n_CH4)=(rad_to_chem(4,1,i,j)*
     &      fact6*tr_mm(n_CH4)-(trm(i,j,1,n_CH4)+
     &      tr3Dsource(i,j,1,nChemistry,n_CH4)*dtsrc))*bydtsrc
          end if
        end do
      end do

! Ox, stratOx, NOx, BrOx and ClOx, have overwriting where P < PltOx hPa:

      !(Interpolate BrOx & ClOx altitude-dependence to model resolution.
      ! Note the PRES2 is all LM layers:)
      CALL LOGPINT(LCOalt,PCOalt,BrOxaltIN,LM,PRES2,BrOxalt,.true.)
      CALL LOGPINT(LCOalt,PCOalt,ClOxaltIN,LM,PRES2,ClOxalt,.true.)

      ! Note L=LS1,LM means it is only allowed in the stratosphere
      ! (constant pressure layers) but can occur in levels above
      ! chemistry. Also since pltOx is independent of topLevelOfChemistry,
      ! you can still do chemistry and overwriting on same layers:
      do L=LS1,LM 
        if(pres2(L) < pltOx)then
          do j=J_0,J_1         
            do i=I_0,IMAXJ(j)
              ! -- Ox --
              tr3Dsource(i,j,L,nOverwrite,n_Ox)=(rad_to_chem(1,L,i,j)*
     &        axyp(i,j)*O3MULT - (trm(i,j,L,n_Ox)+
     &        tr3Dsource(i,j,L,nChemistry,n_Ox)*dtsrc))*bydtsrc
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
              ! -- stratOx --
              tr3Dsource(i,j,L,nOverwrite,n_stratOx)=
     &        (rad_to_chem(1,L,i,j)*axyp(i,j)*O3MULT - (
     &        trm(i,j,L,n_stratOx)
     &        +tr3Dsource(i,j,L,nChemistry,n_stratOx)*dtsrc))*bydtsrc
#endif
              ! -- ClOx --
              tr3Dsource(i,j,L,nOverwrite,n_ClOx)=(1.d-11*ClOxalt(l)
     &        *vol2mass(n_CLOx)*MA(L,i,j)*axyp(i,j) - (
     &        trm(i,j,L,n_ClOx)+tr3Dsource(i,j,L,nChemistry,n_ClOx)
     &        *dtsrc))*bydtsrc    
              ! -- BrOx --
              tr3Dsource(i,j,L,nOverwrite,n_BrOx)=(1.d-11*BrOxalt(l)
     &        *vol2mass(n_BrOx)*MA(L,i,j)*axyp(i,j) - (
     &        trm(i,j,L,n_BrOx)+tr3Dsource(i,j,L,nChemistry,n_BrOx)
     &        *dtsrc))*bydtsrc
              ! -- NOx --
              tr3Dsource(i,j,L,nOverwrite,n_NOx)=(75.d-11 !75=1*300*2.5*.1
     &        *MA(L,i,j)*axyp(i,j)*PIfact(n_NOx)-(trm(i,j,L,n_NOx)+ 
     &        tr3Dsource(i,j,L,nChemistry,n_NOx)*dtsrc))*bydtsrc
            end do ! I 
          end do   ! J
        end if    ! pressure
      end do     ! L


CCCCCCCCCCCCCCCCCC END OVERWRITE SECTION CCCCCCCCCCCCCCCCCCCCCC

c Save new tracer O3 and CH4 fields for use in radiation or elsewhere:
c (radiation code wants atm-cm units):
      do j=J_0,J_1
        if(prnchg)DU_O3(J)=0.d0 ! Drew's diagnostic...
        do i=I_0,imaxj(j) 
          do L=1,LM ! all model layers
            ! Pass Ox to the rad code, except...
            chem_tracer_save(1,L,i,j)=(trm(i,j,L,n_Ox) +
     &      (tr3Dsource(i,j,L,nChemistry,n_Ox) + 
     &      tr3Dsource(i,j,L,nOverwrite,n_Ox))*dtsrc)
     &      *byaxyp(i,j)*byO3MULT
            ! ... if on active chemistry level, pass O3 instead.
            ! (likely, depending on rundeck settings, your Ox above
            ! the top of the chemistry is actually NINT O3 anyway):
            if(L <= topLevelOfChemistry) chem_tracer_save(1,L,i,j)=
     &                          pOx(i,j,L)*chem_tracer_save(1,L,i,j)
            chem_tracer_save(2,L,i,j)=(trm(i,j,L,n_CH4) +
     &      (tr3Dsource(i,j,L,nChemistry,n_CH4) + 
     &      tr3Dsource(i,j,L,nOverwrite,n_CH4))*dtsrc)
     &      *byaxyp(i,j)*avog/(tr_mm(n_CH4)*2.69e20)
            if(prnchg)DU_O3(J)=DU_O3(J)+chem_tracer_save(1,L,i,j)
            ! Above 3D O3 diagnostic in ppbv units is saved (for humans to see).
            ! Here do it in atm-cm units for direct NINT input for rad code.
            taijls(i,j,L,ijlt_O3cmatm)=taijls(i,j,L,ijlt_O3cmatm)+
     &      chem_tracer_save(1,L,i,j)
          end do
#if (defined SHINDELL_STRAT_EXTRA) && (defined ACCMIP_LIKE_DIAGS)
          strato3_tracer_save(1:LM,i,j)=(trm(i,j,1:LM,n_stratOx) +
     &    (tr3Dsource(i,j,1:LM,nChemistry,n_stratOx) +
     &    tr3Dsource(i,j,1:LM,nOverwrite,n_stratOx))*dtsrc)
     &    *byaxyp(i,j)*byO3MULT
#endif
          ! accumulate diag for the column sum of O3 mass hopefully similarly to
          ! how taijn Ox_Total_mass is done. Use O3 for chemistry layers, Ox
          ! tracer above (which is likely anyway actually O3 from NINT intput):
          taijs(i,j,ijs_O3mass)=taijs(i,j,ijs_O3mass)+
     &    sum( pOx(i,j,1:topLevelOfChemistry)*
     &    (trm(i,j,1:topLevelOfChemistry,n_Ox)+
     &    (tr3Dsource(i,j,1:topLevelOfChemistry,nChemistry,n_Ox)+
     &    tr3Dsource(i,j,1:topLevelOfChemistry,nOverwrite,n_Ox))
     &    *dtsrc))*byaxyp(i,j)
          taijs(i,j,ijs_O3mass)=taijs(i,j,ijs_O3mass)+
     &    sum( 
     &    (trm(i,j,topLevelOfChemistry+1:LM,n_Ox)+
     &    (tr3Dsource(i,j,topLevelOfChemistry+1:LM,nChemistry,n_Ox)+
     &    tr3Dsource(i,j,topLevelOfChemistry+1:LM,nOverwrite,n_Ox))
     &    *dtsrc))*byaxyp(i,j)

          ! and the above-chemistry O3 (using the Ox tracer which is probably
          ! O3 from NINT anyway):
          do L=topLevelOfChemistry+1,LM
            taijls(i,j,L,ijlt_O3ppbv)=taijls(i,j,L,ijlt_O3ppbv)+
     &      1.e9*(trm(i,j,L,n_Ox)+(tr3Dsource(i,j,L,nChemistry,n_Ox)+
     &      tr3Dsource(i,j,L,nOverwrite,n_Ox))*dtsrc)*byMA(L,i,j)*
     &      byaxyp(i,j)*mass2vol(n_Ox) ! ppbv
            CALL INC_TAJLS2  ! (V/V air)
     &      (I,J,L,jls_O3vmr,(trm(i,j,L,n_Ox)+
     &      (tr3Dsource(i,j,L,nChemistry,n_Ox)+
     &      tr3Dsource(i,j,L,nOverwrite,n_Ox))*dtsrc)*byMA(L,i,j)*
     &      byaxyp(i,j)*mass2vol(n_Ox))
          end do 

        end do ! i
        if(prnchg)DU_O3(J)=1.d3*DU_O3(J)/IMAXJ(J)
      end do   ! j

#ifdef CACHED_SUBDD
      ! (in future branch combine this accumulation with one for taijs;
      ! here being very conservative to not change taijs by roundoff)
      do j=J_0,J_1
        do i=I_0,imaxj(j)
          o3col(i,j)=0.d0
          o3col(i,j)=o3col(i,j)+
     &    sum( pOx(i,j,1:topLevelOfChemistry)*
     &    (trm(i,j,1:topLevelOfChemistry,n_Ox)+
     &    (tr3Dsource(i,j,1:topLevelOfChemistry,nChemistry,n_Ox)+
     &    tr3Dsource(i,j,1:topLevelOfChemistry,nOverwrite,n_Ox))
     &    *dtsrc))*byaxyp(i,j)
          o3col(i,j)=o3col(i,j)+
     &    sum(
     &    (trm(i,j,topLevelOfChemistry+1:LM,n_Ox)+
     &    (tr3Dsource(i,j,topLevelOfChemistry+1:LM,nChemistry,n_Ox)+
     &    tr3Dsource(i,j,topLevelOfChemistry+1:LM,nOverwrite,n_Ox))
     &    *dtsrc))*byaxyp(i,j)
        end do ! i
      end do ! j
      call find_groups('taijh',grpids,ngroups)
      do igrp=1,ngroups
        subdd => subdd_groups(grpids(igrp))
        diag_loop: do k=1,subdd%ndiags
          select case (subdd%name(k))
          case ('O3col')
            call inc_subdd(subdd,k,o3col) ; cycle diag_loop
          end select
        enddo diag_loop
      enddo ! igroup
#endif

      if(prnchg)then
        call PACK_DATA( grid, DU_O3, DU_O3_glob )
        IF(AM_I_ROOT()) THEN
          write(6,*) 'Ozone column fm -90 to +90'
          write(6,'(46(f4.0,1x))') (DU_O3_glob(J),J=1,JM)
        END IF
      end if

      RETURN

      contains


      subroutine printSS27x2Etc()
      integer :: LPRINT
      if(prnchg .and. J == jprn .and. I == iprn) then
        jay = (J >= J_0 .and. J <= J_1)
        write(out_line,*)'O3pO2 means O3prof from O2 Herz & SRB:'
        call write_parallel(trim(out_line),crit=jay)
        write(out_line,*)
     &  'L, O3pO2, O3pO2*C, OHpptv, HO2pptv, O/O3, NO2/NO, Cl/ClO:'
        call write_parallel(trim(out_line),crit=jay)
        do LPRINT=LS1,topLevelOfChemistry
          if(daylight)then
            ss27x2=2.d0*ss(rj%O2__O_O,LPRINT,i,j)*y(nO2,LPRINT)
     &        *(rr(rrtri%O_O2__O3_M,LPRINT)*y(nO2,LPRINT))
     &        /(rr(rrtri%O_O2__O3_M,LPRINT)*y(nO2,LPRINT)
     &          +rr(rrbi%O_O3__O2_O2,LPRINT)*y(nO3,LPRINT))
          else
            ss27x2=0.d0
          end if
          ss27x2_c=ss27x2*DCOS(SZA*radian) ! prob. no longer wanted
          OHpptv=1.d12*y(nOH,LPRINT)/y(nM,LPRINT)
          HO2pptv=1.d12*y(nHO2,LPRINT)/y(nM,LPRINT)
          ObyO3=y(nO,LPRINT)/y(nO3,LPRINT)
          NO2byNO=y(nNO2,LPRINT)/y(nNO,LPRINT)
          ClbyClO=y(nCl,LPRINT)/y(nClO,LPRINT)
          write(out_line,'(I3,7(1X,E20.5))')
     &    LPRINT,ss27x2,ss27x2_c,OHpptv,HO2pptv,ObyO3,NO2byNO,ClbyClO
          call write_parallel(trim(out_line),crit=jay)
        end do
      end if
      end subroutine printSS27x2Etc


CCCCCCCCCCCCC PRINT SOME CHEMISTRY DIAGNOSTICS CCCCCCCCCCCCCCCC
      subroutine printDaytimeChemistryDiags()
      if(prnchg .and. J == jprn .and. I == iprn) then
       jay = (J >= J_0 .and. J <= J_1) 
       if(lprn <= topLevelOfChemistry) then
         write(out_line,*) ' '
         call write_parallel(trim(out_line),crit=jay)
         write(out_line,*) 'Family ratios at I,J,L: ',i,j,lprn
         call write_parallel(trim(out_line),crit=jay)
         write(out_line,*) 'OH/HO2 = ',y(nOH,lprn)/y(nHO2,lprn)
         call write_parallel(trim(out_line),crit=jay)
         write(out_line,*) 'O/O3 = ',y(nO,lprn)/y(nO3,lprn)
         call write_parallel(trim(out_line),crit=jay)
         write(out_line,*) 'O1D/O3 = ',y(nO1D,lprn)/y(nO3,lprn),
     &    '  J(O1D) = ',ss(rj%O3__O1D_O2,lprn,I,J)
         call write_parallel(trim(out_line),crit=jay)
         write(out_line,*) 'NO/NO2 = ',y(nNO,lprn)/y(nNO2,lprn),
     &    '   J(NO2) = ',ss(rj%NO2__NO_O,lprn,I,J)
         call write_parallel(trim(out_line),crit=jay)
         write(out_line,*) 'conc OH = ',y(nOH,lprn)
         call write_parallel(trim(out_line),crit=jay)
         write(out_line,*) 'Cl,ClO,Cl2O2,OClO,Cl2 = ',y(nCl,lprn),
     &    y(nClO,lprn),y(nCl2O2,lprn),y(nOClO,lprn),y(nCl2,lprn)
         call write_parallel(trim(out_line),crit=jay)
         write(out_line,*) 'Br,BrO = ',y(nBr,lprn),y(nBrO,lprn)
         call write_parallel(trim(out_line),crit=jay)
         write(out_line,*) 'pCl,pClO,pOClO,pBrO = ',pClx(I,J,lprn),
     &    pClOx(I,J,lprn),pOClOx(I,J,lprn),pBrOx(I,J,lprn)
         call write_parallel(trim(out_line),crit=jay)
         write(out_line,*)
     &   'sun, SALBFJ,sza,I,J,Itime= ',albedoToUse,sza,I,J,Itime
         call write_parallel(trim(out_line),crit=jay)
       end if
      end if
      end subroutine printDaytimeChemistryDiags


      subroutine printNightChemistryDiags()
        if(prnchg.and.J == jprn.and.I == iprn.and.L == lprn)then
          jay = (J >= J_0 .and. J <= J_1)
          write(out_line,*)
     &    'dark, SALBFJ,sza,I,J,L,Itime= ',albedoToUse,sza,I,J,L,Itime
          call write_parallel(trim(out_line),crit=jay)
          if(pscX(L))then
            write(out_line,*) 'There are PSCs, T =',ta(L)
            call write_parallel(trim(out_line),crit=jay)
          else
            write(out_line,*) 'There are no PSCs, T =',ta(L)
            call write_parallel(trim(out_line),crit=jay)
          endif
          write(out_line,198) ay(nn_NOx),': ',
     &    changeNOx,' molecules produced; ',
     &    100.d0*(changeNOx)/y(nn_NOx,L),' percent of'
     &    ,y(nn_NOx,L),'(',1.d9*y(nn_NOx,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(nn_HNO3),': ',
     &    changeHNO3,' molecules produced; ',
     &    100.d0*(changeHNO3)/y(nn_HNO3,L),' percent of'
     &    ,y(nn_HNO3,L),'(',1.d9*y(nn_HNO3,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#ifdef TRACERS_HETCHEM
          write(out_line,198) ay(nn_HNO3),': ',
     &    (-krate(i,j,l,1,1)*y(nn_HNO3,l)*dt2),' molecules dest dust ',
     &    (100.d0*(-krate(i,j,l,1,1)*y(nn_HNO3,l)*dt2))/y(nn_HNO3,L),
     &    ' percent of'
     &    ,y(nn_HNO3,L),'(',1.d9*y(nn_HNO3,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#endif
          write(out_line,198) ay(nn_N2O5),': ',
     &    changeN2O5,' net molec produced; ',
     &    100.d0*(changeN2O5)/y(nn_N2O5,L),' percent of'
     &    ,y(nn_N2O5,L),'(',1.d9*y(nn_N2O5,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(nn_N2O5),': ',
     &    gwprodN2O5,' molec prod fm gas;  ',
     &    100.d0*(gwprodN2O5)/y(nn_N2O5,L),' percent of'
     &    ,y(nn_N2O5,L),'(',1.d9*y(nn_N2O5,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(nn_N2O5),': ',
     &    -wprod_sulf,' molec prod fm sulf; ',
     &    -100.d0*(wprod_sulf)/y(nn_N2O5,L),' percent of'
     &    ,y(nn_N2O5,L),'(',1.d9*y(nn_N2O5,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(nn_HCHO),': ',
     &    wprodHCHO,' molecules produced; ',
     &    100.d0*(wprodHCHO)/y(nn_HCHO,L),' percent of'
     &    ,y(nn_HCHO,L),'(',1.d9*y(nn_HCHO,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) 'Aldehyde',': ',
     &    changeAldehyde,' molecules produced; ',
     &    100.d0*(changeAldehyde)/yAldehyde(I,J,L),' percent of'
     &    ,yAldehyde(I,J,L),'(',1.d9*yAldehyde(I,J,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#ifdef TRACERS_dCO
          write(out_line,198) 'd17Oald ',': ',
     &    changed17Oald,' molecules produced; ',
     &    100.d0*(changed17Oald)/yd17Oald(I,J,L),' percent of'
     &    ,yd17Oald(I,J,L),'(',1.d9*yd17Oald(I,J,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) 'd18Oald ',': ',
     &    changed18Oald,' molecules produced; ',
     &    100.d0*(changed18Oald)/yd18Oald(I,J,L),' percent of'
     &    ,yd18Oald(I,J,L),'(',1.d9*yd18Oald(I,J,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) 'd13Cald ',': ',
     &    changed13Cald,' molecules produced; ',
     &    100.d0*(changed13Cald)/yd13Cald(I,J,L),' percent of'
     &    ,yd13Cald(I,J,L),'(',1.d9*yd13Cald(I,J,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#endif  /* TRACERS_dCO */
          write(out_line,198) 'Alkenes ',': ',
     &    changeAlkenes,' molecules produced; ',
     &    100.d0*(changeAlkenes)/y(nn_Alkenes,L),' percent of'
     &    ,y(nn_Alkenes,L),'(',1.d9*y(nn_Alkenes,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#ifdef TRACERS_dCO
          write(out_line,198) 'd13Calke',': ',
     &    changed13Calke,' molecules produced; ',
     &    100.d0*(changed13Calke)/y(nn_d13Calke,L),' percent of'
     &    ,y(nn_d13Calke,L),'(',1.d9*y(nn_d13Calke,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#endif  /* TRACERS_dCO */
#ifdef TRACERS_AEROSOLS_SOA
          write(out_line,198) 'isopp1g ',': ',
     &    changeisopp1g,' molecules produced; ',
     &    100.d0*(changeisopp1g)/y(nn_isopp1g,L),' percent of'
     &    ,y(nn_isopp1g,L),'(',1.d9*y(nn_isopp1g,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) 'isopp2g ',': ',
     &    changeisopp2g,' molecules produced; ',
     &    100.d0*(changeisopp2g)/y(nn_isopp2g,L),' percent of'
     &    ,y(nn_isopp2g,L),'(',1.d9*y(nn_isopp2g,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#ifdef TRACERS_TERP
          write(out_line,198) 'apinp1g ',': ',
     &    changeapinp1g,' molecules produced; ',
     &    100.d0*(changeapinp1g)/y(nn_apinp1g,L),' percent of'
     &    ,y(nn_apinp1g,L),'(',1.d9*y(nn_apinp1g,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) 'apinp2g ',': ',
     &    changeapinp2g,' molecules produced; ',
     &    100.d0*(changeapinp2g)/y(nn_apinp2g,L),' percent of'
     &    ,y(nn_apinp2g,L),'(',1.d9*y(nn_apinp2g,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#endif  /* TRACERS_TERP */
#endif  /* TRACERS_AEROSOLS_SOA */
          write(out_line,198) 'Isoprene',': ',
     &    changeIsoprene,' molecules produced; ',
     &    100.d0*(changeIsoprene)/y(nn_Isoprene,L),' percent of'
     &    ,y(nn_Isoprene,L),'(',1.d9*y(nn_Isoprene,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#ifdef TRACERS_TERP
          write(out_line,198) 'Terpenes',': ',
     &    changeTerpenes,' molecules produced; ',
     &    100.d0*(changeTerpenes)/y(nn_Terpenes,L),' percent of'
     &    ,y(nn_Terpenes,L),'(',1.d9*y(nn_Terpenes,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
#endif  /* TRACERS_TERP */
          write(out_line,198) 'AlkylNit',': ',
     &    changeAlkylNit,' molecules produced; ',
     &    100.d0*(changeAlkylNit)/y(nn_AlkylNit,L),' percent of'
     &    ,y(nn_AlkylNit,L),'(',1.d9*y(nn_AlkylNit,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(nn_ClONO2),': ',
     &    changeClONO2,' molecules produced; ',
     &    100.d0*(changeClONO2)/y(nn_ClONO2,L),' percent of'
     &    ,y(nn_ClONO2,L),'(',1.d9*y(nn_ClONO2,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(nn_ClOx),': ',
     &    changeClOx,' molecules produced; ',
     &    100.d0*(changeClOx)/y(nn_ClOx,L),' percent of'
     &    ,y(nn_ClOx,L),'(',1.d9*y(nn_ClOx,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(nn_HOCl),': ',
     &    changeHOCl,' molecules produced; ',
     &    100.d0*(changeHOCl)/y(nn_HOCl,L),' percent of'
     &    ,y(nn_HOCl,L),'(',1.d9*y(nn_HOCl,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(nn_HCl),': ',
     &    changeHCl,' molecules produced; ',
     &    100.d0*(changeHCl)/y(nn_HCl,L),' percent of'
     &    ,y(nn_HCl,L),'(',1.d9*y(nn_HCl,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(nn_BrONO2),': ',
     &    changeBrONO2,' molecules produced; ',
     &    100.d0*(changeBrONO2)/y(nn_BrONO2,L),' percent of'
     &    ,y(nn_BrONO2,L),'(',1.d9*y(nn_BrONO2,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(nn_BrOx),': ',
     &    changeBrOx,' molecules produced; ',
     &    100.d0*(changeBrOx)/y(nn_BrOx,L),' percent of'
     &    ,y(nn_BrOx,L),'(',1.d9*y(nn_BrOx,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)
          write(out_line,198) ay(nn_HBr),': ',
     &    changeHBr,' molecules produced; ',
     &    100.d0*(changeHBr)/y(nn_HBr,L),' percent of'
     &    ,y(nn_HBr,L),'(',1.d9*y(nn_HBr,L)/y(nM,L),' ppbv)'
          call write_parallel(trim(out_line),crit=jay)

          write(out_line,199) 'NO2, NO3  = ',y(nNO2,L),yNO3(I,J,L)
          call write_parallel(trim(out_line),crit=jay)
          ! diagnose N conservation:
          changeN=changeNOx+2.d0*changeN2O5+changeHNO3+
     &    changeClONO2+changeBrONO2+changeAlkylNit
          write(out_line,197) '===> N conservation: changeN= ',changeN
          call write_parallel(trim(out_line),crit=jay)
        endif
 197    format(a30,e13.3)
 198    format(1x,a8,a2,e13.3,a21,f10.0,a11,2x,e13.3,3x,a1,f12.5,a6)
 199    format(1x,a20,2(2x,e13.3))
        end subroutine printNightChemistryDiags
CCCCCCCCCCCCCCCCCCCC END CHEM DIAG SECT CCCCCCCCCCCCCCCCCCCCCCC


        subroutine checkNighttimeTolerances()
C Make sure nighttime chemistry changes are not too big:
        error=.false.
        if(changeNOx < -1.d15.OR.changeNOx > 1.d15) then
          write(6,*) 'Big chg@ Itime,I,J,L,NOx ',Itime,I,J,L,changeNOx
          write(6,*) 'rlossN,rprodN,ratioN =',rlossN,rprodN,ratioN
          error=.true.
        end if
        if(changeHNO3 < -1.d15.OR.changeHNO3 > 1.d15) then
          write(6,*) 'Big chg@ Itime,I,J,L,HNO3',Itime,I,J,L,changeHNO3
          error=.true.
        end if
        if(changeN2O5 < -1.d15.OR.changeN2O5 > 1.d15) then
          write(6,*) 'Big chg@ Itime,I,J,L,N2O5',Itime,I,J,L,changeN2O5
          error=.true.
        end if
        if(wprodHCHO < -1.d15.OR.wprodHCHO > 1.d15) then
          write(6,*)'Big chg@ Itime,I,J,L,HCHO',Itime,I,J,L,wprodHCHO
          error=.true.
        endif
        if(error)call stop_model('nighttime chem: big changes',255)
        end subroutine checkNighttimeTolerances

      END SUBROUTINE masterchem



      subroutine photo_acetone(I,J,sza)
!@sum calculate photolysis rate for acetone geometrically
!@+ taken from the UK Harwell Model
!@auth Drew Shindell (modelEifications by Greg Faluvegi)
!
      use resolution, only : ls1=>ls1_nominal
      use resolution, only : lm
      use model_com, only: modelEclock
      use atm_com, only: pmid
      use geom, only:  lat2d ! lat is in radians
      use constant, only: twopi,pi,radian,teeny
      use TimeConstants_mod, only: HOURS_PER_DAY, DAYS_PER_YEAR
      use TRCHEM_Shindell_COM, only: topLevelOfChemistry,Jacet

!@var dec declination angle of the earth
!@var lha local hour angle
!@var CC COS(lat)*COS(dec)
!@var SS SIN(lat)*SIN(dec)

      real*8, parameter :: C1=9.269d-7,C2=1.563d0,C3=0.301d0
      real*8 :: dec,lha,CC,SS,sec_func,Jacet0
      real*8, intent(in):: sza ! passed in radians
      integer, intent(in) :: i,j
      integer :: L

      dec=radian*23.455d0*COS( ((modelEclock%getDayOfYear()-173)*twopi)
     * /DAYS_PER_YEAR )
      lha=twopi*real(modelEclock%getHour())/HOURS_PER_DAY

      CC=COS(lat2d(I,J))*COS(dec)
      SS=SIN(lat2d(I,J))*SIN(dec)

      sec_func=1.d0/max(teeny,COS(lha)*CC+SS)
 
      Jacet0=max(0.d0,C1*(COS(sza)*C2)*EXP(-1.*C3*sec_func))
      Jacet(:)=0.d0
      do L=1,min(LS1-1,topLevelOfChemistry)
        Jacet(L)=3.d0*Jacet0/LOG(pmid(L,i,j))
      enddo
      
      return
      end subroutine photo_acetone


      SUBROUTINE Crates(I,J)
!@sum Crates calculate chemical reaction rates for each altitude,
!@+   using JPL 00.  Includes special calculations for pressure
!@+   dependent reactions. Specifically:
!@+   rrbi%CO_OH__HO2_O2, #15 HO2+HO2->H2O2+O2, #16 OH+HNO3->H2O+NO3,
!@+   and reactions #29, and #42.
!@auth Drew Shindell (modelEifications by Greg Faluvegi)

C**** GLOBAL parameters and variables:
      USE RESOLUTION, only  : LM
      USE MODEL_COM, only: Itime,ItimeI
      USE TRACER_COM, only: coupled_chem,trm,nn_N2O5,n_N2O5,n_SO4
      USE TRDIAG_COM, only: jls_N2O5sulf
      USE OldTracer_mod, only: vol2mass
      USE RAD_COM, only  : rad_to_chem
      USE CONSTANT, only : PI, pN2
      USE ATM_COM, only : MA, PMIDL00
      USE TRCHEM_Shindell_COM, only: n_bi,n_tri,n_nst,n_het,ta,ea,rr,pe,
     & cboltz,r1,sb,nst,y,nM,nH2O,ro,sn,which_trop,sulfate,RKBYPIM,dt2,
     & RGAMMASULF,pscX,topLevelOfChemistry,rh,bythick,aero,rrbi,rrhet,
     & boltAvog8byPi

#ifdef TRACERS_AEROSOLS_SOA
      USE TRACER_COM, only: n_isopp1a,n_isopp2a
#ifdef TRACERS_TERP
     &                     ,n_apinp1a,n_apinp2a
#endif  /* TRACERS_TERP */
#ifdef TRACERS_AEROSOLS_Koch
     &                     ,n_SO4
#endif
      USE TRACERS_SOA, only: KpCALC,kpart,kpart_ref,kpart_temp_ref,
     &                       whichsoa,dH_isoprene,dH_apinene
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AMP
      USE TRACER_COM, only: n_M_AKK_SU,n_M_ACC_SU,n_M_DD1_SU
     &                     ,n_M_DS1_SU,n_M_DD2_SU,n_M_DS2_SU
     &                     ,n_M_SSA_SU,n_M_OCC_SU,n_M_BC1_SU
     &                     ,n_M_BC2_SU,n_M_BC3_SU,n_M_DBC_SU
     &                     ,n_M_BOC_SU,n_M_BCS_SU,n_M_MXX_SU
#endif
#ifdef TRACERS_TOMAS
      USE TRACER_COM, only: n_ASO4,nbins
#endif
      USE GEOM, only : lat2d_dg,byaxyp,axyp

      IMPLICIT NONE

C**** Local parameters and variables and arguments:
!@var I,J passed horizontal position indicies
!@var dd,pp,fw,rkp,rk2,rk3M,nb,rrrr,temp dummy "working" variables
!@var L,jj dummy loop variables
!@var byta reciprocal of the local temperature
!@var rkext aerosol extinction from SAGE obs
!@var pscEx NAT PSC surface conc per unit volume (cm^2/cm^3)
!@var beta branching ratio for (HO2+NO) reactions
!@var pcon variable for some pressure conversions
      REAL*8:: byta,dd,pp,fw,rkp,rk2,rk3M,rrrr,temp,beta,pcon,waterPPMV
      real*8 :: associationReaction, activationReaction,pfactor,
     & bypfactor,k0T,k0TM,kinfT,kinfTbyM,RVELN2O5,wprod_sulf,prod_sulf
      INTEGER             :: L,jj,nb
      INTEGER, INTENT(IN) :: I,J
!@var PRES local nominal pressure
!@var LAXt,LAXb lowest and highest levels to have nonzero 
!@+   RAD-code aerosol extinction 
!@var aero array =1 for nonzero rkext, otherwise 0.
      REAL*8, DIMENSION(LM) :: PRES ! = PMIDL00(1:LM). Keeps LM dimension not top of chem
      INTEGER               :: LAXt,LAXb
      real*8, allocatable, dimension(:) :: PSCEX,rkext

      allocate( PSCEX(topLevelOfChemistry) )
      allocate( rkext(topLevelOfChemistry) )

      aero(:)=0
      PRES(1:LM) = PMIDL00(1:LM)
      rkext(:)=0.d0 ! initialize over L
      LAXb=0
      LAXt=0
      if(rad_to_chem(2,1,I,J) /= 0.)call stop_model('kext prob 0',255)
      do L=2,topLevelOfChemistry
        if(rad_to_chem(2,L,I,J) /= 0..and.rad_to_chem(2,L-1,I,J) == 0.)
     &  LAXb=L
        if(rad_to_chem(2,L,I,J) == 0..and.rad_to_chem(2,L-1,I,J) /= 0.)
     &  LAXt=L-1
        if(L==topLevelOfChemistry)then
          if(LAXb > 0 .and. LAXt==0)then   
            if(rad_to_chem(2,L,I,J) /= 0.)then
              LAXt=L
            else
              call stop_model('LAXt failure.',13)     
            end if 
          end if
        end if
      end do
      if(LAXb < 0.or.LAXb > topLevelOfChemistry.or.
     &LAXt < 0.or.LAXt > topLevelOfChemistry) 
     &call stop_model('LAXb or LAXt problem in TRCHEM_master',13)

      do L=1,topLevelOfChemistry            !  ==> BEGIN ALTITUDE LOOP <==
        byta=1.d0/ta(L)
        pcon=y(nM,L)*ta(L)*cboltz/1013.d0
        do jj=1,n_bi+n_nst             ! bimolecular rates start
          rr(jj,L)=pe(jj)*exp(-ea(jj)*byta)
          if (jj==rrbi%O1D_M__O_M) then
!           M is really N2
            rr(jj,L)=rr(jj,L)*pN2
          else if (jj==rrbi%CH4_OH__H2O_CH3O2
#ifdef TRACERS_dCO
     &        .or. jj==rrbi%CH4_OH__H2O_dCH317O2
     &        .or. jj==rrbi%CH4_OH__H2O_dCH318O2
     &        .or. jj==rrbi%CH4_OH__H2O_d13CH3O2
#endif  /* TRACERS_dCO */
     &           ) then
!           based on three-parameters from JPL2011
            rr(jj,L)=rr(jj,L)*(ta(L)**0.667)
          else if (jj==rrbi%CO_OH__HO2_O2
#ifdef TRACERS_dCO
     &        .or. jj==rrbi%dC17O_OH__HO2_O2
     &        .or. jj==rrbi%dC18O_OH__HO2_O2
     &        .or. jj==rrbi%d13CO_OH__HO2_O2
#endif  /* TRACERS_dCO */
     &           ) then
!           based on termolecular reaction from JPL2011
!           (see pages 185-188 and note D1)
            k0TM=y(nM,L)*pe(jj)*((300.d0*byta)**1.4)
            kinfT=1.1d-12*(300.d0*byta)**(-1.3)
            dd=k0TM/kinfT
            pp=0.6d0**(1.d0/(1.d0+(log10(dd))**2.))
            associationReaction=(k0TM/(1.d0+dd))*pp
            k0T=1.5d-13*((300.d0*byta)**(-0.6))
            kinfTbyM=(2.1d9*((300.d0*byta)**(-6.1)))/y(nM,L)
            dd=k0T/kinfTbyM
            pp=0.6d0**(1.d0/(1.d0+(log10(dd))**2.))
            activationReaction=(k0T/(1.d0+dd))*pp
            rr(jj,L)=associationReaction+activationReaction
          else if (jj==rrbi%HO2_HO2__H2O2_O2) then
!           k=(kc+kp)fw, kc=rr
            rkp=2.1d-33*y(nM,L)*exp(920.d0*byta)
            fw=(1.d0+1.4d-21*y(nH2O,L)*exp(2200.d0*byta))
            rr(jj,L)=(rr(jj,L)+rkp)*fw
          else if (jj==rrbi%OH_HNO3__H2O_NO3) then
!           k=[pe*exp(-e(jj)/ta(l))]+k3[M]/(1+k3[M]/k2)
            rk3M=y(nM,l)*6.5d-34*exp(1335.d0*byta)
            rk2=2.7d-17*exp(2199.d0*byta)
            rr(jj,L)=rr(jj,L)+rk3M/(1.d0+(rk3M/rk2))
          else if (jj==rrbi%XO2N_HO2__CH3OOH_O2) then
            rr(jj,L)=rr(rrbi%XO2_HO2__CH3OOH_O2,L)
     &        *rr(rrbi%XO2N_NO__AlkylNit_M,L)
     &        /rr(rrbi%XO2_NO__NO2_M,L)
          else if (jj==rrbi%PAN_M__C2O3_NO2
#ifdef TRACERS_dCO
     &        .or. jj==rrbi%d17OPAN_M__dC217O3_NO2
     &        .or. jj==rrbi%d18OPAN_M__dC218O3_NO2
     &        .or. jj==rrbi%d13CPAN_M__d13C2O3_NO2
#endif  /* TRACERS_dCO */
     &            ) then
!           PAN+M really PAN
            rr(jj,L)=rr(jj,L)/y(nM,L)
          else if (jj==rrbi%ROR_M__Aldehyde_HO2
#ifdef TRACERS_dCO
     &        .or. jj==rrbi%d17OROR_M__d17Oald_HO2
     &        .or. jj==rrbi%d18OROR_M__d18Oald_HO2
     &        .or. jj==rrbi%d13CROR_M__d13Cald_HO2
#endif  /* TRACERS_dCO */
     &            ) then
!           ROR+M really ROR
            rr(jj,L)=rr(jj,L)/y(nM,L)
          else if (jj==rrbi%ROR_M__HO2_M
#ifdef TRACERS_dCO
     &        .or. jj==rrbi%d17OROR_M__HO2_M
     &        .or. jj==rrbi%d18OROR_M__HO2_M
     &        .or. jj==rrbi%d13CROR_M__HO2_M
#endif  /* TRACERS_dCO */
     &            ) then
!           ROR+M really ROR
! WARNING:
! This reaction, together with the other ROR+M one, might have issues.
! Kostas asked Greg who will ask Drew, and together will fix it, if required.
            rr(jj,L)=rr(jj,L)!/y(nM,L)
          else if (jj==rrbi%HO2_NO__OH_NO2
     &        .or. jj==rrbi%HO2_NO__HNO3_M) then
!           calculate branching ratio here Butkovskaya et al J.Phys.Chem 2007
            waterPPMV=1.d6*y(nH2O,L)/y(nM,L)
            if(ta(L)<298.d0 .and. waterPPMV > 100.)then
              beta=(530.d0*byta + 6.4d-4*pcon*760.d0 - 1.73d0)*1.d-2
            else
              beta=0.d0
            endif
            if (jj==rrbi%HO2_NO__HNO3_M) then
                rr(jj,L)=rr(jj,L)*beta
            else if (jj==rrbi%HO2_NO__OH_NO2) then
                rr(jj,L)=rr(jj,L)*(1.d0-beta)
            end if
          end if
        end do                ! bimolecular rates end
                           
        ! here we USED TO tune rr for N2O+O(1D)-->N2+O2 and N2O+O(1D)-->NO+NO
         
        do jj=1,n_tri         ! trimolecular rates start
          rr(n_bi+n_nst+jj,L)=y(nM,L)*ro(jj)*(300.d0*byta)**sn(jj)
          if(r1(jj) .ne. 0.d0)then 
            dd=rr(n_bi+n_nst+jj,L)/(r1(jj)*(300.d0*byta)**sb(jj))
            pp=0.6d0**(1.d0/(1.d0+(log10(dd))**2.))
            rr(n_bi+n_nst+jj,L)=(rr(n_bi+n_nst+jj,L)/(1.d0+dd))*pp
          end if
        end do                ! trimolecular rates end

        do jj=1,n_nst         ! monomolecular rates start
          ! 0.5 for precision,correct following line:
          rrrr=exp(0.5d0*ea(n_bi+jj)*byta)
          rr(n_bi+jj,L)=rr(nst(jj),L)/(rrrr*pe(n_bi+jj)*rrrr*y(nM,l))
        end do              ! monomolecular rates end

c Calculate rates for heterogeneous reactions (Divided by solid
C in Chem1). sticking coefficients from JPL '02:
c       N2O5 + H2O --> 2HNO3          gamma=0.2, 0.0004 (PSC)
c       ClONO2 + H2O --> HOCl + HNO3  gamma=0.8d-2 (aero), 0.004 (PSC)
c       ClONO2 + HCl --> Cl2 + HNO3   gamma=0.2
c       HOCl + HCl --> Cl2 + H2O      gamma=0.1
c       N2O5 + HCl --> ClNO2 + HNO3   gamma=0.003
C
C Aerosols (14-33 km) & PSCs 14-22 km.
C
c Aerosol profiles and latitudinal distribution of extinction 
c coefficients(in km**-1) are from SAGE II data on GISS web site:
        pfactor=axyp(I,J)*MA(L,I,J)/y(nM,L)
        bypfactor=1.d0/pfactor

        if(pres(L) >= 245.d0 .or. pres(L) <= 5.d0)then 
          do jj=n_bi+n_nst+n_tri+1,n_bi+n_nst+n_tri+n_het
            if (jj == rrhet%N2O5_H2O__HNO3_HNO3) cycle
            rr(jj,L)=1.0d-35
          enddo 
          ! Add rxn of N2O5 on sulfate analogous to what is done in darkness:
          if(rh(L)>0.5)then
            rgammasulf = 1.5d-2
          else
            rgammasulf = 5.2d-2 - 2.79d-4*100.d0*rh(L)
            if(ta(L)>290.) rgammasulf=
     &      max(1.d-3,rgammasulf-log10(ta(L)-290.d0)*5.d-2)
          end if
          if (coupled_chem == 1) then
            ! Convert SO4 from mass (kg) to aerosol surface per grid box:
            ! Here there is a factor of 1d-3  that converts kg/m3 to g/cm3
            ! and 1.76d5 is cm2/g from Dentener and Crutzen, 1993.
            ! So 1.d-3*1.76d5=1.76d2, and that value is for a relative
            ! humidity of 0.75 (1/0.75 = 1.33333d0 below). Reciprocal
            ! layer thickness below is in 1/m units:
            sulfate(i,j,L)=0.0
#ifdef TRACERS_AMP
            sulfate(i,j,L)=
     &      trm(i,j,L,n_M_AKK_SU)+trm(i,j,L,n_M_ACC_SU)+
     &      trm(i,j,L,n_M_DD1_SU)+trm(i,j,L,n_M_DS1_SU)+
     &      trm(i,j,L,n_M_DD2_SU)+trm(i,j,L,n_M_DS2_SU)+
     &      trm(i,j,L,n_M_SSA_SU)+trm(i,j,L,n_M_OCC_SU)+
     &      trm(i,j,L,n_M_BC1_SU)+trm(i,j,L,n_M_BC2_SU)+
     &      trm(i,j,L,n_M_BC3_SU)+trm(i,j,L,n_M_DBC_SU)+
     &      trm(i,j,L,n_M_BOC_SU)+trm(i,j,L,n_M_BCS_SU)+
     &      trm(i,j,L,n_M_MXX_SU)
#elif (defined TRACERS_AEROSOLS_Koch)
            sulfate(i,j,L)=trm(i,j,L,n_SO4)
#elif (defined TRACERS_TOMAS)
            do nb=1,nbins
              sulfate(i,j,L)=sulfate(i,j,L)+trm(i,j,L,n_ASO4(nb))
            end do
#endif
            sulfate(i,j,L)=sulfate(i,j,L)*1.76d2*byaxyp(i,j)*bythick(L)
     &      *max(0.1d0,rh(L)*1.33333d0)
            ! just in case loop changes (b/c sulfate is defined to LM):
            if(L>topLevelOfChemistry)sulfate(i,j,L)=0.d0
          end if

          RVELN2O5=SQRT(ta(L)*RKBYPIM)*100.d0
C         Calculate sulfate sink, and cap it at 20% of N2O5:
c         in troposphere loss is rxn on sulfate, in strat rxn w PSC or sulfate
          wprod_sulf=
     &    dt2*sulfate(I,J,L)*y(nn_N2O5,L)*RGAMMASULF*RVELN2O5*0.25d0

          if(wprod_sulf>0.2d0*y(nn_N2O5,L))wprod_sulf=0.2d0*y(nn_N2O5,L)
          prod_sulf=wprod_sulf*pfactor
          CALL INC_TAJLS(I,J,L,jls_N2O5sulf,
     &                   -1.d0*prod_sulf*vol2mass(n_N2O5))
          rr(rrhet%N2O5_H2O__HNO3_HNO3,L)=wprod_sulf/(dt2*y(nn_N2O5,L))

        else  

          if((pres(l) < 245.d0.and.pres(l) > 150.d0) .or. 
     &    LAXb < 1.or.LAXb > topLevelOfChemistry .or.
     &    LAXt < 1.or.LAXt > topLevelOfChemistry)then 
            rkext(L)=0.d0
          else
            if(pres(l) <= 150..and.pres(l) > 31.60)then
              if(l < LAXb) then
                rkext(l)=rad_to_chem(2,LAXb,i,j)
              else if(l > LAXt) then
                rkext(l)=0.33d0*rkext(l-1)
              else
                rkext(l)=rad_to_chem(2,l,i,j)
              end if
            end if
            if(pres(l) <= 31.6d0.and.pres(l) >= 17.8d0)then
              if(l < LAXb) then
                call stop_model('kext problem 1',255)
              else if(l > LAXt) then
                rkext(l)=0.33d0*rkext(l-1)
              else
                rkext(l)=rad_to_chem(2,l,i,j)
              end if
            end if
            if(pres(l) <= 17.8d0.and.pres(l) >= 10.0d0)then
              if(l < LAXb) then
                call stop_model('kext problem 2',255)
              else if(l > LAXt) then
                rkext(l)=0.33d0*rkext(l-1)
              else
                rkext(l)=rad_to_chem(2,l,i,j)
              end if
            end if
            if(pres(l) <= 10.0d0.and.pres(l) >= 4.6d0)then
              if(l < LAXb) then
                call stop_model('kext problem 3',255)
              else if(l > LAXt) then
                rkext(l)=0.33d0*rkext(l-1)
              else
                rkext(l)=rad_to_chem(2,l,i,j)
              end if
            end if
          end if

          ! Convert rkext from optical depth per layer to optical depth
          ! per cm. (Division by thickness added Jan 2019 concurrent
          ! with removal of km-1 --> cm-1 conversion where rkext is
          ! used below):
          rkext(l)=rkext(l)*bythick(l)*1.d-2

          if(rkext(l) /= 0.)aero(l) = 1

          ! NAT PSC surface conc per unit volume (cm^2/cm^3)
          if(pscX(l))then
            pscEx(l)=2.d-6
          else
            pscEx(l)=0.d0
          end if

c         Reaction rrhet%N2O5_H2O__HNO3_HNO3 on sulfate and PSCs:
          temp=sqrt(boltAvog8byPi*ta(l)/108.d0)
          rr(rrhet%N2O5_H2O__HNO3_HNO3,L)=0.5d0*rkext(l)*temp*0.2d0
          if(pres(l) > 31.6d0) rr(rrhet%N2O5_H2O__HNO3_HNO3,L)=
     &      rr(rrhet%N2O5_H2O__HNO3_HNO3,L)+0.25d0*pscEx(l)*temp*4.d-4

c         Reaction rrhet%ClONO2_H2O__HOCl_HNO3 on sulfate and PSCs:
          temp=sqrt(boltAvog8byPi*ta(l)/97.d0)
          rr(rrhet%ClONO2_H2O__HOCl_HNO3,L)=
     &      0.5d0*rkext(l)*temp*0.8d-2
          if(pres(l) > 31.6d0) rr(rrhet%ClONO2_H2O__HOCl_HNO3,L)=
     &      rr(rrhet%ClONO2_H2O__HOCl_HNO3,L)
     &      +0.25d0*pscEx(l)*temp*4.d-3

          if(pres(l) > 31.6d0) then
            rr(rrhet%ClONO2_HCl__Cl_HNO3,L)=0.25d0*pscEx(l)*temp*0.2d0
            rr(rrhet%HOCl_HCl__Cl_H2O,L)=
     &        sqrt(boltAvog8byPi*ta(l)/52.d0)
            rr(rrhet%HOCl_HCl__Cl_H2O,L)=
     &        0.25d0*pscEx(l)*rr(rrhet%HOCl_HCl__Cl_H2O,L)*0.1d0
            rr(rrhet%N2O5_HCl__Cl_HNO3,L)=
     &        sqrt(boltAvog8byPi*ta(l)/108.d0)
            rr(rrhet%N2O5_HCl__Cl_HNO3,L)=
     &        0.25d0*pscEx(l)*rr(rrhet%N2O5_HCl__Cl_HNO3,L)*0.003d0
          end if

        end if

        if(pres(L) < 245.d0 .and. pres(L) > 5.d0)then
          wprod_sulf=dt2*y(nn_N2O5,L)*rr(rrhet%N2O5_H2O__HNO3_HNO3,L)
          prod_sulf=wprod_sulf*pfactor
          CALL INC_TAJLS(I,J,L,jls_N2O5sulf,
     &                   -1.d0*prod_sulf*vol2mass(n_N2O5))
        end if

      end do                  !  ==> END ALTITUDE LOOP <==

#ifdef TRACERS_AEROSOLS_SOA
      ! Kostas had a comment on this do loop that it should be up to LM whether or
      ! not "strat chem is on or off". But I had to change it here because of the
      ! use of ta(L). I will talk with him, but hopefully he just meant he didn't
      ! want it limited to the troposphere?
      do L=1,topLevelOfChemistry
        kpart(L,whichsoa(n_isopp1a))=
     &       KpCALC(dH_isoprene,kpart_ref(whichsoa(n_isopp1a)),ta(L),
     &              kpart_temp_ref(whichsoa(n_isopp1a)))
        kpart(L,whichsoa(n_isopp2a))=
     &       KpCALC(dH_isoprene,kpart_ref(whichsoa(n_isopp2a)),ta(L),
     &              kpart_temp_ref(whichsoa(n_isopp2a)))
#ifdef TRACERS_TERP
        kpart(L,whichsoa(n_apinp1a))=
     &       KpCALC(dH_apinene,kpart_ref(whichsoa(n_apinp1a)),ta(L),
     &              kpart_temp_ref(whichsoa(n_apinp1a)))
        kpart(L,whichsoa(n_apinp2a))=
     &       KpCALC(dH_apinene,kpart_ref(whichsoa(n_apinp2a)),ta(L),
     &              kpart_temp_ref(whichsoa(n_apinp2a)))
#endif  /* TRACERS_TERP */
      end do
#endif  /* TRACERS_AEROSOLS_SOA */

      deallocate( PSCEX )
      deallocate( rkext )
 
      RETURN
      END SUBROUTINE Crates


