#include "rundeck_opts.h"
      SUBROUTINE initTracerGriddedData(is_coldstart)
!@sum init_tracer initializes trace gas attributes
!@calls sync_param, SET_TCON, RDLAND, RDDRYCF
      USE DOMAIN_DECOMP_ATM, only:GRID,getDomainBounds,AM_I_ROOT,
     &     write_parallel,readt8_parallel
      USE RESOLUTION, only : jm,lm
      USE ATM_COM, only: pmidl00
      USE GEOM, only: axyp,byaxyp
      USE ATM_COM, only: MA  ! Air mass of each box (kg/m^2)
      use OldTracer_mod, only: trname
      USE TRACER_COM, only: ntm, tracers, syncProperty
#ifdef TRACERS_ON
      USE TRDIAG_COM
#endif
      USE Dictionary_mod
#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP)  || (defined TRACERS_TOMAS)
      use trdust_drv, only : init_soildust
#endif
#ifdef TRACERS_SPECIAL_Shindell
      USE TRCHEM_Shindell_COM,only:LCOalt,PCOalt,
     &     CH4altINT,CH4altINX,LCH4alt,PCH4alt,
     &     CH4altX,CH4altT,ch4_init_sh,ch4_init_nh,scale_ch4_IC_file,
     &     OxICIN,OxIC,OxICINL,OxICL,
     &     fix_CH4_chemistry,
     &     CH4ICIN,CH4ICX,CH4ICINL,CH4ICL,use_rad_ch4,
     &     COICIN,COIC,COICINL,COICL,Lmax_rad_O3,Lmax_rad_CH4
     &     ,BrOxaltIN,ClOxaltIN,ClONO2altIN,HClaltIN,BrOxalt,
     &     ClOxalt,ClONO2alt,HClalt,N2OICIN,N2OICX,N2OICINL,N2OICL,
     &     CFCICIN,CFCIC,CFCICINL,CFCICL,
     &     use_rad_n2o,use_rad_cfc,cfc_rad95
#endif /* TRACERS_SPECIAL_Shindell */
#ifdef TRACERS_AEROSOLS_SOA
      USE TRACERS_SOA, only: soa_init
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AEROSOLS_VBS
      USE TRACERS_VBS, only: vbs_init
#endif  /* TRACERS_AEROSOLS_VBS */
#if (defined TRACERS_AMP)
      USE AERO_COAG, only : SETUP_KIJ
      USE AERO_SETUP, only: SETUP_CONFIG,SETUP_SPECIES_MAPS,SETUP_DP0,
     &                      SETUP_AERO_MASS_MAP,SETUP_COAG_TENSORS,
     &                      SETUP_EMIS,SETUP_KCI
      USE AERO_NPF, only: SETUP_NPFMASS
      USE AERO_DIAM, only: SETUP_DIAM
#endif
      USE FILEMANAGER, only: openunit,closeunit,nameunit

      implicit none
      logical, intent(in) :: is_coldstart
c
      integer :: l,k,n,kr,m,ns
#ifdef TRACERS_SPECIAL_O18
      real*8 fracls
#endif
#ifdef TRACERS_TOMAS
      integer :: bin
      real*8 :: TOMAS_dens,TOMAS_radius
#endif
#ifdef TRACERS_SPECIAL_Shindell
!@var iu_data unit number
!@var title header read in from file
      integer iu_data,i,j,nq
      character*80 title
      character(len=300) :: out_line
      real*8, dimension(6) :: temp_ghg
      integer :: temp_year
#endif /* TRACERS_SPECIAL_Shindell */

! temp storage for new tracer interfaces
      integer :: values(ntm)
      integer :: val

      INTEGER J_0, J_1, I_0, I_1

C****
C**** Extract useful local domain parameters from "grid"
C****
      call getDomainBounds(grid, J_STRT=J_0,       J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

      do n=1,ntm
#ifdef TRACERS_ON
        select case (trname(n))

        case ('N2O')
#ifdef TRACERS_SPECIAL_Shindell
c***          print*,'HERE!!!!!!!'
          call openunit('N2O_IC',iu_data,.true.,.true.)
          CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),N2OICIN,0)
          call closeunit(iu_data)
          do j=J_0,J_1  ; do i=I_0,I_1
           N2OICINL(:)=N2OICIN(i,j,:) ! now in PPPM
           CALL LOGPINT(LCOalt,PCOalt,N2OICINL,LM,PMIDL00,N2OICL,.true.)
           N2OICX(I,J,:) = N2OICL(:)*MA(:,i,j)*axyp(i,j)
          end do     ; end do
#endif

      case ('CH4')
#ifdef TRACERS_SPECIAL_Shindell
C**** determine initial CH4 distribution if set from rundeck
C**** This is only effective with a complete restart.
          call sync_param("ch4_init_sh",ch4_init_sh)
          call sync_param("ch4_init_nh",ch4_init_nh)
          call sync_param("fix_CH4_chemistry",fix_CH4_chemistry)
          call sync_param("scale_ch4_IC_file",scale_ch4_IC_file)
C         Interpolate CH4 altitude-dependence to model resolution:
          CALL LOGPINT(LCH4alt,PCH4alt,CH4altINT,LM,PMIDL00,CH4altT,
     &         .true.)
          CALL LOGPINT(LCH4alt,PCH4alt,CH4altINX,LM,PMIDL00,CH4altX,
     &         .true.)
          if(fix_CH4_chemistry.eq.-1)then
            call openunit('CH4_IC',iu_data,.true.,.true.)
            CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),
     &           CH4ICIN,0)
            call closeunit(iu_data)
            do j=J_0,J_1  ; do i=I_0,I_1
             CH4ICINL(:)=CH4ICIN(I,J,:)! now in PPPM
             CALL LOGPINT(LCOalt,PCOalt,CH4ICINL,LM,PMIDL00,CH4ICL,
     &            .true.)
             CH4ICX(I,J,:) = CH4ICL(:)*scale_ch4_IC_file*MA(:,i,j)*
     *            axyp(i,j)
            end do     ; end do
          end if
#endif /* TRACERS_SPECIAL_Shindell */

      case ('Ox')
#ifdef TRACERS_SPECIAL_Shindell
          call openunit('Ox_IC',iu_data,.true.,.true.)
          CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),OxICIN,0)
          call closeunit(iu_data)
          do j=J_0,J_1  ; do i=I_0,I_1
           OxICINL(:)=OxICIN(I,J,:)! now in PPPM
           CALL LOGPINT(LCOalt,PCOalt,OxICINL,LM,PMIDL00,OxICL,.true.)
           OxIC(I,J,:) = OxICL(:)*MA(:,i,j)*axyp(i,j)
          end do     ; end do
#endif /* TRACERS_SPECIAL_Shindell */

      case ('CFC')
#ifdef TRACERS_SPECIAL_Shindell
          if(AM_I_ROOT( ))then
C          check on GHG files 1995 value for CFCs:
           call openunit('GHG',iu_data,.false.,.true.)
           do i=1,5; read(iu_data,'(a80)') title; enddo
           temp_year=0
           do while(temp_year <= 1995)
             read(iu_data,*,end=101) temp_year,(temp_ghg(j),j=1,6)
             if(temp_year==1995)then
               temp_ghg(1)=cfc_rad95*0.95d0
               temp_ghg(2)=cfc_rad95*1.05d0
               temp_ghg(3)=(temp_ghg(4)+temp_ghg(5))*1.d-9
               if(temp_ghg(3) < temp_ghg(1) .or.
     &         temp_ghg(3) > temp_ghg(2))then
                 call stop_model('please check on cfc_rad95 2',255)
               endif
             endif
           enddo
 101       continue
           if(temp_year<1995)
     &     call stop_model('please check on cfc_rad95 1',255)
           call closeunit(iu_data)
          endif
C          read the CFC initial conditions:
          call openunit('CFC_IC',iu_data,.true.,.true.)
          CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),CFCICIN,0)
          call closeunit(iu_data)
          do j=J_0,J_1  ; do i=I_0,I_1
           CFCICINL(:)=CFCICIN(I,J,:)! now in PPPM
           CALL LOGPINT(LCOalt,PCOalt,CFCICINL,LM,PMIDL00,CFCICL,.true.)
           CFCIC(I,J,:) = CFCICL(:)*MA(:,i,j)*axyp(i,j)
          end do     ; end do
#endif /* TRACERS_SPECIAL_Shindell */

      case ('CO'
#ifdef TRACERS_dCO
     *     ,'dC17O','dC18O','d13CO'
#endif  /* TRACERS_dCO */
     *     )
#ifdef TRACERS_SPECIAL_Shindell
          call openunit('CO_IC',iu_data,.true.,.true.)
          CALL READT8_PARALLEL(grid,iu_data,NAMEUNIT(iu_data),COICIN,0)
          call closeunit(iu_data)
          do j=J_0,J_1  ; do i=I_0,I_1
           COICINL(:)=COICIN(I,J,:)! now in PPPM
           CALL LOGPINT(LCOalt,PCOalt,COICINL,LM,PMIDL00,COICL,.true.)
           COIC(I,J,:) = COICL(:)*MA(:,i,j)*axyp(i,j)
          end do     ; end do
#endif /* TRACERS_SPECIAL_Shindell */

        end select
#endif
      end do

#ifdef TRACERS_ON
#ifdef TRACERS_AEROSOLS_SOA
      call soa_init
#endif  /* TRACERS_AEROSOLS_SOA */
#ifdef TRACERS_AEROSOLS_VBS
      call vbs_init(ntm)
#endif  /* TRACERS_AEROSOLS_VBS */
C**** Get to_volume_MixRat from rundecks if it exists
      call syncProperty(tracers, "to_volume_MixRat",
     &     set_to_volume_MixRat,to_volume_MixRat())
C**** Get to_conc from rundecks if it exists
      call syncProperty(tracers,"to_conc", set_to_conc, to_conc())

#if (defined TRACERS_DUST) || (defined TRACERS_MINERALS) ||\
    (defined TRACERS_AMP) || (defined TRACERS_TOMAS)
c**** soil dust aerosol initializations
      call init_soildust
#endif

C**** Miscellaneous initialisations

#ifdef TRACERS_DRYDEP
C Read landuse parameters and coefficients for tracer dry deposition:
      CALL RDLAND
      CALL RDDRYCF
#endif
#ifdef BIOGENIC_EMISSIONS
      CALL RDISOPCF
      CALL RDISOBASE
#endif
#ifdef TRACERS_SPECIAL_Shindell
      call cheminit ! **** Initialize the chemistry ****
#endif
#ifdef TRACERS_COSMO
      call init_cosmo
#endif
#endif /* TRACERS_ON */

#ifdef TRACERS_AMP
      CALL SETUP_CONFIG
      CALL SETUP_SPECIES_MAPS
      CALL SETUP_DP0
      CALL SETUP_AERO_MASS_MAP
      CALL SETUP_COAG_TENSORS
      CALL SETUP_DP0
      CALL SETUP_KIJ
      CALL SETUP_EMIS
      CALL SETUP_KCI
      CALL SETUP_NPFMASS
      if(is_coldstart) ! do not overwrite diam during warm starts
     &     CALL SETUP_DIAM
      CALL SETUP_RAD
#endif

      call init_src_dist

      return
      end subroutine initTracerGriddedData
