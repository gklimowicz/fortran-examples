!
! This module provides a mechanism to refer to distributed
! quantities from modelE via alternate names such that the
! gathered "global" quantities can be referred to by their
! usual name.
!
! By isolating the F90 "rename" mechanism is this otherwise
! vacuous module, I hope to minimize the complexity for
! the hybrid_mpi_omp_coupler module.
!
! Tom Clune <Thomas.L.Clune@nasa.gov>
!

module hybrid_mpi_omp_renamer


  ! Variables/entities used from modelE
  USE FLUXES, only: e0_loc      => e0
  USE FLUXES, only: prec_loc    => prec
  USE FLUXES, only: eprec_loc   => eprec
  USE FLUXES, only: evapor_loc  => evapor
  USE FLUXES, only: flowo_loc   => flowo
  USE FLUXES, only: eflowo_loc  => eflowo
  USE FLUXES, only: dmua_loc    => dmua
  USE FLUXES, only: dmva_loc    => dmva

  USE FLUXES, only: erunosi_loc => erunosi
  USE FLUXES, only: runosi_loc  => runosi
  USE FLUXES, only: srunosi_loc => srunosi
  USE FLUXES, only: runpsi_loc  => runpsi
  USE FLUXES, only: srunpsi_loc => srunpsi
  USE FLUXES, only: dmui_loc    => dmui
  USE FLUXES, only: dmvi_loc    => dmvi
  USE FLUXES, only: dmsi_loc    => dmsi
  USE FLUXES, only: dhsi_loc    => dhsi
  USE FLUXES, only: dssi_loc    => dssi

  USE FLUXES, only: gtemp_loc   => gtemp
  USE FLUXES, only: sss_loc     => sss
  USE FLUXES, only: mlhc_loc    => mlhc
  USE FLUXES, only: ogeoza_loc  => ogeoza
  USE FLUXES, only: uosurf_loc  => uosurf
  USE FLUXES, only: vosurf_loc  => vosurf
  USE FLUXES, only: MELTI_loc   => MELTI 
  USE FLUXES, only: EMELTI_loc  => EMELTI
  USE FLUXES, only: SMELTI_loc  => SMELTI

  USE FLUXES, only: gmelt_loc   => gmelt
  USE FLUXES, only: egmelt_loc  => egmelt
  USE FLUXES, only: solar_loc   => solar

  USE SEAICE_COM, only: rsi_loc   => rsi
  USE SEAICE_COM, only: msi_loc   => msi
  USE SEAICE, only:     fsss   ! scalar
  USE SEAICE, only:     tfrez  ! procedure

  USE MODEL_COM, only: focean_loc => focean

  USE MODEL_COM, only: dtsrc,itime,iyear1,nday,jdendofm,aMON,modelEclock

  USE GEOM, only : dxyp
  USE CONSTANT,  only: lhm,shi,shw

  implicit none
  private


  PUBLIC :: e0_loc      
  PUBLIC :: prec_loc    
  PUBLIC :: eprec_loc   
  PUBLIC :: evapor_loc  
  PUBLIC :: flowo_loc   
  PUBLIC :: eflowo_loc  
  PUBLIC :: dmua_loc    
  PUBLIC :: dmva_loc    

  PUBLIC :: erunosi_loc 
  PUBLIC :: runosi_loc  
  PUBLIC :: srunosi_loc 
  PUBLIC :: runpsi_loc  
  PUBLIC :: srunpsi_loc 
  PUBLIC :: dmui_loc    
  PUBLIC :: dmvi_loc    
  PUBLIC :: dmsi_loc    
  PUBLIC :: dhsi_loc    
  PUBLIC :: dssi_loc    

  PUBLIC :: gtemp_loc   
  PUBLIC :: sss_loc     
  PUBLIC :: mlhc_loc    
  PUBLIC :: ogeoza_loc  
  PUBLIC :: uosurf_loc  
  PUBLIC :: vosurf_loc  
  PUBLIC :: MELTI_loc   
  PUBLIC :: EMELTI_loc  
  PUBLIC :: SMELTI_loc  

  PUBLIC :: gmelt_loc   
  PUBLIC :: egmelt_loc  
  PUBLIC :: solar_loc   

  PUBLIC :: rsi_loc   
  PUBLIC :: msi_loc   

  PUBLIC :: fsss   ! scalar parameter
  PUBLIC :: tfrez  ! routine

  PUBLIC :: focean_loc 

  PUBLIC :: dtsrc,itime,iyear1,nday,jdendofm,jyear,jmon,jday,jdate,jhour,aMON

  PUBLIC :: dxyp
  PUBLIC :: lhm,shi,shw

end module hybrid_mpi_omp_renamer
