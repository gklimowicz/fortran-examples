#include "rundeck_opts.h"
!------------------------------------------------------------------------------
module TomasTracersMetadata_mod
!------------------------------------------------------------------------------
!@sum  TomasTracersMetadata_mod encapsulates the TRACERS_TOMAS metadata
!@auth NCCS ASTG
  use sharedTracersMetadata_mod, only: DMS_setspec, SO2_setspec, &
    NH3_setspec, NH4_setspec, H2O2_s_setSpec
  USE CONSTANT, only: pi
  use Dictionary_mod, only: sync_param
  use OldTracer_mod, only: set_emisPerFireByVegType
  use OldTracer_mod, only : set_tr_mm
  use OldTracer_mod, only : set_ntm_power
  use OldTracer_mod, only : set_trpdens
  use OldTracer_mod, only : set_trradius
  use OldTracer_mod, only : set_fq_aer
  use OldTracer_mod, only : set_tr_wd_type
  use OldTracer_mod, only : oldAddTracer
  use OldTracer_mod, only: set_HSTAR
  use OldTracer_mod, only: set_needtrs
  use OldTracer_mod, only: set_trpdens
  use OldTracer_mod, only: set_trradius
  use OldTracer_mod, only: set_tr_wd_TYPE
  use OldTracer_mod, only: set_fq_aer
  use OldTracer_mod, only: set_has_chemistry
  use OldTracer_mod, only: nGAS, nPart
  use TRACER_COM, only: xk, nbins, coupled_chem
  use TRACER_COM, only: n_NH4, n_H2SO4
  use TRACER_COM, only: n_ASO4, n_ANACL, n_AECIL, n_AECOB, &
    n_AOCIL, n_ADUST, n_ANUM, n_AOCOB, n_AH2O, n_SOAgas
  use TRACER_COM, only: set_ntsurfsrc
  use TOMAS_AEROSOL, only : binact10, binact02, fraction10, fraction02
  use RunTimeControls_mod, only: tracers_aerosols_soa
  use RunTimeControls_mod, only: tracers_special_shindell
  use RunTimeControls_mod, only: tracers_drydep
  use Tracer_mod, only: Tracer

  implicit none
  private

  public TOMAS_initMetadata

  integer :: n ! class scoped temporary tracer index

!------------------------------------------------------------------------------
contains
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  subroutine TOMAS_initMetadata(pTracer)
!------------------------------------------------------------------------------
    class (Tracer), pointer :: pTracer
    real*8 :: TOMAS_dens,TOMAS_radius

    CALL initbounds()
#ifdef TOMAS_12_10NM
      call readbinact ("binact10_12.dat",binact10) 
      call readbinact ("binact02_12.dat",binact02) 
      call readfraction("fraction10_12.dat",fraction10) 
      call readfraction("fraction02_12.dat",fraction02) 
#endif
#ifdef TOMAS_12_3NM
      call readbinact ("binact10_12_3nm.dat",binact10) 
      call readbinact ("binact02_12_3nm.dat",binact02) 
      call readfraction("fraction10_12_3nm.dat",fraction10) 
      call readfraction("fraction02_12_3nm.dat",fraction02) 
#endif
    call readmielut           ! aerosol radiation lookup table

    ! For aerosol tracers in TOMAS model, 
    ! fq_aer is determined based on kohler theory. 
    call  TOMAS_H2SO4_setSpec('H2SO4')
    call  DMS_setSpec('DMS')  ! note duplicate with Koch
    call  SO2_setSpec('SO2')
#ifndef TRACERS_AEROSOLS_SOA
    if (.not. tracers_aerosols_soa) &
      call  TOMAS_SOAgas_setSpec('SOAgas')
#endif
    if (.not. tracers_special_shindell .or. coupled_chem.eq.0) then
      call  H2O2_s_setSpec('H2O2_s') ! duplicate with Koch
    endif
    call  NH3_setSpec('NH3')  ! duplicate with nitrate
    call  NH4_setSpec('NH4')  ! duplicate with nitrate

    ! Koch SO4: 1700 kg/m3
    ! Koch BC : 1300 
    ! Koch OC : 1500
    ! Koch DUST : 2500 for  clay and 2650 for silt
    ! Koch SS : 2200 kg/m3

    ! so4 : 1780 kg/m3
    ! ss:  2165 kg/m3
    !bc: 1800 kg/m3 or 2200 kg/m3
    !oc:1400 kg/m3
    !ddust : 2650 kg/m3 
    n_ASO4(:)  = TOMAS_setSpec(TOMAS_ASO4_setSpec,  'ASO4', nbins)
    n_ANACL(:) = TOMAS_setSpec(TOMAS_ANACL_setSpec, 'ANACL',nbins)
    n_AECOB(:) = TOMAS_setSpec(TOMAS_AECOB_setSpec, 'AECOB',nbins)
    n_AECIL(:) = TOMAS_setSpec(TOMAS_AECIL_setSpec, 'AECIL',nbins)
    n_AOCOB(:) = TOMAS_setSpec(TOMAS_AOCOB_setSpec, 'AOCOB',nbins)
    n_AOCIL(:) = TOMAS_setSpec(TOMAS_AOCIL_setSpec, 'AOCIL',nbins)
    n_ADUST(:) = TOMAS_setSpec(TOMAS_ADUST_setSpec, 'ADUST',nbins)
    n_ANUM(:) = TOMAS_setSpec(TOMAS_ANUM_setSpec,   'ANUM', nbins)
    n_AH2O(:)  = TOMAS_setSpec(TOMAS_AH2O_setSpec,  'AH2O', nbins)

!------------------------------------------------------------------------------
  contains
!------------------------------------------------------------------------------

    function TOMAS_setSpec(func, name, nbins) result (indices)
      interface
        integer function func(name, bin)
          character(len=*), intent(in) :: name
          integer, intent(in) :: bin
        end function func
      end interface
      character(len=*), intent(in) :: name
      integer, intent(in) :: nbins
      integer :: indices(nbins)

      integer :: bin
      character(len=len_trim(name) + 4) :: fullName

      do bin = 1, nbins
        if (len_trim(name) == 5) then
          write(fullName,'(a,"_",i2.2)') trim(name), bin
        else
          write(fullName,'(a,"__",i2.2)') trim(name), bin
        end if
        indices(bin) = func(fullName, bin)
      end do

    end function TOMAS_setSpec

    integer function TOMAS_ANUM_setSpec(name, bin) result(n_ANUM)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin

      n = oldAddTracer(name)
      n_ANUM = n  
      TOMAS_dens = 1.5d3
      TOMAS_radius = (sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)**(1./3.)  
      if(bin.le.5) call set_ntm_power(n, 10)
      if(bin.gt.5) call set_ntm_power(n, 8) 

      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, TOMAS_dens)
      call set_trradius(n, TOMAS_radius)
      call set_fq_aer(n, 1.d0)  !not used in wet deposition
      call set_tr_wd_type(n, npart)
    end function TOMAS_ANUM_setSpec

    integer function TOMAS_ASO4_setSpec(name, bin) result(n_ASO4)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin

      n = oldAddTracer(name)
      n_ASO4 = n 
      TOMAS_dens = 1.78d3
      TOMAS_radius = (sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)**(1./3.) 
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 96.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)
    end function TOMAS_ASO4_setSpec

    integer function TOMAS_ANACL_setSpec(name, bin) result(n_ANACL)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin

      n = oldAddTracer(name)
      n_ANACL = n         
      TOMAS_dens = 2.165d3
      TOMAS_radius = (sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)**(1./3.) 
      if(bin.le.10) call set_ntm_power(n, -10)
      if(bin.gt.10) call set_ntm_power(n, -8)
      call set_tr_mm(n, 75.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)
    end function TOMAS_ANACL_setSpec

    integer function TOMAS_AECOB_setSpec(name, bin) result(n_AECOB)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin

      n = oldAddTracer(name)
      n_AECOB = n          
      TOMAS_dens = 1.8d3
      TOMAS_radius = (sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)**(1./3.) 
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)
      call set_has_chemistry(n, .true.)
    end function TOMAS_AECOB_setSpec

    integer function TOMAS_AECIL_setSpec(name, bin) result(n_AECIL)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin

      n = oldAddTracer(name)
      n_AECIL = n        
      TOMAS_dens = 1.8d3
      TOMAS_radius = (sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)**(1./3.) 
      call set_ntm_power(n, -12)
      call set_tr_mm(n, 12.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)
      call set_has_chemistry(n, .true.)
    end function TOMAS_AECIL_setSpec

    integer function TOMAS_AOCOB_setSpec(name, bin) result(n_AOCOB)
      use OldTracer_mod, only: om2oc, set_om2oc
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin
      real*8 :: tmp

      n = oldAddTracer(name)
      n_AOCOB = n 
      TOMAS_dens = 1.4d3
      TOMAS_radius = (sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)**(1./3.) 
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 200.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)        

      if (bin==1) then
        call set_om2oc(n, 1.4d0)
        tmp = om2oc(n_AOCOB)
        call sync_param("OCB_om2oc",tmp)
        call set_om2oc(n_AOCOB, tmp)
      end if
      call set_has_chemistry(n, .true.)
    end function TOMAS_AOCOB_setSpec

    integer function TOMAS_AOCIL_setSpec(name, bin) result(n_AOCIL)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin

      n = oldAddTracer(name)
      n_AOCIL = n  
      TOMAS_dens = 1.4d3
      TOMAS_radius = (sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)**(1./3.) 
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 200.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)        
      call set_has_chemistry(n, .true.)
    end function TOMAS_AOCIL_setSpec

    integer function TOMAS_ADUST_setSpec(name, bin) result(n_ADUST)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin

      n = oldAddTracer(name)
      n_ADUST = n  
      if(bin.le.10) TOMAS_dens= 2.5d3 !clay 
      if(bin.gt.10) TOMAS_dens= 2.65d3 !Silt
      TOMAS_radius = (sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)**(1./3.) 

      if(bin.le.9) call set_ntm_power(n, -11)
      if(bin.gt.9) call set_ntm_power(n, -9)
      call set_tr_mm(n, 1.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)  
    end function TOMAS_ADUST_setSpec

    integer function TOMAS_AH2O_setSpec(name, bin) result(n_AH2O)
      character(len=*), intent(in) :: name
      integer, intent(in) :: bin

      n = oldAddTracer(name)
      n_AH2O = n         
      TOMAS_dens = 1.d3
      TOMAS_radius = (sqrt(xk(bin)*xk(bin+1))/TOMAS_dens/pi/4.*3.)**(1./3.) 
      call set_ntm_power(n, -8)

      call set_tr_mm(n, 18.d+0)
      call set_trpdens(n, TOMAS_dens) !kg/m3 this is sulfate value
      call set_trradius(n, TOMAS_radius) !m
      call set_fq_aer(n, 1.d0   ) !not used in wet deposition
      call set_tr_wd_type(n, npart)  
    end function TOMAS_AH2O_setSpec

    subroutine NH4_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_NH4 = n
      call set_ntm_power(n, -10)
      call set_tr_mm(n, 18.d0)
      call set_trpdens(n, 1.7d3)
      call set_trradius(n, 3.d-7)
      call set_fq_aer(n, 1.0d0   ) !fraction of aerosol that dissolves
      call set_tr_wd_type(n, npart)
      call set_has_chemistry(n, .true.)
    end subroutine NH4_setSpec

    subroutine TOMAS_H2SO4_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_H2SO4 = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 98.d0)
      call set_trpdens(n, 1.78d0)
      call set_fq_aer(n, 1.d0)
      call set_tr_wd_type(n, nGas)
      call set_has_chemistry(n, .true.)
    end subroutine TOMAS_H2SO4_setSpec

    subroutine TOMAS_SOAgas_setSpec(name)
      character(len=*), intent(in) :: name

      n = oldAddTracer(name)
      n_SOAgas = n
      call set_ntm_power(n, -11)
      call set_tr_mm(n, 120.10d0) ! i.e. 10 carbons
      if (tracers_drydep) call set_HSTAR(n,  0.D0)  !no dry dep
      call set_has_chemistry(n, .true.)
    end subroutine TOMAS_SOAgas_setSpec

  end subroutine TOMAS_initMetadata

end module TomasTracersMetadata_mod

