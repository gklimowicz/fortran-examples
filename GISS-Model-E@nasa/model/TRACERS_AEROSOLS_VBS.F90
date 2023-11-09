!===============================================================================
module TRACERS_VBS
!@sum Module for calculating the aerosol concentrations based on the volatility
!@+ basis set (vbs; Robinson et al., 2007). This file is generic enough so that
!@+ it can be run as a box model. It must have NO dependencies on the host model
!@auth Kostas Tsigaridis (ktsigaridis@giss.nasa.gov)
!===============================================================================
!===============================================================================
implicit none
private
public :: vbs_bins, vbs_tracers, vbs_conditions, vbs_init, vbs_calc, &
          ivbs_m2, ivbs_m1, ivbs_m0, ivbs_p1, ivbs_p2, ivbs_p3, ivbs_p4, &
          ivbs_p5, ivbs_p6
!-------------------------------------------------------------------------------
!@param ivbs_m2,ivbs_m1,ivbs_m0,ivbs_p1,ivbs_p2,ivbs_p3,ivbs_p4,ivbs_p5,ivbs_p6
!@+ indices of vbs tracers
integer, parameter :: ivbs_m2=1, ivbs_m1=2, ivbs_m0=3, ivbs_p1=4, &
                      ivbs_p2=5, ivbs_p3=6, ivbs_p4=7, ivbs_p5=8, ivbs_p6=9
!@param vbs_bins Number of bins in the volatility basis set
integer, parameter :: vbs_bins=9
!-------------------------------------------------------------------------------
type vbs_tracers
!@var nbins Number of bins in the volatility basis set, equal to vbs_bins
integer :: nbins=vbs_bins
!@var gas gas-phase concentration of tracer (ug m-3)
!@var aer aerosol-phase concentration of tracer (ug m-3)
real*8, dimension(vbs_bins) :: gas, aer
!@var igas index of gas-phase vbs tracers of the host model
!@var iaer index of aerosol-phase vbs tracers of the host model
integer, dimension(vbs_bins) :: igas, iaer
!@var igasinv inverse igas
!@var iaerinv inverse iaer
integer, allocatable, dimension(:) :: igasinv, iaerinv
!@var chem_prod amount of chemically produced tracer in the gas-phase
!@var chem_loss amount of chemically lost (aged) tracer in the gas-phase
!@var partition amount of tracer that partitioned TO the aerosol phase (FROM, if -)
real*8, dimension(vbs_bins) :: chem_prod, chem_loss, partition
end type vbs_tracers
!-------------------------------------------------------------------------------
type vbs_properties
!@var sat saturation concentration (C*) of the vbs bins (ug m-3). The lower
!@+ the saturation concentration, the less volatile the bin.
!@var dH enthalpy of vaporization of the vbs bins (KJ mol-1). The higher
!@+ the enthalpy of vaporization, the strongest the dependence of saturation
!@+ concentration to temperature.
!@var kOH reaction rate of the aging in the gas-phase. The value in the
!@+ bin i represents the conversion of mass from bin i to bin i-1 (cm3
!@+ molecules-1 s-1)
!@var Tref reference temperature for which sat corresponds to
real*8, dimension(vbs_bins) :: sat, dH, kOH, Tref
end type vbs_properties
!-------------------------------------------------------------------------------
type vbs_conditions
!@var dt timestep (seconds)
!@var OH OH radical concentration (molecules cm-3)
!@var temp temperature
!@var nvoa concentration of non-volatile organic aerosol (ug m-3)
real*8 :: dt, OH, temp, nvoa
end type vbs_conditions
!-------------------------------------------------------------------------------
!@var vbs_tr vbs tracers concentrations (ug m-3), indices and diagnostics
!@+ (ug m-3 timestep-1)
type(vbs_tracers), public :: vbs_tr
!@var vbs_prop vbs bin properties
type(vbs_properties) :: vbs_prop
!@var vbs_cond current contitions of the atmosphere
type(vbs_conditions), public :: vbs_cond
!-------------------------------------------------------------------------------
contains
!===============================================================================
!===============================================================================
subroutine vbs_init(ntm_host)
!===============================================================================
implicit none
!-------------------------------------------------------------------------------
!@var ntm_host number of tracers in the host model
integer, intent(in) :: ntm_host
integer :: bin
!-------------------------------------------------------------------------------
allocate(vbs_tr%igasinv(ntm_host))
allocate(vbs_tr%iaerinv(ntm_host))
vbs_tr%igasinv=0
vbs_tr%iaerinv=0
do bin=1,vbs_bins
  vbs_tr%igasinv(vbs_tr%igas(bin))=bin
  vbs_tr%iaerinv(vbs_tr%iaer(bin))=bin
  vbs_prop%sat(bin)=10.d0**(dble(bin-3)) ! from -2 to +6 with step of 1
enddo
! generic dH/Tref
!vbs_prop%dH(:)=42.d0
!vbs_prop%Tref(:)=298.d0
! Epstein et al., 2010
vbs_prop%dH(:)=-11.d0*dlog10(vbs_prop%sat(:)) + 131.d0
vbs_prop%Tref(:)=300.d0
vbs_prop%kOH(1)=0.d0
vbs_prop%kOH(2:)=1.d-11
!===============================================================================
end subroutine vbs_init
!===============================================================================
!===============================================================================
subroutine vbs_calc(tr,cond)
!===============================================================================
implicit none
!-------------------------------------------------------------------------------
!@var gas input gas-phase vbs concentration (ug m-3)
!@var aer input aerosol-phase vbs concentration (ug m-3)
type(vbs_tracers) :: tr
type(vbs_conditions), intent(in) :: cond
!-------------------------------------------------------------------------------
vbs_tr%gas=tr%gas
vbs_tr%aer=tr%aer
vbs_cond=cond
call vbs_age
call vbs_partition
!===============================================================================
end subroutine vbs_calc
!===============================================================================
!===============================================================================
subroutine vbs_age
!===============================================================================
implicit none
!-------------------------------------------------------------------------------
!@var loss_fraction fraction of the vbs gas-phase tracer that ages
real*8, dimension(vbs_bins) :: loss_fraction
!-------------------------------------------------------------------------------
loss_fraction=1.d0-exp(-vbs_prop%kOH*vbs_cond%OH*vbs_cond%dt)
vbs_tr%chem_loss=-vbs_tr%gas*loss_fraction

vbs_tr%chem_prod(1:vbs_bins-1)=-vbs_tr%chem_loss(2:vbs_bins)
vbs_tr%chem_prod(vbs_bins)=0.d0 ! no chemical production for the most volatile bin

vbs_tr%gas=vbs_tr%gas+vbs_tr%chem_prod+vbs_tr%chem_loss
!===============================================================================
end subroutine vbs_age
!===============================================================================
!===============================================================================
subroutine vbs_partition
!===============================================================================
implicit none
!-------------------------------------------------------------------------------
integer :: iter ! iteration count
real*8, dimension(vbs_bins) :: ksi,vbs_tot ! see Donahue et al., 2006
!@var Ccurr value of C* at current temperature (ug m-3)
real*8, dimension(vbs_bins) :: Ccurr
!@param maxit maximum number of iterations allowed
integer, parameter :: maxit=10000 !!!!!!!!! this is EXTREMELY HIGH!!!!!!!!!
!@var Mo total mass of condensed OA at equilibrium (ug m-3)
!@var Mo_guess total mass of condensed OA while searching for Mo (ug m-3)
real*8 :: Mo, Mo_guess
!-------------------------------------------------------------------------------
vbs_tot=vbs_tr%gas+vbs_tr%aer
if (sum(vbs_tot)+vbs_cond%nvoa == 0.d0) return
call clausius_clapeyron(Ccurr)
Mo=vbs_cond%nvoa
Mo_guess=max(tiny(Mo),Mo)
iter=0
do ! loop indefinitely until a solution is found, or iter > maxit
  iter=iter+1
  ksi=1.d0/(1.d0+Ccurr/Mo_guess)
  Mo=vbs_cond%nvoa+sum(vbs_tot*ksi)
  if (Mo == Mo_guess) then
!    print *,'Solution found after ',iter,' iterations'
!    print *,'ksi=',ksi
!    print *,'temp=',vbs_cond%temp
!    print *,'Ccurr=',Ccurr
    exit
  else
    if (iter > maxit) then
      print *,'iter>maxit, breaking operation. Mo= ',Mo,' Mo_guess= ',Mo_guess
      exit
    endif
  endif
  Mo_guess=Mo
enddo
vbs_tr%partition=vbs_tot*ksi-vbs_tr%aer ! positive means TO the aerosol phase
vbs_tr%gas=vbs_tot*(1.d0-ksi)
vbs_tr%aer=vbs_tot*ksi
!===============================================================================
end subroutine vbs_partition
!===============================================================================
!===============================================================================
subroutine clausius_clapeyron(C)
!===============================================================================
use constant, only: bygasc
implicit none
!-------------------------------------------------------------------------------
real*8, dimension(vbs_bins) :: C
!-------------------------------------------------------------------------------
C=vbs_prop%sat*vbs_prop%Tref/vbs_cond%temp*& ! 1.d3 converts KJ to J
  exp(1.d3*vbs_prop%dH*bygasc*(1.d0/vbs_prop%Tref-1.d0/vbs_cond%temp))
!===============================================================================
end subroutine clausius_clapeyron
!===============================================================================
!===============================================================================
end module TRACERS_VBS
!===============================================================================
