
! Copyright (C) 2010 J. K. Dewhurst, S. Sharma and E. K. U. Gross.
! This file is distributed under the terms of the GNU General Public License.
! See the file COPYING for license details.

module modtddft

!-----------------------------------------!
!     TDDFT linear response variables     !
!-----------------------------------------!
! exchange-correlation kernel type
integer fxctype(3)
! parameters for long-range correction (LRC) kernel
real(8) fxclrc(2)
! number of independent spin components of the f_xc spin tensor
integer nscfxc
! magnetic linear dichroism (MLD) angle between the electric and magnetic fields
real(8) thetamld

!---------------------------------------------!
!     TDDFT real-time evolution variables     !
!---------------------------------------------!
! number of laser pulses defining the time-dependent A-field
integer npulse
! laser pulse parameters: vector amplitude, frequency, phase, chirp rate,
! peak time, full-width at half-maximum, spin components
real(8), allocatable :: pulse(:,:)
! number of A-field ramps
integer nramp
! ramp parameters: vector amplitude, ramp start time, linear, quadratic, cubic,
! quartic coefficients, spin components
real(8), allocatable :: ramp(:,:)
! total simulation time
real(8) tstime
! time step length
real(8) dtimes
! number of time steps
integer ntimes
! starting time step
integer itimes0
! current time step
integer itimes
! time steps
real(8), allocatable :: times(:)
! tdt0 is .true. if the time-evolution starts at t=0
logical tdt0
! phase defining complex direction of time evolution
real(8) tdphi
! number of time steps after which the time-dependent eigenvectors are backed up
integer ntsbackup
! tafieldt is .true. if a time-dependent vector potential is applied
logical tafieldt
! time-dependent total A-field (external plus induced)
real(8), allocatable :: afieldt(:,:)
! induced A-field and its time derivative at the current time step
real(8) afindt(3,0:1)
! induced A-field parameters
real(8) afindpm(0:2)
! if tafindt is .true. then the induced A-field is determined from Maxwell's
! equation and added to the total
logical tafindt
! electric field at current time step
real(8) efieldt(3)
! observables are written to file every ntswrite(1) time steps; this begins
! at or after time step ntswrite(2); writing occurs at the first time step
! irrespective of ntswrite
integer ntswrite(2)
! static charge density
real(8), allocatable :: rhosmt(:,:,:),rhosir(:,:)
! total static charge
real(8) chgstot(3)
! muffin-tin static charge
real(8), allocatable :: chgsmt(:,:)
! gradient of the muffin-tin Coulomb potential of the nucleus and static density
complex(8), allocatable :: gvnsmt(:,:,:)
! the following variables are .true. if the corresponding quantities are to be
! written every ntswrite time steps
logical tdrho1d,tdrho2d,tdrho3d
logical tdmag1d,tdmag2d,tdmag3d
logical tdjr1d,tdjr2d,tdjr3d
logical tddos,tdlsj
! magnitude of complex numbers added to initial eigenvectors
real(8) rndevt0
! number of time steps between force calculations
integer ntsforce
! total atomic forces at each time step
real(8), allocatable :: forcet(:,:,:)

end module

