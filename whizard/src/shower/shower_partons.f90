! WHIZARD 3.1.0 Dec 14 2022
!
! Copyright (C) 1999-2022 by
!     Wolfgang Kilian <kilian@physik.uni-siegen.de>
!     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!     Juergen Reuter <juergen.reuter@desy.de>
!
!     with contributions from
!     cf. main AUTHORS file
!
! WHIZARD is free software; you can redistribute it and/or modify it
! under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2, or (at your option)
! any later version.
!
! WHIZARD is distributed in the hope that it will be useful, but
! WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file has been stripped of most comments.  For documentation, refer
! to the source 'whizard.nw'

module shower_partons

  use kinds, only: default, double
  use constants
  use lorentz
  use particles
  use model_data
  use shower_base
  use rng_base

  implicit none
  private

  public :: parton_t
  public :: parton_pointer_t
  public :: parton_of_particle
  public :: parton_copy
  public :: parton_set_parent
  public :: parton_get_parent
  public :: parton_set_initial
  public :: parton_get_initial
  public :: parton_set_child
  public :: parton_get_child
  public :: P_prt_to_child1
  public :: thetabar
  public :: parton_apply_costheta
  public :: parton_apply_lorentztrafo
  public :: parton_apply_lorentztrafo_recursive
  public :: parton_simulate_stept
  public :: maxzz

  type :: parton_t
     integer :: nr = 0
     integer :: type = 0
     type(shower_settings_t), pointer :: settings => null()
     type(vector4_t) :: momentum = vector4_null
     real(default) :: t  = zero
     real(default) :: mass2  = zero
     real(default) :: scale = zero
     real(default) :: z = zero
     real(default) :: costheta = zero
     real(default) :: x = zero
     logical :: simulated = .false.
     logical :: belongstoFSR = .true.
     logical :: belongstointeraction = .false.
     type(parton_t), pointer :: parent => null ()
     type(parton_t), pointer :: child1 => null ()
     type(parton_t), pointer :: child2 => null ()
     type(parton_t), pointer :: initial => null ()
     integer :: c1 = 0, c2 = 0
     integer :: aux_pt = 0
     integer :: ckkwlabel = 0
     real(default) :: ckkwscale = zero
     integer :: ckkwtype = -1
     integer :: interactionnr = 0
   contains
     procedure :: to_particle => parton_to_particle
     procedure :: to_status => parton_to_status
     procedure :: to_color => parton_to_color
     procedure :: get_costheta => parton_get_costheta
     procedure :: get_costheta_mass => parton_get_costheta_mass
     procedure :: get_costheta_motherfirst => parton_get_costheta_motherfirst
     procedure :: get_beta => parton_get_beta
     procedure :: write => parton_write
     procedure :: is_final => parton_is_final
     procedure :: is_branched => parton_is_branched
     procedure :: set_simulated => parton_set_simulated
     procedure :: is_quark => parton_is_quark
     procedure :: is_squark => parton_is_squark
     procedure :: is_gluon => parton_is_gluon
     procedure :: is_gluino => parton_is_gluino
     procedure :: is_proton => parton_is_proton
     procedure :: is_colored => parton_is_colored
     procedure :: mass => parton_mass
     procedure :: mass_squared => parton_mass_squared
     procedure :: momentum_to_pythia6 => parton_momentum_to_pythia6
     procedure :: generate_ps => parton_generate_ps
     procedure :: generate_ps_ini => parton_generate_ps_ini
     procedure :: next_t_ana => parton_next_t_ana
  end type parton_t

  type :: parton_pointer_t
     type(parton_t), pointer :: p => null ()
  end type parton_pointer_t


  interface
    module function parton_to_particle &
         (parton, model, from_hard_int) result (particle)
      type(particle_t) :: particle
      class(parton_t), intent(in) :: parton
      class(model_data_t), pointer, intent(in) :: model
      logical, intent(in), optional :: from_hard_int
    end function parton_to_particle
    pure module function parton_of_particle (particle, nr) result (parton)
      type(parton_t) :: parton
      type(particle_t), intent(in) :: particle
      integer, intent(in) :: nr
    end function parton_of_particle
    pure module function parton_to_status &
         (parton, from_hard_int) result (status)
      integer :: status
      class(parton_t), intent(in) :: parton
      logical, intent(in), optional :: from_hard_int
    end function parton_to_status
    pure module subroutine parton_to_color (parton, c1, c2, from_hard_int)
      class(parton_t), intent(in) :: parton
      integer, intent(out) :: c1, c2
      logical, intent(in), optional :: from_hard_int
    end subroutine parton_to_color
    module subroutine parton_copy (prt1, prt2)
      type(parton_t), intent(in) :: prt1
      type(parton_t), intent(out) :: prt2
    end subroutine parton_copy
    elemental module function parton_get_costheta (prt) result (costheta)
      class(parton_t), intent(in) :: prt
      real(default) :: costheta
    end function parton_get_costheta
    elemental module function parton_get_costheta_mass (prt) result (costheta)
      class(parton_t), intent(in) :: prt
      real(default) :: costheta
    end function parton_get_costheta_mass
    elemental module function parton_get_costheta_motherfirst &
         (prt) result (costheta)
      class(parton_t), intent(in) :: prt
      real(default) :: costheta
    end function parton_get_costheta_motherfirst
    elemental module function parton_get_beta (prt) result (beta)
      class(parton_t), intent(in) :: prt
      real(default) :: beta
    end function parton_get_beta
    module subroutine parton_write (prt, unit)
      class(parton_t), intent(in) :: prt
      integer, intent(in), optional :: unit
    end subroutine parton_write
    elemental module function parton_is_final (prt) result (is_final)
      class(parton_t), intent(in) :: prt
      logical :: is_final
    end function parton_is_final
    elemental module function parton_is_branched (prt) result (is_branched)
      class(parton_t), intent(in) :: prt
      logical :: is_branched
    end function parton_is_branched
    pure module subroutine parton_set_simulated (prt, sim)
      class(parton_t), intent(inout) :: prt
      logical, intent(in), optional :: sim
    end subroutine parton_set_simulated
    module subroutine parton_set_parent (prt, parent)
      type(parton_t), intent(inout) :: prt
      type(parton_t), intent(in) , target :: parent
    end subroutine parton_set_parent
    module function parton_get_parent (prt) result (parent)
      type(parton_t), intent(in) :: prt
      type(parton_t), pointer :: parent
    end function parton_get_parent
    module subroutine parton_set_initial (prt, initial)
      type(parton_t), intent(inout) :: prt
      type(parton_t), intent(in) , target :: initial
    end subroutine parton_set_initial
    module function parton_get_initial (prt) result (initial)
      type(parton_t), intent(in) :: prt
      type(parton_t), pointer :: initial
    end function parton_get_initial
    module subroutine parton_set_child (prt, child, i)
      type(parton_t), intent(inout) :: prt
      type(parton_t), intent(in), target :: child
      integer, intent(in) ::  i
    end subroutine parton_set_child
    module function parton_get_child (prt, i) result (child)
      type(parton_t), pointer :: child
      type(parton_t), intent(in) :: prt
      integer, intent(in) :: i
    end function parton_get_child
    elemental module function parton_is_quark (prt) result (is_quark)
      class(parton_t), intent(in) ::prt
      logical :: is_quark
    end function parton_is_quark
    elemental module function parton_is_squark (prt) result (is_squark)
      class(parton_t), intent(in) ::prt
      logical :: is_squark
    end function parton_is_squark
    elemental module function parton_is_gluon (prt) result (is_gluon)
      class(parton_t), intent(in) :: prt
      logical :: is_gluon
    end function parton_is_gluon
    elemental module function parton_is_gluino (prt) result (is_gluino)
      class(parton_t), intent(in) :: prt
      logical :: is_gluino
    end function parton_is_gluino
    elemental module function parton_is_proton (prt) result (is_hadron)
      class(parton_t), intent(in) :: prt
      logical :: is_hadron
    end function parton_is_proton
    pure module function parton_is_colored (parton) result (is_colored)
      logical :: is_colored
      class(parton_t), intent(in) :: parton
    end function parton_is_colored
    elemental module function parton_mass (prt) result (mass)
      class(parton_t), intent(in) :: prt
      real(default) :: mass
    end function parton_mass
    elemental module function parton_mass_squared (prt) result (mass_squared)
      class(parton_t), intent(in) :: prt
      real(default) :: mass_squared
    end function parton_mass_squared
    pure module function parton_momentum_to_pythia6 (prt) result (p)
      real(double), dimension(1:5) :: p
      class(parton_t), intent(in) :: prt
    end function parton_momentum_to_pythia6
    module function P_prt_to_child1 (prt) result (retvalue)
      type(parton_t), intent(in) :: prt
      real(default) :: retvalue
    end function P_prt_to_child1
    module function thetabar (prt, recoiler, isr_ang, E3out) result (retvalue)
      type(parton_t), intent(inout) :: prt
      type(parton_t), intent(in) :: recoiler
      real(default), intent(out), optional :: E3out
      logical, intent(in) :: isr_ang
      logical :: retvalue
    end function thetabar
    recursive module subroutine parton_apply_costheta (prt, rng)
      type(parton_t), intent(inout) :: prt
      class(rng_t), intent(inout), allocatable :: rng
    end subroutine parton_apply_costheta
    module subroutine parton_apply_lorentztrafo (prt, L)
      type(parton_t), intent(inout) :: prt
      type(lorentz_transformation_t), intent(in) :: L
    end subroutine parton_apply_lorentztrafo
    recursive module subroutine parton_apply_lorentztrafo_recursive (prt, L)
      type(parton_t), intent(inout) :: prt
      type(lorentz_transformation_t) ,intent(in) :: L
    end subroutine parton_apply_lorentztrafo_recursive
    module subroutine parton_generate_ps (prt, rng)
      class(parton_t), intent(inout) :: prt
      class(rng_t), intent(inout), allocatable :: rng
    end subroutine parton_generate_ps
    module subroutine parton_generate_ps_ini (prt, rng)
      class(parton_t), intent(inout) :: prt
      class(rng_t), intent(inout), allocatable :: rng
    end subroutine parton_generate_ps_ini
    module subroutine parton_next_t_ana (prt, rng)
      class(parton_t), intent(inout) :: prt
      class(rng_t), intent(inout), allocatable :: rng
    end subroutine parton_next_t_ana
    module subroutine parton_simulate_stept &
         (prt, rng, integral, random, gtoqq, lookatsister)
      type(parton_t), intent(inout) :: prt
      class(rng_t), intent(inout), allocatable :: rng
      real(default), intent(inout) :: integral
      real(default), intent(inout) :: random
      integer, intent(out) :: gtoqq
      logical, intent(in), optional :: lookatsister
    end subroutine parton_simulate_stept
    module function maxzz (shat, s, maxz_isr, minenergy_timelike) result (maxz)
      real(default), intent(in) :: shat, s, minenergy_timelike, maxz_isr
      real(default) :: maxz
    end function maxzz
  end interface

end module shower_partons
