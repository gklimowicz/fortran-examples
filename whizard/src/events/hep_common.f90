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

module hep_common

  use kinds, only: default
  use kinds, only: double
  use constants
  use iso_varying_string, string_t => varying_string
  use debug_master, only: debug_on
  use lorentz
  use polarizations
  use model_data
  use particles

  implicit none
  private

  public :: heprup_init
  public :: assure_heprup
  public :: combine_lhef_with_particle_set
  public :: w2p_write_lhef_event
  public :: heprup_get_run_parameters
  public :: heprup_set_lhapdf_id
  public :: heprup_set_process_parameters
  public :: heprup_get_process_parameters
  public :: heprup_write_verbose
  public :: heprup_write_lhef
  public :: heprup_write_ascii
  public :: heprup_read_lhef
  public :: hepeup_init
  public :: hepeup_set_event_parameters
  public :: hepeup_get_event_parameters
  public :: hepeup_set_particle
  public :: hepeup_set_particle_lifetime
  public :: hepeup_set_particle_spin
  public :: hepeup_get_particle
  public :: hepevt_init
  public :: hepevt_set_event_parameters
  public :: hepevt_set_particle
  public :: hepevt_write_verbose
  public :: hepeup_write_verbose
  public :: hepeup_write_lhef
  public :: hepeup_write_lha
  public :: hepevt_write_hepevt
  public :: hepevt_write_ascii
  public :: hepevt_write_athena
  public :: hepevt_write_mokka
  public :: hepeup_read_lhef
  public :: hepeup_from_particle_set
  public :: hepeup_to_particle_set
  public :: hepevt_from_particle_set

  interface hepeup_set_particle_spin
     module procedure hepeup_set_particle_spin_pol
  end interface

  integer, parameter, public :: MAXNUP = 500
  integer, parameter, public :: MAXPUP = 100

  integer, public :: NUP
  integer, public :: IDPRUP
  double precision, public :: XWGTUP
  double precision, public :: SCALUP
  double precision, public :: AQEDUP
  double precision, public :: AQCDUP
  integer, dimension(MAXNUP), public :: IDUP
  integer, dimension(MAXNUP), public :: ISTUP
  integer, dimension(2,MAXNUP), public :: MOTHUP
  integer, dimension(2,MAXNUP), public :: ICOLUP
  double precision, dimension(5,MAXNUP), public :: PUP
  double precision, dimension(MAXNUP), public :: VTIMUP
  double precision, dimension(MAXNUP), public :: SPINUP
  integer, dimension(2), public :: IDBMUP
  double precision, dimension(2), public :: EBMUP
  integer, dimension(2), public :: PDFGUP
  integer, dimension(2), public :: PDFSUP
  integer, public :: IDWTUP
  integer, public :: NPRUP
  double precision, dimension(MAXPUP), public :: XSECUP
  double precision, dimension(MAXPUP), public :: XERRUP
  double precision, dimension(MAXPUP), public :: XMAXUP
  integer, dimension(MAXPUP), public :: LPRUP
  integer, parameter :: NMXHEP = 4000

  integer :: NEVHEP

  integer, public :: NHEP

  integer, dimension(NMXHEP), public :: ISTHEP

  integer, dimension(NMXHEP), public :: IDHEP

  integer, dimension(2, NMXHEP), public :: JMOHEP

  integer, dimension(2, NMXHEP), public :: JDAHEP

  double precision, dimension(5, NMXHEP), public :: PHEP

  double precision, dimension(4, NMXHEP) :: VHEP

  integer, dimension(NMXHEP) :: hepevt_pol

  integer, public :: idruplh

  double precision, public :: eventweightlh

  double precision, public :: alphaqedlh, alphaqcdlh

  double precision, dimension(10), public :: scalelh

  double precision, dimension (3,NMXHEP), public :: spinlh
  integer, dimension (2,NMXHEP), public :: icolorflowlh

  integer :: hepevt_n_out, hepevt_n_remnants

  double precision :: hepevt_weight, hepevt_function_value
  double precision :: hepevt_function_ratio


  common /HEPRUP/ &
       IDBMUP, EBMUP, PDFGUP, PDFSUP, IDWTUP, NPRUP, &
       XSECUP, XERRUP, XMAXUP, LPRUP
  save /HEPRUP/

  common /HEPEUP/ &
       NUP, IDPRUP, XWGTUP, SCALUP, AQEDUP, AQCDUP, &
       IDUP, ISTUP, MOTHUP, ICOLUP, PUP, VTIMUP, SPINUP
  save /HEPEUP/

  common /HEPEVT/ &
       NEVHEP, NHEP, ISTHEP, IDHEP, &
       JMOHEP, JDAHEP, PHEP, VHEP
  save /HEPEVT/

  common /HEPEV4/ &
       eventweightlh, alphaqedlh, alphaqcdlh, scalelh, &
       spinlh, icolorflowlh, idruplh
  save /HEPEV4/


  interface
    module subroutine heprup_init &
         (beam_pdg, beam_energy, n_processes, unweighted, negative_weights)
      integer, dimension(2), intent(in) :: beam_pdg
      real(default), dimension(2), intent(in) :: beam_energy
      integer, intent(in) :: n_processes
      logical, intent(in) :: unweighted
      logical, intent(in) :: negative_weights
    end subroutine heprup_init
    module subroutine assure_heprup (pset)
      type(particle_set_t), intent(in) :: pset
    end subroutine assure_heprup
    module subroutine combine_lhef_with_particle_set &
         (particle_set, u, model_in, model_hadrons)
      type(particle_set_t), intent(inout) :: particle_set
      integer, intent(in) :: u
      class(model_data_t), intent(in), target :: model_in
      class(model_data_t), intent(in), target :: model_hadrons
    end subroutine combine_lhef_with_particle_set
    module subroutine w2p_write_lhef_event (unit)
      integer, intent(in) :: unit
    end subroutine w2p_write_lhef_event
    module subroutine heprup_get_run_parameters &
         (beam_pdg, beam_energy, n_processes, unweighted, negative_weights)
      integer, dimension(2), intent(out), optional :: beam_pdg
      real(default), dimension(2), intent(out), optional :: beam_energy
      integer, intent(out), optional :: n_processes
      logical, intent(out), optional :: unweighted
      logical, intent(out), optional :: negative_weights
    end subroutine heprup_get_run_parameters
    module subroutine heprup_set_lhapdf_id (i_beam, pdf_id)
      integer, intent(in) :: i_beam, pdf_id
    end subroutine heprup_set_lhapdf_id
    module subroutine heprup_set_process_parameters &
         (i, process_id, cross_section, error, max_weight, is_width)
      integer, intent(in) :: i, process_id
      real(default), intent(in), optional :: cross_section, error, max_weight
      logical, intent(in), optional :: is_width
    end subroutine heprup_set_process_parameters
    module subroutine heprup_get_process_parameters  &
         (i, process_id, cross_section, error, max_weight, is_width)
      integer, intent(in) :: i
      integer, intent(out), optional :: process_id
      real(default), intent(out), optional :: cross_section, error, max_weight
      logical, intent(in), optional :: is_width
    end subroutine heprup_get_process_parameters
    module subroutine heprup_write_verbose (unit)
      integer, intent(in), optional :: unit
    end subroutine heprup_write_verbose
    module subroutine heprup_write_lhef (unit)
      integer, intent(in), optional :: unit
    end subroutine heprup_write_lhef
    module subroutine heprup_write_ascii (unit)
      integer, intent(in), optional :: unit
    end subroutine heprup_write_ascii
    module subroutine heprup_read_lhef (u)
      integer, intent(in) :: u
    end subroutine heprup_read_lhef
    module subroutine hepeup_init (n_tot)
      integer, intent(in) :: n_tot
    end subroutine hepeup_init
    module subroutine hepeup_set_event_parameters &
         (proc_id, weight, scale, alpha_qed, alpha_qcd)
      integer, intent(in), optional :: proc_id
      real(default), intent(in), optional :: &
           weight, scale, alpha_qed, alpha_qcd
    end subroutine hepeup_set_event_parameters
    module subroutine hepeup_get_event_parameters &
         (proc_id, weight, scale, alpha_qed, alpha_qcd)
      integer, intent(out), optional :: proc_id
      real(default), intent(out), optional :: &
           weight, scale, alpha_qed, alpha_qcd
    end subroutine hepeup_get_event_parameters
    module subroutine hepeup_set_particle (i, pdg, status, parent, col, p, m2)
      integer, intent(in) :: i
      integer, intent(in) :: pdg, status
      integer, dimension(:), intent(in) :: parent
      type(vector4_t), intent(in) :: p
      integer, dimension(2), intent(in) :: col
      real(default), intent(in) :: m2
    end subroutine hepeup_set_particle
    module subroutine hepeup_set_particle_lifetime (i, lifetime)
      integer, intent(in) :: i
      real(default), intent(in) :: lifetime
    end subroutine hepeup_set_particle_lifetime
    module subroutine hepeup_set_particle_spin_pol (i, p, pol, p_mother)
      integer, intent(in) :: i
      type(vector4_t), intent(in) :: p
      type(polarization_t), intent(in) :: pol
      type(vector4_t), intent(in) :: p_mother
    end subroutine hepeup_set_particle_spin_pol
    module subroutine hepeup_get_particle (i, pdg, status, parent, col, p, m2)
      integer, intent(in) :: i
      integer, intent(out), optional :: pdg, status
      integer, dimension(:), intent(out), optional :: parent
      type(vector4_t), intent(out), optional :: p
      integer, dimension(2), intent(out), optional :: col
      real(default), dimension(5,MAXNUP) :: pup_def
      real(default), intent(out), optional :: m2
    end subroutine hepeup_get_particle
    module subroutine hepevt_init (n_tot, n_out)
      integer, intent(in) :: n_tot, n_out
    end subroutine hepevt_init
    module subroutine hepevt_set_event_parameters &
         (proc_id, weight, function_value, function_ratio, &
         alpha_qcd, alpha_qed, scale, i_evt)
      integer, intent(in), optional :: proc_id
      integer, intent(in), optional :: i_evt
      real(default), intent(in), optional :: weight, function_value, &
         function_ratio, alpha_qcd, alpha_qed, scale
    end subroutine hepevt_set_event_parameters
    module subroutine hepevt_set_particle &
         (i, pdg, status, parent, child, p, m2, hel, vtx, &
         col, pol_status, pol, fill_hepev4)
      integer, intent(in) :: i
      integer, intent(in) :: pdg, status
      integer, dimension(:), intent(in) :: parent
      integer, dimension(:), intent(in) :: child
      logical, intent(in), optional :: fill_hepev4
      type(vector4_t), intent(in) :: p
      real(default), intent(in) :: m2
      integer, dimension(2), intent(in) :: col
      integer, intent(in) :: pol_status
      integer, intent(in) :: hel
      type(polarization_t), intent(in), optional :: pol
      type(vector4_t), intent(in) :: vtx
    end subroutine hepevt_set_particle
    module subroutine hepevt_write_verbose (unit)
      integer, intent(in), optional :: unit
    end subroutine hepevt_write_verbose
    module subroutine hepeup_write_verbose (unit)
      integer, intent(in), optional :: unit
    end subroutine hepeup_write_verbose
    module subroutine hepeup_write_lhef (unit)
      integer, intent(in), optional :: unit
    end subroutine hepeup_write_lhef
    module subroutine hepeup_write_lha (unit)
      integer, intent(in), optional :: unit
    end subroutine hepeup_write_lha
    module subroutine hepevt_write_hepevt (unit)
      integer, intent(in), optional :: unit
    end subroutine hepevt_write_hepevt
    module subroutine hepevt_write_ascii (unit, long)
      integer, intent(in), optional :: unit
      logical, intent(in) :: long
    end subroutine hepevt_write_ascii
    module subroutine hepevt_write_athena (unit)
      integer, intent(in), optional :: unit
    end subroutine hepevt_write_athena
    module subroutine hepevt_write_mokka (unit)
      integer, intent(in), optional :: unit
    end subroutine hepevt_write_mokka
    module subroutine hepeup_read_lhef (u)
      integer, intent(in) :: u
    end subroutine hepeup_read_lhef
    module subroutine hepeup_from_particle_set (pset_in, &
       keep_beams, keep_remnants, tauola_convention)
      type(particle_set_t), intent(in) :: pset_in
      type(particle_set_t), target :: pset
      logical, intent(in), optional :: keep_beams
      logical, intent(in), optional :: keep_remnants
      logical, intent(in), optional :: tauola_convention
    end subroutine hepeup_from_particle_set
    module subroutine hepeup_to_particle_set &
         (particle_set, recover_beams, model, alt_model)
      type(particle_set_t), intent(inout), target :: particle_set
      logical, intent(in), optional :: recover_beams
      class(model_data_t), intent(in), target :: model, alt_model
    end subroutine hepeup_to_particle_set
    module subroutine hepevt_from_particle_set &
         (particle_set, keep_beams, keep_remnants, ensure_order, fill_hepev4)
      type(particle_set_t), intent(in) :: particle_set
      type(particle_set_t), target :: pset_hepevt, pset_tmp
      logical, intent(in), optional :: keep_beams
      logical, intent(in), optional :: keep_remnants
      logical, intent(in), optional :: ensure_order
      logical, intent(in), optional :: fill_hepev4
    end subroutine hepevt_from_particle_set
  end interface

end module hep_common
