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

module sf_lhapdf

  use kinds, only: default
  use iso_varying_string, string_t => varying_string
  use sm_qcd
  use pdg_arrays
  use model_data
  use flavors
  use polarizations
  use sf_base
  use lhapdf !NODEP!

  implicit none
  private

  public :: lhapdf_global_reset
  public :: lhapdf_initialize
  public :: lhapdf_data_t
  public :: lhapdf_t
  public :: alpha_qcd_lhapdf_t

  type :: lhapdf_global_status_t
     private
     logical, dimension(3) :: initialized = .false.
  end type lhapdf_global_status_t

  type, extends (sf_data_t) :: lhapdf_data_t
     private
     type(string_t) :: prefix
     type(string_t) :: file
     type(lhapdf_pdf_t) :: pdf
     integer :: member = 0
     class(model_data_t), pointer :: model => null ()
     type(flavor_t) :: flv_in
     integer :: set = 0
     logical :: invert = .false.
     logical :: photon = .false.
     logical :: has_photon = .false.
     integer :: photon_scheme = 0
     real(default) :: xmin = 0, xmax = 0
     real(default) :: qmin = 0, qmax = 0
     logical, dimension(-6:6) :: mask = .true.
     logical :: mask_photon = .true.
     logical :: hoppet_b_matching = .false.
   contains
       procedure :: init => lhapdf_data_init
       procedure :: write => lhapdf_data_write
       procedure :: get_n_par => lhapdf_data_get_n_par
       procedure :: get_pdg_out => lhapdf_data_get_pdg_out
       procedure :: allocate_sf_int => lhapdf_data_allocate_sf_int
       procedure :: get_pdf_set => lhapdf_data_get_pdf_set
  end type lhapdf_data_t

  type, extends (sf_int_t) :: lhapdf_t
     type(lhapdf_data_t), pointer :: data => null ()
     real(default) :: x = 0
     real(default) :: q = 0
     real(default) :: s = 0
   contains
     procedure :: complete_kinematics => lhapdf_complete_kinematics
     procedure :: recover_x => lhapdf_recover_x
     procedure :: inverse_kinematics => lhapdf_inverse_kinematics
     procedure :: type_string => lhapdf_type_string
     procedure :: write => lhapdf_write
     procedure :: init => lhapdf_init
     procedure :: apply => lhapdf_apply
  end type lhapdf_t

  type, extends (alpha_qcd_t) :: alpha_qcd_lhapdf_t
     type(string_t) :: pdfset_dir
     type(string_t) :: pdfset_file
     integer :: pdfset_member = -1
     type(lhapdf_pdf_t) :: pdf
   contains
     procedure :: write => alpha_qcd_lhapdf_write
     procedure :: get => alpha_qcd_lhapdf_get
     procedure :: init => alpha_qcd_lhapdf_init
     procedure :: get_qmass => alpha_qcd_lhapdf_get_qmass
     procedure :: get_order => alpha_qcd_lhapdf_get_order
  end type alpha_qcd_lhapdf_t


  type(lhapdf_global_status_t), save :: lhapdf_global_status

  interface
     subroutine InitPDFsetM (set, file)
       integer, intent(in) :: set
       character(*), intent(in) :: file
     end subroutine InitPDFsetM
  end interface

  interface
     subroutine InitPDFM (set, mem)
       integer, intent(in) :: set, mem
     end subroutine InitPDFM
  end interface

  interface
     subroutine numberPDFM (set, n_members)
       integer, intent(in) :: set
       integer, intent(out) :: n_members
     end subroutine numberPDFM
  end interface

  interface
     subroutine evolvePDFM (set, x, q, ff)
       integer, intent(in) :: set
       double precision, intent(in) :: x, q
       double precision, dimension(-6:6), intent(out) :: ff
     end subroutine evolvePDFM
  end interface

  interface
     subroutine evolvePDFphotonM (set, x, q, ff, fphot)
       integer, intent(in) :: set
       double precision, intent(in) :: x, q
       double precision, dimension(-6:6), intent(out) :: ff
       double precision, intent(out) :: fphot
     end subroutine evolvePDFphotonM
  end interface

  interface
     subroutine evolvePDFpM (set, x, q, s, scheme, ff)
       integer, intent(in) :: set
       double precision, intent(in) :: x, q, s
       integer, intent(in) :: scheme
       double precision, dimension(-6:6), intent(out) :: ff
     end subroutine evolvePDFpM
  end interface

  interface
     subroutine GetXminM (set, mem, xmin)
       integer, intent(in) :: set, mem
       double precision, intent(out) :: xmin
     end subroutine GetXminM
  end interface

  interface
     subroutine GetXmaxM (set, mem, xmax)
       integer, intent(in) :: set, mem
       double precision, intent(out) :: xmax
     end subroutine GetXmaxM
  end interface

  interface
     subroutine GetQ2minM (set, mem, q2min)
       integer, intent(in) :: set, mem
       double precision, intent(out) :: q2min
     end subroutine GetQ2minM
  end interface

  interface
     subroutine GetQ2maxM (set, mem, q2max)
       integer, intent(in) :: set, mem
       double precision, intent(out) :: q2max
     end subroutine GetQ2maxM
  end interface

  interface
     function has_photon () result(flag)
        logical :: flag
     end function has_photon
  end interface

  interface
     double precision function alphasPDF (Q)
       double precision, intent(in) :: Q
     end function alphasPDF
  end interface


  interface
    module subroutine lhapdf_global_reset ()
    end subroutine lhapdf_global_reset
    module subroutine lhapdf_initialize &
         (set, prefix, file, member, pdf, b_match)
      integer, intent(in) :: set
      type(string_t), intent(inout) :: prefix
      type(string_t), intent(inout) :: file
      type(lhapdf_pdf_t), intent(inout), optional :: pdf
      integer, intent(inout) :: member
      logical, intent(in), optional :: b_match
    end subroutine lhapdf_initialize
    module subroutine lhapdf_complete_kinematics &
         (sf_int, x, xb, f, r, rb, map)
      class(lhapdf_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: r
      real(default), dimension(:), intent(in) :: rb
      logical, intent(in) :: map
    end subroutine lhapdf_complete_kinematics
    module subroutine lhapdf_recover_x (sf_int, x, xb, x_free)
      class(lhapdf_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(out) :: x
      real(default), dimension(:), intent(out) :: xb
      real(default), intent(inout), optional :: x_free
    end subroutine lhapdf_recover_x
    module subroutine lhapdf_inverse_kinematics &
         (sf_int, x, xb, f, r, rb, map, set_momenta)
      class(lhapdf_t), intent(inout) :: sf_int
      real(default), dimension(:), intent(in) :: x
      real(default), dimension(:), intent(in) :: xb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: r
      real(default), dimension(:), intent(out) :: rb
      logical, intent(in) :: map
      logical, intent(in), optional :: set_momenta
    end subroutine lhapdf_inverse_kinematics
    module subroutine lhapdf_data_init &
         (data, model, pdg_in, prefix, file, member, photon_scheme, &
              hoppet_b_matching)
      class(lhapdf_data_t), intent(out) :: data
      class(model_data_t), intent(in), target :: model
      type(pdg_array_t), intent(in) :: pdg_in
      type(string_t), intent(in), optional :: prefix, file
      integer, intent(in), optional :: member
      integer, intent(in), optional :: photon_scheme
      logical, intent(in), optional :: hoppet_b_matching
    end subroutine lhapdf_data_init
    module subroutine lhapdf_data_write (data, unit, verbose)
      class(lhapdf_data_t), intent(in) :: data
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: verbose
    end subroutine lhapdf_data_write
    module function lhapdf_data_get_n_par (data) result (n)
      class(lhapdf_data_t), intent(in) :: data
      integer :: n
    end function lhapdf_data_get_n_par
    module subroutine lhapdf_data_get_pdg_out (data, pdg_out)
      class(lhapdf_data_t), intent(in) :: data
      type(pdg_array_t), dimension(:), intent(inout) :: pdg_out
    end subroutine lhapdf_data_get_pdg_out
    elemental module function lhapdf_data_get_pdf_set (data) result (pdf_set)
      class(lhapdf_data_t), intent(in) :: data
      integer :: pdf_set
    end function lhapdf_data_get_pdf_set
    module function lhapdf_type_string (object) result (string)
      class(lhapdf_t), intent(in) :: object
      type(string_t) :: string
    end function lhapdf_type_string
    module subroutine lhapdf_write (object, unit, testflag)
      class(lhapdf_t), intent(in) :: object
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: testflag
    end subroutine lhapdf_write
    module subroutine lhapdf_init (sf_int, data)
      class(lhapdf_t), intent(out) :: sf_int
      class(sf_data_t), intent(in), target :: data
    end subroutine lhapdf_init
    module subroutine lhapdf_apply &
         (sf_int, scale, negative_sf, rescale, i_sub)
      class(lhapdf_t), intent(inout) :: sf_int
      real(default), intent(in) :: scale
      logical, intent(in), optional :: negative_sf
      class(sf_rescale_t), intent(in), optional :: rescale
      integer, intent(in), optional :: i_sub
    end subroutine lhapdf_apply
    module subroutine alpha_qcd_lhapdf_write (object, unit)
      class(alpha_qcd_lhapdf_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine alpha_qcd_lhapdf_write
    module function alpha_qcd_lhapdf_get (alpha_qcd, scale) result (alpha)
      class(alpha_qcd_lhapdf_t), intent(in) :: alpha_qcd
      real(default), intent(in) :: scale
      real(default) :: alpha
    end function alpha_qcd_lhapdf_get
    module subroutine alpha_qcd_lhapdf_init (alpha_qcd, file, member, path)
      class(alpha_qcd_lhapdf_t), intent(out) :: alpha_qcd
      type(string_t), intent(inout) :: file
      integer, intent(inout) :: member
      type(string_t), intent(inout) :: path
    end subroutine alpha_qcd_lhapdf_init
    module function alpha_qcd_lhapdf_get_qmass (alpha_qcd, i_q) result (mq)
      real(default) :: mq
      class(alpha_qcd_lhapdf_t), intent(in) :: alpha_qcd
      integer, intent(in) :: i_q
    end function alpha_qcd_lhapdf_get_qmass
    module function alpha_qcd_lhapdf_get_order (alpha_qcd) result (order)
      integer :: order
      class(alpha_qcd_lhapdf_t), intent(in) :: alpha_qcd
    end function alpha_qcd_lhapdf_get_order
  end interface

contains

  subroutine lhapdf_data_allocate_sf_int (data, sf_int)
    class(lhapdf_data_t), intent(in) :: data
    class(sf_int_t), intent(inout), allocatable :: sf_int
    allocate (lhapdf_t :: sf_int)
  end subroutine lhapdf_data_allocate_sf_int


end module sf_lhapdf
