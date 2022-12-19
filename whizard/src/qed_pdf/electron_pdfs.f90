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

module electron_pdfs

  use kinds, only: default
  use io_units
  use sm_qed
  
  implicit none
  private

  public :: qed_pdf_t
  public :: elec_asym
  public :: phot_asym
  public :: bar_asym
  public :: rec_num
  public :: recbar
  public :: rechat
  public :: elecbar_asym_p
  public :: photbar_asym_p
  public :: recbar_singlet
  public :: recbar_nonsinglet
  public :: recbar_photon
  public :: rechat_singlet
  public :: rechat_nonsinglet
  public :: rechat_photon
  public :: t_alpha
  public :: endpoint_func_S
  public :: endpoint_func_NS
  public :: endpoint_func_GAM
  public :: elec_pdf

  integer, parameter, public :: EPDF_ELE = 0, EPDF_POS = 1, &
       EPDF_S = 2, EPDF_NS = 3, EPDF_G = 4


  type :: qed_pdf_t     
     private
     integer :: flv = 0
     class(alpha_qed_t), allocatable :: aqed
     real(default) :: mass = 0
     real(default) :: q_max = 0
     real(default) :: alpha = 0
     real(default) :: eps = 0
     real(default), allocatable :: q_in
     integer :: order
     integer :: log_order
     integer :: n_lep
   contains
     procedure :: init => qed_pdf_init
     procedure :: write => qed_pdf_write
     procedure :: set_order => qed_pdf_set_order
     procedure :: evolve_qed_pdf => qed_pdf_evolve_qed_pdf
     procedure :: allocate_aqed => qed_pdf_allocate_aqed
  end type qed_pdf_t


  interface
    module subroutine qed_pdf_init &
         (qed_pdf, mass, alpha, charge, q_max, order, log_order, n_lep)
      class(qed_pdf_t), intent(out) :: qed_pdf
      real(default), intent(in) :: mass, alpha, q_max, charge
      integer, intent(in) :: order, log_order, n_lep
    end subroutine qed_pdf_init
    module subroutine qed_pdf_write (qed_pdf, unit, with_qed)
      class(qed_pdf_t), intent(in) :: qed_pdf
      integer, intent(in), optional :: unit
      logical, intent(in), optional :: with_qed
    end subroutine qed_pdf_write
    module subroutine qed_pdf_set_order (qed_pdf, order)
      class(qed_pdf_t), intent(inout) :: qed_pdf
      integer, intent(in) :: order
    end subroutine qed_pdf_set_order
    module subroutine qed_pdf_evolve_qed_pdf (qed_pdf, x, xb, rb, ff)
      class(qed_pdf_t), intent(inout) :: qed_pdf
      real(default), intent(in) :: x, xb, rb
      real(default), intent(inout) :: ff
    end subroutine qed_pdf_evolve_qed_pdf
    module function elec_asym (epdf, x, scale, alpha, running) result (elec_as)
      type(qed_pdf_t), intent(in) :: epdf
      real(default) :: elec_as
      real(default), intent(in) :: x
      real(default), intent(in) :: scale
      real(default), intent(in) :: alpha
      logical, intent(in) :: running
    end function elec_asym
    module function phot_asym &
         (epdf, x, scale, alpha, nlep, running) result (phot_as)
      type(qed_pdf_t), intent(in) :: epdf
      real(default) :: phot_as
      real(default), intent(in) :: x
      real(default), intent(in) :: scale
      real(default), intent(in) :: alpha
      integer, intent(in) :: nlep
      logical, intent(in) :: running
    end function phot_asym
    module function bar_asym &
         (epdf, flv, x, scale, alpha, running) result (bar_as)
      type(qed_pdf_t), intent(in) :: epdf
      integer, intent(in) :: flv
      real(default) :: bar_as
      real(default), intent(in) :: x
      real(default), intent(in) :: scale
      real(default), intent(in) :: alpha
      logical, intent(in) :: running
    end function bar_asym
    module function rec_num &
         (epdf, flv, x, scale, alpha, running) result (recnum)
      type(qed_pdf_t), intent(in) :: epdf
      integer, intent(in) :: flv
      real(default) :: recnum
      real(default), intent(in) :: x
      real(default), intent(in) :: scale
      real(default), intent(in) :: alpha
      logical, intent(in) :: running
    end function rec_num
    module function recbar &
         (epdf, flv, x, scale, alpha, running) result (bar)
      type(qed_pdf_t), intent(in) :: epdf
      integer, intent(in) :: flv
      real(default) :: bar
      real(default), intent(in) :: x
      real(default), intent(in) :: scale
      real(default), intent(in) :: alpha
      logical, intent(in) :: running
    end function recbar
    module function rechat &
         (epdf, flv, x, scale, alpha, running) result (hat)
      type(qed_pdf_t), intent(in) :: epdf
      integer, intent(in) :: flv
      real(default) :: hat
      real(default), intent(in) :: x
      real(default), intent(in) :: scale
      real(default), intent(in) :: alpha
      logical, intent(in) :: running
    end function rechat
    module subroutine elecbar_asym_p (x, jll_nll, nlep, ln0, order, running)
      real(default), intent(in) :: x
      real(default), dimension(6), intent(out) :: jll_nll
      real(default), intent(in) :: ln0
      integer, intent(in) :: nlep
      logical, dimension(6), intent(in) :: order
      logical, intent(in) :: running
    end subroutine elecbar_asym_p
    module subroutine photbar_asym_p (x, jll_nll, nlep, ln0, order, running)
      real(default), intent(in) :: x
      real(default), dimension(6), intent(out) :: jll_nll
      real(default), intent(in) :: ln0
      integer, intent(in) :: nlep
      logical, dimension(6), intent(in) :: order
      logical, intent(in) :: running
    end subroutine photbar_asym_p
    module subroutine recbar_singlet (x, jll_nll, nlep, ln0, order, running)
      real(default), intent(in) :: x
      real(default), dimension(6), intent(out) :: jll_nll
      real(default), intent(in) :: ln0
      integer, intent(in) :: nlep
      logical, dimension(6), intent(in) :: order
      logical, intent(in) :: running
    end subroutine recbar_singlet
    module subroutine recbar_nonsinglet (x, jll_nll, nlep, ln0, order, running)
      real(default), intent(in) :: x
      real(default), dimension(6), intent(out) :: jll_nll
      real(default), intent(in) :: ln0
      integer, intent(in) :: nlep
      logical, dimension(6), intent(in) :: order
      logical, intent(in) :: running
    end subroutine recbar_nonsinglet
    module subroutine recbar_photon (x, jll_nll, nlep, ln0, order, running)
      real(default), intent(in) :: x
      real(default), dimension(6), intent(out) :: jll_nll
      real(default), intent(in) :: ln0
      integer, intent(in) :: nlep
      logical, dimension(6), intent(in) :: order
      logical, intent(in) :: running
    end subroutine recbar_photon
    module subroutine rechat_singlet (x, jll_nll, nlep, ln0, order, running)
      real(default), intent(in) :: x
      real(default), dimension(6), intent(out) :: jll_nll
      real(default), intent(in) :: ln0
      integer, intent(in) :: nlep
      logical, dimension(6), intent(in) :: order
      logical, intent(in) :: running
    end subroutine rechat_singlet
    module subroutine rechat_nonsinglet (x, jll_nll, nlep, ln0, order, running)
      real(default), intent(in) :: x
      real(default), dimension(6), intent(out) :: jll_nll
      real(default), intent(in) :: ln0
      integer, intent(in) :: nlep
      logical, dimension(6), intent(in) :: order
      logical, intent(in) :: running
    end subroutine rechat_nonsinglet
    module subroutine rechat_photon (x, jll_nll, nlep, ln0, order, running)
      real(default), intent(in) :: x
      real(default), dimension(6), intent(out) :: jll_nll
      real(default), intent(in) :: ln0
      integer, intent(in) :: nlep
      logical, dimension(6), intent(in) :: order
      logical, intent(in) :: running
    end subroutine rechat_photon
    module function t_alpha (epdf, scale) result (t)
      real(default) :: t    
      type(qed_pdf_t), intent(in) :: epdf
      real(default), intent(in) :: scale
    end function t_alpha
    module function endpoint_func_S (x, nlep) result (f)
      real(default), intent(in) :: x
      integer, intent(in) :: nlep
      real(default) :: f
    end function endpoint_func_S
    module function endpoint_func_NS (x) result (f)
      real(default), intent(in) :: x
      real(default) :: f
    end function endpoint_func_NS
    module function endpoint_func_GAM (x) result (f)
      real(default), intent(in) :: x
      real(default) :: f
    end function endpoint_func_GAM
    module function elec_pdf (epdf, flv, x, scale, alpha, &
         running, w_num) result (e_pdf)
      type(qed_pdf_t), intent(in) :: epdf
      integer, intent(in) :: flv
      real(default) :: e_pdf
      real(default), intent(in) :: x
      real(default), intent(in) :: scale
      real(default), intent(in) :: alpha
      logical, intent(in) :: running, w_num
    end function elec_pdf
  end interface

contains

  subroutine qed_pdf_allocate_aqed (qed, order, n_f, n_lep, running)
    class(qed_pdf_t), intent(inout) :: qed
    integer, intent(in) :: order, n_f, n_lep
    logical, intent(in) :: running
    if (running) then
       allocate (alpha_qed_from_scale_t :: qed%aqed)
       select type (aqed => qed%aqed)
       type is (alpha_qed_from_scale_t)
          aqed%order = order
          aqed%nf = n_f
          aqed%nlep = n_lep
       end select
    else
       allocate (alpha_qed_fixed_t :: qed%aqed)
    end if
  end subroutine qed_pdf_allocate_aqed


end module electron_pdfs
