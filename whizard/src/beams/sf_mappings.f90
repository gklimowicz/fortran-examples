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

module sf_mappings

  use kinds, only: default
  use kinds, only: double

  implicit none
  private

  public :: sf_mapping_t
  public :: sf_s_mapping_t
  public :: sf_res_mapping_t
  public :: sf_res_mapping_single_t
  public :: sf_os_mapping_t
  public :: sf_os_mapping_single_t
  public :: sf_ep_mapping_t
  public :: sf_epr_mapping_t
  public :: sf_epo_mapping_t
  public :: sf_ip_mapping_t
  public :: sf_ipr_mapping_t
  public :: sf_ipo_mapping_t
  public :: sf_ei_mapping_t
  public :: sf_eir_mapping_t
  public :: sf_eio_mapping_t
  public :: map_on_shell
  public :: map_on_shell_inverse
  public :: map_on_shell_single
  public :: map_on_shell_single_inverse
  public :: map_power_1
  public :: map_power_inverse_1
  public :: sf_channel_t
  public :: allocate_sf_channels
  public :: any_sf_channel_has_mapping

  integer, parameter :: SFMAP_NONE = 0
  integer, parameter :: SFMAP_SINGLE = 1
  integer, parameter :: SFMAP_MULTI_S = 2
  integer, parameter :: SFMAP_MULTI_RES = 3
  integer, parameter :: SFMAP_MULTI_ONS = 4
  integer, parameter :: SFMAP_MULTI_EP = 5
  integer, parameter :: SFMAP_MULTI_EPR = 6
  integer, parameter :: SFMAP_MULTI_EPO = 7
  integer, parameter :: SFMAP_MULTI_IP = 8
  integer, parameter :: SFMAP_MULTI_IPR = 9
  integer, parameter :: SFMAP_MULTI_IPO = 10
  integer, parameter :: SFMAP_MULTI_EI = 11
  integer, parameter :: SFMAP_MULTI_SRS = 13
  integer, parameter :: SFMAP_MULTI_SON = 14


  type, abstract :: sf_mapping_t
     integer, dimension(:), allocatable :: i
   contains
     procedure (sf_mapping_write), deferred :: write
     procedure :: base_init => sf_mapping_base_init
     procedure :: set_index => sf_mapping_set_index
     procedure :: get_index => sf_mapping_get_index
     procedure :: get_n_dim => sf_mapping_get_n_dim
     procedure (sf_mapping_compute), deferred :: compute
     procedure (sf_mapping_inverse), deferred :: inverse
     procedure :: check => sf_mapping_check
     procedure :: integral => sf_mapping_integral
  end type sf_mapping_t

  type, extends (sf_mapping_t) :: sf_s_mapping_t
     logical :: power_set = .false.
     real(default) :: power = 1
   contains
     procedure :: write => sf_s_mapping_write
     procedure :: init => sf_s_mapping_init
     procedure :: compute => sf_s_mapping_compute
     procedure :: inverse => sf_s_mapping_inverse
  end type sf_s_mapping_t

  type, extends (sf_mapping_t) :: sf_res_mapping_t
     real(default) :: m = 0
     real(default) :: w = 0
   contains
     procedure :: write => sf_res_mapping_write
     procedure :: init => sf_res_mapping_init
     procedure :: compute => sf_res_mapping_compute
     procedure :: inverse => sf_res_mapping_inverse
  end type sf_res_mapping_t

  type, extends (sf_mapping_t) :: sf_res_mapping_single_t
     real(default) :: m = 0
     real(default) :: w = 0
   contains
     procedure :: write => sf_res_mapping_single_write
     procedure :: init => sf_res_mapping_single_init
     procedure :: compute => sf_res_mapping_single_compute
     procedure :: inverse => sf_res_mapping_single_inverse
  end type sf_res_mapping_single_t

  type, extends (sf_mapping_t) :: sf_os_mapping_t
     real(default) :: m = 0
     real(default) :: lm2 = 0
   contains
     procedure :: write => sf_os_mapping_write
     procedure :: init => sf_os_mapping_init
     procedure :: compute => sf_os_mapping_compute
     procedure :: inverse => sf_os_mapping_inverse
  end type sf_os_mapping_t

  type, extends (sf_mapping_t) :: sf_os_mapping_single_t
     real(default) :: m = 0
     real(default) :: lm2 = 0
   contains
     procedure :: write => sf_os_mapping_single_write
     procedure :: init => sf_os_mapping_single_init
     procedure :: compute => sf_os_mapping_single_compute
     procedure :: inverse => sf_os_mapping_single_inverse
  end type sf_os_mapping_single_t

  type, extends (sf_mapping_t) :: sf_ep_mapping_t
     real(default) :: a = 1
   contains
     procedure :: write => sf_ep_mapping_write
     procedure :: init => sf_ep_mapping_init
     procedure :: compute => sf_ep_mapping_compute
     procedure :: inverse => sf_ep_mapping_inverse
  end type sf_ep_mapping_t

  type, extends (sf_mapping_t) :: sf_epr_mapping_t
     real(default) :: a = 1
     real(default) :: m = 0
     real(default) :: w = 0
     logical :: resonance = .true.
   contains
     procedure :: write => sf_epr_mapping_write
     procedure :: init => sf_epr_mapping_init
     procedure :: compute => sf_epr_mapping_compute
     procedure :: inverse => sf_epr_mapping_inverse
  end type sf_epr_mapping_t

  type, extends (sf_mapping_t) :: sf_epo_mapping_t
     real(default) :: a = 1
     real(default) :: m = 0
     real(default) :: lm2 = 0
   contains
     procedure :: write => sf_epo_mapping_write
     procedure :: init => sf_epo_mapping_init
     procedure :: compute => sf_epo_mapping_compute
     procedure :: inverse => sf_epo_mapping_inverse
  end type sf_epo_mapping_t

  type, extends (sf_mapping_t) :: sf_ip_mapping_t
     real(default) :: eps = 0
   contains
     procedure :: write => sf_ip_mapping_write
     procedure :: init => sf_ip_mapping_init
     procedure :: compute => sf_ip_mapping_compute
     procedure :: inverse => sf_ip_mapping_inverse
  end type sf_ip_mapping_t

  type, extends (sf_mapping_t) :: sf_ipr_mapping_t
     real(default) :: eps = 0
     real(default) :: m = 0
     real(default) :: w = 0
     logical :: resonance = .true.
   contains
     procedure :: write => sf_ipr_mapping_write
     procedure :: init => sf_ipr_mapping_init
     procedure :: compute => sf_ipr_mapping_compute
     procedure :: inverse => sf_ipr_mapping_inverse
  end type sf_ipr_mapping_t

  type, extends (sf_mapping_t) :: sf_ipo_mapping_t
     real(default) :: eps = 0
     real(default) :: m = 0
   contains
     procedure :: write => sf_ipo_mapping_write
     procedure :: init => sf_ipo_mapping_init
     procedure :: compute => sf_ipo_mapping_compute
     procedure :: inverse => sf_ipo_mapping_inverse
  end type sf_ipo_mapping_t

  type, extends (sf_mapping_t) :: sf_ei_mapping_t
     type(sf_ep_mapping_t) :: ep
     type(sf_ip_mapping_t) :: ip
   contains
     procedure :: write => sf_ei_mapping_write
     procedure :: init => sf_ei_mapping_init
     procedure :: set_index => sf_ei_mapping_set_index
     procedure :: compute => sf_ei_mapping_compute
     procedure :: inverse => sf_ei_mapping_inverse
  end type sf_ei_mapping_t

  type, extends (sf_mapping_t) :: sf_eir_mapping_t
     type(sf_res_mapping_t) :: res
     type(sf_epr_mapping_t) :: ep
     type(sf_ipr_mapping_t) :: ip
   contains
     procedure :: write => sf_eir_mapping_write
     procedure :: init => sf_eir_mapping_init
     procedure :: set_index => sf_eir_mapping_set_index
     procedure :: compute => sf_eir_mapping_compute
     procedure :: inverse => sf_eir_mapping_inverse
  end type sf_eir_mapping_t

  type, extends (sf_mapping_t) :: sf_eio_mapping_t
     type(sf_os_mapping_t) :: os
     type(sf_epr_mapping_t) :: ep
     type(sf_ipr_mapping_t) :: ip
   contains
     procedure :: write => sf_eio_mapping_write
     procedure :: init => sf_eio_mapping_init
     procedure :: set_index => sf_eio_mapping_set_index
     procedure :: compute => sf_eio_mapping_compute
     procedure :: inverse => sf_eio_mapping_inverse
  end type sf_eio_mapping_t

  type :: sf_channel_t
     integer, dimension(:), allocatable :: map_code
     class(sf_mapping_t), allocatable :: multi_mapping
   contains
     procedure :: write => sf_channel_write
     procedure :: init => sf_channel_init
     generic :: assignment (=) => sf_channel_assign
     procedure :: sf_channel_assign
     procedure :: activate_mapping => sf_channel_activate_mapping
     procedure :: set_s_mapping => sf_channel_set_s_mapping
     procedure :: set_res_mapping => sf_channel_set_res_mapping
     procedure :: set_os_mapping => sf_channel_set_os_mapping
     procedure :: set_ep_mapping => sf_channel_set_ep_mapping
     procedure :: set_epr_mapping => sf_channel_set_epr_mapping
     procedure :: set_epo_mapping => sf_channel_set_epo_mapping
     procedure :: set_ip_mapping => sf_channel_set_ip_mapping
     procedure :: set_ipr_mapping => sf_channel_set_ipr_mapping
     procedure :: set_ipo_mapping => sf_channel_set_ipo_mapping
     procedure :: set_ei_mapping => sf_channel_set_ei_mapping
     procedure :: set_eir_mapping => sf_channel_set_eir_mapping
     procedure :: set_eio_mapping => sf_channel_set_eio_mapping
     procedure :: is_single_mapping => sf_channel_is_single_mapping
     procedure :: is_multi_mapping => sf_channel_is_multi_mapping
     procedure :: get_multi_mapping_n_par => sf_channel_get_multi_mapping_n_par
     procedure :: set_par_index => sf_channel_set_par_index
  end type sf_channel_t


  abstract interface
     subroutine sf_mapping_write (object, unit)
       import
       class(sf_mapping_t), intent(in) :: object
       integer, intent(in), optional :: unit
     end subroutine sf_mapping_write
  end interface

  abstract interface
     subroutine sf_mapping_compute (mapping, r, rb, f, p, pb, x_free)
       import
       class(sf_mapping_t), intent(inout) :: mapping
       real(default), dimension(:), intent(out) :: r, rb
       real(default), intent(out) :: f
       real(default), dimension(:), intent(in) :: p, pb
       real(default), intent(inout), optional :: x_free
     end subroutine sf_mapping_compute
  end interface

  abstract interface
     subroutine sf_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
       import
       class(sf_mapping_t), intent(inout) :: mapping
       real(default), dimension(:), intent(in) :: r, rb
       real(default), intent(out) :: f
       real(default), dimension(:), intent(out) :: p, pb
       real(default), intent(inout), optional :: x_free
     end subroutine sf_mapping_inverse
  end interface


  interface
    module subroutine sf_mapping_base_init (mapping, n_par)
      class(sf_mapping_t), intent(out) :: mapping
      integer, intent(in) :: n_par
    end subroutine sf_mapping_base_init
    module subroutine sf_mapping_set_index (mapping, j, i)
      class(sf_mapping_t), intent(inout) :: mapping
      integer, intent(in) :: j, i
    end subroutine sf_mapping_set_index
    module function sf_mapping_get_index (mapping, j) result (i)
      class(sf_mapping_t), intent(inout) :: mapping
      integer, intent(in) :: j
      integer :: i
    end function sf_mapping_get_index
    module function sf_mapping_get_n_dim (mapping) result (n)
      class(sf_mapping_t), intent(in) :: mapping
      integer :: n
    end function sf_mapping_get_n_dim
    module subroutine sf_mapping_check (mapping, u, p_in, pb_in, fmt_p, fmt_f)
      class(sf_mapping_t), intent(inout) :: mapping
      integer, intent(in) :: u
      real(default), dimension(:), intent(in) :: p_in, pb_in
      character(*), intent(in) :: fmt_p
      character(*), intent(in), optional :: fmt_f
    end subroutine sf_mapping_check
    module function sf_mapping_integral (mapping, n_calls) result (integral)
      class(sf_mapping_t), intent(inout) :: mapping
      integer, intent(in) :: n_calls
      real(default) :: integral
    end function sf_mapping_integral
    module subroutine sf_s_mapping_write (object, unit)
      class(sf_s_mapping_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_s_mapping_write
    module subroutine sf_s_mapping_init (mapping, power)
      class(sf_s_mapping_t), intent(out) :: mapping
      real(default), intent(in), optional :: power
    end subroutine sf_s_mapping_init
    module subroutine sf_s_mapping_compute (mapping, r, rb, f, p, pb, x_free)
      class(sf_s_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_s_mapping_compute
    module subroutine sf_s_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
      class(sf_s_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(in) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_s_mapping_inverse
    module subroutine sf_res_mapping_write (object, unit)
      class(sf_res_mapping_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_res_mapping_write
    module subroutine sf_res_mapping_init (mapping, m, w)
      class(sf_res_mapping_t), intent(out) :: mapping
      real(default), intent(in) :: m, w
    end subroutine sf_res_mapping_init
    module subroutine sf_res_mapping_compute (mapping, r, rb, f, p, pb, x_free)
      class(sf_res_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_res_mapping_compute
    module subroutine sf_res_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
      class(sf_res_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(in) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_res_mapping_inverse
    module subroutine sf_res_mapping_single_write (object, unit)
      class(sf_res_mapping_single_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_res_mapping_single_write
    module subroutine sf_res_mapping_single_init (mapping, m, w)
      class(sf_res_mapping_single_t), intent(out) :: mapping
      real(default), intent(in) :: m, w
    end subroutine sf_res_mapping_single_init
    module subroutine sf_res_mapping_single_compute &
         (mapping, r, rb, f, p, pb, x_free)
      class(sf_res_mapping_single_t), intent(inout) :: mapping
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_res_mapping_single_compute
    module subroutine sf_res_mapping_single_inverse &
         (mapping, r, rb, f, p, pb, x_free)
      class(sf_res_mapping_single_t), intent(inout) :: mapping
      real(default), dimension(:), intent(in) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_res_mapping_single_inverse
    module subroutine sf_os_mapping_write (object, unit)
      class(sf_os_mapping_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_os_mapping_write
    module subroutine sf_os_mapping_init (mapping, m)
      class(sf_os_mapping_t), intent(out) :: mapping
      real(default), intent(in) :: m
    end subroutine sf_os_mapping_init
    module subroutine sf_os_mapping_compute (mapping, r, rb, f, p, pb, x_free)
      class(sf_os_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_os_mapping_compute
    module subroutine sf_os_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
      class(sf_os_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(in) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_os_mapping_inverse
    module subroutine sf_os_mapping_single_write (object, unit)
      class(sf_os_mapping_single_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_os_mapping_single_write
    module subroutine sf_os_mapping_single_init (mapping, m)
      class(sf_os_mapping_single_t), intent(out) :: mapping
      real(default), intent(in) :: m
    end subroutine sf_os_mapping_single_init
    module subroutine sf_os_mapping_single_compute &
         (mapping, r, rb, f, p, pb, x_free)
      class(sf_os_mapping_single_t), intent(inout) :: mapping
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_os_mapping_single_compute
    module subroutine sf_os_mapping_single_inverse &
         (mapping, r, rb, f, p, pb, x_free)
      class(sf_os_mapping_single_t), intent(inout) :: mapping
      real(default), dimension(:), intent(in) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_os_mapping_single_inverse
    module subroutine sf_ep_mapping_write (object, unit)
      class(sf_ep_mapping_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_ep_mapping_write
    module subroutine sf_ep_mapping_init (mapping, a)
      class(sf_ep_mapping_t), intent(out) :: mapping
      real(default), intent(in), optional :: a
    end subroutine sf_ep_mapping_init
    module subroutine sf_ep_mapping_compute (mapping, r, rb, f, p, pb, x_free)
      class(sf_ep_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_ep_mapping_compute
    module subroutine sf_ep_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
      class(sf_ep_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(in) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_ep_mapping_inverse
    module subroutine sf_epr_mapping_write (object, unit)
      class(sf_epr_mapping_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_epr_mapping_write
    module subroutine sf_epr_mapping_init (mapping, a, m, w)
      class(sf_epr_mapping_t), intent(out) :: mapping
      real(default), intent(in) :: a
      real(default), intent(in), optional :: m, w
    end subroutine sf_epr_mapping_init
    module subroutine sf_epr_mapping_compute (mapping, r, rb, f, p, pb, x_free)
      class(sf_epr_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_epr_mapping_compute
    module subroutine sf_epr_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
      class(sf_epr_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(in) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_epr_mapping_inverse
    module subroutine sf_epo_mapping_write (object, unit)
      class(sf_epo_mapping_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_epo_mapping_write
    module subroutine sf_epo_mapping_init (mapping, a, m)
      class(sf_epo_mapping_t), intent(out) :: mapping
      real(default), intent(in) :: a, m
    end subroutine sf_epo_mapping_init
    module subroutine sf_epo_mapping_compute (mapping, r, rb, f, p, pb, x_free)
      class(sf_epo_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_epo_mapping_compute
    module subroutine sf_epo_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
      class(sf_epo_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(in) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: p, pb
      real(default), intent(inout), optional :: x_free
    end  subroutine sf_epo_mapping_inverse
    module subroutine sf_ip_mapping_write (object, unit)
      class(sf_ip_mapping_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_ip_mapping_write
    module subroutine sf_ip_mapping_init (mapping, eps)
      class(sf_ip_mapping_t), intent(out) :: mapping
      real(default), intent(in), optional :: eps
    end subroutine sf_ip_mapping_init
    module subroutine sf_ip_mapping_compute (mapping, r, rb, f, p, pb, x_free)
      class(sf_ip_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_ip_mapping_compute
    module subroutine sf_ip_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
      class(sf_ip_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(in) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_ip_mapping_inverse
    module subroutine sf_ipr_mapping_write (object, unit)
      class(sf_ipr_mapping_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_ipr_mapping_write
    module subroutine sf_ipr_mapping_init (mapping, eps, m, w)
      class(sf_ipr_mapping_t), intent(out) :: mapping
      real(default), intent(in), optional :: eps, m, w
    end subroutine sf_ipr_mapping_init
    module subroutine sf_ipr_mapping_compute (mapping, r, rb, f, p, pb, x_free)
      class(sf_ipr_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_ipr_mapping_compute
    module subroutine sf_ipr_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
      class(sf_ipr_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(in) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_ipr_mapping_inverse
    module subroutine sf_ipo_mapping_write (object, unit)
      class(sf_ipo_mapping_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_ipo_mapping_write
    module subroutine sf_ipo_mapping_init (mapping, eps, m)
      class(sf_ipo_mapping_t), intent(out) :: mapping
      real(default), intent(in), optional :: eps, m
    end subroutine sf_ipo_mapping_init
    module subroutine sf_ipo_mapping_compute (mapping, r, rb, f, p, pb, x_free)
      class(sf_ipo_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_ipo_mapping_compute
    module subroutine sf_ipo_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
      class(sf_ipo_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(in) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_ipo_mapping_inverse
    module subroutine sf_ei_mapping_write (object, unit)
      class(sf_ei_mapping_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_ei_mapping_write
    module subroutine sf_ei_mapping_init (mapping, a, eps)
      class(sf_ei_mapping_t), intent(out) :: mapping
      real(default), intent(in), optional :: a, eps
    end subroutine sf_ei_mapping_init
    module subroutine sf_ei_mapping_set_index (mapping, j, i)
      class(sf_ei_mapping_t), intent(inout) :: mapping
      integer, intent(in) :: j, i
    end subroutine sf_ei_mapping_set_index
    module subroutine sf_ei_mapping_compute (mapping, r, rb, f, p, pb, x_free)
      class(sf_ei_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_ei_mapping_compute
    module subroutine sf_ei_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
      class(sf_ei_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(in) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_ei_mapping_inverse
    module subroutine sf_eir_mapping_write (object, unit)
      class(sf_eir_mapping_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_eir_mapping_write
    module subroutine sf_eir_mapping_init (mapping, a, eps, m, w)
      class(sf_eir_mapping_t), intent(out) :: mapping
      real(default), intent(in) :: a, eps, m, w
    end subroutine sf_eir_mapping_init
    module subroutine sf_eir_mapping_set_index (mapping, j, i)
      class(sf_eir_mapping_t), intent(inout) :: mapping
      integer, intent(in) :: j, i
    end subroutine sf_eir_mapping_set_index
    module subroutine sf_eir_mapping_compute (mapping, r, rb, f, p, pb, x_free)
      class(sf_eir_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_eir_mapping_compute
    module subroutine sf_eir_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
      class(sf_eir_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(in) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_eir_mapping_inverse
    module subroutine sf_eio_mapping_write (object, unit)
      class(sf_eio_mapping_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_eio_mapping_write
    module subroutine sf_eio_mapping_init (mapping, a, eps, m)
      class(sf_eio_mapping_t), intent(out) :: mapping
      real(default), intent(in), optional :: a, eps, m
    end subroutine sf_eio_mapping_init
    module subroutine sf_eio_mapping_set_index (mapping, j, i)
      class(sf_eio_mapping_t), intent(inout) :: mapping
      integer, intent(in) :: j, i
    end subroutine sf_eio_mapping_set_index
    module subroutine sf_eio_mapping_compute (mapping, r, rb, f, p, pb, x_free)
      class(sf_eio_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(out) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(in) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_eio_mapping_compute
    module subroutine sf_eio_mapping_inverse (mapping, r, rb, f, p, pb, x_free)
      class(sf_eio_mapping_t), intent(inout) :: mapping
      real(default), dimension(:), intent(in) :: r, rb
      real(default), intent(out) :: f
      real(default), dimension(:), intent(out) :: p, pb
      real(default), intent(inout), optional :: x_free
    end subroutine sf_eio_mapping_inverse
    module subroutine map_on_shell (r, factor, p, lm2, x_free)
      real(default), dimension(2), intent(out) :: r
      real(default), intent(out) :: factor
      real(default), dimension(2), intent(in) :: p
      real(default), intent(in) :: lm2
      real(default), intent(in), optional :: x_free
    end subroutine map_on_shell
    module subroutine map_on_shell_inverse (r, factor, p, lm2, x_free)
      real(default), dimension(2), intent(in) :: r
      real(default), intent(out) :: factor
      real(default), dimension(2), intent(out) :: p
      real(default), intent(in) :: lm2
      real(default), intent(in), optional :: x_free
    end subroutine map_on_shell_inverse
    module subroutine map_on_shell_single (r, factor, p, lm2, x_free)
      real(default), dimension(1), intent(out) :: r
      real(default), intent(out) :: factor
      real(default), dimension(1), intent(in) :: p
      real(default), intent(in) :: lm2
      real(default), intent(in), optional :: x_free
    end subroutine map_on_shell_single
    module subroutine map_on_shell_single_inverse (r, factor, p, lm2, x_free)
      real(default), dimension(1), intent(in) :: r
      real(default), intent(out) :: factor
      real(default), dimension(1), intent(out) :: p
      real(default), intent(in) :: lm2
      real(default), intent(in), optional :: x_free
    end subroutine map_on_shell_single_inverse
    module subroutine map_power_1 (xb, factor, rb, eps)
      real(default), intent(out) :: xb, factor
      real(default), intent(in) :: rb
      real(default), intent(in) :: eps
    end subroutine map_power_1
    module subroutine map_power_inverse_1 (xb, factor, rb, eps)
      real(default), intent(in) :: xb
      real(default), intent(out) :: rb, factor
      real(default), intent(in) :: eps
    end subroutine map_power_inverse_1
    module subroutine sf_channel_write (object, unit)
      class(sf_channel_t), intent(in) :: object
      integer, intent(in), optional :: unit
    end subroutine sf_channel_write
    module subroutine sf_channel_init (channel, n_strfun)
      class(sf_channel_t), intent(out) :: channel
      integer, intent(in) :: n_strfun
    end subroutine sf_channel_init
    module subroutine sf_channel_assign (copy, original)
      class(sf_channel_t), intent(out) :: copy
      type(sf_channel_t), intent(in) :: original
    end subroutine sf_channel_assign
    module subroutine allocate_sf_channels (channel, n_channel, n_strfun)
      type(sf_channel_t), dimension(:), intent(out), allocatable :: channel
      integer, intent(in) :: n_channel
      integer, intent(in) :: n_strfun
    end subroutine allocate_sf_channels
    module subroutine sf_channel_activate_mapping (channel, i_sf)
      class(sf_channel_t), intent(inout) :: channel
      integer, dimension(:), intent(in) :: i_sf
    end subroutine sf_channel_activate_mapping
    module function sf_channel_is_single_mapping (channel, i_sf) result (flag)
      class(sf_channel_t), intent(in) :: channel
      integer, intent(in) :: i_sf
      logical :: flag
    end function sf_channel_is_single_mapping
    module function sf_channel_is_multi_mapping (channel, i_sf) result (flag)
      class(sf_channel_t), intent(in) :: channel
      integer, intent(in) :: i_sf
      logical :: flag
    end function sf_channel_is_multi_mapping
    module function sf_channel_get_multi_mapping_n_par (channel) result (n_par)
      class(sf_channel_t), intent(in) :: channel
      integer :: n_par
    end function sf_channel_get_multi_mapping_n_par
    module function any_sf_channel_has_mapping (channel) result (flag)
      type(sf_channel_t), dimension(:), intent(in) :: channel
      logical :: flag
    end function any_sf_channel_has_mapping
    module subroutine sf_channel_set_par_index (channel, j, i_par)
      class(sf_channel_t), intent(inout) :: channel
      integer, intent(in) :: j
      integer, intent(in) :: i_par
    end subroutine sf_channel_set_par_index
  end interface

contains

  subroutine sf_channel_set_s_mapping (channel, i_sf, power)
    class(sf_channel_t), intent(inout) :: channel
    integer, dimension(:), intent(in) :: i_sf
    real(default), intent(in), optional :: power
    channel%map_code(i_sf) = SFMAP_MULTI_S
    allocate (sf_s_mapping_t :: channel%multi_mapping)
    select type (mapping => channel%multi_mapping)
    type is (sf_s_mapping_t)
       call mapping%init (power)
    end select
  end subroutine sf_channel_set_s_mapping

  subroutine sf_channel_set_res_mapping (channel, i_sf, m, w, single)
    class(sf_channel_t), intent(inout) :: channel
    integer, dimension(:), intent(in) :: i_sf
    real(default), intent(in) :: m, w
    logical, intent(in) :: single
    if (single) then
       channel%map_code(i_sf) = SFMAP_MULTI_SRS
       allocate (sf_res_mapping_single_t :: channel%multi_mapping)
       select type (mapping => channel%multi_mapping)
       type is (sf_res_mapping_single_t)
          call mapping%init (m, w)
       end select
    else
       channel%map_code(i_sf) = SFMAP_MULTI_RES
       allocate (sf_res_mapping_t :: channel%multi_mapping)
       select type (mapping => channel%multi_mapping)
       type is (sf_res_mapping_t)
          call mapping%init (m, w)
       end select
    end if
  end subroutine sf_channel_set_res_mapping

  subroutine sf_channel_set_os_mapping (channel, i_sf, m, single)
    class(sf_channel_t), intent(inout) :: channel
    integer, dimension(:), intent(in) :: i_sf
    real(default), intent(in) :: m
    logical, intent(in) :: single
    if (single) then
       channel%map_code(i_sf) = SFMAP_MULTI_SON
       allocate (sf_os_mapping_single_t :: channel%multi_mapping)
       select type (mapping => channel%multi_mapping)
       type is (sf_os_mapping_single_t)
          call mapping%init (m)
       end select
    else
       channel%map_code(i_sf) = SFMAP_MULTI_ONS
       allocate (sf_os_mapping_t :: channel%multi_mapping)
       select type (mapping => channel%multi_mapping)
       type is (sf_os_mapping_t)
          call mapping%init (m)
       end select
    end if
  end subroutine sf_channel_set_os_mapping

  subroutine sf_channel_set_ep_mapping (channel, i_sf, a)
    class(sf_channel_t), intent(inout) :: channel
    integer, dimension(:), intent(in) :: i_sf
    real(default), intent(in), optional :: a
    channel%map_code(i_sf) = SFMAP_MULTI_EP
    allocate (sf_ep_mapping_t :: channel%multi_mapping)
    select type (mapping => channel%multi_mapping)
    type is (sf_ep_mapping_t)
       call mapping%init (a = a)
    end select
  end subroutine sf_channel_set_ep_mapping

  subroutine sf_channel_set_epr_mapping (channel, i_sf, a, m, w)
    class(sf_channel_t), intent(inout) :: channel
    integer, dimension(:), intent(in) :: i_sf
    real(default), intent(in) :: a, m, w
    channel%map_code(i_sf) = SFMAP_MULTI_EPR
    allocate (sf_epr_mapping_t :: channel%multi_mapping)
    select type (mapping => channel%multi_mapping)
    type is (sf_epr_mapping_t)
       call mapping%init (a, m, w)
    end select
  end subroutine sf_channel_set_epr_mapping

  subroutine sf_channel_set_epo_mapping (channel, i_sf, a, m)
    class(sf_channel_t), intent(inout) :: channel
    integer, dimension(:), intent(in) :: i_sf
    real(default), intent(in) :: a, m
    channel%map_code(i_sf) = SFMAP_MULTI_EPO
    allocate (sf_epo_mapping_t :: channel%multi_mapping)
    select type (mapping => channel%multi_mapping)
    type is (sf_epo_mapping_t)
       call mapping%init (a, m)
    end select
  end subroutine sf_channel_set_epo_mapping

  subroutine sf_channel_set_ip_mapping (channel, i_sf, eps)
    class(sf_channel_t), intent(inout) :: channel
    integer, dimension(:), intent(in) :: i_sf
    real(default), intent(in), optional :: eps
    channel%map_code(i_sf) = SFMAP_MULTI_IP
    allocate (sf_ip_mapping_t :: channel%multi_mapping)
    select type (mapping => channel%multi_mapping)
    type is (sf_ip_mapping_t)
       call mapping%init (eps)
    end select
  end subroutine sf_channel_set_ip_mapping

  subroutine sf_channel_set_ipr_mapping (channel, i_sf, eps, m, w)
    class(sf_channel_t), intent(inout) :: channel
    integer, dimension(:), intent(in) :: i_sf
    real(default), intent(in), optional :: eps, m, w
    channel%map_code(i_sf) = SFMAP_MULTI_IPR
    allocate (sf_ipr_mapping_t :: channel%multi_mapping)
    select type (mapping => channel%multi_mapping)
    type is (sf_ipr_mapping_t)
       call mapping%init (eps, m, w)
    end select
  end subroutine sf_channel_set_ipr_mapping

  subroutine sf_channel_set_ipo_mapping (channel, i_sf, eps, m)
    class(sf_channel_t), intent(inout) :: channel
    integer, dimension(:), intent(in) :: i_sf
    real(default), intent(in), optional :: eps, m
    channel%map_code(i_sf) = SFMAP_MULTI_IPO
    allocate (sf_ipo_mapping_t :: channel%multi_mapping)
    select type (mapping => channel%multi_mapping)
    type is (sf_ipo_mapping_t)
       call mapping%init (eps, m)
    end select
  end subroutine sf_channel_set_ipo_mapping

  subroutine sf_channel_set_ei_mapping (channel, i_sf, a, eps)
    class(sf_channel_t), intent(inout) :: channel
    integer, dimension(:), intent(in) :: i_sf
    real(default), intent(in), optional :: a, eps
    channel%map_code(i_sf) = SFMAP_MULTI_EI
    allocate (sf_ei_mapping_t :: channel%multi_mapping)
    select type (mapping => channel%multi_mapping)
    type is (sf_ei_mapping_t)
       call mapping%init (a, eps)
    end select
  end subroutine sf_channel_set_ei_mapping

  subroutine sf_channel_set_eir_mapping (channel, i_sf, a, eps, m, w)
    class(sf_channel_t), intent(inout) :: channel
    integer, dimension(:), intent(in) :: i_sf
    real(default), intent(in), optional :: a, eps, m, w
    channel%map_code(i_sf) = SFMAP_MULTI_EI
    allocate (sf_eir_mapping_t :: channel%multi_mapping)
    select type (mapping => channel%multi_mapping)
    type is (sf_eir_mapping_t)
       call mapping%init (a, eps, m, w)
    end select
  end subroutine sf_channel_set_eir_mapping

  subroutine sf_channel_set_eio_mapping (channel, i_sf, a, eps, m)
    class(sf_channel_t), intent(inout) :: channel
    integer, dimension(:), intent(in) :: i_sf
    real(default), intent(in), optional :: a, eps, m
    channel%map_code(i_sf) = SFMAP_MULTI_EI
    allocate (sf_eio_mapping_t :: channel%multi_mapping)
    select type (mapping => channel%multi_mapping)
    type is (sf_eio_mapping_t)
       call mapping%init (a, eps, m)
    end select
  end subroutine sf_channel_set_eio_mapping


end module sf_mappings
