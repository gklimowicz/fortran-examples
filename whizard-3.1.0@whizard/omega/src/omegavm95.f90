!  omegalib.nw --
!
!  Copyright (C) 1999-2022 by
!      Wolfgang Kilian <kilian@physik.uni-siegen.de>
!      Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
!      Juergen Reuter <juergen.reuter@desy.de>
!      with contributions from                                                                                                                                    
!      Fabian Bach <fabian.bach@t-online.de>                                                                                                                 
!      Bijan Chokoufe Nejad <bijan.chokoufe@desy.de>                                                                                                              
!      Christian Speckner <cnspeckn@googlemail.com>     
!
!  WHIZARD is free software; you can redistribute it and/or modify it
!  under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2, or (at your option)
!  any later version.
!
!  WHIZARD is distributed in the hope that it will be useful, but
!  WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module omegavm95
  use kinds, only: default
  use constants
  use iso_varying_string, string_t => varying_string
  use, intrinsic :: iso_fortran_env, only : input_unit, output_unit, error_unit
  use omega95
  use omega95_bispinors, only: bispinor, vectorspinor, veps, pr_grav
  use omega95_bispinors, only: bi_u => u
  use omega95_bispinors, only: bi_v => v
  use omega95_bispinors, only: bi_pr_psi => pr_psi
  use omega_bispinors, only: operator (*), operator (+)
  use omega_color, only: ovm_color_sum, OCF => omega_color_factor
  implicit none
  private
  integer, parameter, public :: stdin = input_unit
  integer, parameter, public :: stdout = output_unit
  integer, parameter, public :: stderr = error_unit
  integer, parameter :: MIN_UNIT = 11, MAX_UNIT = 99
  public :: color_t
  type color_t
  contains
    procedure :: write => color_write
  end type color_t

  public :: col_discrete
  type, extends(color_t) :: col_discrete
    integer :: i
  end type col_discrete

  public :: flavor_t
  type flavor_t
  contains
    procedure :: write => flavor_write
  end type flavor_t

  public :: flv_discrete
  type, extends(flavor_t) :: flv_discrete
    integer :: i
  end type flv_discrete

  public :: helicity_t
  type :: helicity_t
  contains
    procedure :: write => helicity_write
  end type helicity_t

  public :: hel_discrete
  type, extends(helicity_t) :: hel_discrete
    integer :: i
  end type hel_discrete

  public :: hel_trigonometric
  type, extends(helicity_t) :: hel_trigonometric
    real :: theta
  end type hel_trigonometric

  public :: hel_exponential
  type, extends(helicity_t) :: hel_exponential
    real :: phi
  end type hel_exponential

  public :: hel_spherical
  type, extends(helicity_t) :: hel_spherical
    real :: theta, phi
  end type hel_spherical

  integer, parameter :: len_instructions = 8
  integer, parameter :: N_version_lines = 2
  ! Comment lines including the first header description line
  integer, parameter :: N_comments = 6
  ! Actual data lines plus intermediate description lines
  ! 'description \n 1 2 3 \n description \n 3 2 1' would count as 3
  integer, parameter :: N_header_lines = 5
  real(default), parameter, public :: N_ = three

  type :: basic_vm_t
     private
     logical :: verbose
     type(string_t) :: bytecode_file
     integer :: bytecode_fh, out_fh
     integer :: N_instructions, N_levels
     integer :: N_table_lines
     integer, dimension(:, :), allocatable :: instructions
     integer, dimension(:), allocatable :: levels
  end type

  type :: vm_scalar
     logical :: c
     complex(kind=default) :: v
  end type

  type :: vm_spinor
     logical :: c
     type(spinor) :: v
  end type

  type :: vm_conjspinor
     logical :: c
     type(conjspinor) :: v
  end type

  type :: vm_bispinor
     logical :: c
     type(bispinor) :: v
  end type

  type :: vm_vector
     logical :: c
     type(vector) :: v
  end type

  type :: vm_tensor_2
     logical :: c
     type(tensor) :: v
  end type

  type :: vm_tensor_1
     logical :: c
     type(tensor2odd) :: v
  end type

  type :: vm_vectorspinor
     logical :: c
     type(vectorspinor) :: v
  end type

  type, public, extends (basic_vm_t) :: vm_t
     private
     type(string_t) :: version
     type(string_t) :: model
     integer :: N_momenta, N_particles, N_prt_in, N_prt_out, N_amplitudes
     ! helicities = helicity combinations
     integer :: N_helicities, N_col_flows, N_col_indices, N_flavors, N_col_factors

     integer :: N_scalars, N_spinors, N_conjspinors, N_bispinors
     integer :: N_vectors, N_tensors_2, N_tensors_1, N_vectorspinors

     integer :: N_coupl_real, N_coupl_real2, N_coupl_cmplx, N_coupl_cmplx2

     integer, dimension(:, :), allocatable :: table_flavor
     integer, dimension(:, :, :), allocatable :: table_color_flows
     integer, dimension(:, :), allocatable :: table_spin
     logical, dimension(:, :), allocatable :: table_ghost_flags
     type(OCF), dimension(:), allocatable :: table_color_factors
     logical, dimension(:, :), allocatable :: table_flv_col_is_allowed

     real(default), dimension(:), allocatable :: coupl_real
     real(default), dimension(:, :), allocatable :: coupl_real2
     complex(default), dimension(:), allocatable :: coupl_cmplx
     complex(default), dimension(:, :), allocatable :: coupl_cmplx2
     real(default), dimension(:), allocatable :: mass
     real(default), dimension(:), allocatable :: width

     type(momentum), dimension(:), allocatable :: momenta
     complex(default), dimension(:), allocatable :: amplitudes
     complex(default), dimension(:, :, :), allocatable :: table_amplitudes
     class(flavor_t), dimension(:), allocatable :: flavor
     class(color_t), dimension(:), allocatable :: color
     ! gfortran 4.7
     !class(helicity_t), dimension(:), pointer :: helicity => null()
     integer, dimension(:), allocatable :: helicity

     type(vm_scalar), dimension(:), allocatable :: scalars
     type(vm_spinor), dimension(:), allocatable :: spinors
     type(vm_conjspinor), dimension(:), allocatable :: conjspinors
     type(vm_bispinor), dimension(:), allocatable :: bispinors
     type(vm_vector), dimension(:), allocatable :: vectors
     type(vm_tensor_2), dimension(:), allocatable :: tensors_2
     type(vm_tensor_1), dimension(:), allocatable :: tensors_1
     type(vm_vectorspinor), dimension(:), allocatable :: vectorspinors

     logical, dimension(:), allocatable :: hel_is_allowed
     real(default), dimension(:), allocatable :: hel_max_abs
     real(default) :: hel_sum_abs = 0, hel_threshold = 1E10
     integer :: hel_count = 0, hel_cutoff = 100
     integer, dimension(:), allocatable :: hel_map
     integer :: hel_finite
     logical :: cms

     logical :: openmp

  contains
     procedure :: init => vm_init
     procedure :: write => vm_write
     procedure :: reset => vm_reset
     procedure :: run => vm_run
     procedure :: final => vm_final
     procedure :: number_particles_in => vm_number_particles_in
     procedure :: number_particles_out => vm_number_particles_out
     procedure :: number_color_indices => vm_number_color_indices
     procedure :: reset_helicity_selection => vm_reset_helicity_selection
     procedure :: new_event => vm_new_event
     procedure :: color_sum => vm_color_sum
     procedure :: spin_states => vm_spin_states
     procedure :: number_spin_states => vm_number_spin_states
     procedure :: number_color_flows => vm_number_color_flows
     procedure :: flavor_states => vm_flavor_states
     procedure :: number_flavor_states => vm_number_flavor_states
     procedure :: color_flows => vm_color_flows
     procedure :: color_factors => vm_color_factors
     procedure :: number_color_factors => vm_number_color_factors
     procedure :: is_allowed => vm_is_allowed
     procedure :: get_amplitude => vm_get_amplitude
  end type

  integer, parameter :: ovm_ADD_MOMENTA = 1

  integer, parameter :: ovm_LOAD_SCALAR = 10
  integer, parameter :: ovm_LOAD_SPINOR_INC = 11
  integer, parameter :: ovm_LOAD_SPINOR_OUT = 12
  integer, parameter :: ovm_LOAD_CONJSPINOR_INC = 13
  integer, parameter :: ovm_LOAD_CONJSPINOR_OUT = 14
  integer, parameter :: ovm_LOAD_MAJORANA_INC = 15
  integer, parameter :: ovm_LOAD_MAJORANA_OUT = 16
  integer, parameter :: ovm_LOAD_VECTOR_INC = 17
  integer, parameter :: ovm_LOAD_VECTOR_OUT = 18
  integer, parameter :: ovm_LOAD_VECTORSPINOR_INC = 19
  integer, parameter :: ovm_LOAD_VECTORSPINOR_OUT = 20
  integer, parameter :: ovm_LOAD_TENSOR2_INC = 21
  integer, parameter :: ovm_LOAD_TENSOR2_OUT = 22
  integer, parameter :: ovm_LOAD_BRS_SCALAR = 30
  integer, parameter :: ovm_LOAD_BRS_SPINOR_INC = 31
  integer, parameter :: ovm_LOAD_BRS_SPINOR_OUT = 32
  integer, parameter :: ovm_LOAD_BRS_CONJSPINOR_INC = 33
  integer, parameter :: ovm_LOAD_BRS_CONJSPINOR_OUT = 34
  integer, parameter :: ovm_LOAD_BRS_VECTOR_INC = 37
  integer, parameter :: ovm_LOAD_BRS_VECTOR_OUT = 38
  integer, parameter :: ovm_LOAD_MAJORANA_GHOST_INC = 23
  integer, parameter :: ovm_LOAD_MAJORANA_GHOST_OUT = 24
  integer, parameter :: ovm_LOAD_BRS_MAJORANA_INC = 35
  integer, parameter :: ovm_LOAD_BRS_MAJORANA_OUT = 36

  integer, parameter :: ovm_CALC_BRAKET = 2

  integer, parameter :: ovm_FUSE_V_FF = -1
  integer, parameter :: ovm_FUSE_F_VF = -2
  integer, parameter :: ovm_FUSE_F_FV = -3
  integer, parameter :: ovm_FUSE_VA_FF = -4
  integer, parameter :: ovm_FUSE_F_VAF = -5
  integer, parameter :: ovm_FUSE_F_FVA = -6
  integer, parameter :: ovm_FUSE_VA2_FF = -7
  integer, parameter :: ovm_FUSE_F_VA2F = -8
  integer, parameter :: ovm_FUSE_F_FVA2 = -9
  integer, parameter :: ovm_FUSE_A_FF = -10
  integer, parameter :: ovm_FUSE_F_AF = -11
  integer, parameter :: ovm_FUSE_F_FA = -12
  integer, parameter :: ovm_FUSE_VL_FF = -13
  integer, parameter :: ovm_FUSE_F_VLF = -14
  integer, parameter :: ovm_FUSE_F_FVL = -15
  integer, parameter :: ovm_FUSE_VR_FF = -16
  integer, parameter :: ovm_FUSE_F_VRF = -17
  integer, parameter :: ovm_FUSE_F_FVR = -18
  integer, parameter :: ovm_FUSE_VLR_FF = -19
  integer, parameter :: ovm_FUSE_F_VLRF = -20
  integer, parameter :: ovm_FUSE_F_FVLR = -21
  integer, parameter :: ovm_FUSE_SP_FF = -22
  integer, parameter :: ovm_FUSE_F_SPF = -23
  integer, parameter :: ovm_FUSE_F_FSP = -24
  integer, parameter :: ovm_FUSE_S_FF = -25
  integer, parameter :: ovm_FUSE_F_SF = -26
  integer, parameter :: ovm_FUSE_F_FS = -27
  integer, parameter :: ovm_FUSE_P_FF = -28
  integer, parameter :: ovm_FUSE_F_PF = -29
  integer, parameter :: ovm_FUSE_F_FP = -30
  integer, parameter :: ovm_FUSE_SL_FF = -31
  integer, parameter :: ovm_FUSE_F_SLF = -32
  integer, parameter :: ovm_FUSE_F_FSL = -33
  integer, parameter :: ovm_FUSE_SR_FF = -34
  integer, parameter :: ovm_FUSE_F_SRF = -35
  integer, parameter :: ovm_FUSE_F_FSR = -36
  integer, parameter :: ovm_FUSE_SLR_FF = -37
  integer, parameter :: ovm_FUSE_F_SLRF = -38
  integer, parameter :: ovm_FUSE_F_FSLR = -39

  integer, parameter :: ovm_FUSE_G_GG = -40
  integer, parameter :: ovm_FUSE_V_SS = -41
  integer, parameter :: ovm_FUSE_S_VV = -42
  integer, parameter :: ovm_FUSE_S_VS = -43
  integer, parameter :: ovm_FUSE_V_SV = -44
  integer, parameter :: ovm_FUSE_S_SS = -45
  integer, parameter :: ovm_FUSE_S_SVV = -46
  integer, parameter :: ovm_FUSE_V_SSV = -47
  integer, parameter :: ovm_FUSE_S_SSS = -48
  integer, parameter :: ovm_FUSE_V_VVV = -49

  integer, parameter :: ovm_FUSE_S_G2 = -50
  integer, parameter :: ovm_FUSE_G_SG = -51
  integer, parameter :: ovm_FUSE_G_GS = -52
  integer, parameter :: ovm_FUSE_S_G2_SKEW = -53
  integer, parameter :: ovm_FUSE_G_SG_SKEW = -54
  integer, parameter :: ovm_FUSE_G_GS_SKEW = -55
  integer, parameter :: ovm_PROPAGATE_SCALAR = 51
  integer, parameter :: ovm_PROPAGATE_COL_SCALAR = 52
  integer, parameter :: ovm_PROPAGATE_GHOST = 53
  integer, parameter :: ovm_PROPAGATE_SPINOR = 54
  integer, parameter :: ovm_PROPAGATE_CONJSPINOR = 55
  integer, parameter :: ovm_PROPAGATE_MAJORANA = 56
  integer, parameter :: ovm_PROPAGATE_COL_MAJORANA = 57
  integer, parameter :: ovm_PROPAGATE_UNITARITY = 58
  integer, parameter :: ovm_PROPAGATE_COL_UNITARITY = 59
  integer, parameter :: ovm_PROPAGATE_FEYNMAN = 60
  integer, parameter :: ovm_PROPAGATE_COL_FEYNMAN = 61
  integer, parameter :: ovm_PROPAGATE_VECTORSPINOR = 62
  integer, parameter :: ovm_PROPAGATE_TENSOR2 = 63
  integer, parameter :: ovm_PROPAGATE_NONE = 64
contains
  subroutine find_free_unit (u, iostat)
    integer, intent(out) :: u
    integer, intent(out), optional :: iostat
    logical :: exists, is_open
    integer :: i, status
    do i = MIN_UNIT, MAX_UNIT
       inquire (unit = i, exist = exists, opened = is_open, &
            iostat = status)
       if (status == 0) then
          if (exists .and. .not. is_open) then
             u = i
             if (present (iostat)) then
                iostat = 0
             end if
             return
          end if
       end if
    end do
    if (present (iostat)) then
       iostat = -1
    end if
    u = -1
  end subroutine find_free_unit
  subroutine alloc_arrays (vm)
    type(vm_t), intent(inout) :: vm
    integer :: i
    allocate (vm%table_flavor(vm%N_particles, vm%N_flavors))
    allocate (vm%table_color_flows(vm%N_col_indices, vm%N_particles, &
                                   vm%N_col_flows))
    allocate (vm%table_spin(vm%N_particles, vm%N_helicities))
    allocate (vm%table_ghost_flags(vm%N_particles, vm%N_col_flows))
    allocate (vm%table_color_factors(vm%N_col_factors))
    allocate (vm%table_flv_col_is_allowed(vm%N_flavors, vm%N_col_flows))
    allocate (vm%momenta(vm%N_momenta))
    allocate (vm%amplitudes(vm%N_amplitudes))
    allocate (vm%table_amplitudes(vm%N_flavors, vm%N_col_flows, &
                                  vm%N_helicities))
    vm%table_amplitudes = zero
    allocate (vm%scalars(vm%N_scalars))
    allocate (vm%spinors(vm%N_spinors))
    allocate (vm%conjspinors(vm%N_conjspinors))
    allocate (vm%bispinors(vm%N_bispinors))
    allocate (vm%vectors(vm%N_vectors))
    allocate (vm%tensors_2(vm%N_tensors_2))
    allocate (vm%tensors_1(vm%N_tensors_1))
    allocate (vm%vectorspinors(vm%N_vectorspinors))
    allocate (vm%hel_is_allowed(vm%N_helicities))
    vm%hel_is_allowed = .True.
    allocate (vm%hel_max_abs(vm%N_helicities))
    vm%hel_max_abs = 0
    allocate (vm%hel_map(vm%N_helicities))
    vm%hel_map = (/(i, i = 1, vm%N_helicities)/)
    vm%hel_finite = vm%N_helicities
  end subroutine alloc_arrays

  subroutine vm_init (vm, bytecode_file, version, model, &
      coupl_real, coupl_real2, coupl_cmplx, coupl_cmplx2, &
      mass, width, verbose, out_fh, openmp)
    class(vm_t), intent(out) :: vm
    type(string_t), intent(in) :: bytecode_file
    type(string_t), intent(in) :: version
    type(string_t), intent(in) :: model
    real(default), dimension(:), optional, intent(in) :: coupl_real
    real(default), dimension(:, :), optional, intent(in) :: coupl_real2
    complex(default), dimension(:), optional, intent(in) :: coupl_cmplx
    complex(default), dimension(:, :), optional, intent(in) :: coupl_cmplx2
    real(default), dimension(:), optional, intent(in) :: mass
    real(default), dimension(:), optional, intent(in) :: width
    logical, optional, intent(in) :: verbose
    integer, optional, intent(in) :: out_fh
    logical, optional, intent(in) :: openmp
    vm%bytecode_file = bytecode_file
    vm%version = version
    vm%model = model
    if (present (coupl_real)) then
       allocate (vm%coupl_real (size (coupl_real)), source=coupl_real)
    end if
    if (present (coupl_real2)) then
       allocate (vm%coupl_real2 (2, size (coupl_real2, 2)), source=coupl_real2)
    end if
    if (present (coupl_cmplx)) then
       allocate (vm%coupl_cmplx (size (coupl_cmplx)), source=coupl_cmplx)
    end if
    if (present (coupl_cmplx2)) then
       allocate (vm%coupl_cmplx2 (2, size (coupl_cmplx2, 2)), &
                 source=coupl_cmplx2)
    end if
    if (present (mass)) then
       allocate (vm%mass(size(mass)), source=mass)
    end if
    if (present (width)) then
       allocate (vm%width(size (width)), source=width)
    end if
    if (present (openmp)) then
       vm%openmp = openmp
    else
       vm%openmp = .false.
    end if
    vm%cms = .false.

    call basic_init (vm, verbose, out_fh)
  end subroutine vm_init

  subroutine vm_reset (vm, &
      coupl_real, coupl_real2, coupl_cmplx, coupl_cmplx2, &
      mass, width, verbose, out_fh)
    class(vm_t), intent(inout) :: vm
    real(default), dimension(:), optional, intent(in) :: coupl_real
    real(default), dimension(:, :), optional, intent(in) :: coupl_real2
    complex(default), dimension(:), optional, intent(in) :: coupl_cmplx
    complex(default), dimension(:, :), optional, intent(in) :: coupl_cmplx2
    real(default), dimension(:), optional, intent(in) :: mass
    real(default), dimension(:), optional, intent(in) :: width
    logical, optional, intent(in) :: verbose
    integer, optional, intent(in) :: out_fh
    if (present (coupl_real)) then
       vm%coupl_real = coupl_real
    end if
    if (present (coupl_real2)) then
       vm%coupl_real2 = coupl_real2
    end if
    if (present (coupl_cmplx)) then
       vm%coupl_cmplx = coupl_cmplx
    end if
    if (present (coupl_cmplx2)) then
       vm%coupl_cmplx2 = coupl_cmplx2
    end if
    if (present (mass)) then
       vm%mass = mass
    end if
    if (present (width)) then
       vm%width = width
    end if
    if (present (verbose)) then
       vm%verbose = verbose
    end if
    if (present (out_fh)) then
       vm%out_fh = out_fh
    end if
  end subroutine vm_reset

  subroutine vm_write (vm)
    class(vm_t), intent(in) :: vm
    integer :: i, j, k
    call basic_write (vm)
    write(vm%out_fh, *) 'table_flavor              = ', vm%table_flavor
    write(vm%out_fh, *) 'table_color_flows         = ', vm%table_color_flows
    write(vm%out_fh, *) 'table_spin                = ', vm%table_spin
    write(vm%out_fh, *) 'table_ghost_flags         = ', vm%table_ghost_flags
    write(vm%out_fh, *) 'table_color_factors       = '
    do i = 1, size(vm%table_color_factors)
       write(vm%out_fh, *)  vm%table_color_factors(i)%i1, &
            vm%table_color_factors(i)%i2, &
            vm%table_color_factors(i)%factor
    end do

    write(vm%out_fh, *) 'table_flv_col_is_allowed  = ', &
                      vm%table_flv_col_is_allowed
    do i = 1, vm%N_flavors
       do j = 1, vm%N_col_flows
          do k = 1, vm%N_helicities
             write(vm%out_fh, *) 'table_amplitudes(f,c,h), f, c, h = ', vm%table_amplitudes(i,j,k), i, j, k
          end do
       end do
    end do
    if (allocated(vm%coupl_real)) then
       write(vm%out_fh, *) 'coupl_real          = ', vm%coupl_real
    end if
    if (allocated(vm%coupl_real2)) then
       write(vm%out_fh, *) 'coupl_real2         = ', vm%coupl_real2
    end if
    if (allocated(vm%coupl_cmplx)) then
       write(vm%out_fh, *) 'coupl_cmplx         = ', vm%coupl_cmplx
    end if
    if (allocated(vm%coupl_cmplx2)) then
       write(vm%out_fh, *) 'coupl_cmplx2        = ', vm%coupl_cmplx2
    end if
    write(vm%out_fh, *) 'mass                = ', vm%mass
    write(vm%out_fh, *) 'width               = ', vm%width
    write(vm%out_fh, *) 'momenta             = ', vm%momenta
    ! gfortran 4.7
    !do i = 1, size(vm%flavor)
       !call vm%flavor(i)%write (vm%out_fh)
    !end do
    !do i = 1, size(vm%color)
       !call vm%color(i)%write (vm%out_fh)
    !end do
    !do i = 1, size(vm%helicity)
       !call vm%helicity(i)%write (vm%out_fh)
    !end do
    write(vm%out_fh, *) 'helicity            = ', vm%helicity
    write(vm%out_fh, *) 'amplitudes       = ', vm%amplitudes
    write(vm%out_fh, *) 'scalars       = ', vm%scalars
    write(vm%out_fh, *) 'spinors       = ', vm%spinors
    write(vm%out_fh, *) 'conjspinors   = ', vm%conjspinors
    write(vm%out_fh, *) 'bispinors     = ', vm%bispinors
    write(vm%out_fh, *) 'vectors       = ', vm%vectors
    write(vm%out_fh, *) 'tensors_2     = ', vm%tensors_2
    write(vm%out_fh, *) 'tensors_1     = ', vm%tensors_1
    !!! !!! !!! Regression with ifort 16.0.0
    !!! write(vm%out_fh, *) 'vectorspinors = ', vm%vectorspinors
    write(vm%out_fh, *) 'N_momenta       = ', vm%N_momenta
    write(vm%out_fh, *) 'N_particles     = ', vm%N_particles
    write(vm%out_fh, *) 'N_prt_in        = ', vm%N_prt_in
    write(vm%out_fh, *) 'N_prt_out       = ', vm%N_prt_out
    write(vm%out_fh, *) 'N_amplitudes    = ', vm%N_amplitudes
    write(vm%out_fh, *) 'N_helicities    = ', vm%N_helicities
    write(vm%out_fh, *) 'N_col_flows     = ', vm%N_col_flows
    write(vm%out_fh, *) 'N_col_indices   = ', vm%N_col_indices
    write(vm%out_fh, *) 'N_flavors       = ', vm%N_flavors
    write(vm%out_fh, *) 'N_col_factors   = ', vm%N_col_factors
    write(vm%out_fh, *) 'N_scalars       = ', vm%N_scalars
    write(vm%out_fh, *) 'N_spinors       = ', vm%N_spinors
    write(vm%out_fh, *) 'N_conjspinors   = ', vm%N_conjspinors
    write(vm%out_fh, *) 'N_bispinors     = ', vm%N_bispinors
    write(vm%out_fh, *) 'N_vectors       = ', vm%N_vectors
    write(vm%out_fh, *) 'N_tensors_2     = ', vm%N_tensors_2
    write(vm%out_fh, *) 'N_tensors_1     = ', vm%N_tensors_1
    write(vm%out_fh, *) 'N_vectorspinors = ', vm%N_vectorspinors
    write(vm%out_fh, *) 'Overall size of VM: '
    ! GNU extension
    ! write(vm%out_fh, *) 'sizeof(wavefunctions) = ', &
    !   sizeof(vm%scalars) + sizeof(vm%spinors) + sizeof(vm%conjspinors) + &
    !   sizeof(vm%bispinors) + sizeof(vm%vectors) + sizeof(vm%tensors_2) + &
    !   sizeof(vm%tensors_1) +  sizeof(vm%vectorspinors)
    ! write(vm%out_fh, *) 'sizeof(mometa) = ', sizeof(vm%momenta)
    ! write(vm%out_fh, *) 'sizeof(amplitudes) = ', sizeof(vm%amplitudes)
    ! write(vm%out_fh, *) 'sizeof(tables) = ', &
    !   sizeof(vm%table_amplitudes) + sizeof(vm%table_spin) + &
    !   sizeof(vm%table_flavor) + sizeof(vm%table_flv_col_is_allowed) + &
    !   sizeof(vm%table_color_flows) + sizeof(vm%table_color_factors) + &
    !   sizeof(vm%table_ghost_flags)
  end subroutine vm_write

  subroutine vm_final (vm)
    class(vm_t), intent(inout) :: vm
    deallocate (vm%table_flavor)
    deallocate (vm%table_color_flows)
    deallocate (vm%table_spin)
    deallocate (vm%table_ghost_flags)
    deallocate (vm%table_color_factors)
    deallocate (vm%table_flv_col_is_allowed)
    if (allocated (vm%coupl_real)) then
       deallocate (vm%coupl_real)
    end if
    if (allocated (vm%coupl_real2)) then
       deallocate (vm%coupl_real2)
    end if
    if (allocated (vm%coupl_cmplx)) then
       deallocate (vm%coupl_cmplx)
    end if
    if (allocated (vm%coupl_cmplx2)) then
       deallocate (vm%coupl_cmplx2)
    end if
    if (allocated (vm%mass)) then
       deallocate (vm%mass)
    end if
    if (allocated (vm%width)) then
       deallocate (vm%width)
    end if
    deallocate (vm%momenta)
    deallocate (vm%flavor)
    deallocate (vm%color)
    deallocate (vm%helicity)
    deallocate (vm%amplitudes)
    deallocate (vm%table_amplitudes)
    deallocate (vm%scalars)
    deallocate (vm%spinors)
    deallocate (vm%conjspinors)
    deallocate (vm%bispinors)
    deallocate (vm%vectors)
    deallocate (vm%tensors_2)
    deallocate (vm%tensors_1)
    deallocate (vm%vectorspinors)
  end subroutine vm_final

  subroutine vm_run (vm, mom, flavor, color, helicity)
    class(vm_t), intent(inout) :: vm
    real(default), dimension(0:3, *), intent(in) :: mom
    class(flavor_t), dimension(:), optional, intent(in) :: flavor
    class(color_t), dimension(:), optional, intent(in) :: color
    ! gfortran 4.7
    !class(helicity_t), dimension(:), optional, target, intent(in) :: helicity
    integer, dimension(:), optional, intent(in) :: helicity
    integer :: i, h, hi
    do i = 1, vm%N_particles
      if (i <= vm%N_prt_in) then
        vm%momenta(i) = - mom(:, i)          ! incoming, crossing symmetry
      else
        vm%momenta(i) = mom(:, i)            ! outgoing
      end if
    end do
    if (present (flavor)) then
       allocate(vm%flavor(size(flavor)), source=flavor)
    else
       if (.not. (allocated (vm%flavor))) then
          allocate(flv_discrete::vm%flavor(vm%N_particles))
       end if
    end if
    if (present (color)) then
       allocate(vm%color(size(color)), source=color)
    else
       if (.not. (allocated (vm%color))) then
          allocate(col_discrete::vm%color(vm%N_col_flows))
       end if
    end if
    ! gfortran 4.7
    if (present (helicity)) then
       !vm%helicity => helicity
       vm%helicity = helicity
       call vm_run_one_helicity (vm, 1)
    else
      !if (.not. (associated (vm%helicity))) then
         !allocate(hel_discrete::vm%helicity(vm%N_particles))
      !end if
      if (.not. (allocated (vm%helicity))) then
         allocate(vm%helicity(vm%N_particles))
      end if
      if (vm%hel_finite == 0) return
      do hi = 1, vm%hel_finite
         h = vm%hel_map(hi)
         !<Work around [[gfortran 4.7 Bug 56731]] Implementation>>
         vm%helicity = vm%table_spin(:,h)
         call vm_run_one_helicity (vm, h)
      end do
    end if
  end subroutine vm_run

  subroutine vm_run_one_helicity (vm, h)
    class(vm_t), intent(inout) :: vm
    integer, intent(in) :: h
    integer :: f, c, i
    vm%amplitudes = zero
    if (vm%N_levels > 0) then
       call null_all_wfs (vm)
       call iterate_instructions (vm)
    end if
    i = 1
    do c = 1, vm%N_col_flows
       do f = 1, vm%N_flavors
          if (vm%table_flv_col_is_allowed(f,c)) then
             vm%table_amplitudes(f,c,h) = vm%amplitudes(i)
             i = i + 1
          end if
       end do
    end do
  end subroutine

  subroutine null_all_wfs (vm)
    type(vm_t), intent(inout) :: vm
    integer :: i, j
    vm%scalars%c = .False.
    vm%scalars%v = zero
    vm%spinors%c = .False.
    vm%conjspinors%c = .False.
    vm%bispinors%c = .False.
    vm%vectorspinors%c = .False.
    do i = 1, 4
       vm%spinors%v%a(i) = zero
       vm%conjspinors%v%a(i) = zero
       vm%bispinors%v%a(i) = zero
       do j = 1, 4
          vm%vectorspinors%v%psi(i)%a(j) = zero
       end do
    end do
    vm%vectors%c = .False.
    vm%vectors%v%t = zero
    vm%tensors_1%c = .False.
    vm%tensors_2%c = .False.
    do i = 1, 3
       vm%vectors%v%x(i) = zero
       vm%tensors_1%v%e(i) = zero
       vm%tensors_1%v%b(i) = zero
       do j = 1, 3
          vm%tensors_2%v%t(i,j) = zero
       end do
    end do
  end subroutine

  subroutine load_header (vm, IO)
    type(vm_t), intent(inout) :: vm
    integer, intent(inout) :: IO
    integer, dimension(len_instructions) :: line
    read(vm%bytecode_fh, fmt = *, iostat = IO) line
    vm%N_momenta = line(1)
    vm%N_particles = line(2)
    vm%N_prt_in = line(3)
    vm%N_prt_out = line(4)
    vm%N_amplitudes = line(5)
    vm%N_helicities = line(6)
    vm%N_col_flows = line(7)
    if (vm%N_momenta == 0) then
       vm%N_col_indices = 2
    else
       vm%N_col_indices = line(8)
    end if
    read(vm%bytecode_fh, fmt = *, iostat = IO)
    read(vm%bytecode_fh, fmt = *, iostat = IO) line
    vm%N_flavors = line(1)
    vm%N_col_factors = line(2)
    vm%N_scalars = line(3)
    vm%N_spinors = line(4)
    vm%N_conjspinors = line(5)
    vm%N_bispinors = line(6)
    vm%N_vectors = line(7)
    vm%N_tensors_2 = line(8)
    read(vm%bytecode_fh, fmt = *, iostat = IO)
    read(vm%bytecode_fh, fmt = *, iostat = IO) line
    vm%N_tensors_1 = line(1)
    vm%N_vectorspinors = line(2)
    ! Add 1 for seperating label lines like 'Another table'
    vm%N_table_lines = vm%N_helicities + 1 + vm%N_flavors + 1 + vm%N_col_flows &
      + 1 + vm%N_col_flows + 1 + vm%N_col_factors + 1 + vm%N_col_flows
  end subroutine load_header

  subroutine read_tables (vm, IO)
    type(vm_t), intent(inout) :: vm
    integer, intent(inout) :: IO
    integer :: i
    integer, dimension(2) :: tmpcf
    integer, dimension(3) :: tmpfactor
    integer, dimension(vm%N_flavors) :: tmpF
    integer, dimension(vm%N_particles) :: tmpP
    real(default) :: factor
    do i = 1, vm%N_helicities
      read(vm%bytecode_fh, fmt = *, iostat = IO) vm%table_spin(:, i)
    end do

    read(vm%bytecode_fh, fmt = *, iostat = IO)
    do i = 1, vm%N_flavors
      read(vm%bytecode_fh, fmt = *, iostat = IO) vm%table_flavor(:, i)
    end do

    read(vm%bytecode_fh, fmt = *, iostat = IO)
    do i = 1, vm%N_col_flows
      read(vm%bytecode_fh, fmt = *, iostat = IO) vm%table_color_flows(:, :, i)
    end do

    read(vm%bytecode_fh, fmt = *, iostat = IO)
    do i = 1, vm%N_col_flows
      read(vm%bytecode_fh, fmt = *, iostat = IO) tmpP
      vm%table_ghost_flags(:, i) = int_to_log(tmpP)
    end do

    read(vm%bytecode_fh, fmt = *, iostat = IO)
    do i = 1, vm%N_col_factors
      read(vm%bytecode_fh, fmt = '(2I9)', iostat = IO, advance='no') tmpcf
      factor = zero
      do
        read(vm%bytecode_fh, fmt = '(3I9)', iostat = IO, advance='no', EOR=10) tmpfactor
        factor = factor + color_factor(tmpfactor(1), tmpfactor(2), tmpfactor(3))
      end do
      10 vm%table_color_factors(i) = OCF(tmpcf(1), tmpcf(2), factor)
    end do

    read(vm%bytecode_fh, fmt = *, iostat = IO)
    do i = 1, vm%N_col_flows
      read(vm%bytecode_fh, fmt = *, iostat = IO) tmpF
      vm%table_flv_col_is_allowed(:, i) = int_to_log(tmpF)
    end do
  end subroutine read_tables

  subroutine extended_version_check (vm, IO)
    type(vm_t), intent(in) :: vm
    integer, intent(inout) :: IO
    character(256) :: buffer
    read(vm%bytecode_fh, fmt = "(A)", iostat = IO) buffer
    if (vm%version /= buffer) then
      print *, "Warning: Bytecode has been generated with an older O'Mega version."
    else
      if (vm%verbose) then
         write (vm%out_fh, fmt = *) "Bytecode version fits."
      end if
    end if
  end subroutine extended_version_check

  subroutine basic_init (vm, verbose, out_fh)
    type(vm_t), intent(inout) :: vm
    logical, optional, intent(in) :: verbose
    integer, optional, intent(in) :: out_fh
    if (present (verbose)) then
       vm%verbose = verbose
    else
       vm%verbose = .true.
    end if
    if (present (out_fh)) then
       vm%out_fh = out_fh
     else
       vm%out_fh = stdout
    end if
    call set_stream (vm)
    call alloc_and_count (vm)
    if (vm%N_levels > 0) then
       call read_bytecode (vm)
       call sanity_check (vm)
    end if
    close (vm%bytecode_fh)
  end subroutine basic_init

  subroutine basic_write (vm)
    type(vm_t), intent(in) :: vm
    integer :: i
    write (vm%out_fh, *) '=====> VM ', char(vm%version), ' <====='
    write (vm%out_fh, *) 'verbose          =    ', vm%verbose
    write (vm%out_fh, *) 'bytecode_file    =    ', char (vm%bytecode_file)
    write (vm%out_fh, *) 'N_instructions   =    ', vm%N_instructions
    write (vm%out_fh, *) 'N_levels         =    ', vm%N_levels
    write (vm%out_fh, *) 'instructions     =    '
    do i = 1, vm%N_instructions
       write (vm%out_fh, *) vm%instructions(:, i)
    end do
    write (vm%out_fh, *) 'levels           =    ', vm%levels
  end subroutine basic_write

  subroutine alloc_and_count (vm)
    type(vm_t), intent(inout) :: vm
    integer, dimension(len_instructions) :: line
    character(256) :: buffer
    integer :: i, IO
    read(vm%bytecode_fh, fmt = "(A)", iostat = IO) buffer
    if (vm%model /= buffer) then
      print *, "Warning: Bytecode has been generated with an older O'Mega version."
    else
      if (vm%verbose) then
         write (vm%out_fh, fmt = *) "Using the model: "
         write (vm%out_fh, fmt = *) char(vm%model)
      end if
    end if
    call extended_version_check (vm, IO)
    if (vm%verbose) then
       write (vm%out_fh, fmt = *) "Trying to allocate."
    end if
    do i = 1, N_comments
      read(vm%bytecode_fh, fmt = *, iostat = IO)
    end do
    call load_header (vm, IO)
    call alloc_arrays (vm)
    if (vm%N_momenta /= 0) then
       do i = 1, vm%N_table_lines + 1
         read(vm%bytecode_fh, fmt = *, iostat = IO)
       end do
       vm%N_instructions = 0
       vm%N_levels = 0
       do
         read(vm%bytecode_fh, fmt = *, end = 42) line
         if (line(1) /= 0) then
           vm%N_instructions = vm%N_instructions + 1
         else
           vm%N_levels = vm%N_levels + 1
         end if
       end do
       42 rewind(vm%bytecode_fh, iostat = IO)
       allocate (vm%instructions(len_instructions, vm%N_instructions))
       allocate (vm%levels(vm%N_levels))
       if (IO /= 0) then
         print *, "Error: vm.alloc : Couldn't load bytecode!"
         stop 1
       end if
    end if
  end subroutine alloc_and_count

  subroutine read_bytecode (vm)
    type(vm_t), intent(inout) :: vm
    integer, dimension(len_instructions) :: line
    integer :: i, j, IO
    ! Jump over version number, comments, header and first table description
    do i = 1, N_version_lines + N_comments + N_header_lines + 1
      read (vm%bytecode_fh, fmt = *, iostat = IO)
    end do
    call read_tables (vm, IO)
    read (vm%bytecode_fh, fmt = *, iostat = IO)
    i = 0; j = 0
    do
      read (vm%bytecode_fh, fmt = *, iostat = IO) line
      if (IO /= 0) exit
      if (line(1) == 0) then
        if (j <= vm%N_levels) then
          j = j + 1
          vm%levels(j) = i                 ! last index of a level is saved
        else
          print *, 'Error: vm.read_bytecode: File has more levels than anticipated!'
          stop 1
        end if
      else
        if (i <= vm%N_instructions) then
          i = i + 1                        ! A valid instruction line
          vm%instructions(:, i) = line
        else
          print *, 'Error: vm.read_bytecode: File is larger than anticipated!'
          stop 1
        end if
      end if
    end do
  end subroutine read_bytecode

  subroutine iterate_instructions (vm)
    type(vm_t), intent(inout) :: vm
    integer :: i, j
    if (vm%openmp) then
       !$omp parallel
       do j = 1, vm%N_levels - 1
          !$omp do schedule (static)
          do i = vm%levels (j) + 1, vm%levels (j + 1)
             call decode (vm, i)
          end do
          !$omp end do
       end do
       !$omp end parallel
    else
       do j = 1, vm%N_levels - 1
          do i = vm%levels (j) + 1, vm%levels (j + 1)
             call decode (vm, i)
          end do
       end do
    end if
  end subroutine iterate_instructions

  subroutine set_stream (vm)
    type(vm_t), intent(inout) :: vm
    integer :: IO
    call find_free_unit (vm%bytecode_fh, IO)
    open (vm%bytecode_fh, file = char (vm%bytecode_file), form = 'formatted', &
      access = 'sequential', status = 'old', position = 'rewind', iostat = IO, &
      action = 'read')
    if (IO /= 0) then
      print *, "Error: vm.set_stream: Bytecode file '", char(vm%bytecode_file), &
               "' not found!"
      stop 1
    end if
  end subroutine set_stream

  subroutine sanity_check (vm)
    type(vm_t), intent(in) :: vm
    if (vm%levels(1) /= 0) then
       print *, "Error: vm.vm_init: levels(1) != 0"
       stop 1
    end if
    if (vm%levels(vm%N_levels) /= vm%N_instructions) then
       print *, "Error: vm.vm_init: levels(N_levels) != N_instructions"
       stop 1
    end if
    if (vm%verbose) then
      write(vm%out_fh, *) "vm passed sanity check. Starting calculation."
    end if
  end subroutine sanity_check

    ! pure & ! if no warnings
    subroutine decode (vm, instruction_index)
      type(vm_t), intent(inout) :: vm
      integer, intent(in) :: instruction_index
      integer, dimension(len_instructions) :: i, curr
      complex(default) :: braket
      integer :: tmp
      real(default) :: w
      i = vm%instructions (:, instruction_index)
      select case (i(1))
      case ( : -1)       ! Jump over subinstructions

      case (ovm_ADD_MOMENTA)
         vm%momenta(i(4)) = vm%momenta(i(5)) + vm%momenta(i(6))
         if (i(7) > 0) then
            vm%momenta(i(4)) = vm%momenta(i(4)) + vm%momenta(i(7))
         end if

      case (ovm_LOAD_SCALAR)
        vm%scalars(i(4))%v = one
        vm%scalars(i(4))%c = .True.

      case (ovm_LOAD_SPINOR_INC)
         call load_spinor(vm%spinors(i(4)), - vm%momenta(i(5)), vm%mass(i(2)), &
                          vm%helicity(i(5)), ovm_LOAD_SPINOR_INC)

      case (ovm_LOAD_SPINOR_OUT)
         call load_spinor(vm%spinors(i(4)), vm%momenta(i(5)), vm%mass(i(2)), &
                          vm%helicity(i(5)), ovm_LOAD_SPINOR_OUT)

      case (ovm_LOAD_CONJSPINOR_INC)
         call load_conjspinor(vm%conjspinors(i(4)), - vm%momenta(i(5)), &
           vm%mass(i(2)), vm%helicity(i(5)), ovm_LOAD_CONJSPINOR_INC)

      case (ovm_LOAD_CONJSPINOR_OUT)
         call load_conjspinor(vm%conjspinors(i(4)), vm%momenta(i(5)), &
           vm%mass(i(2)), vm%helicity(i(5)), ovm_LOAD_CONJSPINOR_OUT)

      case (ovm_LOAD_MAJORANA_INC)
         call load_bispinor(vm%bispinors(i(4)), - vm%momenta(i(5)), &
           vm%mass(i(2)), vm%helicity(i(5)), ovm_LOAD_MAJORANA_INC)

      case (ovm_LOAD_MAJORANA_OUT)
         call load_bispinor(vm%bispinors(i(4)), vm%momenta(i(5)), vm%mass(i(2)), &
                          vm%helicity(i(5)), ovm_LOAD_MAJORANA_OUT)

      case (ovm_LOAD_VECTOR_INC)
         call load_vector(vm%vectors(i(4)), - vm%momenta(i(5)), vm%mass(i(2)), &
                          vm%helicity(i(5)), ovm_LOAD_VECTOR_INC)

      case (ovm_LOAD_VECTOR_OUT)
         call load_vector(vm%vectors(i(4)), vm%momenta(i(5)), vm%mass(i(2)), &
                          vm%helicity(i(5)), ovm_LOAD_VECTOR_OUT)

      case (ovm_LOAD_VECTORSPINOR_INC)
        !select type (h => vm%helicity(i(5)))
        !type is (hel_discrete)
           !vm%vectorspinors(i(4))%v = veps(vm%mass(i(2)), - vm%momenta(i(5)), &
                                           !h%i)
        !end select
        vm%vectorspinors(i(4))%v = veps(vm%mass(i(2)), - vm%momenta(i(5)), &
                                        vm%helicity(i(5)))
        vm%vectorspinors(i(4))%c = .True.

      case (ovm_LOAD_VECTORSPINOR_OUT)
        !select type (h => vm%helicity(i(5)))
        !type is (hel_discrete)
           !vm%vectorspinors(i(4))%v = veps(vm%mass(i(2)), vm%momenta(i(5)), &
                                           !h%i)
        !end select
        vm%vectorspinors(i(4))%v = veps(vm%mass(i(2)), vm%momenta(i(5)), &
                                        vm%helicity(i(5)))
        vm%vectorspinors(i(4))%c = .True.

      case (ovm_LOAD_TENSOR2_INC)
        !select type (h => vm%helicity(i(5)))
        !type is (hel_discrete)
           !vm%tensors_2(i(4))%v = eps2(vm%mass(i(2)), - vm%momenta(i(5)), &
                                       !h%i)
        !end select
        vm%tensors_2(i(4))%c = .True.

      case (ovm_LOAD_TENSOR2_OUT)
        !select type (h => vm%helicity(i(5)))
        !type is (hel_discrete)
           !vm%tensors_2(i(4))%v = eps2(vm%mass(i(2)), vm%momenta(i(5)), h%i)
        !end select
        vm%tensors_2(i(4))%c = .True.

      case (ovm_LOAD_BRS_SCALAR)
        vm%scalars(i(4))%v = (0, -1) * (vm%momenta(i(5)) * vm%momenta(i(5)) - &
                                        vm%mass(i(2))**2)
        vm%scalars(i(4))%c = .True.

      case (ovm_LOAD_BRS_SPINOR_INC)
        print *, 'not implemented'
        stop 1
      case (ovm_LOAD_BRS_SPINOR_OUT)
        print *, 'not implemented'
        stop 1
      case (ovm_LOAD_BRS_CONJSPINOR_INC)
        print *, 'not implemented'
        stop 1
      case (ovm_LOAD_BRS_CONJSPINOR_OUT)
        print *, 'not implemented'
        stop 1
      case (ovm_LOAD_BRS_VECTOR_INC)
        print *, 'not implemented'
        stop 1
      case (ovm_LOAD_BRS_VECTOR_OUT)
        print *, 'not implemented'
        stop 1
      case (ovm_LOAD_MAJORANA_GHOST_INC)
        print *, 'not implemented'
        stop 1
      case (ovm_LOAD_MAJORANA_GHOST_OUT)
        print *, 'not implemented'
        stop 1
      case (ovm_LOAD_BRS_MAJORANA_INC)
        print *, 'not implemented'
        stop 1
      case (ovm_LOAD_BRS_MAJORANA_OUT)
        print *, 'not implemented'
        stop 1

      case (ovm_CALC_BRAKET)
         
         tmp = instruction_index + 1
         do
           if (tmp > vm%N_instructions) exit
           curr = vm%instructions(:, tmp)
           if (curr(1) >= 0) exit                        ! End of fusions
           select case (curr(1))
           case (ovm_FUSE_V_FF, ovm_FUSE_VL_FF, ovm_FUSE_VR_FF)
             braket = vm%vectors(curr(4))%v * vec_ff(vm, curr)

           case (ovm_FUSE_F_VF, ovm_FUSE_F_VLF, ovm_FUSE_F_VRF)
             braket = vm%conjspinors(curr(4))%v * ferm_vf(vm, curr)

           case (ovm_FUSE_F_FV, ovm_FUSE_F_FVL, ovm_FUSE_F_FVR)
             braket = ferm_fv(vm, curr) * vm%spinors(curr(4))%v

           case (ovm_FUSE_VA_FF)
             braket = vm%vectors(curr(4))%v * vec_ff2(vm, curr)

           case (ovm_FUSE_F_VAF)
             braket = vm%conjspinors(curr(4))%v * ferm_vf2(vm, curr)

           case (ovm_FUSE_F_FVA)
             braket = ferm_fv2(vm, curr) * vm%spinors(curr(4))%v

           case (ovm_FUSE_S_FF, ovm_FUSE_SP_FF)
             braket = vm%scalars(curr(4))%v * scal_ff(vm, curr)

           case (ovm_FUSE_F_SF, ovm_FUSE_F_SPF)
             braket = vm%conjspinors(curr(4))%v * ferm_sf(vm, curr)

           case (ovm_FUSE_F_FS, ovm_FUSE_F_FSP)
             braket = ferm_fs(vm, curr) * vm%spinors(curr(4))%v

           case (ovm_FUSE_G_GG)
             braket = vm%vectors(curr(4))%v * &
               g_gg(sgn_coupl_cmplx(vm, curr(2)), &
                    vm%vectors(curr(5))%v, vm%momenta(curr(6)), &
                    vm%vectors(curr(7))%v, vm%momenta(curr(8)))

           case (ovm_FUSE_S_VV)
             braket = vm%scalars(curr(4))%v * sgn_coupl_cmplx(vm, curr(2)) * &
                      (vm%vectors(curr(5))%v * vm%vectors(curr(6))%v)

           case (ovm_FUSE_V_SS)
             braket = vm%vectors(curr(4))%v * &
                      v_ss(sgn_coupl_cmplx(vm, curr(2)), vm%scalars(curr(5))%v, vm%momenta(curr(6)), &
                                  vm%scalars(curr(7))%v, vm%momenta(curr(8)))

           case (ovm_FUSE_S_G2, ovm_FUSE_S_G2_SKEW)
             braket = vm%scalars(curr(4))%v * scal_g2(vm, curr)

           case (ovm_FUSE_G_SG, ovm_FUSE_G_GS, ovm_FUSE_G_SG_SKEW, ovm_FUSE_G_GS_SKEW)
             braket = vm%vectors(curr(4))%v * gauge_sg(vm, curr)

           case (ovm_FUSE_S_VS)
             braket = vm%scalars(curr(4))%v * &
               s_vs(sgn_coupl_cmplx(vm, curr(2)), &
                    vm%vectors(curr(5))%v, vm%momenta(curr(6)), &
                    vm%scalars(curr(7))%v, vm%momenta(curr(8)))

           case (ovm_FUSE_V_SV)
             braket = (vm%vectors(curr(4))%v * vm%vectors(curr(6))%v) * &
                      (sgn_coupl_cmplx(vm, curr(2)) * vm%scalars(curr(5))%v)

           case (ovm_FUSE_S_SS)
             braket = vm%scalars(curr(4))%v * &
               sgn_coupl_cmplx(vm, curr(2)) * &
               (vm%scalars(curr(5))%v * vm%scalars(curr(6))%v)

           case (ovm_FUSE_S_SSS)
             braket = vm%scalars(curr(4))%v * &
               sgn_coupl_cmplx(vm, curr(2)) * &
               (vm%scalars(curr(5))%v * vm%scalars(curr(6))%v * &
                vm%scalars(curr(7))%v)

           case (ovm_FUSE_S_SVV)
             braket = vm%scalars(curr(4))%v * &
               sgn_coupl_cmplx(vm, curr(2)) * &
               vm%scalars(curr(5))%v * (vm%vectors(curr(6))%v * &
                                        vm%vectors(curr(7))%v)

           case (ovm_FUSE_V_SSV)
             braket = vm%vectors(curr(4))%v * &
               (sgn_coupl_cmplx(vm, curr(2)) * vm%scalars(curr(5))%v * &
                vm%scalars(curr(6))%v) * vm%vectors(curr(7))%v

           case (ovm_FUSE_V_VVV)
             braket = sgn_coupl_cmplx(vm, curr(2)) * &
               (vm%vectors(curr(5))%v * vm%vectors(curr(6))%v) * &
               (vm%vectors(curr(4))%v * vm%vectors(curr(7))%v)

           case default
             print *, 'Braket', curr(1), 'not implemented'
             stop 1

           end select
           vm%amplitudes(i(4)) = vm%amplitudes(i(4)) + curr(3) * braket
           tmp = tmp + 1
         end do

         vm%amplitudes(i(4)) = vm%amplitudes(i(4)) * i(2)
         if (i(5) > 1) then
           vm%amplitudes(i(4)) = vm%amplitudes(i(4)) * &         ! Symmetry factor
                                 (one / sqrt(real(i(5), kind=default)))
         end if

      
      case (ovm_PROPAGATE_SCALAR : ovm_PROPAGATE_NONE)
         tmp = instruction_index + 1
         do
           curr = vm%instructions(:,tmp)
           if (curr(1) >= 0) exit                        ! End of fusions
           select case (curr(1))
           case (ovm_FUSE_V_FF, ovm_FUSE_VL_FF, ovm_FUSE_VR_FF)
             vm%vectors(curr(4))%v = vm%vectors(curr(4))%v + curr(3) * &
                                     vec_ff(vm, curr)

           case (ovm_FUSE_F_VF, ovm_FUSE_F_VLF, ovm_FUSE_F_VRF)
             vm%spinors(curr(4))%v = vm%spinors(curr(4))%v + curr(3) * &
                                     ferm_vf(vm, curr)

           case (ovm_FUSE_F_FV, ovm_FUSE_F_FVL, ovm_FUSE_F_FVR)
             vm%conjspinors(curr(4))%v = vm%conjspinors(curr(4))%v + curr(3) * &
                                         ferm_fv(vm, curr)

           case (ovm_FUSE_VA_FF)
             vm%vectors(curr(4))%v = vm%vectors(curr(4))%v + curr(3) * &
                                     vec_ff2(vm, curr)

           case (ovm_FUSE_F_VAF)
             vm%spinors(curr(4))%v = vm%spinors(curr(4))%v + curr(3) * &
                                     ferm_vf2(vm, curr)

           case (ovm_FUSE_F_FVA)
             vm%conjspinors(curr(4))%v = vm%conjspinors(curr(4))%v + curr(3) * &
                                         ferm_fv2(vm, curr)

           case (ovm_FUSE_S_FF, ovm_FUSE_SP_FF)
             vm%scalars(curr(4))%v = vm%scalars(curr(4))%v + curr(3) * &
                                     scal_ff(vm, curr)

           case (ovm_FUSE_F_SF, ovm_FUSE_F_SPF)
             vm%spinors(curr(4))%v = vm%spinors(curr(4))%v + curr(3) * &
                                     ferm_sf(vm, curr)

           case (ovm_FUSE_F_FS, ovm_FUSE_F_FSP)
             vm%conjspinors(curr(4))%v = vm%conjspinors(curr(4))%v + curr(3) * &
                                         ferm_fs(vm, curr)

           case (ovm_FUSE_G_GG)
             vm%vectors(curr(4))%v = vm%vectors(curr(4))%v + curr(3) * &
               g_gg(sgn_coupl_cmplx(vm, curr(2)), vm%vectors(curr(5))%v, &
                    vm%momenta(curr(6)), vm%vectors(curr(7))%v, &
                    vm%momenta(curr(8)))

           case (ovm_FUSE_S_VV)
             vm%scalars(curr(4))%v = vm%scalars(curr(4))%v + curr(3) * &
               sgn_coupl_cmplx(vm, curr(2)) * &
               (vm%vectors(curr(5))%v * vm%vectors(curr(6))%v)

           case (ovm_FUSE_V_SS)
             vm%vectors(curr(4))%v = vm%vectors(curr(4))%v + curr(3) * &
                      v_ss(sgn_coupl_cmplx(vm, curr(2)), vm%scalars(curr(5))%v, vm%momenta(curr(6)), &
                                  vm%scalars(curr(7))%v, vm%momenta(curr(8)))


           case (ovm_FUSE_S_G2, ovm_FUSE_S_G2_SKEW)
             vm%scalars(curr(4))%v = vm%scalars(curr(4))%v + &
                                     scal_g2(vm, curr) * curr(3)

           case (ovm_FUSE_G_SG, ovm_FUSE_G_GS, ovm_FUSE_G_SG_SKEW, ovm_FUSE_G_GS_SKEW)
             vm%vectors(curr(4))%v = vm%vectors(curr(4))%v + &
                                     gauge_sg(vm, curr) * curr(3)

           case (ovm_FUSE_S_VS)
             vm%scalars(curr(4))%v = vm%scalars(curr(4))%v + &
               s_vs(sgn_coupl_cmplx(vm, curr(2)), &
                    vm%vectors(curr(5))%v, vm%momenta(curr(6)), &
                    vm%scalars(curr(7))%v, vm%momenta(curr(8))) * curr(3)

           case (ovm_FUSE_V_SV)
             vm%vectors(curr(4))%v = vm%vectors(curr(4))%v + &
               vm%vectors(curr(6))%v * &
               (sgn_coupl_cmplx(vm, curr(2)) * vm%scalars(curr(5))%v * curr(3))

           case (ovm_FUSE_S_SS)
             vm%scalars(curr(4))%v = vm%scalars(curr(4))%v + &
               sgn_coupl_cmplx(vm, curr(2)) * &
               (vm%scalars(curr(5))%v * vm%scalars(curr(6))%v) * curr(3)

           case (ovm_FUSE_S_SSS)
             vm%scalars(curr(4))%v = vm%scalars(curr(4))%v + &
               sgn_coupl_cmplx(vm, curr(2)) * &
               (vm%scalars(curr(5))%v * vm%scalars(curr(6))%v * &
                vm%scalars(curr(7))%v) * curr(3)

           case (ovm_FUSE_S_SVV)
             vm%scalars(curr(4))%v = vm%scalars(curr(4))%v + &
               sgn_coupl_cmplx(vm, curr(2)) * &
               vm%scalars(curr(5))%v * (vm%vectors(curr(6))%v * &
                                        vm%vectors(curr(7))%v) * curr(3)

           case (ovm_FUSE_V_SSV)
             vm%vectors(curr(4))%v = vm%vectors(curr(4))%v + &
               (sgn_coupl_cmplx(vm, curr(2)) * vm%scalars(curr(5))%v * &
                vm%scalars(curr(6))%v) * vm%vectors(curr(7))%v * curr(3)

           case (ovm_FUSE_V_VVV)
             vm%vectors(curr(4))%v = vm%vectors(curr(4))%v + &
               (sgn_coupl_cmplx(vm, curr(2)) * (vm%vectors(curr(5))%v * &
                vm%vectors(curr(6))%v)) * curr(3) * vm%vectors(curr(7))%v

           case default
             print *, 'Fusion', curr(1), 'not implemented'
             stop 1

           end select
           tmp = tmp + 1
         end do

         select case (i(3))
         case (0)
           w = zero

         case (1)
           w = vm%width(i(2))
           vm%cms = .false.

         case (2)
           w = wd_tl(vm%momenta(i(5)), vm%width(i(2)))

         case (3)
           w = vm%width(i(2))
           vm%cms = .true.

         case (4)
           w = wd_run(vm%momenta(i(5)), vm%mass(i(2)), vm%width(i(2)))

         case default
            print *, 'not implemented'
            stop 1

         end select

         select case (i(1))
         case (ovm_PROPAGATE_SCALAR)
           vm%scalars(i(4))%v = pr_phi(vm%momenta(i(5)), vm%mass(i(2)), &
                w, vm%scalars(i(4))%v)
           vm%scalars(i(4))%c = .True.

         case (ovm_PROPAGATE_COL_SCALAR)
           vm%scalars(i(4))%v = - one / N_ * pr_phi(vm%momenta(i(5)), &
                vm%mass(i(2)), w, vm%scalars(i(4))%v)
           vm%scalars(i(4))%c = .True.

         case (ovm_PROPAGATE_GHOST)
           vm%scalars(i(4))%v = imago * pr_phi(vm%momenta(i(5)), vm%mass(i(2)), &
                w, vm%scalars(i(4))%v)
           vm%scalars(i(4))%c = .True.

         case (ovm_PROPAGATE_SPINOR)
           vm%spinors(i(4))%v = pr_psi(vm%momenta(i(5)), vm%mass(i(2)), &
                w, vm%cms, vm%spinors(i(4))%v)
           vm%spinors(i(4))%c = .True.

         case (ovm_PROPAGATE_CONJSPINOR)
           vm%conjspinors(i(4))%v = pr_psibar(vm%momenta(i(5)), vm%mass(i(2)), &
                w, vm%cms, vm%conjspinors(i(4))%v)
           vm%conjspinors(i(4))%c = .True.

         case (ovm_PROPAGATE_MAJORANA)
           vm%bispinors(i(4))%v = bi_pr_psi(vm%momenta(i(5)), vm%mass(i(2)), &
                w, vm%cms, vm%bispinors(i(4))%v)
           vm%bispinors(i(4))%c = .True.

         case (ovm_PROPAGATE_COL_MAJORANA)
           vm%bispinors(i(4))%v = (- one / N_) * &
                bi_pr_psi(vm%momenta(i(5)), vm%mass(i(2)), &
                w, vm%cms, vm%bispinors(i(4))%v)
           vm%bispinors(i(4))%c = .True.

         case (ovm_PROPAGATE_UNITARITY)
           vm%vectors(i(4))%v = pr_unitarity(vm%momenta(i(5)), vm%mass(i(2)), &
                w, vm%cms, vm%vectors(i(4))%v)
           vm%vectors(i(4))%c = .True.

         case (ovm_PROPAGATE_COL_UNITARITY)
           vm%vectors(i(4))%v = - one / N_ * pr_unitarity(vm%momenta(i(5)), &
                vm%mass(i(2)), w, vm%cms, vm%vectors(i(4))%v)
           vm%vectors(i(4))%c = .True.

         case (ovm_PROPAGATE_FEYNMAN)
           vm%vectors(i(4))%v = pr_feynman(vm%momenta(i(5)), vm%vectors(i(4))%v)
           vm%vectors(i(4))%c = .True.

         case (ovm_PROPAGATE_COL_FEYNMAN)
           vm%vectors(i(4))%v = - one / N_ * &
                pr_feynman(vm%momenta(i(5)), vm%vectors(i(4))%v)
           vm%vectors(i(4))%c = .True.

         case (ovm_PROPAGATE_VECTORSPINOR)
           vm%vectorspinors(i(4))%v = pr_grav(vm%momenta(i(5)), vm%mass(i(2)), &
                w, vm%vectorspinors(i(4))%v)
           vm%vectorspinors(i(4))%c = .True.

         case (ovm_PROPAGATE_TENSOR2)
           vm%tensors_2(i(4))%v = pr_tensor(vm%momenta(i(5)), vm%mass(i(2)), &
                w, vm%tensors_2(i(4))%v)
           vm%tensors_2(i(4))%c = .True.

         case (ovm_PROPAGATE_NONE)
         ! This will not work with color MC. Appropriate type%c has to be set to
         ! .True.

         end select

      case (0)
        print *, 'Error: Levelbreak put in decode! Line:', &
                 instruction_index
        stop 1
      case default
        print *, "Error: Decode has case not catched! Line: ", &
                 instruction_index
        stop 1
      end select
    end subroutine decode

  subroutine load_bispinor(wf, p, m, h, opcode)
    type(vm_bispinor), intent(out) :: wf
    type(momentum), intent(in) :: p
    real(default), intent(in) :: m
    !class(helicity_t), intent(in) :: h
    integer, intent(in) :: h
    integer, intent(in) :: opcode
    procedure(bi_u), pointer :: load_wf
    
    select case (opcode)
    case (ovm_LOAD_MAJORANA_INC)
       load_wf => bi_u
    case (ovm_LOAD_MAJORANA_OUT)
       load_wf => bi_v
    case default
       load_wf => null()
    end select
    !select type (h)
    !type is (hel_trigonometric)
       !wf%v = (cos (h%theta) * load_wf (m, p, + 1) + &
               !sin (h%theta) * load_wf (m, p, - 1)) * sqrt2
    !type is (hel_exponential)
       !wf%v = exp (+ imago * h%phi) * load_wf (m, p, + 1) + &
              !exp (- imago * h%phi) * load_wf (m, p, - 1)
    !type is (hel_spherical)
       !wf%v = (exp (+ imago * h%phi) * cos (h%theta) * load_wf (m, p, + 1) + &
               !exp (- imago * h%phi) * sin (h%theta) * load_wf (m, p, - 1)) * &
              !sqrt2
    !type is(hel_discrete)
         !wf%v = load_wf (m, p, h%i)
    !end select
    wf%v = load_wf (m, p, h)
    wf%c = .True.
  end subroutine load_bispinor

  subroutine load_spinor(wf, p, m, h, opcode)
    type(vm_spinor), intent(out) :: wf
    type(momentum), intent(in) :: p
    real(default), intent(in) :: m
    !class(helicity_t), intent(in) :: h
    integer, intent(in) :: h
    integer, intent(in) :: opcode
    procedure(u), pointer :: load_wf
    
    select case (opcode)
    case (ovm_LOAD_SPINOR_INC)
       load_wf => u
    case (ovm_LOAD_SPINOR_OUT)
       load_wf => v
    case default
       load_wf => null()
    end select
    !select type (h)
    !type is (hel_trigonometric)
       !wf%v = (cos (h%theta) * load_wf (m, p, + 1) + &
               !sin (h%theta) * load_wf (m, p, - 1)) * sqrt2
    !type is (hel_exponential)
       !wf%v = exp (+ imago * h%phi) * load_wf (m, p, + 1) + &
              !exp (- imago * h%phi) * load_wf (m, p, - 1)
    !type is (hel_spherical)
       !wf%v = (exp (+ imago * h%phi) * cos (h%theta) * load_wf (m, p, + 1) + &
               !exp (- imago * h%phi) * sin (h%theta) * load_wf (m, p, - 1)) * &
              !sqrt2
    !type is(hel_discrete)
         !wf%v = load_wf (m, p, h%i)
    !end select
    wf%v = load_wf (m, p, h)
    wf%c = .True.
  end subroutine load_spinor

  subroutine load_conjspinor(wf, p, m, h, opcode)
    type(vm_conjspinor), intent(out) :: wf
    type(momentum), intent(in) :: p
    real(default), intent(in) :: m
    !class(helicity_t), intent(in) :: h
    integer, intent(in) :: h
    integer, intent(in) :: opcode
    procedure(ubar), pointer :: load_wf
    
    select case (opcode)
    case (ovm_LOAD_CONJSPINOR_INC)
       load_wf => vbar
    case (ovm_LOAD_CONJSPINOR_OUT)
       load_wf => ubar
    case default
       load_wf => null()
    end select
    !select type (h)
    !type is (hel_trigonometric)
       !wf%v = (cos (h%theta) * load_wf (m, p, + 1) + &
               !sin (h%theta) * load_wf (m, p, - 1)) * sqrt2
    !type is (hel_exponential)
       !wf%v = exp (+ imago * h%phi) * load_wf (m, p, + 1) + &
              !exp (- imago * h%phi) * load_wf (m, p, - 1)
    !type is (hel_spherical)
       !wf%v = (exp (+ imago * h%phi) * cos (h%theta) * load_wf (m, p, + 1) + &
               !exp (- imago * h%phi) * sin (h%theta) * load_wf (m, p, - 1)) * &
              !sqrt2
    !type is(hel_discrete)
         !wf%v = load_wf (m, p, h%i)
    !end select
    wf%v = load_wf (m, p, h)
    wf%c = .True.
  end subroutine load_conjspinor

  subroutine load_vector(wf, p, m, h, opcode)
    type(vm_vector), intent(out) :: wf
    type(momentum), intent(in) :: p
    real(default), intent(in) :: m
    !class(helicity_t), intent(in) :: h
    integer, intent(in) :: h
    integer, intent(in) :: opcode
    procedure(eps), pointer :: load_wf
    
    load_wf => eps
    !select type (h)
    !type is (hel_trigonometric)
       !wf%v = (cos (h%theta) * load_wf (m, p, + 1) + &
               !sin (h%theta) * load_wf (m, p, - 1)) * sqrt2
    !type is (hel_exponential)
       !wf%v = exp (+ imago * h%phi) * load_wf (m, p, + 1) + &
              !exp (- imago * h%phi) * load_wf (m, p, - 1)
    !type is (hel_spherical)
       !wf%v = (exp (+ imago * h%phi) * cos (h%theta) * load_wf (m, p, + 1) + &
               !exp (- imago * h%phi) * sin (h%theta) * load_wf (m, p, - 1)) * &
              !sqrt2
    !type is(hel_discrete)
         !wf%v = load_wf (m, p, h%i)
    !end select
    wf%v = load_wf (m, p, h)
    wf%c = .True.
    if (opcode == ovm_LOAD_VECTOR_OUT) then
       wf%v = conjg(wf%v)
    end if
  end subroutine load_vector

  function ferm_vf(vm, curr) result (x)
    type(spinor) :: x
    class(vm_t), intent(in) :: vm
    integer, dimension(:), intent(in) :: curr
    procedure(f_vf), pointer :: load_wf
    select case (curr(1))
    case (ovm_FUSE_F_VF)
       load_wf => f_vf
    case (ovm_FUSE_F_VLF)
       load_wf => f_vlf
    case (ovm_FUSE_F_VRF)
       load_wf => f_vrf
    case default
       load_wf => null()
    end select
    x = load_wf(sgn_coupl_cmplx(vm, curr(2)), vm%vectors(curr(5))%v, vm%spinors(curr(6))%v)
  end function ferm_vf

  function ferm_vf2(vm, curr) result (x)
    type(spinor) :: x
    class(vm_t), intent(in) :: vm
    integer, dimension(:), intent(in) :: curr
    procedure(f_vaf), pointer :: load_wf
    select case (curr(1))
    case (ovm_FUSE_F_VAF)
       load_wf => f_vaf
    case default
       load_wf => null()
    end select
    x = f_vaf(sgn_coupl_cmplx2(vm, curr(2), 1), sgn_coupl_cmplx2(vm, curr(2), 2), vm%vectors(curr(5))%v, vm%spinors(curr(6))%v)
  end function ferm_vf2

  function ferm_sf(vm, curr) result (x)
    type(spinor) :: x
    class(vm_t), intent(in) :: vm
    integer, dimension(:), intent(in) :: curr
    select case (curr(1))
    case (ovm_FUSE_F_SF)
       x = f_sf(sgn_coupl_cmplx(vm, curr(2)), vm%scalars(curr(5))%v, vm%spinors(curr(6))%v)
    case (ovm_FUSE_F_SPF)
       x = f_spf(sgn_coupl_cmplx2(vm, curr(2), 1), sgn_coupl_cmplx2(vm, curr(2), 2), vm%scalars(curr(5))%v, vm%spinors(curr(6))%v)
    case default
    end select
  end function ferm_sf

  function ferm_fv(vm, curr) result (x)
    type(conjspinor) :: x
    class(vm_t), intent(in) :: vm
    integer, dimension(:), intent(in) :: curr
    procedure(f_fv), pointer :: load_wf
    select case (curr(1))
    case (ovm_FUSE_F_FV)
       load_wf => f_fv
    case (ovm_FUSE_F_FVL)
       load_wf => f_fvl
    case (ovm_FUSE_F_FVR)
       load_wf => f_fvr
    case default
       load_wf => null()
    end select
    x = load_wf(sgn_coupl_cmplx(vm, curr(2)), vm%conjspinors(curr(5))%v, vm%vectors(curr(6))%v)
  end function ferm_fv

  function ferm_fv2(vm, curr) result (x)
    type(conjspinor) :: x
    class(vm_t), intent(in) :: vm
    integer, dimension(:), intent(in) :: curr
    procedure(f_fva), pointer :: load_wf
    select case (curr(1))
    case (ovm_FUSE_F_FVA)
       load_wf => f_fva
    case default
       load_wf => null()
    end select
    x = f_fva(sgn_coupl_cmplx2(vm, curr(2), 1), sgn_coupl_cmplx2(vm, curr(2), 2), &
              vm%conjspinors(curr(5))%v, vm%vectors(curr(6))%v)
  end function ferm_fv2

  function ferm_fs(vm, curr) result (x)
    type(conjspinor) :: x
    class(vm_t), intent(in) :: vm
    integer, dimension(:), intent(in) :: curr
    procedure(f_fs), pointer :: load_wf
    select case (curr(1))
    case (ovm_FUSE_F_FS)
       x = f_fs(sgn_coupl_cmplx(vm, curr(2)), vm%conjspinors(curr(5))%v, vm%scalars(curr(6))%v)
    case (ovm_FUSE_F_FSP)
       x = f_fsp(sgn_coupl_cmplx2(vm, curr(2), 1), sgn_coupl_cmplx2(vm, curr(2), 2), &
            vm%conjspinors(curr(5))%v, vm%scalars(curr(6))%v)
    case default
       x%a = zero
    end select
  end function ferm_fs

  function vec_ff(vm, curr) result (x)
    type(vector) :: x
    class(vm_t), intent(in) :: vm
    integer, dimension(:), intent(in) :: curr
    procedure(v_ff), pointer :: load_wf
    select case (curr(1))
    case (ovm_FUSE_V_FF)
       load_wf => v_ff
    case (ovm_FUSE_VL_FF)
       load_wf => vl_ff
    case (ovm_FUSE_VR_FF)
       load_wf => vr_ff
    case default
       load_wf => null()
    end select
    x = load_wf(sgn_coupl_cmplx(vm, curr(2)), vm%conjspinors(curr(5))%v, vm%spinors(curr(6))%v)
  end function vec_ff

  function vec_ff2(vm, curr) result (x)
    type(vector) :: x
    class(vm_t), intent(in) :: vm
    integer, dimension(:), intent(in) :: curr
    procedure(va_ff), pointer :: load_wf
    select case (curr(1))
    case (ovm_FUSE_VA_FF)
       load_wf => va_ff
    case default
       load_wf => null()
    end select
    x = load_wf(sgn_coupl_cmplx2(vm, curr(2), 1), sgn_coupl_cmplx2(vm, curr(2), 2), &
                vm%conjspinors(curr(5))%v, vm%spinors(curr(6))%v)
  end function vec_ff2

  function scal_ff(vm, curr) result (x)
    complex(default) :: x
    class(vm_t), intent(in) :: vm
    integer, dimension(:), intent(in) :: curr
    select case (curr(1))
    case (ovm_FUSE_S_FF)
       x = s_ff(sgn_coupl_cmplx(vm, curr(2)), &
            vm%conjspinors(curr(5))%v, vm%spinors(curr(6))%v)
    case (ovm_FUSE_SP_FF)
       x = sp_ff(sgn_coupl_cmplx2(vm, curr(2), 1), sgn_coupl_cmplx2(vm, curr(2), 2), &
            vm%conjspinors(curr(5))%v, vm%spinors(curr(6))%v)
    case default
       x = zero
    end select
  end function scal_ff

  function scal_g2(vm, curr) result (x)
    complex(default) :: x
    class(vm_t), intent(in) :: vm
    integer, dimension(:), intent(in) :: curr
    select case (curr(1))
    case (ovm_FUSE_S_G2)
       x = sgn_coupl_cmplx(vm, curr(2)) * ((vm%momenta(curr(6)) * vm%vectors(curr(7))%v) * &
                    (vm%momenta(curr(8)) * vm%vectors(curr(5))%v) - &
                    (vm%momenta(curr(6)) * vm%momenta(curr(8))) * &
                    (vm%vectors(curr(7))%v * vm%vectors(curr(5))%v))
    case (ovm_FUSE_S_G2_SKEW)
       x = - phi_vv(sgn_coupl_cmplx(vm, curr(2)), vm%momenta(curr(6)), vm%momenta(curr(8)), &
                    vm%vectors(curr(5))%v, vm%vectors(curr(7))%v)
    case default
       x = zero
    end select
  end function scal_g2

  pure function gauge_sg(vm, curr) result (x)
    type(vector) :: x
    class(vm_t), intent(in) :: vm
    integer, dimension(:), intent(in) :: curr
    select case (curr(1))
    case (ovm_FUSE_G_SG)
       x = sgn_coupl_cmplx(vm, curr(2)) * vm%scalars(curr(5))%v * ( &
            -((vm%momenta(curr(6)) + vm%momenta(curr(8))) * &
               vm%vectors(curr(7))%v) * vm%momenta(curr(8)) - &
            (-(vm%momenta(curr(6)) + vm%momenta(curr(8))) * &
               vm%momenta(curr(8))) * vm%vectors(curr(7))%v)
    case (ovm_FUSE_G_GS)
       x = sgn_coupl_cmplx(vm, curr(2)) * vm%scalars(curr(5))%v * ( &
            -((vm%momenta(curr(6)) + vm%momenta(curr(8))) * &
               vm%vectors(curr(7))%v) * vm%momenta(curr(8)) - &
            (-(vm%momenta(curr(6)) + vm%momenta(curr(8))) * &
               vm%momenta(curr(8))) * vm%vectors(curr(7))%v)
    case (ovm_FUSE_G_SG_SKEW)
       x = - v_phiv(sgn_coupl_cmplx(vm, curr(2)), vm%scalars(curr(5))%v, vm%momenta(curr(6)), &
                    vm%momenta(curr(8)), vm%vectors(curr(7))%v)
    case (ovm_FUSE_G_GS_SKEW)
       x = - v_phiv(sgn_coupl_cmplx(vm, curr(2)), vm%scalars(curr(7))%v, vm%momenta(curr(6)), &
                    vm%momenta(curr(8)), vm%vectors(curr(5))%v)
    case default
       x = [zero, zero, zero, zero]
    end select
  end function gauge_sg

  elemental function sgn_coupl_cmplx(vm, j) result (s)
    class(vm_t), intent(in) :: vm
    integer, intent(in) :: j
    complex(default) :: s
    s = isign(1, j) * vm%coupl_cmplx(abs(j))
  end function sgn_coupl_cmplx

  elemental function sgn_coupl_cmplx2(vm, j, i) result (s)
    class(vm_t), intent(in) :: vm
    integer, intent(in) :: j, i
    complex(default) :: s
    if (i == 1) then
       s = isign(1, j) * vm%coupl_cmplx2(i, abs(j))
     else
       s = isign(1, j) * vm%coupl_cmplx2(i, abs(j))
    end if
  end function sgn_coupl_cmplx2

  elemental function int_to_log(i) result(yorn)
    integer, intent(in) :: i
    logical :: yorn
    if (i /= 0) then
      yorn = .true.
    else
      yorn = .false.
    end if
  end function

  elemental function color_factor(num, den, pwr) result (cf)
    integer, intent(in) :: num, den, pwr
    real(kind=default) :: cf
    if (pwr == 0) then
      cf = (one * num) / den
    else
      cf = (one * num) / den * (N_**pwr)
    end if
  end function color_factor

  elemental function vm_number_particles_in (vm) result (n)
    class(vm_t), intent(in) :: vm
    integer :: n
    n = vm%N_prt_in
  end function vm_number_particles_in

  elemental function vm_number_particles_out (vm) result (n)
    class(vm_t), intent(in) :: vm
    integer :: n
    n = vm%N_prt_out
  end function vm_number_particles_out

  elemental function vm_number_spin_states (vm) result (n)
    class(vm_t), intent(in) :: vm
    integer :: n
    n = vm%N_helicities
  end function vm_number_spin_states

  pure subroutine vm_spin_states (vm, a)
    class(vm_t), intent(in) :: vm
    integer, dimension(:,:), intent(out) :: a
    a = vm%table_spin
  end subroutine vm_spin_states

  elemental function vm_number_flavor_states (vm) result (n)
    class(vm_t), intent(in) :: vm
    integer :: n
    n = vm%N_flavors
  end function vm_number_flavor_states

  pure subroutine vm_flavor_states (vm, a)
    class(vm_t), intent(in) :: vm
    integer, dimension(:,:), intent(out) :: a
    a = vm%table_flavor
  end subroutine vm_flavor_states

  elemental function vm_number_color_indices (vm) result (n)
    class(vm_t), intent(in) :: vm
    integer :: n
    n = vm%N_col_indices
  end function vm_number_color_indices

  elemental function vm_number_color_flows (vm) result (n)
    class(vm_t), intent(in) :: vm
    integer :: n
    n = vm%N_col_flows
  end function vm_number_color_flows

  pure subroutine vm_color_flows (vm, a, g)
    class(vm_t), intent(in) :: vm
    integer, dimension(:,:,:), intent(out) :: a
    logical, dimension(:,:), intent(out) :: g
    a = vm%table_color_flows
    g = vm%table_ghost_flags
  end subroutine vm_color_flows

  elemental function vm_number_color_factors (vm) result (n)
    class(vm_t), intent(in) :: vm
    integer :: n
    n = vm%N_col_factors
  end function vm_number_color_factors

  pure subroutine vm_color_factors (vm, cf)
    class(vm_t), intent(in) :: vm
    type(OCF), dimension(:), intent(out) :: cf
    cf = vm%table_color_factors
  end subroutine vm_color_factors

  ! pure & ! pure unless OpenMp
  function vm_color_sum (vm, flv, hel) result (amp2)
    class(vm_t), intent(in) :: vm
    integer, intent(in) :: flv, hel
    real(default) :: amp2
    amp2 = ovm_color_sum (flv, hel, vm%table_amplitudes, vm%table_color_factors)
  end function vm_color_sum

  subroutine vm_new_event (vm, p)
    class(vm_t), intent(inout) :: vm
    real(default), dimension(0:3,*), intent(in) :: p
    logical :: mask_dirty
    integer :: hel
    call vm%run (p)
    if ((vm%hel_threshold .gt. 0) .and. (vm%hel_count .le. vm%hel_cutoff)) then
       call omega_update_helicity_selection (vm%hel_count, vm%table_amplitudes, &
         vm%hel_max_abs, vm%hel_sum_abs, vm%hel_is_allowed, vm%hel_threshold, &
         vm%hel_cutoff, mask_dirty)
       if (mask_dirty) then
          vm%hel_finite = 0
          do hel = 1, vm%N_helicities
             if (vm%hel_is_allowed(hel)) then
                vm%hel_finite = vm%hel_finite + 1
                vm%hel_map(vm%hel_finite) = hel
             end if
          end do
      end if
    end if
  end subroutine vm_new_event

  pure subroutine vm_reset_helicity_selection (vm, threshold, cutoff)
    class(vm_t), intent(inout) :: vm
    real(kind=default), intent(in) :: threshold
    integer, intent(in) :: cutoff
    integer :: i
    vm%hel_is_allowed = .True.
    vm%hel_max_abs = 0
    vm%hel_sum_abs = 0
    vm%hel_count = 0
    vm%hel_threshold = threshold
    vm%hel_cutoff = cutoff
    vm%hel_map = (/(i, i = 1, vm%N_helicities)/)
    vm%hel_finite = vm%N_helicities
  end subroutine vm_reset_helicity_selection

  pure function vm_is_allowed (vm, flv, hel, col) result (yorn)
    class(vm_t), intent(in) :: vm
    logical :: yorn
    integer, intent(in) :: flv, hel, col
    yorn = vm%table_flv_col_is_allowed(flv,col) .and. vm%hel_is_allowed(hel)
  end function vm_is_allowed

  pure function vm_get_amplitude (vm, flv, hel, col) result (amp_result)
    class(vm_t), intent(in) :: vm
    complex(kind=default) :: amp_result
    integer, intent(in) :: flv, hel, col
    amp_result = vm%table_amplitudes(flv, col, hel)
  end function vm_get_amplitude

  subroutine color_write (color, fh)
    class(color_t), intent(in) :: color
    integer, intent(in) :: fh
    select type(color)
    type is (col_discrete)
       write(fh, *) 'color_discrete%i           = ', color%i
    end select
  end subroutine color_write

  subroutine helicity_write (helicity, fh)
    class(helicity_t), intent(in) :: helicity
    integer, intent(in) :: fh
    select type(helicity)
    type is (hel_discrete)
       write(fh, *) 'helicity_discrete%i           = ', helicity%i
    type is (hel_trigonometric)
       write(fh, *) 'helicity_trigonometric%theta  = ', helicity%theta
    type is (hel_exponential)
       write(fh, *) 'helicity_exponential%phi      = ', helicity%phi
    type is (hel_spherical)
       write(fh, *) 'helicity_spherical%phi        = ', helicity%phi
       write(fh, *) 'helicity_spherical%theta      = ', helicity%theta
    end select
  end subroutine helicity_write

  subroutine flavor_write (flavor, fh)
    class(flavor_t), intent(in) :: flavor
    integer, intent(in) :: fh
    select type(flavor)
    type is (flv_discrete)
       write(fh, *) 'flavor_discrete%i           = ', flavor%i
    end select
  end subroutine flavor_write

end module omegavm95
