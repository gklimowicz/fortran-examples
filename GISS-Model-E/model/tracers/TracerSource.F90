module TracerSource_mod
  use TimeConstants_mod, only: INT_HOURS_PER_DAY
  use constant, only : undef
  implicit none
  private

  public :: TracerSource
  public :: TracerSource3D
  public :: N_MAX_SECT


!@param n_max_sect maximum number of sectors for emissions altering
  integer, parameter :: N_MAX_SECT = 10

  type TracerSource
    integer :: num_tr_sectors = 0! number of sectors for a particular tracer and source
    integer :: tr_sect_index(N_MAX_SECT)  = 0 ! array hold the sector index for given tracer/source
    logical :: applyDiurnalCycle = .false. ! whether or not to apply the factors in the diurnalCycle array
    character(len=10) :: tr_sect_name(N_MAX_SECT) = ' '! array hold the sector name for given tracer/source
    real*8 :: diurnalCycle(INT_HOURS_PER_DAY) = undef ! an optional diurnal cycle to be read from a file
  end type TracerSource

  type, extends(TracerSource) :: TracerSource3D
  end type TracerSource3D

end module TracerSource_mod
