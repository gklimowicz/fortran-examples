module Hidden_mod
  private
  public :: Hidden_type

  type Hidden_type
    integer :: mpi_comm
    integer :: numProcesses
    integer :: numAllProcesses
    integer :: mpi_tag
    
    logical :: hasSouthPole
    logical :: hasNorthPole
    logical :: hasEquator
    
    logical :: periodicBC
    integer :: log_unit
    
!@var lookup_pet index of PET for a given J
    integer, dimension(:), pointer :: lookup_pet
    
  end type Hidden_type

end module Hidden_mod
      
