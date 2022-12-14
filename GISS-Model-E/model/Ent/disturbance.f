      module disturbance

      use ent_types

      implicit none


      contains
      !*********************************************************************

      subroutine fire_frequency_cell(dtsec,time, entcell)
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: time  !Greenwich Mean Time
      type(entcelltype) :: entcell

      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------


      end subroutine fire_frequency_cell

      !*********************************************************************
      subroutine calc_cell_disturbance_rates(dtsec,time,entcell)
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: time  !Greenwich Mean Time
      type(entcelltype) :: entcell
      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------

      end subroutine calc_cell_disturbance_rates


      !*********************************************************************

      end module disturbance
