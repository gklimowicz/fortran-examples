      module reproduction
!@sum Routines to calculate reproduction rates in an entcell.

      use ent_types

      implicit none


      contains
      !*********************************************************************

      subroutine reproduction_calc(dtsec, time, pptr)
      real*8 :: dtsec           !dt in seconds
      type(timestruct) :: time  !Greenwich Mean Time
      type(patch),pointer :: pptr

      !---------------------------------------------------------------
      !              FILL IN CODE                                    
      !---------------------------------------------------------------
      !Not needed in GISS replication test.

      end subroutine reproduction_calc

      !*********************************************************************
      
      end module reproduction
