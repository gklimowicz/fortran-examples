! WHIZARD <<Version>> <<Date>>

! dummy replacement routine for F77 LoopTools wrapper
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! call routine C0 to calculate the 3-point integral J0
      double complex function J0_LoopTools (p2, k2, q2, m2)
      double complex p2, k2, q2, m2
      write (0, "(A)")  "***********************************************
     +*****************"
      write (0, "(A)")  "*** LoopTools: Error: library not linked, WHIZA
     +RD terminates ***"
      write (0, "(A)")  "***********************************************
     +*****************"
      end
