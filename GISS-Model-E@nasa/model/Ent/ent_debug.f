      module ent_debug_mod
      implicit none

      integer, parameter :: NPFT=16
      integer, parameter :: SIZE_ENT_DEBUG=NPFT*(16+11 +1)

      type ent_debug
        sequence
        real*8 :: vf(NPFT)      ! 1
        real*8 :: Anet(NPFT)    ! 2
        real*8 :: Atot(NPFT)    ! 3
        real*8 :: Rd(NPFT)      ! 4
        real*8 :: GCANOPY(NPFT) ! 5
        real*8 :: TRANS_SW(NPFT)! 6
        real*8 :: LAI(NPFT)     ! 7
        real*8 :: Resp_fol(NPFT)! 8
        real*8 :: Resp_sw(NPFT) ! 9
        real*8 :: Resp_lab(NPFT) ! 10 
        real*8 :: Resp_root(NPFT)! 11
        real*8 :: Resp_maint(NPFT) ! 12
        real*8 :: Resp_growth_1(NPFT) ! 13
        real*8 :: Resp_growth(NPFT) ! 14
        real*8 :: GPP(NPFT)     ! 15
        real*8 :: R_auto(NPFT)  !16

        real*8 :: total(NPFT)  ! 17
        real*8 :: C_lab(NPFT)  ! 18
        real*8 :: C_fol(NPFT)  ! 19
        real*8 :: C_sw(NPFT)  ! 20
        real*8 :: C_hw(NPFT)  ! 21
        real*8 :: C_froot(NPFT)  ! 22
        real*8 :: C_croot(NPFT)  ! 23
        real*8 :: C_soil(NPFT)  ! 24
        real*8 :: Resp_soil(NPFT)  ! 25
        real*8 :: phenofactor(NPFT) ! 26
        real*8 :: betad(NPFT) ! 27

        ! cell vars................. 28
        real*8 :: soiltemp_10d ! 1
        real*8 :: airtemp_10d  ! 2
        real*8 :: par_10d      ! 3
        real*8 :: gdd          ! 4
        real*8 :: ncd          ! 5
        real*8 :: sgdd         ! 6
        real*8 :: ld           ! 7
        real*8 :: ld0          ! 8
        real*8 :: foo(NPFT-8)

       end type ent_debug

      type(ent_debug), target :: ent_d
      !real*8 :: ent_dl(SIZE_ENT_DEBUG)
      !equivalence(ent_d,ent_dl)

      contains

      subroutine get_ent_debug_ptr(ptr)
      use ISO_C_BINDING
      real*8, pointer :: ptr(:)
      !---
      type (C_PTR) :: cptr

      cptr = c_loc(ent_d)
      call c_f_pointer( cptr, ptr, (/SIZE_ENT_DEBUG/) )

      end subroutine get_ent_debug_ptr

      end module ent_debug_mod
