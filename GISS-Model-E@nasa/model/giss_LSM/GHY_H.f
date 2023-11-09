!@sum header file for TerraE Global Land Model
      module ghy_h
!@sum module with main parameters for GHY
      integer, parameter, public :: ngm=6, imt=5, nlsn=3
      integer, parameter, public :: LS_NFRAC=3

      type ghy_tr_str
!@var trpr flux of tracers in precipitation (?/m^2 s)
!@var trdd flux of tracers as dry deposit (?/m^2 s)
!@var trirrig flux of tracers in irrigation (?/m^2 s)
!@var tr_w amount of tracers in the soil (?)
!@var tr_wsn amount of tracers in the snow (?)
        integer :: ntg  ! number of GHY tracers
        integer, pointer :: ntixw(:) ! indices of water tracers
        logical, pointer :: is_water(:)
        ! in
        real*8, pointer :: trpr(:), trdd(:), trirrig(:),  tr_surf(:)
        ! inout
        real*8, pointer :: tr_w(:,:,:), tr_wsn(:,:,:)
        ! out
        real*8, pointer :: atr_evap(:),atr_rnff(:),atr_g(:)
      end type ghy_tr_str

      end module ghy_h
