      subroutine OAtoAA(oW,oA,aA)
      use GEOM, only : ima=>IM, jma=>JM
      use OCEAN, only : imo=>IM,jmo=>JM
      use oceanr_dim, only : oGRID
      use DOMAIN_DECOMP_ATM, only: agrid=>grid,atm_unpack=>unpack_data
      use regrid_com
      implicit none
      character*80 :: TITLE,name
      real*8, intent(inout) ::
     &     oW(imo,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),
     &     oA(imo,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO), 
     &     aA(agrid%i_strt_halo:agrid%i_stop_halo,
     &     agrid%j_strt_halo:agrid%j_stop_halo)

      real*8, allocatable :: aAglob(:,:,:)
      real*4, allocatable :: aAglob4(:,:,:)
      real*8, allocatable :: aArea(:,:,:)
      real*8 :: missing

      allocate( aAglob(ima,jma,agrid%ntiles),
     &     aAglob4(ima,jma,agrid%ntiles),
     &     aArea(ima,jma,agrid%ntiles) )

      missing=-1.e30

      call repr_regrid_wt(xO2A,oW,missing,oA,aAglob,aArea)   
      
c      if (am_i_root()) then
c         name="regridO2A_repr"
c         write(*,*) name
c         open(20, FILE=name,FORM='unformatted', STATUS='unknown')
c         aAglob4=aAglob
c         TITLE="OTOA reproducible regridding"
c         write(unit=20) TITLE, aAglob4
c         close(20)
c      endif
      
      call atm_unpack(agrid,aAglob,aA)
      
      deallocate(aAglob4,aAglob,aArea)

      end subroutine OAtoAA
