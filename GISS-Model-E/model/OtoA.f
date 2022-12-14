      subroutine OAtoAA_WO(oA,aA)
!@sum Regridding from ocean to atmosphere for several types of grid. Different cases:
!@+   -ocean and atmosphere use latlon grids with same resolution, 
!@+   -ocean and atmo use latlon grids with different resolution,
!@+   -ocean uses latlon grid, atmosphere uses cubed-sphere grid
!@auth Gary Russell & Denis Gueyffier
      use GEOM, only : ima=>IM, jma=>JM
#ifndef CUBED_SPHERE
     &     ,adlatm=>DLATM
#endif
      use OCEAN, only : imo=>IM, jmo=>JM,odlatm=>DLATM,
     &     odxyp=>DXYPO,ofocean=>FOCEAN
      use DOMAIN_DECOMP_1D , only : ocean_pack=>pack_data,
     &     ocean_unpack=>unpack_data, am_i_root  ! ocean always uses 1d decomposition
      use DOMAIN_DECOMP_ATM, only: agrid=>grid,atm_pack=>pack_data,
     &     atm_unpack=>unpack_data  ! atm uses either 1d or 2d decomposition
      use oceanr_dim, only : ogrid
      use regrid_com, only : xO2A
      implicit none
      real*8, intent(inout) ::
     &     oA(imo,oGRID%J_STRT_HALO:oGRID%J_STOP_HALO),    
     &     aA(agrid%i_strt_halo:agrid%i_stop_halo,
     &     agrid%j_strt_halo:agrid%j_stop_halo)
      real*8, allocatable :: ofocean_loc(:,:)
      real*8 :: oAglob(imo,jmo),aAglob(ima,jma,agrid%ntiles)     ! ntiles = 1 : latlon grid, 6 : cubed-sphere
      real*8 :: aArea(ima,jma,agrid%ntiles) ! ntiles = 1 : latlon grid, 6 : cubed-sphere
      real*8 :: DENOM
      Integer*4 J1O,JNO, J1A,JNA,J1P,JNP, I,J
C**** 
      
      if (agrid%ntiles .eq. 1) then
         
         if (IMA==IMO .and. JMA==JMO)  GoTo 100
         if (IMA==144 .and. JMA== 90 .and. IMO==288 .and. JMO==180)
     *        GoTo 200
C**** 
C**** Non-standard regular latitude-longitude resolution
C**** 
         call ocean_pack (oGRID,oA,oAGLOB)
         if (AM_I_ROOT())  Then
#ifndef CUBED_SPHERE
            Call HNTR80 (IMO,JMO,0d0,oDLATM, IMA,JMA,0d0,aDLATM, 0d0)
#endif
            Call HNTR8P (oFOCEAN,oAGLOB,aAGLOB) 
         endIf
         call atm_unpack (aGRID,aAGLOB,aA)
         return
C**** 
C**** IMA=IMO and JMA=JMO
C**** 
 100     aA(:,:) = oA(:,:)
         return
C**** 
C**** IMA=144, JMA=90, IMO=288 and JMO=180
C**** 
 200     J1A = aGRID%J_STRT  ;  JNA = aGRID%J_STOP
         J1P = Max (J1A,2)   ;  JNP = Min (JNA,JMA-1)
         do J=J1P,JNP
            do I=1,IMA
               DENOM = oDXYP(J*2-1)*(oFOCEAN(I*2-1,J*2-1)
     +              +oFOCEAN(I*2,J*2-1))
     +              + oDXYP(J*2)  *(oFOCEAN(I*2-1,J*2)  
     +              +oFOCEAN(I*2,J*2))
               aA(I,J) = 0.
               If (DENOM .eq. 0.)  cycle
               aA(I,J) = (oDXYP(J*2-1)*(oFOCEAN(I*2-1,J*2-1)
     *              *oA(I*2-1,J*2-1) +
     +              oFOCEAN(I*2  ,J*2-1)*oA(I*2  ,J*2-1))
     +              + oDXYP(J*2)  *(oFOCEAN(I*2-1,J*2)
     *              *oA(I*2-1,J*2) +
     +              oFOCEAN(I*2  ,J*2)*oA(I*2,J*2))) / DENOM
            enddo
         enddo
C**** South Pole cell
         If (J1A==1)  Then
            DENOM = oDXYP(1)*IMO*oFOCEAN(1,1) 
     +           + oDXYP(2)*Sum(oFOCEAN(:,2))
            aA(I,J) = 0
            If (DENOM == 0)
     *           aA(1,1) = (oDXYP(1)*IMO*oFOCEAN(1,1)*oA(1,1) +
     +           oDXYP(2)*Sum(oFOCEAN(:,2)*oA(:,2))) / DENOM
            aA(2:,1) = aA(1,1)    
         endif
C**** North Pole cell
         if (JNA==JMA)  Then
            DENOM = oDXYP(JMO)*IMO*oFOCEAN(1,JMO) +
     +           oDXYP(JMO-1)*Sum(oFOCEAN(:,JMO-1))
            aA(1,JMA) = 0
#ifndef CUBED_SPHERE
            if (DENOM == 0)
     *           aA(1,JMA) = (oDXYP(JMO)*IMO
     *           *oFOCEAN(1,JMO)*oA(1,JMO) +
     +           oDXYP(JMO-1)*Sum(oFOCEAN(1,JMO-1)*oA(:,JMO-1))) *
     *           aDXYP(JMA)
#endif
            aA(2:,JMA) = aA(1,JMA)  
         endif
         return
                     
C**** 
C**** ATMO = 'CUBE-SPHERE'
C**** 
      else if (agrid%ntiles .eq. 6) then
         allocate(
     &        ofocean_loc(imo, oGRID%J_STRT_HALO:oGRID%J_STOP_HALO)
     &        )

         call ocean_unpack(ogrid,oFOCEAN,ofocean_loc)

c*       regridding local ocean array to global atm. array
         call parallel_regrid_wt(xO2A,oA,aAglob,aArea,ofocean_loc) 

         call atm_unpack(aGrid,aAglob,aA)

         deallocate(ofocean_loc)
      endif
      
      end subroutine OAtoAA_WO
