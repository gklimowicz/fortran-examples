#include "rundeck_opts.h"

      subroutine OCN_mesosc(kappam3d)

      USE MODEL_COM,  only : nstep=>itime,itimei,nday
     . ,iyear1,aMON,dtsrc,xlabel,lrunid
      use JulianCalendar_mod, only: jdendofm
      USE CONSTANT,   only : grav,omega
      USE OCEANR_DIM, only : ogrid
      USE OCEANRES,   only : idm=>imo,jdm=>jmo,kdm=>lmo,dzo
      USE OFLUXES,    only : oRSI,oAPRESS
      USE OCEAN,      only : ZOE=>ZE,g0m,s0m,mo,dxypo,focean,lmm
     .                      ,oLON_DG,oLAT_DG,uo,vo,sinpo,im,dxpo,dypo
      USE KPP_COM,    only : kpl

      USE GM_COM, only: RHOX, RHOY
#ifdef OCN_Mesoscales
      USE ODIAG, only: oij=>oij_loc,oijl=>oijl_loc
     .            ,ij_eke,ij_rd,ijl_ueddy,ijl_veddy,ijl_n2
      USE OCEAN,      only : auvel,avvel
#endif

      USE DOMAIN_DECOMP_1D, only: AM_I_ROOT
     ., HALO_UPDATE, NORTH, SOUTH
      use TimerPackage_mod


      implicit none

      integer, parameter :: itest=97,jtest=128
      integer i,j,k,l,ndepi
      integer i_0,i_1,j_0,j_1

      Real*8,External   :: VOLGSP,TEMGSP,TEMGS
      Real*8,External   :: VOLGS
      real*8  g,s,pres
      real*8   p1d(kdm+1),dp1d(kdm),temp1d(kdm),saln1d(kdm)
      real*8   rho_water,amld_cgs,Rd
      real*8   K0,Ustar(kdm),Vstar(kdm)
      real*8   z_cm(kdm+1),dens_cgs(kdm)
     .        ,uvel_cgs(kdm),vvel_cgs(kdm)
     .        ,drhodz_cgs(kdm),coriol,n2(kdm)
     .        ,drhodx_cgs(kdm),drhody_cgs(kdm)
      real*8  drhopdz_cgs(kdm),densulave_cgs(kdm)
      REAL*8,ALLOCATABLE :: presi(:),densu(:),densl(:),temp3d(:,:,:)
      REAL*8, DIMENSION(idm,ogrid%j_strt_halo:ogrid%j_stop_halo,kdm) ::
     .   kappam3d
      REAL*8 wta
      integer ip1,im1

      logical vrbos
c
      kappam3d=0. 

      if (nstep.eq.0) return


      i_0=ogrid%I_STRT
      i_1=ogrid%I_STOP
      j_0=ogrid%J_STRT
      j_1=ogrid%J_STOP

      CALL HALO_UPDATE(ogrid,
     *                 RHOX(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
     *                 FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(ogrid,
     *                 RHOY(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
     *                 FROM=SOUTH+NORTH)
      CALL HALO_UPDATE(ogrid,
     *                 vo(:,ogrid%j_strt_halo:ogrid%j_stop_halo,:),
     *                 FROM=SOUTH+NORTH)

      do 1000 j=j_0,j_1
      do 1000 i=i_0,i_1
      IF(FOCEAN(I,J).gt.0.) THEN

      im1=i-1
      IF(i.eq.1) im1=idm

      vrbos=.false.
c     vrbos=.true.
      if (i.eq.itest.and.j.eq.jtest) vrbos=.true.

! max depth cell
      ndepi=lmm(i,j)

C-- TONY - 03/07/11
C-- Calculate an interface pressure array
      ALLOCATE(presi(0:lmm(i,j)),densu(lmm(i,j)),densl(lmm(i,j)))
      presi(0) = 0.
      DO k=1,lmm(i,j)
        presi(k) = presi(k-1) + MO(i,j,k)*GRAV
      ENDDO
C-- Calculate a pot. density referenced to the upper and lower interfaces
      DO k=1,lmm(i,j)
        g=G0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
        s=S0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
        densu(k) = 1d0/VOLGSP(g,s,presi(k-1)) 
        densl(k) = 1d0/VOLGSP(g,s,presi(k))
c       densulave_cgs(k) = 0.5*(densu(k)+densl(k))*1.d-3
c       densu(k) =1d0/volgs(g,s)
      ENDDO
C--

!change units
       pres = oAPRESS(i,j)    !surface atm. pressure
       do k=1,lmm(i,j)
         pres=pres+MO(I,J,k)*GRAV*.5
         g=G0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
         s=S0M(I,J,k)/(MO(I,J,k)*DXYPO(J))
         temp1d(k)=TEMGSP(g,s,pres)    !in situ   temperature
         saln1d(k)=s*1000.             !convert to psu (eg. ocean mean salinity=35psu)
         rho_water = 1d0/VOLGSP(g,s,pres)
         dp1d(k)=MO(I,J,K)/rho_water   !local thickenss of each layer in meters
         !add missing part of density to get to the bottom of the layer
         !now pres is at the bottom of the layer
         pres=pres+MO(I,J,k)*GRAV*.5
      enddo

C-- convert in CGS
      p1d(1)=0.
      do k=2,kdm+1
      p1d(k)=p1d(k-1)+dp1d(k-1)
      enddo
      do k=1,kdm
       z_cm(k) = p1d(k)*100.d0                !depth in cm
       dens_cgs(k) = MO(i,j,k)/max(1.d0,dp1d(k)) *1.d-3 
       uvel_cgs(k) = uo(i,j,k) * 100.d0
       vvel_cgs(k) = vo(i,j,k) * 100.d0
       IF(uo(i,j,k).ne.0..AND.uo(im1,j,k).ne.0.)
     *   uvel_cgs(k) = 0.5*(uo(i,j,k)+uo(im1,j,k)) * 100.d0
       IF(vo(i,j,k).ne.0..AND.vo(i,j-1,k).ne.0.)
     *   vvel_cgs(k) = 0.5*(vo(i,j,k)+vo(i,j-1,k)) * 100.d0
      enddo

c     IF(j.eq.90) WRITE(*,*)'uvel:',i,uvel_cgs(1)

      z_cm(kdm+1) = p1d(kdm+1)*100.d0                !depth in cm
      coriol = 2d0*omega*sinpo(j)      !evaluated at tracer points
!drhodz_cgs
      do k=1,kdm-1
         if (dens_cgs(k+1).ne.0.d0.and.dens_cgs(k).ne.0.d0) then
           drhodz_cgs(k) = (dens_cgs(k+1) - dens_cgs(k))          ! at (i,j,k+1/2)   interface
c    .                 /(max(1.d0,z_cm(k+1)+z_cm(k)))/2.
     .                 /(max(1.d0,z_cm(k+1)-z_cm(k)))
         elseif (dens_cgs(k+1).eq.0.d0.
     *       and.dens_cgs(k).ne.0.d0.
     *       and.dens_cgs(max(1,k-1)).ne.0.d0) then
           drhodz_cgs(k)=drhodz_cgs(max(1,k-1))
         else
           drhodz_cgs(k) = 1.d0
         endif
      enddo
      do k=1,kdm
         drhodx_cgs(k)=rhox(i,j,k)*1.d-3/100.d0
         drhody_cgs(k)=rhoy(i,j,k)*1.d-3/100.d0
         IF(rhox(i,j,k).ne.0..AND.rhox(im1,j,k).ne.0.)
     *     drhodx_cgs(k)=0.5*(rhox(i,j,k)+rhox(im1,j,k))*1.d-3/100.d0
         IF(rhoy(i,j,k).ne.0..AND.rhoy(i,j-1,k).ne.0.)
     *     drhody_cgs(k)=0.5*(rhoy(i,j,k)+rhoy(i,j-1,k))*1.d-3/100.d0
      enddo
C-- TONY - 03/07/11
C-- Calculate a new drhodz using the potential densities
      drhopdz_cgs = 1.d30
      DO k=1,lmm(i,j)-1
        drhopdz_cgs(k)=(1.d-3)*(densu(k+1)-densl(k))/(z_cm(k+1)-z_cm(k))
      ENDDO
      drhopdz_cgs(lmm(i,j)) = drhopdz_cgs(max(1,lmm(i,j)-1))
C
! compute Brunt-Vaisala
      do k=1,kdm-1
c       if (drhodz_cgs(k).ne.1d30) then
c       n2(k) = GRAV* drhodz_cgs(k)*100.d0     !no minus sign here, because drhodz defined as rho(k+1)-rho(k)
c    .            /(max(1.d0,dens_cgs(k+1)+dens_cgs(k))/2.)   !at (i,j,k+1/2)   interface
C-- TONY - 03/07/11 - Updated n2 using the pot density instead of the in situ one
        if (drhopdz_cgs(k+1).ne.1d30) then
        n2(k) = GRAV*100.d0 * drhopdz_cgs(k)
     *         /(max(1.d0,(1.d-3)*(densu(k+1)+densl(k)))/2.)
        else
        n2(k) = 1.d30
        endif
        IF(k.eq.lmm(i,j)) THEN
        IF(lmm(i,j).eq.2) THEN
          n2(k)=n2(k-1)
        ELSE
          n2(k)=n2(k-1)+((n2(k-1)-n2(k-2))/(z_cm(k-1)-z_cm(k-2)))
     *          *(z_cm(k)-z_cm(k-1))
        ENDIF
        ENDIF
      enddo
      DEALLOCATE(presi,densu,densl)
! mixed layer depth
      amld_cgs = 0.d0
      do k=1,kpl(i,j)
      if (dp1d(k).ne.1.d30) then
          amld_cgs = amld_cgs + dp1d(k)*100.d0      ! at (i,j,k) the middle of last layer in MLD
      endif
      enddo
      if (amld_cgs.gt.p1d(ndepi)*100.d0) amld_cgs=p1d(ndepi)*100.d0

#ifdef OCN_Mesoscales
C-- time-averaging of u,v
      wta=exp(-1./240.)
      auvel(i,j,:) = wta*auvel(i,j,:) + (1.-wta)*uvel_cgs(:)
      avvel(i,j,:) = wta*avvel(i,j,:) + (1.-wta)*vvel_cgs(:)

      uvel_cgs=auvel(i,j,:)
      vvel_cgs=avvel(i,j,:)
#endif

      call mesoscales1d(kdm,ndepi,
     .      dens_cgs,uvel_cgs,vvel_cgs,
     .      n2,drhodx_cgs,drhody_cgs,drhopdz_cgs,coriol,amld_cgs
     .     ,Rd,K0,Ustar,Vstar,i,j,kappam3d(i,j,:))


#ifdef OCN_Mesoscales
      if (vrbos) then
      write(*,'(a,2i6)')'MESOSCALES:',nstep
      endif


       OIJ(I,J,IJ_eke)  = OIJ(I,J,IJ_eke) + K0      ! eddy kinetic energy, ocean
       OIJ(I,J,IJ_rd )  = OIJ(I,J,IJ_rd ) + Rd      ! Rossby radius of deformation
C
       DO k=1,kdm
          OIJL(I,J,k,IJL_n2   )= OIJL(I,J,k,IJL_n2   ) + n2(k)    ! brunt vaisala squared
          OIJL(I,J,k,IJL_ueddy)= OIJL(I,J,k,IJL_ueddy) + ustar(k) ! ustar, eddy induced velocity (Canuto)
          OIJL(I,J,k,IJL_veddy)= OIJL(I,J,k,IJL_veddy) + vstar(k) ! vstar, eddy induced velocity (Canuto)
       ENDDO
#endif

      endif    !focean

 1000 continue

C--   Convert kappam3d from CGS to MKS
      kappam3d=kappam3d*1.d-4
C--   Maximum of 15000 m2/s
C--   Minimum value of 1 m2/s
      DO k=1,kdm
       DO j=j_0,j_1
         DO i=i_0,i_1
           kappam3d(i,j,k)=min(kappam3d(i,j,k),15000.d0)
           kappam3d(i,j,k)=max(kappam3d(i,j,k),1.d0)
         ENDDO
       ENDDO
      ENDDO
C--


      end subroutine OCN_mesosc

      subroutine mesoscales1d(km,ndepi,
     .      rhoi,UI,VI,n2,dxrho,dyrho,dzrho,f,ml
     .     ,rdm,K02,Ustar,Vstar,
     .     igrid,jgrid,KAPPAMZa) 

        USE MODEL_COM,  only : nstep=>itime
        USE ODIAG, only : zoc
        USE OCEAN_DYN, Only : DH

        IMPLICIT NONE

        real*8, parameter :: a02=0.03
        INTEGER k,km,kmli,igrid,jgrid,ndepi
        REAL*8 z_cm(km),zma(km),zh(km)
        REAL*8 Ustar(km),Vstar(km)
        REAL*8 KAPPAM,KAPPAMZa(km)
        REAL*8 SIGMAT,f,ml,yt,rd,rdm,pi
        PARAMETER(SIGMAT=1.)
        INTEGER kkpmax2,kkpint,ktap,kkpint2,ktap2
        REAL*8 DZLXB1AVE,DZLYB1AVE,kpmax2
        REAL*8 UI(km),VI(km)
        REAL*8 UB1AVE,VB1AVE,UL,VL
        REAL*8 rhoi(km)
        REAL*8 frd,n2(km),lr
        REAL*8 K02,CK,K02A,K02D,CA,CD,INTGAMMAA,INTGAMMAD,INTGAMMA,
     *    INTGAMMA1,INTEXP,CD1,INTGAMMAS
        REAL*8 dxrho(km),dyrho(km),dzrho(km)
        REAL*8, ALLOCATABLE :: n2t(:),n2t1(:),z(:),zm(:),b1(:),b1t(:),
     *   U(:),V(:),rho(:),n2tr(:),
     *   LX(:),LY(:),DZLX(:),DZLY(:),
     *   UINT(:),VINT(:),UT(:),VT(:),
     *   dxb(:),dyb(:),KP(:),KINT(:),
     *   zint(:),UREV(:),VREV(:),KPINT(:),
     *   b1rev(:),DZLXREV(:),DZLYREV(:),
     *   UTINT(:),VTINT(:),DZU(:),DZV(:),
     *   F1X(:),F2X(:),F1Y(:),F2Y(:),KP2(:),KINT2(:),KP2INT(:),
     *   KAPPAMZ(:),GAMMA(:),n2sm(:),Tap(:),Adrag(:),DZU2(:)
        INTEGER im,kmi
        REAL*8 Nm
        REAL*8 zw(km),zwt(km),zt(km)
        REAL*8 k02count
        REAL*8, ALLOCATABLE :: ud(:),vd(:),MK(:)
 
      pi=4.*atan(1.)

      zh(1)=DH(igrid,jgrid,1)/2.
      zma(1)=DH(igrid,jgrid,1)
      DO k=2,km-1
        zh(k)=zma(k-1)+DH(igrid,jgrid,k)/2.
        zma(k)=zma(k-1)+DH(igrid,jgrid,k)
      ENDDO
      zh(km)=zoc(km)

      z_cm=zh*100.

      if (ndepi.eq.0) then
        K02=0.
        KAPPAMZa=0.
        goto 10
      endif

      do k=1,km
         zw(k) = z_cm(k)
         zt(k) =-z_cm(k)
      enddo
      zwt(1)=0.
      do k=2,km
        zwt(k) = zw(k-1)
      enddo

      kmi = ndepi

      kmli=km+1   !outside bounds
      DO k=1,kmi-1
        IF(zw(k).le.ml.and.ml.lt.zw(k+1)) kmli=k+1
      ENDDO
C
      kmli=min(kmli,kmi)
C
      ALLOCATE(n2t(kmi),n2t1(kmi),z(kmi),zm(kmi),b1(kmi),b1t(kmi),
     * U(kmi),V(kmi),rho(kmi),n2tr(kmi),
     * LX(kmi),LY(kmi),DZLX(kmi),DZLY(kmi),
     * UINT(kmi),VINT(kmi),UT(kmi),VT(kmi),
     * dxb(kmi),dyb(kmi),KP(kmi),
     * zint(kmi),UREV(kmi),VREV(kmi),KPINT(kmi),
     * b1rev(kmi),DZLXREV(kmi),DZLYREV(kmi),
     * UTINT(kmi),VTINT(kmi),DZU(kmi),DZV(kmi),
     * F1X(kmi),F1Y(kmi),F2X(kmi),F2Y(kmi),KP2(kmi),KP2INT(kmi))
C-- M1: Calculate B1, frd and rd
      n2t(1)=n2(1)
C     n2t(1)=max(n2t(1),1.d-8)
      DO k=2,kmi
        n2t(k)=n2(k-1)
C       n2t(k)=max(n2t(k),1.d-8)
      ENDDO
      DO k=1,kmi-1
        n2tr(k)=0.5*(n2(k)+n2(k+1))
      ENDDO
      n2tr(kmi)=n2(kmi-1)
C
      z=-zw(1:kmi)
      zm=-z
C     n2t=((7.d-3)**2)*exp(2.*z/100000.)
c     n2t1=n2t
      n2t1=n2(:)
C
c     CALL baroclin1st(zm,n2t1,kmi,b1,im,Nm,frd)
      CALL baroclin1st(zm,n2(:),kmi,b1,im,Nm,frd,igrid,jgrid)
      rd=dabs(frd/f)
c     rdm=rd
      DO k=1,kmi-1
        b1t(k)=0.5*(b1(k)+b1(k+1))
      ENDDO
      b1t(kmi)=b1(kmi)
      b1=b1t
      U=UI(1:kmi)
      V=VI(1:kmi)
      rho=rhoi(1:kmi)
c
      CALL d1sym(kmi,z,U,DZU)
      CALL d1sym(kmi,z,V,DZV)
c
      ALLOCATE(n2sm(kmi))
      n2sm = n2(1:kmi)
      DO k=1,50
      CALL z121_mesosc(n2sm,kmi,kmi)
      ENDDO
      do k=1,kmi
c        if (n2tr(k).ne.0.d0)then 
c               LX(k)=-(f/n2tr(k))*DZV(k)
c               LY(k)=+(f/n2tr(k))*DZU(k)
         if (n2sm(k).ne.0.d0)then
                LX(k)=-(f/n2sm(k))*DZV(k)
                LY(k)=+(f/n2sm(k))*DZU(k)
c               LX(k)=-dxrho(k)/(-n2tr(k)*rho(k)/981.)
c               LY(k)=-dyrho(k)/(-n2tr(k)*rho(k)/981.)
         else
                LX(k)=0.d0
                LY(k)=0.d0
         endif
C-- TONY 03/08/11 -- limiter for isopycnal slopes
c        LX(k)=min(1.d-3,LX(k))
c        LY(k)=min(1.d-3,LY(k))
      enddo
C
      UREV=0.
      VREV=0.
      DO k=1,kmi
        zint(k)=zt(kmi-k+1)
        UREV(k)=U(kmi-k+1)
        VREV(k)=V(kmi-k+1)
        b1rev(k)=b1(kmi-k+1)
      ENDDO
C
      z=zt(1:kmi)
      CALL B1AVERAGE(UREV,zint,b1rev,kmi,UB1AVE,kmli+0)
      CALL B1AVERAGE(VREV,zint,b1rev,kmi,VB1AVE,kmli+0)
C
      CALL d1sym(kmi,z,LX,DZLX)
      CALL d1sym(kmi,z,LY,DZLY)

      DO k=1,kmi
        DZLXREV(k)=DZLX(kmi-k+1)
        DZLYREV(k)=DZLY(kmi-k+1)
      ENDDO

      CALL B1AVERAGE(DZLXREV,zint,b1rev,kmi,DZLXB1AVE,kmli+0)
      CALL B1AVERAGE(DZLYREV,zint,b1rev,kmi,DZLYB1AVE,kmli+0)
C--
C-- Calculate UT and VT
      DO k=1,kmi
        CALL dqtfg(zint(k:kmi),UREV(k:kmi),UINT(k:kmi),kmi-k+1)
        CALL dqtfg(zint(k:kmi),VREV(k:kmi),VINT(k:kmi),kmi-k+1)
        if (zint(k).ne.0.d0) then
          UTINT(k)=-(1./zint(k))*UINT(kmi)
          VTINT(k)=-(1./zint(k))*VINT(kmi)
        else
          UTINT(k)=0.d0
          VTINT(k)=0.d0
        endif
      ENDDO
      UTINT(kmi)=0.
      VTINT(kmi)=0.
      DO k=1,kmi
        UT(k)=UTINT(kmi-k+1)
        VT(k)=VTINT(kmi-k+1)
      ENDDO
c
C-- Calulate the Rhines Scale and replace rd by lr between +/-30 degrees
      lr=dsqrt(dsqrt(UT(kmli)**2+VT(kmli)**2)/(2.*2.1*1.d-13))
      rd=min(rd,lr)
      IF(jgrid.ge.70.AND.jgrid.le.110.AND.rd.eq.0.) rd=lr
      rdm=rd
 
C-- Calculate F1
      if (rd.ne.0.d0) then                !Natassa
         F1X = DZLXB1AVE + (1.+1./SIGMAT)*(1./f)*(1./rd**2)*(-VB1AVE)
      else
         F1X = 0.d0
      endif

      if (rd.ne.0.d0) then                !Natassa
          F1Y = DZLYB1AVE + (1.+1./SIGMAT)*(1./f)*(1./rd**2)*(+UB1AVE)
      else
          F1Y = 0.d0
      endif

C-- Calculate F2
      if (rd.ne.0.d0) then                !Natassa
      F2X = (1./f)*(1./rd**2)*((1./SIGMAT)*(+VT)+V)
      F2Y = (1./f)*(1./rd**2)*((1./SIGMAT)*(-UT)-U)
      else
      F2X = 0.d0
      F2Y = 0.d0
      endif

C-- Calculate DXB and DYB
      dxb=-981.*dxrho(1:kmi)/rho(1:kmi)
      dyb=-981.*dyrho(1:kmi)/rho(1:kmi)

C-- Calculate Tapering function T(z)
      ktap=1
      ALLOCATE(Tap(kmi))
      Tap=0.
      DO k=1,kmi
        IF(((LX(k)**2+LY(k)**2)*n2tr(k)).ne.0.)
     *  Tap(k)=-z(k)*((F1X(k)+F2X(k))*dxb(k)+(F1Y(k)+F2Y(k))*dyb(k))
     *      /((LX(k)**2+LY(k)**2)*n2tr(k))
        IF(Tap(k).gt.1.) THEN
          ktap=k
          GO TO 116
        ENDIF
      ENDDO
 116  CONTINUE
      DEALLOCATE(Tap)

C-- Calculate GAMMA
      ALLOCATE(GAMMA(kmi))
      GAMMA=(a02+(dabs(b1))**2)/(1.+a02)

C-- Calculate new K
      CK=27.
      KP2=sqrt(GAMMA)*z*((F1X+F2X)*dxb+(F1Y+F2Y)*dyb)
c     KP2=-z*((F1X+F2X)*sx+(F1Y+F2Y)*sy)*n2tr
      DO k=1,kmi
        KP2INT(k)=KP2(kmi-k+1)
      ENDDO
c     kkpmax2=1
c     kpmax2=0.
c     DO k=1,kmi
c       IF(dabs(KP2(k)).gt.kpmax2) THEN
c         kpmax2=dabs(KP2(k))
c         kkpmax2=k
c       ENDIF
c     ENDDO
C
c     kkpint=max(1,kkpmax2)
c     kkpint=max(1,kmli)
c     kkpint=max(1,max(kmli,kkpmax2))

c     ktap2=max(ktap,kmli)
c     kkpint=max(1,ktap2)
      kkpint=max(1,ktap)
c     IF(ktap.le.kmli) kkpint=1

      K02D=0.
      ALLOCATE(KINT2(kkpint))
      CALL dqtfg(zint(kmi-(kkpint-1):kmi),KP2INT(kmi-(kkpint-1):kmi),
     *           KINT2,kkpint)
      if (rd.ne.0d0.and.z(kkpint).ne.0.d0) then
        K02D=-(rd**2)*KINT2(kkpint)
      else
        K02D=0.d0
      endif
      DEALLOCATE(KINT2)

      IF(K02D.lt.0.) K02D=0.
c     IF(ktap.eq.kmi.OR.kmli.eq.kmi) K02D=0.

      KAPPAMZa=0.

C-- TONY -  1/20/2012 - additional term
      K02A=0.
      KP2=0.
      KP2INT=0.
      KP2=sqrt(GAMMA)*(LX**2+LY**2)*n2tr
      DO k=1,kmi
        KP2INT(k)=KP2(kmi-k+1)
      ENDDO
      kkpint=max(1,ktap)
      kkpint2=max(1,kmi-kkpint)
      IF(ktap.eq.1) kkpint2=max(1,kmi-kmli)
      ALLOCATE(KINT2(kkpint2))
      KINT2=0.
      CALL dqtfg(zint(1:kkpint2),KP2INT(1:kkpint2),
     *           KINT2,kkpint2)
      IF(ktap.ne.kmi.and.kmi.ne.1)
     *  K02A = (rd**2)*KINT2(kkpint2)
      DEALLOCATE(KINT2)

      IF(ktap.eq.kmi) K02A=0.
c     IF(K02D.eq.0.) K02A=0.
      IF(K02A.lt.0.) K02A=0.

C-- K02=K02/int(GAMMA**3/2) between -H and 0
C--- CA=((3/2)*Ko)^3/2 with Ko=4-> CA=14.7
      CA=14.7
C--- CD=((3/2)*Ko)^3/2 with Ko=4 -> CD=14.7
      CD=14.7
 
      INTGAMMAA=0.
      INTGAMMAD=0.

      KP2=0.
      KP2INT=0.
      KP2=GAMMA**(3./2.)
      DO k=1,kmi
        KP2INT(k)=KP2(kmi-k+1)
      ENDDO

      ALLOCATE(KINT2(kkpint2))
      KINT2=0.
      CALL dqtfg(zint(1:kkpint2),KP2INT(1:kkpint2),
     *           KINT2,kkpint2)
      INTGAMMAA=KINT2(kkpint2)
      DEALLOCATE(KINT2)

      ALLOCATE(KINT2(kkpint))
      KINT2=0.
      CALL dqtfg(zint(kmi-(kkpint-1):kmi),KP2INT(kmi-(kkpint-1):kmi),
     *           KINT2,kkpint)
      INTGAMMAD=KINT2(kkpint)
      DEALLOCATE(KINT2)

c?    IF(K02A.EQ.0.) INTGAMMAA=0.
c?    IF(K02D.EQ.0.) INTGAMMAD=0.

      INTGAMMA=INTGAMMAA/CA+INTGAMMAD/CD
      K02=0.

      INTGAMMA1=0.
      KP2=0.
      KP2INT=0.
      KP2=GAMMA
      DO k=1,kmi
        KP2INT(k)=KP2(kmi-k+1)
      ENDDO
      ALLOCATE(KINT2(kmi))
      KINT2=0.
      CALL dqtfg(zint,KP2INT,KINT2,kmi)
      INTGAMMA1=KINT2(kmi)
      DEALLOCATE(KINT2)

      INTEXP=0.
      CD1=0.
      KP2=0.
      KP2INT=0.
      KP2=GAMMA*exp(-((z+abs(z(kmi)))**2)/(2.*4000.**2))
      DO k=1,kmi
        KP2INT(k)=KP2(kmi-k+1)
      ENDDO
      ALLOCATE(KINT2(kkpint2))
      KINT2=0.
      CALL dqtfg(zint(1:kkpint2),KP2INT(1:kkpint2),
     *           KINT2,kkpint2)
      INTEXP=KINT2(kkpint2)
      DEALLOCATE(KINT2)
      CD1=INTEXP*sqrt(2./pi)*(abs(z(kmi)/4000.))

      IF(INTGAMMA1.ne.0.) INTGAMMAS=CD1*INTGAMMA/INTGAMMA1+INTGAMMAD/CD
      IF(INTGAMMAS.ne.0.) K02=(K02A+K02D)/INTGAMMAS

      IF(GAMMA(kmi).eq.1.) K02=0.

      IF (K02.le.0.) THEN
c       K02=0.d0
        K02=1.
        k02count=1.
      ELSE
        k02count=0.
      ENDIF
 
C---  Calculate udrift (ud,vd)
C---- Calculate ud0,vd0
      ALLOCATE(ud(kmi),vd(kmi))
      ud=0.
      vd=0.
      ud = UB1AVE - 0.5 * f * (rd**2) * (-DZLYB1AVE)
      vd = VB1AVE - 0.5 * f * (rd**2) * (+DZLXB1AVE)
C-- Calculate factor M
      ALLOCATE(MK(kmi))
      MK=0.
c     IF(K02.ne.0.) MK=1./(1.+0.75*((U-ud)**2+(V-vd)**2)/K02)
      DO k=1,kmi
      IF(K02*GAMMA(k).ne.0.)
     *MK(k)=1./(1.+0.5*((U(k)-ud(k))**2+(V(k)-vd(k))**2)/(K02*GAMMA(k)))
      ENDDO

C-- Calculate KAPPAM(Z)
      ALLOCATE(KAPPAMZ(kmi))
      KAPPAMZ=rd*dsqrt(K02*GAMMA)*MK
c     KAPPAMZ=rd*dsqrt(K02*GAMMA)
      KAPPAMZa=0.
      KAPPAMZa(1:kmi)=KAPPAMZ
      DEALLOCATE(KAPPAMZ)

      Ustar=0.
      Vstar=0.
c     Ustar(1:kmi)=KAPPAMZa

c     rdm=k02count

      DEALLOCATE(ud,vd,MK)
      DEALLOCATE(GAMMA)
      DEALLOCATE(n2sm)
      DEALLOCATE(n2t,n2t1,z,zm,b1,b1t,
     *  U,V,rho,n2tr,
     *  LX,LY,DZLX,DZLY,
     *  UINT,VINT,UT,VT,
     *  dxb,dyb,KP,
     *  zint,UREV,VREV,KPINT,
     *  b1rev,DZLXREV,DZLYREV,
     *  UTINT,VTINT,DZU,DZV,
     *  F1X,F1Y,F2X,F2Y,KP2,KP2INT)
 
 10   CONTINUE

 
      END  subroutine mesoscales1d


C*****
      SUBROUTINE B1AVERAGE(F,Z,B1,K,FAVE,KML)
      INTEGER K,KML
      REAL*8 F(K),Z(K),B1(K),FAVE,FB1AVE(K-KML+1),
     *  B1AVE(K-KML+1),a02
C
      a02=0.03
c     a02=0.5
C
      CALL dqtfg(Z(1:K-KML+1),F(1:K-KML+1)
     *           *dsqrt((a02+B1(1:K-KML+1)**2)/(1.+a02)),
     *           FB1AVE,K-KML+1)     
      CALL dqtfg(Z(1:K-KML+1),dsqrt((a02+B1(1:K-KML+1)**2)/(1.+a02)),
     *           B1AVE,K-KML+1)
C     WRITE(*,*),"B1AVE",FB1AVE/B1AVE
      FAVE=0.
      IF(B1AVE(K-KML+1).NE.0.)
     *FAVE=FB1AVE(K-KML+1)/B1AVE(K-KML+1)
C
      RETURN
      END
C*****
        SUBROUTINE d1sym(n,x,y,dyodx)
C080626 Input n,x(n),y(n)               Output dyodx(n)
C       Calculates the numerical ordinary derivative of y with respect to x, dyodx,
C       using the symmetrical approach of considering the nearest neighbor points
C       on both sides and assuming the derivative is changing linearly.
C       The ratio of differences on each side thus give the midpoint derivatives
C       on the respective sides and the averages of the midpoint derivatives
C       weighted each by the other midpoint's  distance give the derivative at the point.
C       Note for equal spacing same as ratio of neighbor differences with each other.
C       At top and bottom revert to ratio of differences with lone nearest neighbor.
        IMPLICIT NONE
        INTEGER n
        REAL*8 x(n),y(n),dyodx(n)
        REAL*8 dyodxm,dyodxp,dxm,dxp,dym,dyp

        INTEGER i


        DO i=1,n
           IF(i.EQ.1) THEN
              dyodx(i) = ((y(i+1)-y(i))/(x(i+1)-x(i)))
           ELSE IF(i.EQ.n) THEN
c             dyodx(i) = ((y(i)-y(i-1))/(x(i)-x(i-1)))
C-- TONY - 03/09/11 - updated the bottom calculation with a 3 points derivative
              IF(n.eq.2) THEN
                dyodx(i) = 0.
              ELSE
                dxm = x(i-1)-x(i-2)
                dxp = x(i)-x(i-1)
                dym = y(i-1)-y(i-2)
                dyp = y(i)-y(i-1)
                dyodxm = dym/dxm
                dyodxp = dyp/dxp
                dyodx(i) = (dxm*dyodxp + dxp*dyodxm)/(dxm+dxp)
              ENDIF
           ELSE
              dxm = x(i)-x(i-1)
              dxp = x(i+1)-x(i)
              dym = y(i)-y(i-1)
              dyp = y(i+1)-y(i)
              dyodxm = dym/dxm
              dyodxp = dyp/dxp
              dyodx(i) = (dxm*dyodxp + dxp*dyodxm)/(dxm+dxp)
           END IF
        END DO

        RETURN
        END
C*****
C080625-30AH Subroutine to calculate the first baroclinic mode in CD WKB approximation
C       from an input N^2 profile, artificially treating N^2 as zero where it is less than
C       zero since its square and fourth roots appear in the WKB approximation formula,
C       B_1(z)=[N(z)/N_max]^1/2 cos[(f r_d)^-1 \Integral_z^z_max N(z) dz (where z_max was
C       called "z_0" by CD and take (f r_d) = {\Integral -H^z_max N(z) dz \over pi},
C       but with B_1(z) kept at 1 above z_max as formula wasn't intended for that range,
C       corrected and adapted from Cheng's 080522 program that calls dqtfg , named B1b.f:
!@sum B1b.f based on B1a.f but only reads form one data file
!@+   and ignors the last line in the data file
!@    5-22-2008

        SUBROUTINE baroclin1st(z,n2,m,ba1,im,Nm,frd,igrid,jgrid)
C       Inputs:
C               z(m)    !Depth                                                  [m]
C               n2(m)   !Square of Brunt Vaisala frequency                      [s^-2]
C
C       Outputs:
C               ba1(m)  !First baroclinic mode
C               im      !Index of maximum N on column
C               Nm      !Maximum N on column                                    [s^-2]
C               frd     !(f r_d)                                                [m/s]
C
C
C       Internal:
C               lifout  !Logical switch set .TRUE. iff diagnostic output here

      implicit none
C080625AH
        integer igrid,jgrid
        INTEGER m       ! Number of ocean levels
        REAL*8 ba1(m)   ! Array for first baroclinic mode as calculated by this routine
        LOGICAL lifout
        PARAMETER(lifout=.TRUE.)
      real*8 n2,z,n,ni,pi,a,nm,arg,B1
      dimension n2(m)   !Array for N^2
      dimension z(m),n(m),ni(m) !Arrays for depth, SQRT(N^2), and z integral of N^2
        REAL*8 frd      !"f r_d" calculated in this routine
C******AH
C080625AH Work variables introduced.
        INTEGER im,i
C******AH
C080630AH Depth integrations are only done from maximum N level, "the pycnocline depth",
C       to the bottom so points above the depth of Maximum N are excluded unlike in B1b.f.
C       Introduce new arrays for depth and SQRT(N^2) which start at "pycnocline depth".
        REAL*8 zsubpyc(m),nsubpyc(m)
        INTEGER isubpyc(m),msubpyc
C******AH

      pi=acos(-1.d0)
      a=.1


C080625AH Canuto's WKB approximation formula contains (N^2)^1/2 and (N^2)^(1/4)
C       and therefore CANNOT HANDLE NEGATIVE N^2. Artificially set N^2=0 .
        DO i=1,m
           n2(i)=MAX(0.D0,n2(i))
        END DO
C******AH


      do i=1,m
         n(i)=sqrt(n2(i))
      end do

      ! find Nm = N_max
      Nm=N(1)
C080625AH Write out maximum N and depth when lifout=.TRUE. .
      im=1
      do i=2,m
         IF(N(i).GT.Nm) im=i
         Nm=max(N(i),Nm)
      end do
      IF(lifout) THEN
c       IF(ik.eq.177.and.jk.eq.163) THEN
c       write(*,*) "Nm=",Nm
c       WRITE(*,*) "im=",im
c       ENDIF
      END IF
C******AH

C080630AH Fill z & N arrays which go from max. N level down and integrate N on same range.
        DO i=im,m
           isubpyc(i)=i+1-im
           zsubpyc(isubpyc(i))=z(i)
           nsubpyc(isubpyc(i))=n(i)
        END DO
        msubpyc=m+1-im

      call dqtfg(zsubpyc,nsubpyc,ni,msubpyc)
C******AH

      do i=1,m ! z in meters
C80630AH Integrate N only over subpycnocline to find "f r_d". Keep B_1=1 above pycnocline.
         frd=ni(msubpyc)/pi
         IF(i.LT.im) THEN
           B1=1.D0
         ELSE
           if(frd.ne.0.d0) then           !Natassa
           arg=1./frd*ni(isubpyc(i))
           B1=(nsubpyc(isubpyc(i))/nm)**.5 * cos(arg)
           else
           B1=1.D0
           endif
         END IF
C******AH
C080625AH Store first baroclinic mode values in an array.
         ba1(i)=B1
C        Gam=1./(1+a) * (a+B1**2)
C        DMW=1./4.+3./4.*exp(-abs(z(i))/500.)
C        write(23,'(9e14.6)') z(i),Gam,sqrt(Gam)
C &        ,n2(i)/n2(1),DMW
C******AH
c        IF(ik.eq.177.and.jk.eq.163.and.i.eq.m) THEN
c          WRITE(*,*)'b1:ba1=',i,ba1(i)
c          WRITE(*,*)'b1:nsubpyc(isubpyc(i))=',nsubpyc(isubpyc(i))
c          WRITE(*,*)'b1:nm=',nm
c          WRITE(*,*)'b1:cos(arg)=',cos(arg)
c          WRITE(*,*)'b1:arg=',arg
c          WRITE(*,*)'b1:frd=',frd
c          WRITE(*,*)'b1:ni(isubpyc(i))=',ni(isubpyc(i))
c        ENDIF
         end do

c     IF(im.eq.m) WRITE(*,*)'b1:im,m,frd=',ik,jk,im,m,frd

      end
C*****
      subroutine dqtfg(x,y,z,n)
      !@sum integrate y from x(1) to x(i) and store it in z(i)
      implicit none
      real*8 sum1,sum2
      integer n,i
      real*8 x(n),y(n),z(n)
      sum2=0.
      if(n-1)4,3,1
    1 do 2 i=2,n
      sum1=sum2
      sum2=sum2+.5d0*(x(i)-x(i-1))*(y(i)+y(i-1))
    2 z(i-1)=sum1
    3 z(n)=sum2
    4 return
      end

C*********
      subroutine z121_mesosc (v0,kmtj,km)
!@sum z121 Apply 121 smoothing in k to 2-d array V(k=1,km)
!@+   top (0) value is used as a dummy
!@+   bottom (km+1) value is set to input value from above.
      IMPLICIT NONE
      REAL*8, PARAMETER :: p5=5d-1, p25=2.5d-1
      INTEGER, INTENT (IN) :: kmtj,km
      REAL*8, INTENT (INOUT) :: V0(km)  ! 2-D array to be smoothed
      REAL*8 V(0:km+1)
      INTEGER K
      REAL*8 tmp
      V(1:km)=V0
      V(0)      =  p25 * V(1)
      V(kmtj+1) =        V(kmtj)

      do k=1,kmtj
         tmp      =  V(k)
         V(k)   =  V(0)  + p5 * V(k) + p25 * V(k+1)
         V(0)   =  p25 * tmp
      end do
      V0=V(1:km)
      return
      end subroutine z121_mesosc
