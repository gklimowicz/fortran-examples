#include "rundeck_opts.h"
!    **************************************************************
!@sum   SETTOMAS_LEV                                              
!    **************************************************************
!@+    This is to compute aerosol optical depth for each component. 
!@+  2=external mixing (=icomp-2) radiation calls  |
!@+  1=internal mixing (but AECOB is externally mixed) (ANUM_01) radiation call
!@+  2 is only available now.
!@auth  Yunha Lee, May 2006
C                        
                       
      subroutine SETTOMAS_LEV(i,j,l)

      USE domain_decomp_atm,ONLY: am_i_root

      USE RESOLUTION,  only: lm
      USE MODEL_COM,   only: itime,itimeI
      use OldTracer_mod, only: trName, TRPDENS
      USE TRACER_COM,  only: TRM,NBINS,n_ASO4,n_ANUM,xk
      USE RADPAR,      only: aesqex,aesqsc,aesqcb !Diagnostics
      USE ATM_COM, only : t            ! potential temperature (C)
     $                     ,q            ! saturated pressure
      ! aerosol radiative properties from lookup table
      USE TOMAS_AEROSOL, only : TOMAS_qext,TOMAS_qsca,TOMAS_gsca,icomp
      USE CONSTANT,   only : pi,lhe
      USE ATM_COM,   only: pmid,pk   ! midpoint pressure in hPa (mb)
      USE GEOM,        only: BYDXYP ! inverse area of gridbox [m-2]

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: i,j,l

C DEFINE VARIABLES***************************************************
      integer nwave,w,k,c,trnum
      parameter (nwave=6)
      integer wavelength(nwave) !just median value between two boundaries
C                         -------------NIR------------      VIS
C                     L=    1     2     3     4     5        6
C             WavA (nm)=  2200  1500  1250   860   770      300
C             WavB (nm)=  4000  2200  1500  1250   860      770
      data wavelength/3100.,1850.,1375.,1055.,815.,535./
      real rho_h2o,rho_spec(icomp-2) !density [kg/m3]
      real*8 density(icomp-2)
      real*8 mh2o,m_spec(icomp-2)    !mass [kg]
      real*8 rhe,qsat,temp
      real*8 reffwet,mp,mtot,vwetp,ntot
      real minX,maxX,stepX,x1
      real*8 vol_h2o(icomp-2),vol_spec(icomp-2),num_spec(icomp-2)
      integer iwave,ix1,irefre1,ilogrefim1
      real refre_spec(nwave,icomp-2),refim_spec(nwave,icomp-2)
      real refre_h2o(nwave),refim_h2o(nwave)
      real minrefre,maxrefre,steprefre,minLogrefim,maxLogrefim,refmed
     &     ,stepLogrefim,refre,refim   
      real*8 burden1

      real*8 waterso4, waternacl, waterocil,wrmass !for wateruptake function
      external waterso4, waternacl, waterocil
      real*8 help ! dummy variable for radiative forcing diagnostics

      real*8 qext,qsca,gsca
c$$$      data n_h2o/1.333e0/       ! Hale and Querry (1973)
c$$$      data k_h2o/1.96e-9/
c$$$      data n_bc/1.950e0/        ! Bond and Bergstrom (2006)
c$$$      data k_bc/0.790e0/
c$$$      data n_oc/1.530e0/        ! Kinne et al. (2003)
c$$$      data k_oc/6.00e-3/ 
c$$$      data n_dust/1.56e0/       ! Torres et al. (2002)
c$$$      data k_dust/0.006e0/ 
c$$$      data n_so4/1.43e0/        ! Hess et al (1998)
c$$$      data k_so4/1.00e-8/  
c$$$      data n_nacl/1.49e0/       ! Shettle and Fenn (1979) 
c$$$      data k_nacl/1.e-6/     

      data refre_spec/1.46099,1.48313,1.49719,1.50793,1.52000,1.52815, !Sulfate
     &     1.46390, 1.45000, 1.45000, 1.45000,1.45000,1.45000, ! Seasalt
     &  2.0, 1.98, 1.92,1.86,1.85,1.85,  !BC (ECOB)
     &  2.0, 1.98, 1.92,1.86,1.85,1.85,  !BC (ECIL) 
     &  1.46099,1.48313,1.49719,1.50805,1.52000,1.52775, ! OC (OCOB)
     &  1.46099,1.48313,1.49719,1.50805,1.52000,1.52775, ! OC (OCIL)
     &  1.47978,1.50719,1.51608,1.52998,1.54000, 1.56448 / !DUST 

!IMPORTANT - 2.15 for BC is changed to 2.(upper limit in TOMAS lookup table) 
!Perhaps can use extrapolation (like MATRIX model)

      data refim_spec/0.0764233,0.000516502,1.98240e-05,1.64469e-06,
     &     1.00000e-07, 1.00000e-07, !Sulfate
     &    0.00571719, 1.0e-09, 1.0e-09, 1.0e-09, 1.0e-09, 1.0e-09,  ! Seasalt
     &  1.05,0.79,0.75, 0.73,0.71,0.71,   !BC (ECOB)
     &  1.05,0.79,0.75, 0.73,0.71,0.71,  !BC (ECIL) 
     &  0.0761930,0.004700,0.004700,0.00480693,0.005400,0.0144927, ! OC (OCOB)
     & 0.0761930,0.004700,0.004700,0.00480693,0.005400,0.0144927, ! OC (OCIL)
     & 0.0211233,0.00584169,0.00378434,0.00178703,0.000800,0.00221463  / !DUST 

!Due to the mininum refim lookup table, set 0.0 to 1.e-9

      data refre_h2o/1.26304,1.31148,1.32283,1.32774,1.33059,1.33447/ !H2O 

      data refim_h2o/0.0774872, 0.000347758, 0.000115835, 3.67435e-06,
     &     1.58222e-07,3.91074e-08/

      data rho_spec/1780.,2165.,1800.,1800.,1400.,1400.,2650./

      data rho_h2o/1000./

C*********************************************************************

!      if (itime.ne.itimeI) then 
      temp = pk(l,i,j)*t(i,j,l) !should be in [K]
      rhe =100.d0* MIN(1.,q(i,j,l)/QSAT(temp,lhe,pmid(l,i,j))) ! rhe [0-100]
      if (rhe .gt. 99.d0) rhe=99.d0
      if (rhe .lt. 1.d0) rhe=1.d0

      do K=1,NBINS 
        mtot=0.d0
        do c=1,icomp-2
          trnum=n_aso4(1)-1+k+nbins*(c-1)
          m_spec(c)=trm(i,j,l,trnum) !aerosol mass in a size bin [kg]
          
          if(c.eq.1) then  !Sulfate water uptake
            m_spec(c)=trm(i,j,l,trnum)*1.1875 !so4 => (NH4)2SO4
            wrmass=waterso4(rhe)
          elseif (c.eq.2)then !Sea-salt water uptake
            wrmass=waternacl(rhe)
          elseif(c.eq.6)then 
            wrmass=waterocil(rhe) !OC water uptake 
          else
            wrmass=1.d0
          endif
          
          mh2o=m_spec(c)*(wrmass-1.d0)  !aerosol-water mass[kg]
          vol_h2o(c)=mh2o/rho_h2o
          mtot=m_spec(c)+mtot   !dry total mass in a size bin [kg]
        enddo
        
!Negligible number and mass ==> zero AOD
        if(trm(i,j,l,n_anum(1)-1+k).lt.1.d-5
     &       .or.mtot.le.0.) cycle !skip this bin

        mp=mtot/trm(i,j,l,n_anum(1)-1+k) ! dry particle diameter [kg/a particle]
        
        if(mp.gt.1000.*xk(nbins+1).or.mp.lt.xk(1)/10.)then
          print*,'mp is out of range',mp,k,rhe,c,trnum,n_aso4(1)
        endif
        
        ntot=0.d0
        do c=1,icomp-2
      
          num_spec(c)=m_spec(c)/mp ! number of each species [#]
          vol_spec(c)=m_spec(c)/rho_spec(c) !each species volume [m3]

          if(m_spec(c).gt.0.and.vol_spec(c).gt.0.) then
            vwetp=(vol_spec(c)+vol_h2o(c))/num_spec(c) !one wet particle's vol [m3]
            reffwet =(3./4.*vwetp/pi)**(1./3.) !m
            density(c) = (rho_spec(c)*vol_spec(c)+rho_h2o*vol_h2o(c))/
     &             (vol_spec(c)+vol_h2o(c)) !average density
          else
            reffwet=5.0e-10        ! Dp=1nm 
            density=0.d0
          endif

!dmw: skip calculation if density is tiny, 0, or negative
          if (density(c) .le. 1.d-10) cycle

          
C     BHMIE RADIATIVE PROPERTIES LOOKUP TABLE
          refmed = 1.0          ! refractive index of medium(air)-its imag. refr.index is close to zero. 
          minX=0.001            !nm lower bound
          maxX=1500.            !nm upper bound
          stepX=(sqrt(2.))**(1./3.)
          minRefRe=1.
          maxRefRe=2.
          stepRefRe=0.01
          minLogRefIm=-9.
          maxLogRefIm=0.
          stepLogRefIm=0.1
C     THIS IS FOR EXTERNALLY MIXED PARTICLE'S REFRACTIVE INDEX (EXTERNAL MIXING)
          do w=1, nwave

            if(vol_spec(c).gt.0.)then
              
             refre=(vol_spec(c)*refre_spec(w,c)+vol_h2o(c)*refre_h2o(w))
     &             /(vol_spec(c)+vol_h2o(c))
              
             refim=(vol_spec(c)*refim_spec(w,c)+vol_h2o(c)*refim_h2o(w))
     &             /(vol_spec(c)+vol_h2o(c))
              
            else
              refre=minrefre
              refim=1.e-9
            endif

C     Determine size parameter
            x1=pi*2.*reffwet*10**(9)/wavelength(w)
            if (x1.lt.1.e-3)then ! OUT OF RANGE
              print*, 'x1 is too small', x1
              x1=1.e-3
            endif
            if (x1.gt.1500.)then ! OUT OF RANGE
              print*, 'x1 is too big', x1
              x1=1500.
            endif

!     TO MATCH WITH  BHMIE LOOKUP TABLE
            ix1=1+int(log(x1/minX)/log(stepX))   

! temporal fix for refractive index
            if(refre.gt.maxRefRe) refre=maxrefre
            if(refre.lt.minRefRe) refre=minrefre
            if(log10(refim).gt.maxLogRefIm) refim=1.
            if(log10(refim).lt.minLogRefIm) refim=1.e-9

            irefre1=1+int((refre-minrefre)/steprefre)                        
            iLogRefIm1=1+int((log10(refim)-minLogRefIm)
     &           /stepLogRefIm) 

!     To read extinction efficiency and other properties    
            if(ix1.gt.124.or.irefre1.gt.101.
     &           or.ilogrefim1.gt.91) then
              print*,'bad lookup table',
     &           x1,refre,refim 
              CALL STOP_MODEL('bad rad lookup table',255) 
            endif
            
            qext=TOMAS_qext(ix1,irefre1,ilogrefim1) !extinction efficiency [no unit]
            qsca=TOMAS_qsca(ix1,irefre1,ilogrefim1) !scattering efficiency
            gsca=TOMAS_gsca(ix1,irefre1,ilogrefim1) !asymmetry parameter
            burden1 = m_spec(c)*bydxyp(j) ! [kg/m2] 

!            print*, 'rfwet', reffwet
!            print*, 'density', density
!            print*, 'qext', qext
!            print*, 'burden1',burden1
            aesqex(L,w,c)=(0.75d0/reffwet/density(c)    
     &           *qext*burden1)+aesqex(L,w,c)
            
            help = (aesqcb(L,w,c)*aesqsc(L,w,c))
     &           +(gsca*0.75d0/reffwet/density(c) 
     &           *qsca*burden1)
            
            aesqsc(L,w,c)=(0.75d0/reffwet/density(c)    
     &           *qsca*burden1)+aesqsc(L,w,c)
            
            aesqcb(L,w,c)=help/(aesqsc(L,w,c)+1.d-15)

          enddo                 !w:wavelength   
        enddo                   !c: chem components

      enddo                     !K=1,NBINS

!     endif !for timeI
      return
      end SUBROUTINE SETTOMAS_LEV

                    
!    **************************************************************
!@sum   SETTOMAS                                              
!    **************************************************************
!@+    This subroutine computes total aerosol optical depth (all comp) 
!@+    Currently, it assumes external mixing state
!@+    and no absorption in the longwave length. 
!@+    Need new subroutine for internal-mixing case
!@auth  Yunha Lee (modified from the existing modelE code)
C
C ************************************************************   
                       
      subroutine SETTOMAS(TOMAS_EXT,TOMAS_SCT,TOMAS_GCB,TOMAS_TRBALK)

      USE domain_decomp_atm,ONLY: am_i_root
      USE RESOLUTION,  only: lm
      USE MODEL_COM,   only: itime,itimeI
      USE RADPAR,      only: aesqex,aesqsc,aesqcb,FSTOPX,FTTOPX !Diagnostics

      USE TOMAS_AEROSOL, only : icomp
  
      IMPLICIT NONE

      ! Arguments: 
      INTEGER l,nc,w
      REAL*8 HELP, TOMAS_TAB(LM,33,ICOMP-2)   
      ! Arguments: Optical Parameters dimension(lm,wavelength)
      REAL*8, INTENT(OUT) :: TOMAS_EXT(LM,6)     ! Extinction AOD, SW
      REAL*8, INTENT(OUT) :: TOMAS_SCT(LM,6)     ! Scattering AOD, SW
      REAL*8, INTENT(OUT) :: TOMAS_GCB(LM,6)     ! Asymmetry Factor, SW
      REAL*8, INTENT(OUT) :: TOMAS_TRBALK(LM,33) ! Thermal absorption AOD?, LW

!=============================================================================
      TOMAS_EXT(:,:)=0.d0
      TOMAS_SCT(:,:)=0.d0
      TOMAS_GCB(:,:)=0.d0
      TOMAS_TRBALK(:,:)=0.d0
      TOMAS_TAB(:,:,:)=0.d0 ! zero for now
      
      if (itime.ne.itimeI) then 
        do L = 1,LM             !radiation has 3 extra levels on the top - aerosol are zero
          
          do nc=1,icomp-2
            do w=1,6
c     SW
              TOMAS_EXT(l,w) = TOMAS_EXT(l,w)+ aesqex(l,w,nc)*FSTOPX(nc)
              HELP = (TOMAS_SCT(l,w)*TOMAS_GCB(l,w))
     &             +(aesqcb(l,w,nc)*aesqsc(l,w,nc)*FSTOPX(nc))
              
              TOMAS_SCT(l,w) = TOMAS_SCT(l,w)+ aesqsc(l,w,nc)*FSTOPX(nc)
              TOMAS_GCB(l,w) = HELP/(TOMAS_SCT(l,w)+1.d-15)

            enddo
c     LW
              TOMAS_TRBALK(l,:) = TOMAS_TRBALK(l,:)
     &           +0.d0*TOMAS_TAB(l,:,nc)*FSTOPX(nc) !no absorption for longwave
          enddo
        enddo
      endif
     
      return
      end SUBROUTINE SETTOMAS
      
!    **************************************************************
!@sum   readmielut                                              
!    **************************************************************
!@+    This subroutine reads aerosol optical properties from the 
!@+    lookup table, which computed them by size parameter and chemical 
!@+    components.   
!@auth  Yunha Lee
C
C ************************************************************ 

      subroutine readmielut
c to read a mie lookup table for aerosol optical properties
      USE TOMAS_AEROSOL, only : TOMAS_qext, TOMAS_qsca,TOMAS_qabs,
     &     TOMAS_gsca!,TOMAS_qback

      character*80 infile
      integer innum, ii, jj, kk
      parameter (innum=61)     
      infile='bhmie_new.dat'

 1    format(e11.4,e11.4,e11.4,e11.4,e11.4)
      open(unit=innum,file=infile,FORM='FORMATTED',STATUS='OLD')
      do ii=1,124                !nX
         do jj=1,101            ! nrefre
            do kk=1,91          ! nlogrefim
               read(innum,1)
     &              TOMAS_qext(ii,jj,kk),TOMAS_qsca(ii,jj,kk),
     &              TOMAS_qabs(ii,jj,kk),TOMAS_gsca(ii,jj,kk),
     &              dummy
!TOMAS_qback is not used. Set to dummy variable!. 
            enddo
         enddo
      enddo
      close(innum)
      return
      end subroutine readmielut
      
