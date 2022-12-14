      program fix_ocean_rsf
C**** take ocean rsf files in modelE format and ensure that sea level is
C**** near zero with same pot. T and S
C**** compile with:
C****   gmake aux RUN=Exyz
C****   cd ../aux ; gmake fix_ocean_rsf.o
C**** f90 -o fix_ocean_rsf fix_ocean_rsf.o ../decks/Exyz_bin/Exyz.a
C****    -O2 -64 -mips4 -static -OPT:reorg_comm=off -w2 -listing -mp
C**** Note that since it uses modules and routines from the model, it
C**** must be compiled after the model
      USE DOMAIN_DECOMP_1D, ONLY : init_decomp
      USE CONSTANT, only : grav,rhows,bygrav
      USE MODEL_COM, only : im,jm,lm,p,xlabel,iowrite,ptop,irerun
     *     ,iowrite_mon,focean,nday,itime,itimei,itimee,itime0,iyear1
      USE GEOM, only : dxyp,geom_b
      USE SOMTQ_COM
      USE GHYCOM, only : snowe,tearth,wearth,aiearth,snoage,wbare,wvege
     *     ,htbare,htvege,snowbv,ngm,evap_max_ij,fr_sat_ij,qg_ij
      use veg_com, only : cint,qfol,cnc_ij
      USE SEAICE_COM, only : rsi,msi,hsi,snowi,ssi
      USE SEAICE, only : ace1i,xsi
      USE FILEMANAGER
      USE OCEAN, only : im,jm,lmo,g0m,gzmo,szmo,mo,hocean,lmm,lmu,lmv
     *     ,dypo,dxvo,opress,ogeoz,imaxj,ratoc,hatmo,ze,zerat,ze1
     *     ,GXMO,GYMO,SXMO,SYMO,S0M,ogeoz_sv
      USE OCFUNC, only : vgsp,tgsp,hgsp,agsp,bgsp,cgs
      USE OCEAN_DYN, only : vbar,mmi
      USE FLUXES, only : ogeoza
      IMPLICIT NONE
      CHARACTER infile*60, outfile*60,title*80
      INTEGER IARGC,I,J,L,N,ioerr,iu_TOPO,iu_OFTAB
      REAL*8 TAUX,X
      REAL*8 MSI1,tfo,MSL,AREA,PHIE,RATIO,MSL1,MSL2
      REAL*8, DIMENSION(IM,JM,LMO) :: UM0=0,VM0=0,MN,SN,GN
      INTEGER NS,mnow,ItimeX

      IF (IARGC().lt.1) THEN
        PRINT*,"Fix ocean rsf files to have 0m mean sea level"
        PRINT*," Needs TOPO set"
        PRINT*,"fix_ocean_rsf <rsf_file>"
        PRINT*,"Output is in <rsf_file>.adj"
        STOP
      END IF

!AOO calls to init routines for dynamically allocated arrays:part 2 of 3
      call init_decomp(grid,im,jm)
      call alloc_drv()
!
      CALL GETARG(1,infile)

      call io_rsf(trim(infile),ItimeX,irerun,ioerr)

C**** read in FLAKE/FOCEAN data
      call openunit ("TOPO", iu_TOPO, .true., .true.)
      CALL READT (iu_TOPO,0,IM*JM,FOCEAN,1) ! ocean fraction
      CALL READT (iu_TOPO,0,IM*JM,HATMO ,4) ! Atmo. Topography
      CALL READT (iu_TOPO,0,IM*JM,HOCEAN,1) ! Ocean depths
      call closeunit (iu_TOPO)

C**** get geometry
      call geom_b
      call geomo

      DO 110 L=0,LMO
  110 ZE(L) = ZE1*(ZERAT**L-1d0)/(ZERAT-1d0)

C**** Calculate LMM and modify HOCEAN
      DO 170 J=1,JM
      DO 170 I=1,IM
      LMM(I,J) =0
      IF(FOCEAN(I,J).LE.0.)  GO TO 170
      DO 150 L=1,LMO-1
  150 IF(HATMO(I,J)+HOCEAN(I,J) .le. 5d-1*(ZE(L)+ZE(L+1)))  GO TO 160
C     L=LMO
  160 LMM(I,J)=L
      HOCEAN(I,J) = -HATMO(I,J) + ZE(L)
  170 CONTINUE

C**** Read in table function for specific volume
      CALL openunit("OFTAB",iu_OFTAB,.TRUE.,.TRUE.)
      READ  (iu_OFTAB) TITLE,VGSP
      WRITE (6,*) 'Read from unit ',iu_OFTAB,': ',TITLE
      call closeunit(iu_OFTAB)

C**** Calculate mean sea level
C**** Calculate pressure anomaly at ocean surface (and scale for areas)
      DO J=1,JM
      DO I=1,IMAXJ(J)
        IF(FOCEAN(I,J).GT.0.) THEN
          OPRESS(I,J) = RATOC(J)*(100.*(P(I,J)+PTOP-1013.25d0)+
     *     RSI(I,J)*(SNOWI(I,J)+ACE1I+MSI(I,J))*GRAV)
        END IF
      END DO
      END DO

C**** calculate ocean density
      CALL OVTOM (MMI,UM0,VM0)
      CALL OPGF0
C****
C**** Calculate geopotential PHI (m**2/s**2)
C****
      DO J=1,JM
      DO I=1,IMAXJ(J)
        IF(FOCEAN(I,J).GT.0.) THEN
C**** Calculate geopotential by integrating from the bottom up
          PHIE = -HOCEAN(I,J)*GRAV
          DO L=LMM(I,J),1,-1
            PHIE = PHIE + VBAR(I,J,L)*MO(I,J,L)*GRAV
          END DO
          OGEOZA(I,J) = PHIE
C**** add in sea ice
          OGEOZA(I,J)=OGEOZA(I,J)*BYGRAV+
     *         RSI(I,J)*(MSI(I,J)+SNOWI(I,J)+ACE1I)/RHOWS
        END IF
      END DO
      END DO

      OGEOZA(2:IM, 1) = OGEOZA(1,1)
      OGEOZA(2:IM,JM) = OGEOZA(1,JM)

C**** calculate area weighted mean
      MSL=0.
      MSL1=0. ; MSL2=0.
      AREA=0.
      DO J=1,JM
      DO I=1,IM
        IF(FOCEAN(I,J).GT.0.) THEN
          MSL=MSL+ OGEOZA(I,J)*DXYP(J)
          MSL1=MSL1+ (OGEOZ(I,J)*BYGRAV+
     *         RSI(I,J)*(MSI(I,J)+SNOWI(I,J)+ACE1I)/RHOWS)*DXYP(J)
          MSL2=MSL2+ (OGEOZ_SV(I,J)*BYGRAV+
     *         RSI(I,J)*(MSI(I,J)+SNOWI(I,J)+ACE1I)/RHOWS)*DXYP(J)
          AREA=AREA+DXYP(J)
        END IF
      END DO
      END DO
      MSL=MSL/AREA
      print*,"Mean sea level",MSL,MSL1/AREA,MSL2/AREA

C**** Adjust mass to offset mean sea level change
      DO J=1,JM
      DO I=1,IM
        IF(FOCEAN(I,J).GT.0.) THEN
          RATIO=(1.-MSL/HOCEAN(I,J))
          MO(I,J,:)=MO(I,J,:)*RATIO
          S0M(I,J,:)=S0M(I,J,:)*RATIO
          G0M(I,J,:)=G0M(I,J,:)*RATIO
          GXMO(I,J,:)=GXMO(I,J,:)*RATIO
          GYMO(I,J,:)=GYMO(I,J,:)*RATIO
          GZMO(I,J,:)=GZMO(I,J,:)*RATIO
          SXMO(I,J,:)=SXMO(I,J,:)*RATIO
          SYMO(I,J,:)=SYMO(I,J,:)*RATIO
          SZMO(I,J,:)=SZMO(I,J,:)*RATIO
#ifdef TRACERS_OCEAN
          TRMO(I,J,:,:) = TRMO(I,J,:,:)*RATIO
          TXMO(I,J,:,:) = TXMO(I,J,:,:)*RATIO
          TYMO(I,J,:,:) = TYMO(I,J,:,:)*RATIO
          TZMO(I,J,:,:) = TZMO(I,J,:,:)*RATIO
#endif
        END IF
      END DO
      END DO

C**** check calculation
C**** calculate ocean density
      CALL OVTOM (MMI,UM0,VM0)
      CALL OPGF0
C****
C**** Calculate geopotential PHI (m**2/s**2)
C****
      DO J=1,JM
      DO I=1,IMAXJ(J)
        IF(FOCEAN(I,J).GT.0.) THEN
C**** Calculate geopotential by integrating from the bottom up
          PHIE = -HOCEAN(I,J)*GRAV
          DO L=LMM(I,J),1,-1
            PHIE = PHIE + VBAR(I,J,L)*MO(I,J,L)*GRAV
          END DO
          OGEOZA(I,J) = PHIE
C**** save new values
          OGEOZ(I,J)=OGEOZA(I,J)
          OGEOZ_SV(I,J)=OGEOZA(I,J)
C**** add in sea ice
          OGEOZA(I,J)=OGEOZA(I,J)*BYGRAV+
     *         RSI(I,J)*(MSI(I,J)+SNOWI(I,J)+ACE1I)/RHOWS
        END IF
      END DO
      END DO

      OGEOZA(2:IM, 1) = OGEOZA(1,1)
      OGEOZA(2:IM,JM) = OGEOZA(1,JM)

C**** calculate area weighted mean
      MSL=0.
      AREA=0.
      DO J=1,JM
      DO I=1,IM
        IF(FOCEAN(I,J).GT.0.) THEN
          MSL=MSL+ OGEOZA(I,J)*DXYP(J)
          AREA=AREA+DXYP(J)
        END IF
      END DO
      END DO
      MSL=MSL/AREA
      print*,"New mean sea level",MSL

      outfile=trim(infile)//".adj"
      call io_rsf(trim(outfile),ItimeX,iowrite_mon,ioerr)
      print*,ioerr

      print*,"New rsf file written out to ",trim(outfile)
      stop
 800  print*,"Error reading in AIC file "
      stop
 810  print*,"Error reading in AIC file - end reached"
      stop
      end
