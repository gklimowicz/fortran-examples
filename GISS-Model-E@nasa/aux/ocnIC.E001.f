      program ocnIC
!@sum ocnIC integrates ocean heat fluxes through 1 year and
!@+   writes an output disk file containing ALL ocean data
!@+   The values are obtained by integrating in time and
!@+   applying subroutine OSTRUC.
C****
C**** Output:
C****       TOCEAN(1) = mixed layer temperature
C****       TOCEAN(2) = mean temperature from mixed layer to ann max
C****       TOCEAN(3) = ocean temperature at annual maximum mixed layer
C****         Z1O = current mixed layer depth (at end of year)
C****
C**** Input: OSST = climatological ocean data
C****        SICE = sea ice data
C****        OCNML = mixed layer depth
C****        MLMAX = annual maximal mixed layer depths
C****        TOPO = topography
C****        SNOW = daily snow amounts (from vertflux)
C****
!AOO use statements added for domain_decomp and dynamics to pull in
!AOO dynamically allocated arrays: part 1 of 3
      use domain_decomp_1d, only : init_app, init_grid,grid, finish_app
!AOO                        end of part 1 of 3
      USE STATIC_OCEAN
      USE SEAICE, only : ace1i,ac2oim
      USE SEAICE_COM, only : snowi,msi,ssi
      USE FLUXES, only : sss
      USE FILEMANAGER
      USE GEOM
      implicit none
      integer i, j, k, m, last_day, kday, jday0, IH,
     *     months, month, iu_TOPO, iu_MLMAX, iu_SNOW, iu_OCNOUT
      REAL*4 month_day(12),focn4(im,jm),z4(im,jm)
      real*8 z1ox(im,jm),z12o_max
      CHARACTER*80 TITLE
      data month_day /31,28,31,30,31,30,31,31,30,31,30,31/
      off_line = .true. ! skip unneeded parts in OCLIM
!AOO calls to init routines for dynamically allocated arrays:part 2 of 3
      call init_app()
      call init_grid(grid,im,jm,lm)
      call alloc_drv()
!AOO end of part 2 of 3
      call getarg(1,title)
      read (title,*) months
      call getarg(2,title)
      read (title,*) z12o_max
C****
C**** Read in FOCEAN - ocean fraction
C****
      call openunit("TOPO",iu_TOPO,.true.,.true.)
      CALL READT_PARALLEL(grid,iu_TOPO,"TOPO",focean,1) ! Ocean fraction
      call closeunit(iu_TOPO)
C*
      fland = 1.- focean
      focn4(:,:) = focean(1:im,1:jm)
C*
C**** Read in aux. sea-ice file
C*
      call openunit("SICE",iu_SICE,.true.,.true.)
      CALL READT_PARALLEL(grid,iu_SICE,"SICE",DM,1)
C****
C**** Calculate spherical geometry
C****
      call GEOM_B
C**** set up unit numbers for ocean climatologies
      call openunit("OSST",iu_OSST,.true.,.true.)
C**** Set up unit number of mixed layer depth climatogies
      call openunit("OCNML",iu_OCNML,.true.,.true.)
C**** find and limit ocean ann max mix layer depths
      z12o = 0.
      do m=1,12
        CALL READT (iu_OCNML,0,IM*JM,z1ox,1)
        do j=1,jm
        do i=1,im
ccc       z12o(i,j)=min( z12o_max , max(z12o(i,j),z1ox(i,j)) )
ccc   the above line could substitute for next 3 lines: same results ?
          if (focean(i,j).gt.0. .and. z1ox(i,j).gt.z12o_max)
     *        z1ox(i,j)=z12o_max
          if (z1ox(i,j).gt.z12o(i,j)) z12o(i,j)=z1ox(i,j)
        end do
        end do
      end do
      rewind iu_OCNML
      write(6,*) 'Mixed Layer Depths limited to',z12o_max
C**** open snow file
      call openunit("SNOW",iu_SNOW,.true.,.true.)
C**** define sea surface salinity (needed for OCLIM)
      sss(:,:)=sss0
C**** initialise sea ice mass etc.
      do j=1,jm
        do i=1,im
          msi(i,j)=ac2oim
          ssi(1:2,i,j)=ssi0*ace1i*xsi(1:2)
          ssi(3:4,i,j)=ssi0*msi(i,j)*xsi(3:4)
        end do
      end do
C****
C**** Loop over days of the year
C****
      jday = JDendOfM(months-1)
      do m = months, months+11
        month = 1 + mod(m-1,12)
        last_day = month_day(month)
        do kday = 1,last_day
C*
C**** Interpolate daily ocean data from the monthly climatology
C*
          kocean = 0
          jmon = month
          jdate = kday

C***  Read in ocean ice snow data
          READ(iu_SNOW) TITLE,SNOWI
          WRITE (6,*) TITLE

          CALL OCLIM (.true.)

C***  Read in the ocean mixed layer depth data
C***  and interpolate them for the current day
          kocean = 1
          jday = jday + 1
          CALL OCLIM (.true.)

          IF(m.eq.months .and. kday.eq.1) THEN
C**** Initialize TOCEAN(2) and TOCEAN(3) on Day 1
            DO J = 1,JM
              DO I = 1,IM
                TOCEAN(2:3,I,J) = TOCEAN(1,I,J)
              END DO
            END DO
          ELSE
C**** Restructure the ocean temperature profile on subsequent days
            CALL OSTRUC(.false.)
          END IF
        end do
        print*,"Z1O,Z12O,TOCEAN(1:3)",Z1O(71,23),Z12O(71,23),TOCEAN(1:3
     *       ,71,23)
      end do
      CLOSE(iu_SNOW)
C****
C**** Write Z1O, TOCEAN to a disk file
C****
      print*,"OCNOUT: TOC(1:3),Z1O"
      print*,"   Arc",TOCEAN(:,1,45),Z1O(1,45)
      print*,"   Equ",TOCEAN(:,1,23),Z1O(1,23)
      call openunit("OCNOUT",iu_OCNOUT,.true.,.false.)
      WRITE (iu_OCNOUT) TOCEAN,Z1O
      call closeunit(iu_OCNOUT)
      WRITE (6,940)
C****
C**** PRODUCE MAPS OF OCEAN DATA ON LAST DAY
C****
      jday0=1
      IH=24*(JDAY0-1)
      TITLE = '  TGO = Ocean Temperature of Mixed Layer (C)'
      z4 = TOCEAN(1,1:im,1:jm)
      CALL MAP1(IM,JM,IH,TITLE,z4                 ,focn4,1.,0.,0)
      TITLE = ' TG2O = OCEAN TEMPERATURE OF SECOND LAYER (C)'
      z4 = TOCEAN(2,1:im,1:jm)
      CALL MAP1(IM,JM,IH,TITLE,z4                 ,focn4,1.,0.,0)
      TITLE = 'TG12O = OCEAN TEMPERATURE AT BOTTOM OF SECOND LAYER (C)'
      z4 = TOCEAN(3,1:im,1:jm)
      CALL MAP1(IM,JM,IH,TITLE,z4                 ,focn4,1.,0.,0)
      TITLE = '  Z1O = MIXED LAYER DEPTH (M)'
      z4 = Z1O(1:im,1:jm)
      CALL MAP1(IM,JM,IH,TITLE,z4       ,focn4,1.,0.,0)
      TITLE = ' Z12O = DEPTH OF BOTTOM OF SECOND LAYER (M)'
      z4 = Z12O(1:im,1:jm)
      CALL MAP1(IM,JM,IH,TITLE,z4        ,focn4,1.,0.,0)

!AOO not sure if this is needed, but just in case ... part 3 of 3
      call finish_app()
!AOO end of part 3 of 3
      STOP
C****
  940 FORMAT ('0Z1O, TG2O and TG12O written on unit 2,',
     *  ' DEC31.M25OD')
      END

