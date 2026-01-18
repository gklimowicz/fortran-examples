#include "rundeck_opts.h"
      MODULE COSMO_SOURCES

      USE TRACER_COM
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds

      SAVE
!@var be7_ and be10_src_3d Source functions for 7Be & 10Be (kg tracer)/( (kg air)/m^2 s )
      real*8, allocatable, dimension(:,:,:) :: be7_src_3d, be10_src_3d
!@var be7_src_param global multiplier of be7_src_3d to match obs
      real*8 :: be7_src_param=1    !default value
      real*8, allocatable, dimension(:,:) :: BE7W_acc, BE7D_acc 
      INTEGER :: variable_phi

      END MODULE COSMO_SOURCES


      SUBROUTINE init_cosmo
!      IMPLICIT NONE
      USE COSMO_SOURCES, only : be7_src_3d, be10_src_3d,
     $     BE7W_acc, BE7D_acc, variable_phi
      USE Dictionary_mod, only : sync_param
      USE TRACER_COM
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      INTEGER :: J_1H, J_0H, I_0H, I_1H

      call sync_param("variable_phi",variable_phi)

      call getDomainBounds(grid,J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      I_0H = grid%I_STRT_HALO
      I_1H = grid%I_STOP_HALO
      
      ALLOCATE(be7_src_3d (I_0H:I_1H, J_0H:J_1H, 1:lm))
      ALLOCATE(be10_src_3d (I_0H:I_1H, J_0H:J_1H, 1:lm))
      ALLOCATE(BE7W_acc (I_0H:I_1H, J_0H:J_1H))
      ALLOCATE(BE7D_acc (I_0H:I_1H, J_0H:J_1H))
 
      if (variable_phi .eq. 0) call read_Be_source_noAlpha
      print*, "variable_phi = ", variable_phi

      if (variable_phi .eq. 1) call read_Be_source
      print*, "variable_phi = ", variable_phi
         
      if (variable_phi .eq. 2) call update_annual_phi
      print*, "variable_phi = ", variable_phi
      
      if (variable_phi .eq. 3) call update_daily_phi   
      print*, "variable_phi = ", variable_phi

      END SUBROUTINE init_cosmo
      
      SUBROUTINE read_Be_source_noAlpha
!@sum reads in cosmogenic Be7 source from appropriate "old" (no alpha particles) versions 
!@sum of the Beer production files
!@auth C Salyk
      USE CONSTANT, only : byavog
      USE COSMO_SOURCES, only: be7_src_3d, be10_src_3d
      USE TRACER_COM
      USE GEOM, only: axyp
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds 
      USE FILEMANAGER, only: openunit,closeunit
      IMPLICIT NONE
      integer, parameter :: layers=23
      real*8, dimension(jm,lm) :: ibe, ibm
      integer iuc, j,i,l
C**** constants used in file to make numbers neater
      real*8 :: tfacti, tfact2
      real*8 :: be7_src_param=1 
      INTEGER :: J_1, J_0, I_0, I_1

      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Open source file
      call openunit('BE7_COSMO_NO_ALPHA', iuc, .false., .true.)
      read(iuc,*) tfacti
      read(iuc,*) ibe
C**** ibe has units atoms/g/s
      call closeunit(iuc)

C**** convert from atoms/g/s to (kg tracer)/ (kg air/m^2) /s
      do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
        be7_src_3d(i,j,l)=ibe(j,l)*axyp(i,j)*(tr_mm(n_Be7)*tfacti
     &    *byavog)
      end do ; end do ; end do

C**** multiply by air mass to put in the right units
      do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
         be10_src_3d(i,j,l) = 0.52d0 * be7_src_param * be7_src_3d(i,j,l)
     $        * tr_mm(n_Be10)/tr_mm(n_Be7)
          end do; end do; end do

         print*, "be7_src_param = ", be7_src_param
         print*, "tr_mm(n_Be10)/tr_mm(n_Be7) = ", tr_mm(n_Be10)
     $        /tr_mm(n_Be7)
         print*, "be7_src_3d(10,24,1) = ",be7_src_3d(10,24,1) 
         print*, "be10_src_3d(10,24,1) = ",be10_src_3d(10,24,1)


      END SUBROUTINE read_Be_source_noAlpha



      SUBROUTINE read_Be_source
!@sum reads in cosmogenic Be7 source from appropriate file
!@auth C Salyk
      USE CONSTANT, only : byavog
      USE COSMO_SOURCES, only: be7_src_3d, be10_src_3d
      USE TRACER_COM
      USE GEOM, only: axyp
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE FILEMANAGER, only: openunit,closeunit
      IMPLICIT NONE
      integer, parameter :: layers=23
      real*8, dimension(jm,lm) :: ibe, ibm
      real*8, dimension(jm,lm) :: ibe_10
      integer iuc, j,i,l
C**** constants used in file to make numbers neater
      real*8 :: tfacti, tfact2, tfacti_10
      real*8 :: be7_src_param=1 
      INTEGER :: J_1, J_0, I_0,I_1

!      call getDomainBounds(grid, J_STRT_HALO=J_0H, J_STOP_HALO=J_1H)
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

C**** Open source file
      print*, "layers = ", layers
      print*, "about to open Be7_cosmo"
      call openunit('BE7_COSMO', iuc, .false., .true.)
      print*, "opened Be7_cosmo file"
      read(iuc,*) tfacti
      print*, "read tfacti"
      print*, "tfacti = ", tfacti
      print*, "reading ibe"
      read(iuc,*) ibe
      print*, "read ibe"
C**** ibe has units atoms/g/s
      call closeunit(iuc)
      print*, "closed be7 cosmo file"
      
C**** convert from atoms/g/s to (kg tracer)/ (kg air/m^2) /s
      print*, "converting"
      do l=1,lm; do j=J_0,J_1 ; do i=I_0,I_1
        be7_src_3d(i,j,l)=ibe(j,l)*axyp(i,j)*(tr_mm(n_Be7)*tfacti
     &    *byavog)
      end do ; end do ; end do

C     repeat for Be10:
      print*, "about to open Be10_cosmo"
      call openunit('BE10_COSMO', iuc, .false., .true.)
      print*, "opened Be10_cosmo file"
      read(iuc,*) tfacti_10
      print*, "read tfacti_10"
      read(iuc,*) ibe_10
      print*, "read ibe_10"
C**** ibe has units atoms/g/s
      call closeunit(iuc)
      print*, "closed be10 cosmo file"
      
C**** convert from atoms/g/s to (kg tracer) (kg air/m^2) /s
      print*, "converting"
      do l=1,lm; do j=J_0,J_1 ; do i=I_0,I_1
        be10_src_3d(i,j,l)=ibe_10(j,l)*axyp(i,j)*(tr_mm(n_Be10)
     *       *tfacti_10*byavog)

      end do ; end do ; end do
      print*, "finished converting"

      
      END SUBROUTINE read_Be_source

      
      SUBROUTINE update_annual_phi
!@sum interpolates betw. diff. Be7 production values for each year to get correct Be7 
!@production corresponding to phi value for a given year.
!@auth C Field

      USE FILEMANAGER, only: openunit,closeunit
      USE GEOM, only: axyp
      use model_com, only: modelEclock
      USE CONSTANT, only : byavog
      USE TRACER_COM
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE COSMO_SOURCES, only : be7_src_3d, be10_src_3d
     
      IMPLICIT NONE
      character title*80
      integer, parameter :: npress=23, nlat=46, nphi=21
      integer, parameter :: nyrs=1002
      real :: year(nyrs), phi_record(nyrs), phi_yr, phi_list(nphi), phi
      real*8 :: be7_prod(1:nphi,1:nlat,1:npress)
      real*8 :: be10_prod(1:nphi,1:nlat,1:npress)
      real*8 :: slope_10(1:nlat, 1:npress)
      real*8 :: slope_7(1:nlat, 1:npress)
      real*8 :: tfacti_10, tfacti_7, delta_phi
      real*8 :: be7_src_param=1 
      real*8, allocatable, dimension(:,:) :: ibe_10, ibe_7
      integer :: iuc, i, j, k, l, n, J_1, J_0
 
      call getDomainBounds(grid, J_STRT=J_0, J_STOP=J_1)
      ALLOCATE(ibe_10 (J_0:J_1, 1:lm))
      ALLOCATE(ibe_7 (J_0:J_1, 1:lm))
      

      print*, "reading file for annual Be10 prod"
      call openunit('BE10_ANN_PROD', iuc, .false., .true.)
      print*, "opened"
      read(iuc,*) tfacti_10
      print*, "read tfacti_10"
      print*, "tfacti_10 = ", tfacti_10
      read(iuc,*) be10_prod
      print*, "be10 prod: ", be10_prod(12,12,12)
 20   call closeunit(iuc)
      
      print*, "reading file for annual Be7 prod"
      call openunit('BE7_ANN_PROD', iuc, .false., .true.)
      print*, "opened"
      read(iuc,*) tfacti_7
      print*, "read tfacti_7"
      print*, "tfacti_7 = ", tfacti_7
      read(iuc,*) be7_prod
      print*, "be7 prod: ", be7_prod(12,12,12)
 30   call closeunit(iuc)

      print*, "reading file for annual phi values"
      call openunit('PHI_ANN', iuc, .false., .true.)
      do i = 1,nyrs
         read(iuc,*) year(i), phi_record(i)
c        print*, "phi_record (1) = ", phi_record(1)
      end do
      call closeunit(iuc)  

      phi_yr = modelEclock%getYear() + 0.5
      print*, "phi_yr = ", phi_yr
      do i = 1,nyrs
         if (phi_yr .eq. year(i)) then
            phi = phi_record(i)
            print*, "found phi_record: ", phi
         end if
      end do
      
      phi_list(1) = 0.
      do n = 2,nphi
         phi_list(n) = phi_list(n-1) + 50.
         print*, "Phi = ", phi_list(n)
      end do


      do i = 1,nphi-1
         if ((phi_list(i) .le. phi) .and.(phi_list(i+1) .gt. phi)) then
            print*, "phi values: ", phi, phi_list(i), phi_list(i+1)
            delta_phi = phi - phi_list(i)
            print*, "delta_phi = ", delta_phi
            print*, "10: ", be10_prod(i+1,12,12), be10_prod(i,12,12)
            print*, "7: ", be7_prod(i+1,12,12), be7_prod(i,12,12)

            print*, "be10 =", (be10_prod(i+1,12,12)-be10_prod(i,12,12))
            print*, "diff phi =", ((phi_list(i+1))-(phi_list(i)))

            slope_10=(be10_prod(i+1,:,:)-be10_prod(i,:,:))/((phi_list(i
     $           +1))-(phi_list(i)))
            print*, "slope 10 = ", slope_10(12,12)
           
            slope_7=(be7_prod(i+1,:,:)-be7_prod(i,:,:))/((phi_list(i+1))
     $           -(phi_list(i)))
            print*, "slope 7 = ", slope_7(12,12)
            
            ibe_10(J_0:J_1,:) = be10_prod(i,J_0:J_1,:)
     *           +(slope_10(:,:)*delta_phi)
            ibe_7(J_0:J_1,:) = be7_prod(i,J_0:J_1,:)
     *           +(slope_7(:,:)*delta_phi)
            print*, "ibes: ", ibe_10(J_0:J_1,12), ibe_7(J_0:J_1,12)
            print*, "J_0 = ", J_0
            
         else if (phi_list(nphi) .le. phi) then
            print*, "phi values: ", phi, phi_list(nphi)
            delta_phi = 0.
            
            ibe_10(J_0:J_1,:) = be10_prod(nphi,J_0:J_1,:)
            ibe_7(J_0:J_1,:) = be7_prod(nphi,J_0:J_1,:)
            print*, "ibes: ", ibe_10(J_0,12), ibe_7(J_0,12)
            print*, 'J_0 = ', J_0
            
            
         end if
      end do
      
C**** convert from atoms/g/s to (kg tracer)/ (kg air/m^2) /s
      print*, "converting units for Be10 and Be7"
      do l=1,lm; do j=J_0,J_1; do i=I_0,I_1
        be10_src_3d(i,j,l)=ibe_10(j,l)*axyp(i,j)*(tr_mm(n_Be10)
     *       *tfacti_10*byavog)
         
        be7_src_3d(i,j,l)=ibe_7(j,l)*axyp(i,j)*(tr_mm(n_Be7)*tfacti_7
     $       *byavog)
      end do ; end do ; end do
      
      print*, "be7_src_param = ", be7_src_param
      print*, "tr_mm(n_Be10)/tr_mm(n_Be7) = ", tr_mm(n_Be10)
     $     /tr_mm(n_Be7)
      print*, "be7_src_3d(J_0,12,12) = ",be7_src_3d(J_0,12,12) 
      print*, "be10_src_3d(J_0,12,12) = ",be10_src_3d(J_0,12,12)
      print*, "J_0 = ", J_0
      
      END SUBROUTINE update_annual_phi


 
      SUBROUTINE update_daily_phi
!@sum to be used with doing model runs for Usoskin experiments, Jan. 2005 - Feb. 2005
!@sum reads in phi timeseries, Be7 production values and geomagnetic (Pc) values.
!@sum Production values are in atoms/g/s.
!@auth C Field
      USE FILEMANAGER, only: openunit,closeunit
      USE GEOM, only: axyp
      use model_com, only: modelEclock
      USE MODEL_COM, only : itime
      USE CONSTANT, only : byavog
      USE TRACER_COM
      USE DOMAIN_DECOMP_ATM, only : GRID, getDomainBounds
      USE COSMO_SOURCES, only : be7_src_3d

      IMPLICIT NONE
      character title*80
      integer, parameter :: npress=23, npc=41, nphi=31
      integer, parameter :: nmonths=2, ndays=59
      real :: day(ndays), daily_phi(ndays)
      real :: phi_day
      real :: pc_list(npc), phi_list(nphi)
      real*8 :: be7_prod(1:npress,1:npc,1:nphi)
      real*8 :: be7_scr(1:npress,1:npc)
      real*8 pc_month
      real*4, dimension(im,jm,nmonths) :: pc_table
      real*8, dimension(lm) :: new_prod, new_prod1, new_prod2, scr_prod
      real*8 :: be7_src_param=1 
      real :: phi_hi, phi_low, pc_hi, pc_low
      real :: delta_pc, delta_phi, delta_phi_1, delta_pc_3
      real :: slope(npress), slope_1(npress)
      real :: slope_2(npress), slope_3(npress), slope_4(npress) 
      integer :: iuc, i, j, k, l, m, n, iphi, ipc
      INTEGER :: J_1, J_0, I_0, I_1
      integer :: year, month, day, hour, date

      call modelEclock%get(month=month, dayOfYear=day)

      call getDomainBounds(grid,J_STRT=J_0, J_STOP=J_1)
      I_0 = grid%I_STRT
      I_1 = grid%I_STOP

! files needed:
!     BE_DAILY_PROD = Be23_TEST_all_UsoskinLinear_v3.dat = contains production values for all vertical layers, all pc values and all phi values (units are atoms/g/s).
!     BE_PC_VALUES = Pc_output_JF_Data_TEST.dat = contains pc values for each I/J box and each month from 1955 to 2006 

      print*, "reading files for daily phi"
      call openunit('BE7_DAILY_PROD', iuc, .false., .true.)
      print*, "opened"
      read(iuc,*,end=50) (((be7_prod(i,j,k), i=1,npress), j=1,npc), k=1
     $     ,nphi)
 !     print*, be7_prod(1,1,1)
 50   call closeunit(iuc)

      print*, "reading files for SCR production"
      call openunit('BE7_SCR_PROD', iuc, .false., .true.)
      print*, "opened"
      read(iuc,*,end=60) ((be7_scr(i,j), i=1,npress), j=1,npc)
 !     print*, be7_scr(1,1)
 60   call closeunit(iuc)
      
      print*, "reading"
      call openunit('PC_LOOKUP', iuc, .true., .true.)
      do k=1,nmonths
         read(iuc) title,pc_table(:,:,k)
         print*,"read ",title
      end do
!      print*, "pc_table (1,1,1) = ", pc_table(1,1,1) 
      call closeunit(iuc)  
      
       print*, "reading PC list"
      call openunit('PC_LIST', iuc, .false., .true.)
         read(iuc,*) (pc_list(i), i=1,npc)
!         print*, "pc_list (1) = ", pc_list(1) 
      call closeunit(iuc)

       print*, "reading Phi list"
      call openunit('PHI_LIST', iuc, .false., .true.)
         read(iuc,*) (phi_list(i), i=1,nphi)
         do j = 1,nphi
        print*, "phi_list = ", phi_list(j)
        end do

      call closeunit(iuc)

      print*, "reading PHI_BY_DAY"
      call openunit('PHI_BY_DAY', iuc, .false., .true.)
      do i = 1,ndays
         read(iuc,*) day(i), daily_phi(i)
!         print*, "daily_phi (1) = ", daily_phi(1)
      end do
      call closeunit(iuc)  

      

c     for each day, calculate Be7 production based on that month's pc
c     value and that day's phi value      

      if (jday > ndays) then
        phi_day = daily_phi(ndays) ! ? default 
      else
        phi_day = daily_phi(jday)
      end if
        print*, "phi_day = ", phi_day

c        phi_day = daily_phi(1)

      do m = 2,nphi
         if ((phi_day .lt. phi_list(m)) .and. (phi_day .ge. phi_list(m-1
     $        ))) then
            phi_low = phi_list(m-1)
            phi_hi = phi_list(m)
            iphi = m-1
            exit 
         end if
      end do

c      phi_low = phi_list(14)
c      phi_hi = phi_list(15)
c      iphi = 14

!      print*, "phi values: ", phi_low, phi_hi, iphi
      

      do j = J_0,J_1
         do i = I_0,I_1
            pc_month =  pc_table(i,j,jmon)
            if ((i .eq. 10) .and. (j .eq. J_1)) then
               print*, "pc_month = ", pc_month, jmon,i,j
            end if
            
            do n = 2,npc
               if ((pc_month .lt. pc_list(n)) .and. (pc_month
     $              .ge. pc_list(n-1))) then
                  pc_low = pc_list(n-1) 
                  pc_hi = pc_list(n) 
                  ipc = n-1
                  exit
               end if
            end do
            if ((i .eq. 10) .and. (j .eq. J_1)) then
               print*, "pc values (n. pole): ", pc_month, pc_low, pc_hi,
     $              ipc
            end if
            

c     there are 4 possible cases:
c     1. no interpolation needed
c     2. interpolation needed over pc values
c     3. interpolation needed over phi values
c     4. interpolation needed over both pc and phi
c     all four scenarios are accounted for in the following code:

c     first interpolate between the 2 phis at pc_lo:
               delta_phi_1 = phi_day - phi_low
               
               slope_1 = (be7_prod(:,ipc,iphi+1) - be7_prod(:,ipc,iphi)
     $              )/(phi_hi - phi_low)  
               
               new_prod1(:) = be7_prod(:,ipc,iphi) + (slope_1
     $              *delta_phi_1)
               if ((i .eq. 1) .and. (j .eq. 1)) then
                  print*, "new_prod1 = ", new_prod1(1)
               end if

c     next interpolate betw. the 2 phis at pc_hi:
               slope_2 = (be7_prod(:,ipc+1,iphi+1) - be7_prod(:,ipc+1
     $              ,iphi))/(phi_hi - phi_low)
               
               new_prod2(:) = be7_prod(:,ipc+1,iphi) + (slope_2
     $              *delta_phi_1)
               if ((i .eq. 1) .and. (j .eq. 1)) then
                  print*, "new_prod2 = ", new_prod2(1)
               end if
               

c     now interpolate betw. the 2 interpolated values that you've just
c     created:
               delta_pc_3 = pc_month - pc_low
               
               slope_3 = (new_prod2(:) - new_prod1(:))/(pc_hi -
     $              pc_low)
               
               new_prod(:) = new_prod1(:) + (slope_3*delta_pc_3)
               if ((i .eq. 1) .and. (j .eq. 1)) then
                  print*, "new_prod = ", new_prod(1)
               end if

!     calculate Be7 production from SCR if JDAY = Jan.20:
               if (jday .eq. 20) then
                  slope_4 = (be7_scr(:,ipc+1) - be7_scr(:,ipc))/(pc_hi -
     $                 pc_low)
                  scr_prod(:) = be7_scr(:,ipc) + (slope_4*delta_pc_3)
                  
                  new_prod = new_prod + scr_prod
               end if
            

C**** convert from atoms/g/s to (kg tracer)/ (kg air/m^2) /s
            do k=1,npress
!               print*, "converting atoms/g/s to kg/kg"
              be7_src_3d(i,j,k)=new_prod(k)*axyp(i,j)*(tr_mm(n_Be7)
     $              *byavog)
              if ((i .eq. 10) .and. (j .eq. 46) .and. (k .eq. 1)) then
                 print*, "DXYP = ", axyp(i,j)
                 print*, "tr_mm(n_Be7) = ", tr_mm(n_Be7)
                 print*, " "
              end if
              if ((i .eq. 37) .and. (j .eq. J_1)) then
                 print*, "itime = ", itime, "be7_src_3d (n. pole) = ",
     $                be7_src_3d(37,J_1,k)
              end if
           end do 
            
         end do
      end do
   

      print*, "closing update_daily_phi"

      END SUBROUTINE update_daily_phi
      









      



