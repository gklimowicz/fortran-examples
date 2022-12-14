! (c) British Crown Copyright 2008, the Met Office.
! All rights reserved.
! $Revision: 88 $, $Date: 2013-11-13 07:08:38 -0700 (Wed, 13 Nov 2013) $
! $URL: http://cfmip-obs-sim.googlecode.com/svn/stable/v1.4.0/cosp_lidar.F90 $
!
! Redistribution and use in source and binary forms, with or without modification, are permitted
! provided that the following conditions are met:
!
!     * Redistributions of source code must retain the above copyright notice, this list
!       of conditions and the following disclaimer.
!     * Redistributions in binary form must reproduce the above copyright notice, this list
!       of conditions and the following disclaimer in the documentation and/or other materials
!       provided with the distribution.
!     * Neither the name of the Met Office nor the names of its contributors may be used
!       to endorse or promote products derived from this software without specific prior written
!       permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
! IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
! FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
! DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
! IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
! OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

!
! History:
! Jul 2007 - A. Bodas-Salcedo - Initial version
! Oct 2008 - S. Bony          - Instructions "Call for large-scale cloud" removed  -> sgx%frac_out is used instead.
!                               Call lidar_simulator changed (lsca, gbx%cca and depol removed;
!                               frac_out changed in sgx%frac_out)
! Jun 2011 - G. Cesana        - Added betaperp_tot argument
!
!
MODULE MOD_COSP_LIDAR
  USE MOD_COSP_CONSTANTS
  USE MOD_COSP_TYPES
  IMPLICIT NONE

CONTAINS


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!------------------- SUBROUTINE COSP_LIDAR ------------------------
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_LIDAR(gbx,sgx,sghydro,y)

  ! Arguments
  type(cosp_gridbox),intent(in) :: gbx  ! Gridbox info
  type(cosp_subgrid),intent(in) :: sgx  ! Subgrid info
  type(cosp_sghydro),intent(in) :: sghydro  ! Subgrid info for hydrometeors
  type(cosp_sglidar),intent(inout) :: y ! Subgrid output

  ! Local variables
  integer :: i
  real :: presf(sgx%Npoints, sgx%Nlevels + 1)
  real,dimension(sgx%Npoints, sgx%Nlevels) :: lsca,mr_ll,mr_li,mr_cl,mr_ci
  real,dimension(sgx%Npoints, sgx%Nlevels) :: beta_tot,tau_tot
  real,dimension(sgx%Npoints, sgx%Nlevels) :: betaperp_tot
  real,dimension(sgx%Npoints, PARASOL_NREFL)  :: refle

  ! Debugging
  ! WRITE(0,FMT='(20x, "in COSP_LIDAR")')
  ! write(0,fmt='(/"gbx%T(1) Min = ",2E10.3," Max = ",E10.3," NZ = ",I10)') &
  !   MINVAL(gbx%T(:,1)), MINVAL(gbx%T(:,1),mask=gbx%T(:,1).gt.0), &
  !   MAXVAL(gbx%T(:,1)), COUNT(mask=gbx%T(:,1).gt.0)
  ! write(0,fmt='(/"gbx%T(Nlevels) Min = ",2E10.3," Max = ",E10.3," NZ = ",I10)') &
  !   MINVAL(gbx%T(:,sgx%Nlevels)), MINVAL(gbx%T(:,sgx%Nlevels),mask=gbx%T(:,sgx%Nlevels).gt.0), &
  !   MAXVAL(gbx%T(:,sgx%Nlevels)), COUNT(mask=gbx%T(:,sgx%Nlevels).gt.0)

  ! Fix?
  y%temp_tot(:,:) = gbx%T(:,:)

  presf(:,1:sgx%Nlevels) = gbx%ph
  presf(:,sgx%Nlevels + 1) = 0.0
  lsca = gbx%tca-gbx%cca
  do i=1,sgx%Ncolumns
      ! Temporary arrays for simulator call
      mr_ll(:,:) = sghydro%mr_hydro(:,i,:,I_LSCLIQ)
      mr_li(:,:) = sghydro%mr_hydro(:,i,:,I_LSCICE)
      mr_cl(:,:) = sghydro%mr_hydro(:,i,:,I_CVCLIQ)
      mr_ci(:,:) = sghydro%mr_hydro(:,i,:,I_CVCICE)
      ! Debugging
      !WRITE(0,FMT='(20x, "call lidar_simulator()", I10)') i

      call lidar_simulator(sgx%Npoints, sgx%Nlevels, 4, PARASOL_NREFL, LIDAR_UNDEF  &
                 , gbx%p, presf, gbx%T, mr_ll, mr_li, mr_cl, mr_ci &
                 , gbx%Reff(:,:,I_LSCLIQ), gbx%Reff(:,:,I_LSCICE) &
                 , gbx%Reff(:,:,I_CVCLIQ), gbx%Reff(:,:,I_CVCICE) &
                 , gbx%lidar_ice_type, y%beta_mol, beta_tot &
                 , betaperp_tot, tau_tot, refle )

      y%betaperp_tot(:,i,:) = betaperp_tot(:,:)
      y%beta_tot(:,i,:) = beta_tot(:,:)
      y%tau_tot(:,i,:)  = tau_tot(:,:)
      y%refl(:,i,:)     = refle(:,:)
  enddo

! ! Debugging
!   WRITE(0,FMT='(/"results COSP_LIDAR()")')
!   write(0,fmt='(/"gbx%T(1) Min = ",2E10.3," Max = ",E10.3," NZ = ",I10)') &
!     MINVAL(gbx%T(:,1)), MINVAL(gbx%T(:,1),mask=gbx%T(:,1).gt.0), &
!     MAXVAL(gbx%T(:,1)), COUNT(mask=gbx%T(:,1).gt.0)
!   write(0,fmt='(/"gbx%T(Nlevels) Min = ",2E10.3," Max = ",E10.3," NZ = ",I10)') &
!     MINVAL(gbx%T(:,sgx%Nlevels)), MINVAL(gbx%T(:,sgx%Nlevels),mask=gbx%T(:,sgx%Nlevels).gt.0), &
!     MAXVAL(gbx%T(:,sgx%Nlevels)), COUNT(mask=gbx%T(:,sgx%Nlevels).gt.0)
END SUBROUTINE COSP_LIDAR

END MODULE MOD_COSP_LIDAR
