! /*
! 
! Copyright (C) 2008-2022 Michele Martone
! 
! This file is part of librsb.
! 
! librsb is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published
! by the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
! 
! librsb is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
! FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
! License for more details.
! 
! You should have received a copy of the GNU Lesser General Public
! License along with librsb; see the file COPYING.
! If not, see <http://www.gnu.org/licenses/>.
! 
! */
!
! Parallel Sparse BLAS fortran interface testing code
!
!
!> @cond INNERDOC
      PROGRAM main

      USE psb_base_mod
      USE psb_mvsv_tester
      IMPLICIT NONE
      INTEGER :: res,passed=0,failed=0,fi=0
      INTEGER            :: ictxt, iam=-1, np=-1
      CHARACTER(LEN=psb_fidasize_) :: afmt
      CALL psb_init(ictxt)
      CALL psb_info(ictxt,iam,np)
      IF(iam<0)THEN
            GOTO 9999
      ENDIF
      DO fi=1,2
      IF(fi.EQ.1)afmt=psb_csr_afmt_
      IF(fi.EQ.2)afmt=psb_coo_afmt_
      CALL       ts_sg_de_usmv_2_n_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_t_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_c_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_n_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_t_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_c_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_n_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_t_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_c_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_n_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_t_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_c_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_n_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_t_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_c_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_n_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_t_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_c_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_n_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_t_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_c_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_n_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_t_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       ts_sg_de_usmv_2_c_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_n_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_t_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_c_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_n_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_t_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_c_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_n_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_t_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_c_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_n_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_t_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_c_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_n_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_t_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_c_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_n_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_t_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_c_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_n_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_t_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_c_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_n_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_t_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       td_sg_de_usmv_2_c_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_n_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_t_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_c_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_n_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_t_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_c_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_n_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_t_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_c_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_n_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_t_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_c_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_n_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_t_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_c_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_n_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_t_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_c_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_n_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_t_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_c_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_n_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_t_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tc_sg_de_usmv_2_c_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_n_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_t_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_c_ap3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_n_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_t_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_c_ap3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_n_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_t_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_c_ap1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_n_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_t_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_c_ap1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_n_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_t_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_c_anr1_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_n_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_t_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_c_anr1_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_n_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_t_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_c_anr3_bp1_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_n_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_t_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      CALL       tz_sg_de_usmv_2_c_anr3_bnr0_ix1_iy1(errval,afmt,ictxt)
      IF(errval.NE.0)failed=failed+1
      IF(errval.EQ.0)passed=passed+1
      errval=0
      
      ENDDO
9999      CONTINUE
      PRINT *,"PASSED:",passed
      PRINT *,"FAILED:",failed
      CALL psb_exit(ictxt)
      END PROGRAM
!> @endcond



