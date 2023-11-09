#include "rundeck_opts.h"

      module hycom_dim_glob

#if !(defined CUBED_SPHERE) && !defined(STANDALONE_HYCOM)
      use hycom_dim, only : iia,jja
#endif

      use hycom_dim, only : idm, jdm, ms, kdm, ntrcr,
     &     iio,jjo,ii,jj,jchunk,kk,ii1,J_0,J_1,
     &     ogrid,
     &     ip_loc => ip,
     &     iu_loc => iu,
     &     iv_loc => iv,
     &     iq_loc => iq,
     .     ifp_loc => ifp,
     &     ilp_loc => ilp,
     &     isp_loc => isp,
     &     jfp => jfp,
     &     jlp => jlp,
     &     jsp => jsp,
     .     ifq_loc => ifq,
     &     ilq_loc => ilq,
     &     isq_loc => isq,
     &     jfq => jfq,
     &     jlq => jlq,
     &     jsq => jsq,
     .     ifu_loc => ifu,
     &     ilu_loc => ilu,
     &     isu_loc => isu,
     &     jfu => jfu,
     &     jlu => jlu,
     &     jsu => jsu,
     &     ifv_loc => ifv,
     &     ilv_loc => ilv,
     &     isv_loc => isv,
     &     jfv => jfv,
     &     jlv => jlv,
     &     jsv => jsv,
     &     msk_loc => msk
ccc
cddd     &     ip => ip,
cddd     &     iu => iu,
cddd     &     iv => iv,
cddd     &     iq => iq,
cddd     .     ifp => ifp,
cddd     &     ilp => ilp,
cddd     &     isp => isp,
cddd     &     jfp => jfp,
cddd     &     jlp => jlp,
cddd     &     jsp => jsp,
cddd     .     ifq => ifq,
cddd     &     ilq => ilq,
cddd     &     isq => isq,
cddd     &     jfq => jfq,
cddd     &     jlq => jlq,
cddd     &     jsq => jsq,
cddd     .     ifu => ifu,
cddd     &     ilu => ilu,
cddd     &     isu => isu,
cddd     &     jfu => jfu,
cddd     &     jlu => jlu,
cddd     &     jsu => jsu,
cddd     .     ifv => ifv,
cddd     &     ilv => ilv,
cddd     &     isv => isv,
cddd     &     jfv => jfv,
cddd     &     jlv => jlv,
cddd     &     jsv => jsv,
cddd     .     msk => msk
ccc
cddd     &     ip ,
cddd     &     iu ,
cddd     &     iv ,
cddd     &     iq ,
cddd     .     ifp,
cddd     &     ilp,
cddd     &     isp,
cddd     &     jfp,
cddd     &     jlp,
cddd     &     jsp,
cddd     .     ifq,
cddd     &     ilq,
cddd     &     isq,
cddd     &     jfq,
cddd     &     jlq,
cddd     &     jsq,
cddd     .     ifu,
cddd     &     ilu,
cddd     &     isu,
cddd     &     jfu,
cddd     &     jlu,
cddd     &     jsu,
cddd     .     ifv,
cddd     &     ilv,
cddd     &     isv,
cddd     &     jfv,
cddd     &     jlv,
cddd     &     jsv,
cddd     .     msk

      implicit none
      private

#if !(defined CUBED_SPHERE) && !defined(STANDALONE_HYCOM)
      public iia,jja
#endif

      public idm, jdm, ms, kdm, ntrcr
      public iio,jjo,ii,jj,jchunk,kk,ii1
      public gather_hycom_dim, alloc_hycom_dim_glob
      public ip
      public iu
      public iv
      public iq
      public ifp
      public ilp
      public isp
      public jfp
      public jlp
      public jsp
      public ifq
      public ilq
      public isq
      public jfq
      public jlq
      public jsq
      public ifu
      public ilu
      public isu
      public jfu
      public jlu
      public jsu
      public ifv
      public ilv
      public isv
      public jfv
      public jlv
      public jsv
      public msk


      integer, allocatable :: ip(:,:),iu(:,:),iv(:,:),iq(:,:),
     .ifp(:,:),ilp(:,:),isp(:),!jfp(:,:),jlp(:,:),jsp(:),
     .ifq(:,:),ilq(:,:),isq(:),!jfq(:,:),jlq(:,:),jsq(:),
     .ifu(:,:),ilu(:,:),isu(:),!jfu(:,:),jlu(:,:),jsu(:),
     .ifv(:,:),ilv(:,:),isv(:),!jfv(:,:),jlv(:,:),jsv(:),
     .msk(:,:)

      contains

      subroutine gather_hycom_dim
      USE DOMAIN_DECOMP_1D, ONLY: PACK_DATA,PACK_DATAJ

      !write(0,*) "ok ",__FILE__,__LINE__

      !write(803,*) ip_loc(:,J_0:J_1)
      !write(0,*) "ok ",__FILE__,__LINE__
      !write(0,*) shape(ip_loc), shape(ip)
      call pack_data( ogrid, ip_loc,  ip )
      !write(0,*) "ok ",__FILE__,__LINE__
      !write(804,*) ip(:,1:jdm)

      call pack_data( ogrid, iu_loc,  iu )
      call pack_data( ogrid, iv_loc,  iv )
      call pack_data( ogrid, iq_loc,  iq )

      call pack_dataj( ogrid, ifp_loc,  ifp )
      call pack_dataj( ogrid, ilp_loc,  ilp )
      call pack_data( ogrid, isp_loc,  isp )
      !call broadcast ( ogrid,  jfp, jfp_loc )
      !call broadcast ( ogrid,  jlp, jlp_loc )
      !call broadcast ( ogrid,  jsp, jsp_loc )
      call pack_dataj( ogrid, ifq_loc,  ifq )
      call pack_dataj( ogrid, ilq_loc,  ilq )
      call pack_data( ogrid, isq_loc,  isq )
      !call broadcast ( ogrid,  jfq, jfq_loc )
      !call broadcast ( ogrid,  jlq, jlq_loc )
      !call broadcast ( ogrid,  jsq, jsq_loc )
      call pack_dataj( ogrid, ifu_loc,  ifu )
      call pack_dataj( ogrid, ilu_loc,  ilu )
      call pack_data( ogrid, isu_loc,  isu )
      !call broadcast ( ogrid,  jfu, jfu_loc )
      !call broadcast ( ogrid,  jlu, jlu_loc )
      !call broadcast ( ogrid,  jsu, jsu_loc )
      call pack_dataj( ogrid, ifv_loc,  ifv )
      call pack_dataj( ogrid, ilv_loc,  ilv )
      call pack_data( ogrid, isv_loc,  isv )
      !call broadcast ( ogrid,  jfv, jfv_loc )
      !call broadcast ( ogrid,  jlv, jlv_loc )
      !call broadcast ( ogrid,  jsv, jsv_loc )
      call pack_data( ogrid, msk_loc,  msk )

cddd       ip= ip_loc
cddd       iu= iu_loc
cddd       iv= iv_loc
cddd       iq= iq_loc
cddd
cddd        ifp= ifp_loc
cddd        ilp= ilp_loc
cddd       isp= isp_loc
cddd        ifq= ifq_loc
cddd        ilq= ilq_loc
cddd       isq= isq_loc
cddd        ifu= ifu_loc
cddd        ilu= ilu_loc
cddd       isu= isu_loc
cddd        ifv= ifv_loc
cddd        ilv= ilv_loc
cddd       isv= isv_loc
cddd
cddd       msk= msk_loc


      end subroutine gather_hycom_dim


      subroutine alloc_hycom_dim_glob

      allocate(
     . ip(idm,jdm),iu(idm,jdm),iv(idm,jdm),iq(idm,jdm),
     .ifp(jdm,ms),ilp(jdm,ms),isp(jdm),!jfp(idm,ms),jlp(idm,ms),jsp(idm),
     .ifq(jdm,ms),ilq(jdm,ms),isq(jdm),!jfq(idm,ms),jlq(idm,ms),jsq(idm),
     .ifu(jdm,ms),ilu(jdm,ms),isu(jdm),!jfu(idm,ms),jlu(idm,ms),jsu(idm),
     .ifv(jdm,ms),ilv(jdm,ms),isv(jdm),!jfv(idm,ms),jlv(idm,ms),jsv(idm),
     .msk(idm,jdm) )

      end subroutine alloc_hycom_dim_glob



      end module hycom_dim_glob
