#include "rundeck_opts.h"
      module hycom_arrays_glob_renamer

      USE HYCOM_ARRAYS, only :
     . u_loc => u,v_loc => v
     .,dp_loc => dp,dpold_loc => dpold
     .,dpu_loc => dpu,dpv_loc => dpv
     .,p_loc => p
     .,pu_loc => pu,pv_loc => pv
     .,latij_loc => latij,lonij_loc => lonij
     .,corio_loc => corio
     .,potvor_loc => potvor
     .,temp_loc => temp
     .,saln_loc => saln
     .,th3d_loc => th3d
     .,thstar_loc => thstar
     .,wgtkap_loc => wgtkap
     .,psikk_loc => psikk
     .,thkk_loc => thkk
     .,dpmixl_loc => dpmixl
     .,srfhgt_loc => srfhgt
c
     .,montg_loc => montg
     .,defor1_loc => defor1,defor2_loc => defor2
     .,ubavg_loc => ubavg,vbavg_loc => vbavg
     .,pbavg_loc => pbavg
     .,ubrhs_loc => ubrhs,vbrhs_loc => vbrhs
     .,utotm_loc => utotm,vtotm_loc => vtotm
     .,utotn_loc => utotn,vtotn_loc => vtotn
     .,uflux_loc => uflux,vflux_loc => vflux
     .,uflux1_loc => uflux1,vflux1_loc => vflux1
     .,uflux2_loc => uflux2,vflux2_loc => vflux2
     .,uflux3_loc => uflux3,vflux3_loc => vflux3
     .,uflx_loc => uflx,vflx_loc => vflx
     .,bolusu_loc => bolusu,bolusv_loc => bolusv
c
     .,uav_loc => uav, vav_loc =>   vav
     .,dpuav_loc => dpuav,dpvav_loc => dpvav
     .,temav_loc => temav,salav_loc => salav
     .,th3av_loc => th3av, dpav_loc =>  dpav
     .,pbavav_loc => pbavav,sfhtav_loc => sfhtav
     .,uflxav_loc => uflxav,vflxav_loc => vflxav
     .,ufxavp_loc => ufxavp,vfxavp_loc => vfxavp
     .,diaflx_loc => diaflx
     .,salflav_loc => salflav,brineav_loc => brineav
     .,eminpav_loc => eminpav
     .,surflav_loc => surflav
     .,tauxav_loc => tauxav
     .,tauyav_loc => tauyav
     .,ufxcum_loc => ufxcum,vfxcum_loc => vfxcum,dpinit_loc => dpinit
     .,dpmxav_loc => dpmxav,oiceav_loc => oiceav
     .,util1_loc => util1, util2_loc => util2
     .,util3_loc => util3, util4_loc => util4
c
     .,scpx_loc => scpx,scpy_loc => scpy
     .,scux_loc => scux,scuy_loc => scuy
     .,scvx_loc => scvx,scvy_loc => scvy
     .,scqx_loc => scqx,scqy_loc => scqy
     .,scu2_loc => scu2,scv2_loc => scv2
     .,scp2_loc => scp2,scq2_loc => scq2
     .,scuxi_loc => scuxi,scvyi_loc => scvyi
     .,scp2i_loc => scp2i,scq2i_loc => scq2i
c
     .,pgfx_loc => pgfx,pgfy_loc => pgfy
     .,gradx_loc => gradx,grady_loc => grady
     .,depthu_loc => depthu,depthv_loc => depthv
     .,pvtrop_loc => pvtrop
     .,depths_loc => depths
     .,drag_loc => drag
     .,glue_loc => glue
     .,zone_loc => zone
     .,dampu_loc => dampu,dampv_loc => dampv
c
     .,uja_loc => uja,ujb_loc => ujb
     .,via_loc => via,vib_loc => vib
     .,pbot_loc => pbot
     .,tracer_loc => tracer
     .,diadff_loc => diadff
     .,tprime_loc => tprime
     .,sgain_loc => sgain
     .,surflx_loc => surflx
     .,salflx_loc => salflx
     .,sflxcum_loc => sflxcum
     .,hflxcum_loc => hflxcum
c    .,thkice_loc => thkice
c    .,covice_loc => covice
c    .,temice_loc => temice
     .,odhsi_loc => odhsi
     .,odmsi_loc => odmsi
     .,omlhc_loc => omlhc
     .,dmfz_loc => dmfz
c
     &,klist_loc => klist
     &,ijlist_loc => ijlist
c
     .,taux_loc => taux
     .,tauy_loc => tauy
c    .,wndspd_loc => wndspd
c    .,airtmp_loc => airtmp
c    .,vapmix_loc => vapmix
c    .,oprec_loc => oprec
c    .,oevap_loc => oevap
     .,oemnp_loc => oemnp
     .,oicemlt_loc => oicemlt
     .,oflxa2o_loc => oflxa2o,oice_loc => oice
     .,ustar_loc => ustar
     .,ustarb_loc => ustarb
     .,osalt_loc => osalt
c
     .,freshw_loc => freshw
     .,diafor_loc => diafor
     .,diag1_loc => diag1
     .,diag2_loc => diag2
     .,diag3_loc => diag3
     .,diag4_loc => diag4

      implicit none

      private

      public u_loc
      public v_loc
      public dp_loc
      public dpold_loc
      public dpu_loc
      public dpv_loc
      public p_loc
      public pu_loc
      public pv_loc
      public latij_loc
      public lonij_loc
      public corio_loc
      public potvor_loc
      public temp_loc
      public saln_loc
      public th3d_loc
      public thstar_loc
      public wgtkap_loc
      public psikk_loc
      public thkk_loc
      public dpmixl_loc
      public srfhgt_loc
      public montg_loc
      public defor1_loc
      public defor2_loc
      public ubavg_loc
      public vbavg_loc
      public pbavg_loc
      public ubrhs_loc
      public vbrhs_loc
      public utotm_loc
      public vtotm_loc
      public utotn_loc
      public vtotn_loc
      public uflux_loc
      public vflux_loc
      public uflux1_loc
      public vflux1_loc
      public uflux2_loc
      public vflux2_loc
      public uflux3_loc
      public vflux3_loc
      public uflx_loc
      public vflx_loc
      public bolusu_loc
      public bolusv_loc
      public uav_loc
      public vav_loc
      public dpuav_loc
      public dpvav_loc
      public temav_loc
      public salav_loc
      public th3av_loc
      public dpav_loc
      public pbavav_loc
      public sfhtav_loc
      public uflxav_loc
      public vflxav_loc
      public ufxavp_loc
      public vfxavp_loc
      public diaflx_loc
      public salflav_loc
      public brineav_loc
      public eminpav_loc
      public surflav_loc
      public tauxav_loc
      public tauyav_loc
      public ufxcum_loc
      public vfxcum_loc
      public dpinit_loc
      public dpmxav_loc
      public oiceav_loc
      public util1_loc
      public util2_loc
      public util3_loc
      public util4_loc
      public scpx_loc
      public scpy_loc
      public scux_loc
      public scuy_loc
      public scvx_loc
      public scvy_loc
      public scqx_loc
      public scqy_loc
      public scu2_loc
      public scv2_loc
      public scp2_loc
      public scq2_loc
      public scuxi_loc
      public scvyi_loc
      public scp2i_loc
      public scq2i_loc
      public pgfx_loc
      public pgfy_loc
      public gradx_loc
      public grady_loc
      public depthu_loc
      public depthv_loc
      public pvtrop_loc
      public depths_loc
      public drag_loc
      public glue_loc
      public zone_loc
      public dampu_loc
      public dampv_loc
      public uja_loc
      public ujb_loc
      public via_loc
      public vib_loc
      public pbot_loc
      public tracer_loc
      public diadff_loc
      public tprime_loc
      public sgain_loc
      public surflx_loc
      public salflx_loc
      public sflxcum_loc
      public hflxcum_loc
      public odhsi_loc
      public odmsi_loc
      public omlhc_loc
      public dmfz_loc
      public taux_loc
      public tauy_loc
      public oemnp_loc
      public oicemlt_loc
      public oflxa2o_loc
      public oice_loc
      public ustar_loc
      public ustarb_loc
      public osalt_loc
      public freshw_loc
      public diafor_loc
      public klist_loc
      public ijlist_loc
      public diag1_loc,diag2_loc,diag3_loc,diag4_loc

      end module hycom_arrays_glob_renamer


