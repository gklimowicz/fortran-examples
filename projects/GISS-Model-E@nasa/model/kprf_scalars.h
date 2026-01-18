      real
     . dp00,cb,tmljmp,rigr,ribc,rinfty,ricr,bldmin,bldmax,cekman,
     . cmonob,difm0,difs0,difmiw,difsiw,dsfmax,rrho0,cs,cstar,cv,c11,
     . vonk,zmin,zmax,umin,umax,epsilon,deltaz,deltau,vtc,cg,dp0enh,
     . qdif0,qdifiw

      integer hblflg

      logical locsig,bblkpp,shinst,dbdiff,nonloc,latdiw,botdiw,difsmo
     .       ,mxlkpp,mxlgis

      common /kprf_real/
     . dp00,cb,tmljmp,rigr,ribc,rinfty,ricr,bldmin,bldmax,cekman,
     . cmonob,difm0,difs0,difmiw,difsiw,dsfmax,rrho0,cs,cstar,cv,c11,
     . vonk,zmin,zmax,umin,umax,epsilon,deltaz,deltau,vtc,cg,dp0enh,
     . qdif0,qdifiw

      common /kprf_integ/hblflg
c
      common /kprf_char/
     .  locsig,bblkpp,shinst,dbdiff,nonloc,latdiw,botdiw,difsmo
     . ,mxlkpp,mxlgis
c
c --- nasa giss variables
c
      integer nextrtbl0,ifexpabstable,nextrtbl1,
     &        nextrtbl,nposapprox,mt0,mt,ntbl,
     &        mt_ra_r,n_theta_r_oct,nbig

      common/gissi1/
     &        nextrtbl0,ifexpabstable,nextrtbl1,
     &        nextrtbl,nposapprox,mt0,mt,ntbl,
     &        mt_ra_r,n_theta_r_oct,nbig
      save  /gissi1/
c
      real    deltheta_r,pidbl,rri

      common/gissr1/
     &        deltheta_r,pidbl,rri
      save  /gissr1/
c
      integer ifback,ifsali,ifepson2,ifrafgmax,
     &        ifsalback,ifchengcon,ifunreal,idefmld,
     &        ifpolartablewrite,ifbg_theta_interp

      common/gissi3/
     &        ifback,ifsali,ifepson2,ifrafgmax,
     &        ifsalback,ifchengcon,ifunreal,idefmld,
     &        ifpolartablewrite,ifbg_theta_interp
      save  /gissi3/
c
      real   back_ph_0,adjust_gargett,back_k_0,back_del_0,back_s2,
     &       ri0,ebase,epson2_ref,
     &       eps_bot0,scale_bot,      !for bottom-enhanced
     &       eplatidepmin,wave_30,    !and latitude dependent mixing
     &       deltemld,delrhmld,
     &       back_sm2,v_back0,t_back0,
     &       s_back0,ri_internal,backfrac,backfact,ako,tpvot0,sgmt,
     &       tptot0,tpcot0,ttot0,tcot0,tctot0,tpvot,tptot,tpcot,
     &       ttot,tcot,tctot,back_l_0

      common/gissr3/
     &       back_ph_0,adjust_gargett,back_k_0,back_del_0,back_s2,
     &       ri0,ebase,epson2_ref,
     &       eps_bot0,scale_bot,      !for bottom-enhanced
     &       eplatidepmin,wave_30,    !and latitude dependent mixing
     &       deltemld,delrhmld,
     &       back_sm2,v_back0,t_back0,
     &       s_back0,ri_internal,backfrac,backfact,ako,tpvot0,sgmt,
     &       tptot0,tpcot0,ttot0,tcot0,tctot0,tpvot,tptot,tpcot,
     &       ttot,tcot,tctot,back_l_0
      save  /gissr3/
c
