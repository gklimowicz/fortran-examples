      module kprf_arrays_loc_renamer

      USE KPRF_ARRAYS, only :
     & klstsv => klstsv_loc,
     & jerlov => jerlov_loc,
     & t1sav => t1sav_loc,
     & s1sav => s1sav_loc,
     & tmlb => tmlb_loc,
     & smlb => smlb_loc,
     & hekman => hekman_loc,
     & hmonob => hmonob_loc,
     & dpbl => dpbl_loc,
     & dpbbl => dpbbl_loc,
     & dpmold => dpmold_loc,
     & tmix => tmix_loc,
     & smix => smix_loc,
     & thmix => thmix_loc,
     & umix => umix_loc,
     & vmix => vmix_loc,
     & betard, betabl, redfac,
     & akpar => akpar_loc,
     & zgrid => zgrid_loc,
     & vcty => vcty_loc,
     & difs => difs_loc,
     & dift => dift_loc,
     & ghats => ghats_loc,
     & buoflx => buoflx_loc,
     & bhtflx => bhtflx_loc,
     & mixflx => mixflx_loc,
     & sswflx => sswflx_loc,
     & lookup,
     & irimax,
     & nb,
     & ribtbl,
     & ridb,
     & slq2b
     &,dri
     &,smb
     &,shb
     &,ssb
     &,back_ra_r
     &,sisamax
     &,ra_rmax
     &,c_y_r0
     &,sm_r1
     &,sh_r1
     &,ss_r1
     &,slq2_r1
     &,b1,visc_cbu_limit,diff_cbt_limit
     &,theta_rcrp,theta_rcrn

      implicit none

      private

      public
     & klstsv
     &,jerlov
     &,t1sav
     &,s1sav
     &,tmlb
     &,smlb
     &,hekman
     &,hmonob
     &,dpbl
     &,dpbbl
     &,dpmold
     &,tmix
     &,smix
     &,thmix
     &,umix
     &,vmix
     &,betard, betabl, redfac
     &,akpar
     &,zgrid
     &,vcty
     &,difs
     &,dift
     &,ghats
     &,buoflx
     &,bhtflx
     &,mixflx
     &,sswflx
     &,lookup 
     &,irimax 
     &,nb 
     &,ribtbl 
     &,ridb 
     &,slq2b
     &,dri
     &,smb
     &,shb
     &,ssb
     &,back_ra_r
     &,sisamax
     &,ra_rmax
     &,c_y_r0
     &,sm_r1
     &,sh_r1
     &,ss_r1
     &,slq2_r1
     &,b1,visc_cbu_limit,diff_cbt_limit
     &,theta_rcrp,theta_rcrn

      end module kprf_arrays_loc_renamer
