! Program to collect AJL and AGC outputs into the traditional order
! in the PRT file.
      program fuse_ajlagc
      implicit none
      character(len=256) :: ajlfile,agcfile,outfile,custfile
      include 'netcdf.inc'
      integer :: status,ajl_fid,agc_fid,xxx_fid,ofid,ivid,nvars,vid_
      integer :: ifile,lstr,cnt,maxcnt,ivar,iunit,nargs,ios
      character(len=64) :: vname,vname_
      character(len=64), dimension(:), allocatable :: names
      logical :: custom_order,ex

      nargs = iargc()
      if(nargs.ne.3 .and. nargs.ne.4) then
        write(6,*)
     &       'usage: fuse_ajlagc ajlfile agcfile outfile [custfile]'
        write(6,*) 'Lines in custfile should be of the form'
        write(6,*) 'ajl:varname'
        write(6,*) 'or'
        write(6,*) 'agc:varname'
        stop
      endif
      call getarg(1,ajlfile)
      call getarg(2,agcfile)
      call getarg(3,outfile)
      custom_order = nargs==4

      call handle_err(nf_open(trim(ajlfile),nf_nowrite,ajl_fid),
     &     'error opening '//trim(ajlfile))
      call handle_err(nf_open(trim(agcfile),nf_nowrite,agc_fid),
     &     'error opening '//trim(agcfile))

      maxcnt = 1024
      allocate(names(maxcnt))
      cnt = 0

      if(custom_order) then
        call getarg(4,custfile)
        inquire(file=trim(custfile),exist=ex)
        if(.not.ex) then
          write(6,*) 'nonexistent file '//trim(custfile)
          stop
        endif
        open(iunit,file=trim(custfile),form='formatted',iostat=ios)
        if(ios.ne.0) then
          write(6,*) 'error opening '//trim(custfile)
          stop
        endif
        call read_custom_order(iunit,maxcnt,cnt,names)
        close(iunit)
      else
        call get_default_order(maxcnt,cnt,names)
      endif

      status = nf_create(trim(outfile),nf_clobber,ofid)

c copy global attributes
      vname = 'nf_global'
      call copy_varinfo(ajl_fid,ofid,vname)

c copy auxiliary variables from both files
      do ifile=1,2
        if(ifile.eq.1) then
          xxx_fid = ajl_fid
        else
          xxx_fid = agc_fid
        endif
        status = nf_inq_nvars(xxx_fid,nvars)
        do ivid=1,nvars
          vname = ''
          status = nf_inq_varname(xxx_fid,ivid,vname)
          lstr = len_trim(vname)
          if(lstr.ge.7) then
            if(vname(lstr-5:lstr).eq.'_hemis') cycle
            if(vname(lstr-5:lstr).eq.'_vmean') cycle
          endif
          if(nf_inq_varid(
     &         xxx_fid,trim(vname)//'_hemis',vid_).eq.nf_noerr)
     &         cycle
          call copy_varinfo(xxx_fid,ofid,vname)
        enddo
      enddo

c copy metadata for requested variables
      do ivar=1,cnt
        vname = names(ivar)
        if(vname(1:3).eq.'ajl') then
          xxx_fid = ajl_fid
          vname = vname(5:len_trim(vname))
        elseif(vname(1:3).eq.'agc') then
          xxx_fid = agc_fid
          vname = vname(5:len_trim(vname))
        else
          write(6,*) 'unrecognized category ',vname(1:3)
          cycle
        endif
        call copy_varinfo(xxx_fid,ofid,vname)
        vname_ = trim(vname)//'_hemis'
        call copy_varinfo(xxx_fid,ofid,vname_)
        vname_ = trim(vname)//'_vmean'
        call copy_varinfo(xxx_fid,ofid,vname_)
      enddo

      status = nf_enddef(ofid)

c copy all output variables
      call copy_shared_vars(ajl_fid,ofid)
      call copy_shared_vars(agc_fid,ofid)

      status = nf_close(ajl_fid)
      status = nf_close(agc_fid)
      status = nf_close(ofid)

      end program fuse_ajlagc

      subroutine get_default_order(maxcnt,cnt,names)
      implicit none
      integer :: maxcnt,cnt
      character(len=64), dimension(maxcnt) :: names

#define NEXT_ cnt=cnt+1;names(cnt)
! incrementer not working properly under gfortran
!#define NEXT_ names(next())

      cnt = 0

      NEXT_ = 'agc:npts_avg1'
      NEXT_ = 'agc:dp_cp2'
      NEXT_ = 'agc:stdev_dp'

      NEXT_ = 'ajl:tx'
      NEXT_ = 'ajl:tx_radonly'
      NEXT_ = 'ajl:height'
      NEXT_ = 'ajl:q'
      NEXT_ = 'ajl:rh'

      NEXT_ = 'ajl:cldh2o'
      NEXT_ = 'ajl:cldwtr'
      NEXT_ = 'ajl:cldice'

      NEXT_ = 'agc:u'
      NEXT_ = 'agc:v'
      NEXT_ = 'agc:vstar'

      NEXT_ = 'agc:psi_cp'
      NEXT_ = 'agc:psi_tem'

      NEXT_ = 'agc:wstar'
      NEXT_ = 'agc:vvel'

      NEXT_ = 'agc:stand_eddy_ke'
      NEXT_ = 'agc:eddy_ke'
      NEXT_ = 'agc:tot_ke'

      NEXT_ = 'agc:pot_temp'
      NEXT_ = 'agc:pot_vort'

      NEXT_ = 'agc:nt_sheat_eddy'
      NEXT_ = 'agc:nt_dse_stand_eddy'
      NEXT_ = 'agc:nt_dse_eddy'
      NEXT_ = 'agc:tot_nt_dse'

      NEXT_ = 'agc:nt_lh_eddy'
      NEXT_ = 'agc:jl_tot_nt_lh'
      NEXT_ = 'agc:nt_lh_stand_eddy'

      NEXT_ = 'agc:nt_lh_e'
      NEXT_ = 'agc:tot_nt_lh'
      NEXT_ = 'agc:nt_se_eddy'

      NEXT_ = 'agc:tot_nt_se'

      NEXT_ = 'agc:tot_nt_ke'

      NEXT_ = 'agc:nt_u_stand_eddy'
      NEXT_ = 'agc:nt_u_eddy'
      NEXT_ = 'agc:tot_nt_u'

      NEXT_ = 'agc:dyn_conv_dse'
      NEXT_ = 'agc:dyn_conv_eddy_geop'

      NEXT_ = 'agc:baroc_eddy_ke_gen'

      NEXT_ = 'agc:p2k_eddy_pgf'
      NEXT_ = 'agc:vt_geopot_eddy'

      NEXT_ = 'agc:vt_dse_e'
      NEXT_ = 'agc:tot_vt_dse'

      NEXT_ = 'agc:vt_lh_eddy1'
      NEXT_ = 'agc:jl_tot_vt_lh'
      NEXT_ = 'agc:vt_lh_eddy'
      NEXT_ = 'agc:tot_vt_lh'

      NEXT_ = 'agc:vt_se_eddy'
      NEXT_ = 'agc:tot_vt_se'

      NEXT_ = 'agc:tot_vt_ke'

      NEXT_ = 'agc:vt_u_eddy'
      NEXT_ = 'agc:tot_vt_u'

      NEXT_ = 'agc:vt_pv'
      NEXT_ = 'agc:vt_pv_eddy'

      NEXT_ = 'agc:nt_eddy_qgpv'
      NEXT_ = 'agc:del_qgpv'

      NEXT_ = 'agc:refr_ind_wave1'
      NEXT_ = 'agc:refr_ind_wave2'
      NEXT_ = 'agc:refr_ind_wave3'
      NEXT_ = 'agc:refr_ind_wave6'
      NEXT_ = 'agc:refr_ind_wave9'

      NEXT_ = 'agc:tot_dudt'

      NEXT_ = 'agc:dudt_sum1'
      NEXT_ = 'agc:dudt_meanadv'
      NEXT_ = 'agc:dudt_eddycnv'
      NEXT_ = 'agc:dudt_trnsadv'
      NEXT_ = 'agc:dudt_epflxdiv'
      NEXT_ = 'agc:dudt_fderr1'
      NEXT_ = 'agc:dudt_fderr2'

      NEXT_ = 'agc:dudt_mean_advec'
      NEXT_ = 'agc:dudt_eddy_conv'
      NEXT_ = 'agc:dudt_advec_tem'
      NEXT_ = 'agc:dudt_epdiv'

      NEXT_ = 'ajl:dudt_mtndrg'
      NEXT_ = 'ajl:dudt_dfmdrg'
      NEXT_ = 'ajl:dudt_shrdrg'

      NEXT_ = 'ajl:dudt_mcdrgpm10'
      NEXT_ = 'ajl:dudt_mcdrgpm40'
      NEXT_ = 'ajl:dudt_mcdrgpm20'

      NEXT_ = 'ajl:dudt_sumdrg'

      NEXT_ = 'ajl:dudt_sdiff'
      NEXT_ = 'ajl:dudt_vdiff'

      NEXT_ = 'ajl:dudt_sdrag'

      NEXT_ = 'agc:dtempdt'

      NEXT_ = 'agc:dtempdt_mean_advec'

      NEXT_ = 'agc:dtempdt_eddy_conv'

      NEXT_ = 'agc:dtempdt_advec_tem'

      NEXT_ = 'ajl:dtempdt_sdrag'

      NEXT_ = 'ajl:dtempdt_dynamics'

      NEXT_ = 'ajl:mc_mflx'
      NEXT_ = 'ajl:mc_dflx'

      NEXT_ = 'ajl:srad_heat'
      NEXT_ = 'ajl:srad_heat_radonly'
      NEXT_ = 'ajl:trad_cool'
      NEXT_ = 'ajl:trad_cool_radonly'
      NEXT_ = 'ajl:rad_cool'
      NEXT_ = 'ajl:rad_cool_radonly'

      NEXT_ = 'ajl:totcld'
      NEXT_ = 'ajl:sscld'
      NEXT_ = 'ajl:mccld'

      NEXT_ = 'ajl:rhe'

      NEXT_ = 'ajl:wcld'
      NEXT_ = 'ajl:icld'

      NEXT_ = 'ajl:wcod'
      NEXT_ = 'ajl:icod'

      NEXT_ = 'ajl:wcsiz'
      NEXT_ = 'ajl:icsiz'

      NEXT_ = 'ajl:tke'

      NEXT_ = 'ajl:lscond_heat'
      NEXT_ = 'ajl:turb_heat'

      NEXT_ = 'ajl:turb_lat'
      NEXT_ = 'ajl:moist_lat'

      NEXT_ = 'ajl:tot_ht_mc'

      NEXT_ = 'ajl:tot_ht_deepmc'
      NEXT_ = 'ajl:tot_ht_shlwmc'

      NEXT_ = 'ajl:tot_dry_mc'

      NEXT_ = 'ajl:csizmc'
      NEXT_ = 'ajl:csizss'

      NEXT_ = 'ajl:cnumwm'
      NEXT_ = 'ajl:cnumws'
      NEXT_ = 'ajl:cnumim'
      NEXT_ = 'ajl:cnumis'

      NEXT_ = 'agc:avail_pe'

      NEXT_ = 'ajl:dudt_dc'
      NEXT_ = 'ajl:dudt_mc'

      NEXT_ = 'ajl:u_epac'
      NEXT_ = 'ajl:v_epac'
      NEXT_ = 'ajl:vvel_epac'

      NEXT_ = 'ajl:u_wpac'
      NEXT_ = 'ajl:v_wpac'
      NEXT_ = 'ajl:vvel_wpac'

      NEXT_ = 'agc:epflx_north'
      NEXT_ = 'agc:epflx_vert'
      NEXT_ = 'agc:epflx_div'

      NEXT_ = 'agc:phi_amp_wave1'
      NEXT_ = 'agc:phi_amp_wave2'
      NEXT_ = 'agc:phi_amp_wave3'
      NEXT_ = 'agc:phi_amp_wave4'

      NEXT_ = 'agc:phi_phase_wave1'
      NEXT_ = 'agc:phi_phase_wave2'
      NEXT_ = 'agc:phi_phase_wave3'
      NEXT_ = 'agc:phi_phase_wave4'

      contains

      integer function next()
      cnt = cnt + 1
      if(cnt.gt.maxcnt) stop 'cnt > maxcnt'
      next = cnt
      end function next

      end subroutine get_default_order

      subroutine read_custom_order(iunit,maxcnt,cnt,names)
      implicit none
      integer :: iunit,maxcnt,cnt
      character(len=64), dimension(maxcnt) :: names
      character(len=64) :: name_
      integer :: ios
      cnt = 0
      do
        name_ = ''
        read(iunit,*,iostat=ios) name_
        if(ios.ne.0) exit
        if(len_trim(name_).gt.0) then
          if(name_(1:3).ne.'ajl' .and. name_(1:3).ne.'agc') then
            write(6,*) trim(name_)
            stop 'bad request in read_custom_order'
          endif
          cnt = cnt + 1
          names(cnt) = name_
        endif
      enddo
      end subroutine read_custom_order

      subroutine copy_varinfo(ifid,ofid,vname)
      implicit none
      include 'netcdf.inc'
      integer :: ifid,ofid
      character(len=64) :: vname
      integer :: status,ndims,n,dimsiz
      character(len=40) :: dimname,att_name,vname_
      integer :: vtype,natts,ivid,ovid
      integer, dimension(7) :: idimids,odimids

      ! first check if variable is already defined
      if(nf_inq_varid(ofid,trim(vname),ovid).eq.nf_noerr) return

      if(trim(vname).eq.'nf_global') then
        status = nf_inq_varnatts(ifid,nf_global,natts)
        ivid = nf_global
        ovid = nf_global
      else
        status = nf_inq_varid(ifid,trim(vname),ivid)
        if(status.ne.nf_noerr) return
        status = nf_inq_var(ifid,ivid,vname_,vtype,ndims,idimids,natts)
c copy dimension info as needed
        do n=1,ndims
          status = nf_inq_dimlen(ifid,idimids(n),dimsiz)
          status = nf_inq_dimname(ifid,idimids(n),dimname)
          if(nf_inq_dimid(ofid,trim(dimname),odimids(n)).ne.nf_noerr)
     &         status = nf_def_dim(ofid,dimname,dimsiz,odimids(n))
        enddo
c define variable
        status = nf_def_var(ofid,trim(vname),vtype,ndims,odimids,ovid)
      endif
c copy attributes
      do n=1,natts
        status = nf_inq_attname(ifid,ivid,n,att_name)
        status = nf_copy_att(ifid,ivid,att_name,ofid,ovid)
      enddo
      return
      end subroutine copy_varinfo
