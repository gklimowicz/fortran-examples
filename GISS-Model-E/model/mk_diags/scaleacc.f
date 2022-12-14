!@sum scaleacc is a generic scaling routine for modelE acc-files.
!@+   See conventions.txt for documentation.
!@auth M. Kelley
      subroutine scaleacc(fids,accname,accfile)
      use regrid_to_latlon_mod, only : setup_remap,regrid_4d
      implicit none
      include 'netcdf.inc'
      integer :: fids(2)              ! input file IDs
      character(len=*) :: accname     ! name of acc-array to scale
      character(len=200) :: accfile    ! name of acc-file
c
      character(len=200) :: ofile_base ! basename of output file
      character(len=20) :: dcat,dcat_cdl
      character(len=20) :: dcat_list(100)
      integer :: ndcats,icat
      real*4, dimension(:), allocatable :: scale_acc
      integer, dimension(:), allocatable :: ia_acc,denom_acc
      character(len=30), dimension(:), allocatable :: sname_acc
      character(len=30) :: vname
      real*4, dimension(:), allocatable :: xout,xout_,xden
      real*4, dimension(:), allocatable ::
     &     accarr,accarr_hemis,accarr_vmean
      real*4, dimension(:,:,:,:,:), allocatable :: xcs
      real*4, dimension(:,:,:,:), allocatable :: xll
      real*4, dimension(:), allocatable :: xout_hemis,xden_hemis
      real*4, dimension(:), allocatable :: xout_vmean,xden_vmean
      real*4, dimension(:), allocatable :: xout_all,xout_hemis_all
     &     ,xout_vmean_all
      integer :: idacc(12)
      integer :: n,k,kd,kacc,arrsize,arrsize_out,ndims,ndimsh,sdim,jdim,
     &     kk
      integer, dimension(7) :: srt,cnt,dimids,
     &     accsizes,shpout,hemi_sizes,vmean_sizes
      integer, dimension(7) :: cnt_hemis,cnt_vmean
      integer :: fid,status,ofid,accid,varid,varid_vmean,varid_hemis,
     &     accid_hemis,jdimid,accid_vmean,nvars,ndims_out
      real*4, parameter :: undef=-1.e30
      character(len=132) :: xlabel
      character(len=100) :: fromto
      logical :: do_hemis,do_vmean
      integer :: remap_fid,jmdid,nl,nk,imlon,jmlat,im,order
      logical :: remap_output
      integer :: tile_dim_out,d1,d2,d3
      integer :: slices_remaining,slices_total,arrsize_out_all
      integer :: n1,n2,n1_hemis,n2_hemis,n1_vmean,n2_vmean

      fid = fids(1)
      remap_fid = fids(2)
      remap_output = .false.
      if(remap_fid.ne.-99) remap_output = .true.

      xlabel=''; fromto=''
      status = nf_get_att_text(fid,nf_global,'xlabel',xlabel)
      status = nf_get_att_text(fid,nf_global,'fromto',fromto)

      if(trim(accname).eq.'all') then ! all diag. categories requested
        status = nf_inq_nvars(fid,nvars)
        ndcats = 0
        do varid=1,nvars
          vname=''
          status = nf_inq_varname(fid,varid,vname)
          if(index(vname,'_latlon').gt.0) cycle
          if(vname(1:4).eq.'cdl_') then
            ndcats = ndcats + 1
            dcat_list(ndcats) = vname(5:len_trim(vname))
          endif
        enddo
      else ! comma-separated list or single category
        ndcats = 1
        n = 1
        do k=1,len_trim(accname)-1
          if(accname(k:k).eq.',') then
            if(k.eq.n) stop 'parse error'
            dcat_list(ndcats) = accname(n:k-1)
            n = k + 1
            ndcats = ndcats + 1
          endif
          dcat_list(ndcats) = accname(n:len_trim(accname))
        enddo
      endif

      if(remap_output .and. ndcats.gt.1) then
        write(6,*) 'error: remapped output is currently available'
        write(6,*) 'only for single-category requests (not via "all")'
        stop
      endif

c
c Loop over the requested diagnostics categories
c
      dcat_loop: do icat=1,ndcats

      dcat = trim(dcat_list(icat))

      dcat_cdl = dcat

      call handle_err(nf_inq_varid(fid,trim(dcat),accid),
     &     'finding '//trim(dcat)//' in input file')
      call get_vdimsizes(fid,trim(dcat),ndims,accsizes)
      srt(:) = 1
      cnt(1:ndims) = accsizes(1:ndims)
      ndims_out = ndims-1
      tile_dim_out = ndims_out
      if(.not.remap_output) then
        status = nf_get_att_int(fid,accid,'tile_dim_out',tile_dim_out)
      endif
      write(6,*) 'processing ',trim(dcat)

c
c Find the size of the dimension along which to split the data
c
      status = nf_get_att_int(fid,accid,'split_dim',sdim)
      if(status.eq.nf_noerr) then
        kacc = accsizes(sdim)
        cnt(sdim) = 1
      else
        sdim = 7 ! necessary?
        kacc = 1
      endif
c
c Calculate the size of one input array
c
      arrsize = product(cnt(1:ndims))
      arrsize_out = arrsize
#ifdef HIMEM
      allocate(accarr(arrsize*kacc))
      status = nf_get_var_real(fid,accid,accarr)
#endif
c
c Setup some regridding info
c
      if(remap_output) then
        call setup_remap(remap_fid,imlon,jmlat)
        call get_dimsize(fid,'im',im)
        dcat_cdl = trim(dcat)//'_latlon'
        status = nf_inq_dimid(fid,'jm',jmdid)
        status = nf_inq_vardimid(fid,accid,dimids)
        do n=1,ndims
          if(dimids(n).eq.jmdid) jdim=n
        enddo
        nl = 1
        if(jdim.gt.2) nl = product(accsizes(1:jdim-2))
        if(sdim.lt.jdim-1) nl = nl/kacc
        nk = 1
        if(jdim.lt.ndims-1) nk = product(accsizes(jdim+1:ndims-1))
        if(sdim.gt.jdim) nk = nk/kacc
        arrsize_out = imlon*jmlat*nl*nk
        allocate(xcs(nl,im,im,nk,6),xll(nl,imlon,jmlat,nk))
      endif

      status = nf_inq_varid(fid,'hemis_'//trim(dcat),accid_hemis)
      do_hemis = status.eq.nf_noerr
      if(do_hemis) then
        call get_vdimsizes(fid,'hemis_'//trim(dcat),ndimsh,hemi_sizes)
        call get_dimsize(fid,'shnhgm',k)
        if(k.ne.3) stop 'bad size for shnhgm dimension'
      endif
      status = nf_inq_varid(fid,'vmean_'//trim(dcat),accid_vmean)
      do_vmean = status.eq.nf_noerr
     &     .and. (index(dcat,'ajl').gt.0 .or. index(dcat,'agc').gt.0)

c
c Allocate space for metadata and output arrays
c
      allocate(ia_acc(kacc),denom_acc(kacc),scale_acc(kacc))
      allocate(sname_acc(kacc))
      allocate(xout(arrsize_out),xden(arrsize_out))
      if(tile_dim_out.ne.ndims_out) then
        allocate(xout_(arrsize_out))
        shpout(1:sdim-1) = cnt(1:sdim-1)
        shpout(sdim:ndims-1) = cnt(sdim+1:ndims)
        d1 = product(shpout(1:tile_dim_out-1))
        d3 = 6
        d2 = arrsize_out/(d1*d3)
      endif
c
c Some setup for hemispheric/global means.
c We assume only that shnhgm and the latitude dimension
c have the same index.
c
      if(do_hemis) then
        status = nf_inq_vardimid(fid,accid_hemis,dimids)
        status = nf_inq_dimid(fid,'shnhgm',jdimid)
        do jdim=1,ndimsh
          if(jdimid.eq.dimids(jdim)) exit
        enddo
        cnt_hemis(1:ndimsh) = hemi_sizes(1:ndimsh)
        cnt_hemis(sdim) = 1
        arrsize = product(cnt_hemis(1:ndimsh))
        allocate(xout_hemis(arrsize),xden_hemis(arrsize))
#ifdef HIMEM
        allocate(accarr_hemis(arrsize*kacc))
        status = nf_get_var_real(fid,accid_hemis,accarr_hemis)
#endif
      endif

c
c Setup for vertical means (currently only for jl-type arrays).
c
      do_vmean = do_vmean .and. jdim.eq.1
      if(do_vmean) then
        vmean_sizes(1:3) = (/ accsizes(1)+3, 1, kacc /)
        allocate(xout_vmean(accsizes(1)+3),
     &           xden_vmean(accsizes(1)+3))
        cnt_vmean(:) = cnt(:)
        cnt_vmean(1) = accsizes(1)+3
        cnt_vmean(2) = 1
#ifdef HIMEM
        allocate(accarr_vmean(cnt_vmean(1)*kacc))
        status = nf_get_var_real(fid,accid_vmean,accarr_vmean)
#endif
      endif

c
c Define the output file
c
      k = index(accfile,'.acc')
      kk = index(accfile,'.subdd')
      if(k.gt.0) then
        kk = k+4
      elseif(kk.gt.0) then
        k = kk
        kk = kk+6
      else
        stop 'expecting *.acc*nc or *.subdd*.nc input file'
      endif
      ofile_base = accfile(1:k)//trim(dcat)//
     &     accfile(kk:index(accfile,'.nc')-1)
      k = index(ofile_base,'/',back=.true.)
      if(k.gt.0) ofile_base=ofile_base(k+1:len(ofile_base))
      call parse_cdl(fid,dcat_cdl,ofile_base,xlabel,fromto,
     &     remap_output,remap_fid)

c
c Read acc metadata needed for scaling
c
      call get_var_int(fid,'idacc',idacc)
      call get_var_real(fid,'scale_'//trim(dcat),scale_acc)
      status = nf_inq_varid(fid,'denom_'//trim(dcat),varid)
      if(status.eq.nf_noerr) then ! this acc array needs denom info
        call get_var_int(fid,'denom_'//trim(dcat),denom_acc)
      else
        denom_acc = 0
      endif
      call get_var_text(fid,'sname_'//trim(dcat),sname_acc)
      status = nf_inq_varid(fid,'ia_'//trim(dcat),varid)
      if(status.eq.nf_noerr) then ! this acc array has idacc-info
        call get_var_int(fid,'ia_'//trim(dcat),ia_acc)
      else
        ia_acc(:) = 1
        status = nf_inq_varid(fid,'ntime_'//trim(dcat),varid)
        if(status.eq.nf_noerr) then
          ! this acc array has a custom counter, store in idacc(1)
          call get_var_int(fid,'ntime_'//trim(dcat),idacc(1))
        else
          idacc(1) = 1 ! not a traditional accumulation (e.g. min, max)
        endif
      endif

c
c open output file
c
      status = nf_open(trim(ofile_base)//'.nc',nf_write,ofid)

c
c copy coordinate and other info from the acc file to the output file
c
      call copy_shared_vars(fid,ofid)
      if(remap_output) call copy_shared_vars(remap_fid,ofid)

c
c loop over outputs
c
      slices_remaining = 0
      do k=1,kacc
        if(slices_remaining.eq.0) then
          status = nf_inq_varid(ofid,trim(sname_acc(k)),varid)
          if(status.ne.nf_noerr) cycle ! this output was not requested
          call get_varsize(ofid,trim(sname_acc(k)),arrsize_out_all)
          arrsize_out_all = max(arrsize_out_all,arrsize_out)
          slices_total = arrsize_out_all/arrsize_out
          if(arrsize_out_all.ne.arrsize_out*slices_total) then
            write(6,*) 'scaleacc: size mismatch'
            stop
          endif
          slices_remaining = slices_total
          allocate(xout_all(arrsize_out_all))
          n1 = 1
          if(do_hemis) then
            status = nf_inq_varid(ofid,
     &           trim(sname_acc(k))//'_hemis',varid_hemis)
            allocate(xout_hemis_all(size(xout_hemis)*slices_total))
            n1_hemis = 1
          endif
          if(do_vmean) then
            status = nf_inq_varid(ofid,trim(sname_acc(k))//'_vmean',
     &           varid_vmean)
            if(status.eq.nf_noerr) then
              allocate(
     &             xout_vmean_all(size(xout_vmean)*slices_total))
            else
              varid_vmean = -99
            endif
            n1_vmean = 1
          endif
        endif
        slices_remaining = slices_remaining - 1

c
c scale this field
c
        srt(sdim) = k
        if(remap_output) then ! accumulation needs regridding first
#ifdef HIMEM
          call get_slice_real(accarr,ndims,accsizes,sdim,k,xcs)
#else
          status = nf_get_vara_real(fid,accid,srt,cnt,xcs)
#endif
          order = 1
          status = nf_get_att_int(ofid,varid,'regrid_order',order)
          call regrid_4d(xcs,xll,nl,nk,order)
          xout = reshape(xll,shape(xout))
        else
#ifdef HIMEM
          call get_slice_real(accarr,ndims,accsizes,sdim,k,xout)
#else
          status = nf_get_vara_real(fid,accid,srt,cnt,xout)
#endif
        endif
        xout = xout*scale_acc(k)/idacc(ia_acc(k))
        kd = denom_acc(k)
        if(kd.gt.0) then
          srt(sdim) = kd
          if(remap_output) then
#ifdef HIMEM
            call get_slice_real(accarr,ndims,accsizes,sdim,kd,xcs)
#else
            status = nf_get_vara_real(fid,accid,srt,cnt,xcs)
#endif
            call regrid_4d(xcs,xll,nl,nk,order)
            xden = reshape(xll,shape(xden))
          else
#ifdef HIMEM
            call get_slice_real(accarr,ndims,accsizes,sdim,kd,xden)
#else
            status = nf_get_vara_real(fid,accid,srt,cnt,xden)
#endif
          endif
          where(xden.ne.0.)
            xout = xout*idacc(ia_acc(kd))/xden
          elsewhere
            xout = undef
          end where
        endif

        if(tile_dim_out.ne.ndims_out) then
          xout_(:) = xout(:)
          call swap_outer_dims(xout_,xout,d1,d2,d3)
        endif
c
c write this field to the output file
c
        n2 = n1 + arrsize_out - 1
        xout_all(n1:n2) = xout
        n1 = n2 + 1
        if(slices_remaining.eq.0) then
          status = nf_put_var_real(ofid,varid,xout_all)
          deallocate(xout_all)
        endif

c
c scale/write the hemispheric/global means of this field if present
c
        if(do_hemis) then
          srt(sdim) = k
#ifdef HIMEM
          call get_slice_real(accarr_hemis,ndimsh,
     &         hemi_sizes,sdim,k,xout_hemis)
#else
          status = nf_get_vara_real(fid,accid_hemis,srt,cnt_hemis,
     &         xout_hemis)
#endif
          xout_hemis = xout_hemis*scale_acc(k)/idacc(ia_acc(k))
          if(kd.gt.0) then
            srt(sdim) = kd
#ifdef HIMEM
            call get_slice_real(accarr_hemis,ndimsh,
     &           hemi_sizes,sdim,kd,xden_hemis)
#else
            status = nf_get_vara_real(fid,accid_hemis,srt,cnt_hemis,
     &           xden_hemis)
#endif
            where(xden_hemis.ne.0.)
              xout_hemis = xout_hemis*idacc(ia_acc(kd))/xden_hemis
            elsewhere
              xout_hemis = undef
            end where
          endif

          n2_hemis = n1_hemis + size(xout_hemis) - 1
          xout_hemis_all(n1_hemis:n2_hemis) = xout_hemis
          n1_hemis = n2_hemis + 1
          if(slices_remaining.eq.0) then
            status = nf_put_var_real(ofid,varid_hemis,xout_hemis_all)
            deallocate(xout_hemis_all)
          endif
        endif

c
c scale/write the vertical means of this field if present
c
        if(do_vmean .and. varid_vmean.gt.0) then
          srt(sdim) = k
#ifdef HIMEM
          call get_slice_real(accarr_vmean,3,
     &         vmean_sizes,sdim,k,xout_vmean)
#else
          status = nf_get_vara_real(fid,accid_vmean,srt,cnt_vmean,
     &         xout_vmean)
#endif
          xout_vmean = xout_vmean*scale_acc(k)/idacc(ia_acc(k))
          if(kd.gt.0) then
            srt(sdim) = kd
#ifdef HIMEM
            call get_slice_real(accarr_vmean,3,
     &           vmean_sizes,sdim,kd,xden_vmean)
#else
            status = nf_get_vara_real(fid,accid_vmean,srt,cnt_vmean,
     &           xden_vmean)
#endif
            where(xden_vmean.ne.0.)
              xout_vmean = xout_vmean*idacc(ia_acc(kd))/xden_vmean
            elsewhere
              xout_vmean = undef
            end where
          endif

          n2_vmean = n1_vmean + size(xout_vmean) - 1
          xout_vmean_all(n1_vmean:n2_vmean) = xout_vmean
          n1_vmean = n2_vmean + 1
          if(slices_remaining.eq.0) then
            status = nf_put_var_real(ofid,varid_vmean,xout_vmean_all)
            deallocate(xout_vmean_all)
          endif

        endif
      enddo

      status = nf_close(ofid)

c
c Deallocate workspace arrays
c
      deallocate(ia_acc,denom_acc,scale_acc)
      deallocate(sname_acc)
      deallocate(xout,xden)
      if(tile_dim_out.ne.ndims_out) deallocate(xout_)
#ifdef HIMEM
      deallocate(accarr)
      if(allocated(accarr_hemis)) deallocate(accarr_hemis)
      if(allocated(accarr_vmean)) deallocate(accarr_vmean)
#endif

      if(allocated(xout_hemis)) deallocate(xout_hemis)
      if(allocated(xden_hemis)) deallocate(xden_hemis)

      if(allocated(xout_vmean)) deallocate(xout_vmean)
      if(allocated(xden_vmean)) deallocate(xden_vmean)

      enddo dcat_loop ! end loop over diagnostics categories

      return
      end subroutine scaleacc

      subroutine parse_cdl(fid,dcat,ofile_base,xlabel,fromto,
c Fow now, using the ncgen utility to parse the metadata
     &     define_lonlat,lonlat_defs_fid)
      implicit none
      integer :: fid,lonlat_defs_fid
      character(len=*) :: dcat,ofile_base,xlabel,fromto
      logical :: define_lonlat
      character(len=3) :: llstr
      character(len=100), dimension(:), allocatable :: cdl
      character(len=200) :: cdlfile
      integer :: k,kcdl,kgw,kend,imlon,jmlat
      
      if(define_lonlat) then
        call get_dimsize(lonlat_defs_fid,'lon',imlon)
        call get_dimsize(lonlat_defs_fid,'lat',jmlat)
      endif

      call get_dimsize(fid,'kcdl_'//trim(dcat),kcdl)      
      allocate(cdl(kcdl))
      cdl = ''
      call get_var_text(fid,'cdl_'//trim(dcat),cdl)
      cdl(1) = 'netcdf '//trim(ofile_base)//' { '
      cdlfile = trim(ofile_base)//'.cdl'
      open(10,file=cdlfile)
      if(define_lonlat) then ! fill in appropriate dim sizes
        do k=1,kcdl
          if(index(cdl(k),' lon = ').gt.0) then
            write(llstr,'(i3)') imlon
            cdl(k) = ' lon = '//trim(llstr)//';'
            exit
          endif
        enddo
        do k=1,kcdl
          if(index(cdl(k),' lat = ').gt.0) then
            write(llstr,'(i3)') jmlat
            cdl(k) = ' lat = '//trim(llstr)//';'
            exit
          endif
        enddo
      endif
      do k=kcdl,1,-1
        if(index(cdl(k),'}').gt.0) then
          kend = k; exit
        endif
      enddo
      kgw = kend
      do k=kend,1,-1
        if(index(cdl(k),'data:').gt.0) then
          kgw = k; exit
        endif
      enddo
      do k=1,kend
        if(k.eq.kgw) then
          write(10,*) '// global attributes:'
          write(10,'(a)') '    :xlabel = "'//trim(xlabel)//'" ;'
          if(len_trim(fromto).gt.0) then
            write(10,'(a)') '    :fromto = "'//fromto//'" ;'
          endif
        endif
        write(10,'(a)') cdl(k)
      enddo
      close(10)
      deallocate(cdl)
      call system('ncgen -b -x '//trim(cdlfile))
      call system('rm -f '//trim(cdlfile))
      return
      end subroutine parse_cdl

#ifdef HIMEM
      subroutine get_slice_real(arr,nd,dsizes,sd,s,arrout)
      implicit none
      real :: arr,arrout
      integer :: nd,dsizes(7),sd,s
      integer :: nl,nk
      nl = 1; nk = 1
      if(sd.gt.1) nl = product(dsizes(1:min(nd,sd-1)))
      if(sd.lt.nd) nk = product(dsizes(sd+1:nd))
      call copy_slice_real(arr,nl,dsizes(sd),nk,s,arrout)
      return
      end subroutine get_slice_real

      subroutine copy_slice_real(arr,nl,ns,nk,s,arrout)
      implicit none
      integer :: nl,ns,nk,s
      real :: arr(nl,ns,nk),arrout(nl,nk)
      arrout(1:nl,1:nk) = arr(1:nl,s,1:nk)
      return
      end subroutine copy_slice_real
#endif /* HIMEM */

      subroutine swap_outer_dims(a,b,s1,s2,s3)
      implicit none
      integer :: s1,s2,s3
      real*4 a(s1,s2,s3)
      real*4 b(s1,s3,s2)
      integer :: k
      do k=1,s2
        b(:,:,k) = a(:,k,:)
      enddo
      end subroutine swap_outer_dims
