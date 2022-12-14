      type(dist_grid) :: grid
      integer :: fid
      logical, intent(in), optional :: r4_on_disk
      logical, intent(in), optional :: with_record_dim
      character(len=*), intent(in), optional :: defby
      character(len=*) :: varinfo
      integer :: dtypex,nd,shp(7),rc,vid
      logical :: with_record_dim_
      integer*8 :: attlen8
      dtypex = dtype
      if(present(r4_on_disk)) then
        if(r4_on_disk) dtypex = nf_float
      endif
      with_record_dim_ = .false.
      if(present(with_record_dim)) with_record_dim_ = with_record_dim
      nd = size(shape(arr))
      shp(1:nd) = shape(arr)
      call broadcast(shp) ! in case size is zero on non-root PEs 
      call define_var(fid,dtypex,varinfo,nd,shp,
     &     grid%npx,grid%npy,grid%ntiles,rc,vid,grid%am_i_globalroot,
     &     with_record_dim_)
      call stoprc(rc,0)
      if(present(defby)) then
        attlen8 = len_trim(defby)
        rc = nfmpi_put_att_text(fid,vid,'defined_by',attlen8,defby)
      endif
