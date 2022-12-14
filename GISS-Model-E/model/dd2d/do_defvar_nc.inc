      type(dist_grid) :: grid
      integer :: fid
      logical, intent(in), optional :: r4_on_disk
      character(len=*), intent(in), optional :: defby
      character(len=*) :: varinfo
      integer :: dtypex,nd,shp(7),rc,vid
      if(grid%am_i_globalroot) then
      dtypex = dtype
      if(present(r4_on_disk)) then
        if(r4_on_disk) dtypex = nf_float
      endif
      nd = size(shape(arr))
      shp(1:nd) = shape(arr)
      call define_var(fid,dtypex,varinfo,nd,shp,
     &     grid%npx,grid%npy,grid%ntiles,rc,vid)
      endif
      call stoprc(rc,0)
      if(grid%am_i_globalroot) then
      if(present(defby)) then
        rc = nf_put_att_text(fid,vid,'defined_by',len_trim(defby),
     &       defby)
      endif
      endif
