      program ll2cs
!@auth D,Gueyffier
!@ver 1.0
!@sum    ll2cs is a generic tool to convert latlon input files to the cubed sphere
!@+  it is used together with companion script remap.pl and regridding parameter file
!@+
!@usage  to compile: cd aux; gmake regrid 
!@+  to run : ./remap.pl -par parameter_file -in input_file -out output_file
!@+  where parameter_file is e.g. regrid.par
!@+  and input_file is e.g. top_index_360x180.ij.ext
!@+  output_file is the resulting cubed sphere file
!@+
!@+  The parameter file contains the following entries
!@+  filedir      =  path to input file (e.g. /discover/nobackup/projects/giss/prod_input_files/)
!@+  regridfile   =  remapping file containing interpolation weights (e.g. remap360-180C90-90.nc)
!@+  imsource     =  im resolution of source (latlon) grid
!@+  jmsource     =  jm resolution of source (latlon) grid
!@+  ntilessource =  number of faces of source (latlon) grid (=1)
!@+  imtarget     =  im resolution of target (cubed sphere) grid
!@+  jmtarget     =  jm resolution of target (cubed sphere) grid
!@+  ntilestarget =  number of faces of target (cubed sphere) grid (=6)
!@+  format       =  'ij' : each record contains title, data(imsource,jmsource)
!@+                  'ijl': each record contains title, data(imsource,jmsource,levels)
!@+                  'lij': each record contains title, data(levels,imsouce,jmsource)
!@+  nfields      =  sets the number of fields per record
!@+                         i.e.  title, data_1(im,jm), data_2(im,jm)..data_nfields(im,jm)
!@+  levels (optional) = for 'ijl' or 'lij' arrays only, sets the number of "l" levels 
!@+  title        =  defines if each record contains a title or not
!@+  maintile     =  defines if the first record is a title with no data
!@+
!@+   Program is standalone and runs in serial mode only 
!@+
!@+   The latlon input file has to conform to the GISS file format:
!@+   char*80 title, real*4 data
!@+
!@+   if the parameter real is set to real8 then the program expects
!@+   char*80 title, real*8 data
!@+
!@auth Denis Gueyffier
      use regrid, only : init_regrid, x_2grids, do_regrid
      implicit none
      character*150 :: filesource,filetarget,regridfile,cims,cjms,cnts,
     &     cimt,cjmt,cntt,cformat,cfields,fsource,ftarget,ctitle,
     &     clevels,cmaintitle,creal
      character*80 :: title
      real*4, allocatable :: sij4(:,:,:,:) , tij4(:,:,:,:)   ! ij arrays
      real*8, allocatable :: sij(:,:,:,:)  , tij(:,:,:,:)
      real*4, allocatable :: sijl4(:,:,:,:) , tijl4(:,:,:,:)  !ijl arrays
      real*8, allocatable :: sijl(:,:,:,:)  , tijl(:,:,:,:)
      real*4, allocatable :: slij4(:,:,:,:) , tlij4(:,:,:,:)  !lij arrays
      real*8, allocatable :: slij(:,:,:,:)  , tlij(:,:,:,:)
      integer :: ims,jms,nts,imt,jmt,ntt,nfields,nlevels,n,
     &     nargs,maxrec,iuin,iuout
      type (x_2grids) :: x2grids

      nargs = IARGC()
      IF(nargs.lt.15) write(*,*) "ll2cs needs 15 arguments";

      call getarg(1,filesource)
      call getarg(2,filetarget)
      call getarg(3,regridfile)
      call getarg(4,cims)
      call getarg(5,cjms)
      call getarg(6,cnts)
      call getarg(7,cimt)
      call getarg(8,cjmt)
      call getarg(9,cntt)
      call getarg(10,cformat)
      call getarg(11,cfields)
      call getarg(12,ctitle)
      call getarg(13,clevels)
      call getarg(14,cmaintitle)
      call getarg(15,creal)
      read(cims,'(I4)') ims
      read(cjms,'(I4)') jms
      read(cnts,'(I4)') nts
      read(cimt,'(I4)') imt
      read(cjmt,'(I4)') jmt
      read(cntt,'(I4)') ntt
      read(cfields,'(I4)') nfields
      read(clevels,'(I4)') nlevels

      write(*,*) filesource,filetarget,regridfile,cims,cjms,cnts,
     &     cimt,cjmt,cntt,cformat,cfields,ctitle,clevels,cmaintitle,
     &     creal

c***  Initialize exchange grid
      call init_regrid(x2grids,ims,jms,nts,imt,jmt,ntt,regridfile)

c***  Read data (must be in GISS format)
      fsource=trim(filesource)
      ftarget=trim(filetarget)

      iuin=200
      iuout=300
      maxrec=0

      open(iuin,FILE=fsource,FORM='unformatted', STATUS='old')
      open(iuout,FILE=ftarget,FORM='unformatted', STATUS='unknown')

      cformat = trim(cformat)
      ctitle = trim(ctitle)
      cmaintitle=trim(cmaintitle)

      if (cmaintitle .eq. 'yes' .or. cmaintitle .eq. 'YES'
     &        .or. cmaintitle .eq. 'y' .or. cmaintitle .eq. 'Y') then
         read(unit=iuin) title
         write(unit=iuout) title
      endif

      if (cformat .eq. 'ij' .or. cformat .eq. 'IJ' 
     &     .or. cformat .eq. 'giss' .or. cformat .eq. 'GISS') then

         allocate(sij(ims,jms,nts,nfields),sij4(ims,jms,nts,nfields))
         allocate(tij(imt,jmt,ntt,nfields),tij4(imt,jmt,ntt,nfields))         
         
         if (ctitle .eq. 'yes' .or. ctitle .eq. 'YES'
     &        .or. ctitle .eq. 'y' .or. ctitle .eq. 'Y') then
            do 
               if (creal.ne.'real8') then
                  read(unit=iuin,end=30) title,sij4
                  sij=sij4
                  do n=1,nfields
                     call do_regrid(x2grids,sij(:,:,:,n),tij(:,:,:,n))
                  enddo
                  tij4=tij
                  write(unit=iuout) title,tij4
               else
                  read(unit=iuin,end=30) title,sij
                  do n=1,nfields
                     call do_regrid(x2grids,sij(:,:,:,n),tij(:,:,:,n))
                  enddo
                  write(unit=iuout) title,tij
               endif

               maxrec=maxrec+1
            enddo
 30      continue

         else if (ctitle .eq. 'no' .or. ctitle .eq. 'NO'
     &        .or. ctitle .eq. 'n' .or. ctitle .eq. 'N') then
            do 
               if (creal.ne.'real8') then
                  read(unit=iuin,end=40) sij4
                  sij=sij4
                  do n=1,nfields
                     call do_regrid(x2grids,sij(:,:,:,n),tij(:,:,:,n))
                  enddo
                  tij4=tij
                  write(unit=iuout) tij4
               else
                  read(unit=iuin,end=40) sij
                  do n=1,nfields
                     call do_regrid(x2grids,sij(:,:,:,n),tij(:,:,:,n))
                  enddo
                  write(unit=iuout) tij
               endif

               maxrec=maxrec+1
            enddo
 40      continue
         endif

         deallocate(sij4, sij) 
         deallocate(tij4, tij)

         
      else if (cformat .eq. 'ijl' .or. cformat .eq. 'IJL') then
         write(*,*) "IJL"            
         allocate(sijl(ims,jms,nlevels,nts),sijl4(ims,jms,nlevels,nts))
         allocate(tijl(imt,jmt,nlevels,ntt),tijl4(imt,jmt,nlevels,ntt))         
         
         if (ctitle .eq. 'yes' .or. ctitle .eq. 'YES'
     &        .or. ctitle .eq. 'y' .or. ctitle .eq. 'Y') then
            do 
               if (creal.ne.'real8') then
                  read(unit=iuin,end=50) title,sijl4
                  sijl=sijl4
                  do n=1,nlevels
                     call do_regrid(x2grids,
     &                    sijl(:,:,n,:),tijl(:,:,n,:))
                  enddo
                  tijl4=tijl
                  write(unit=iuout) title,tijl4
               else
                  read(unit=iuin,end=50) title,sijl
                  do n=1,nlevels
                     call do_regrid(x2grids,
     &                    sijl(:,:,n,:),tijl(:,:,n,:))
                  enddo
                  write(unit=iuout) title,tijl
               endif
               maxrec=maxrec+1
            enddo
 50      continue
         else if (ctitle .eq. 'no' .or. ctitle .eq. 'NO'
     &        .or. ctitle .eq. 'n' .or. ctitle .eq. 'N') then
c            write(*,*) "NOTITLE IJL"
c            open(777,FILE='tc90',FORM='unformatted', STATUS='unknown')
            do 
               if (creal.ne.'real8') then
                  read(unit=iuin,end=60) sijl4
                  sijl=sijl4
                  do n=1,nlevels
                     call do_regrid(x2grids,
     &                    sijl(:,:,n,:),tijl(:,:,n,:))
                  enddo
                  tijl4=tijl
                  write(unit=iuout) tijl4
               else
                  read(unit=iuin,end=60) sijl
                  do n=1,nlevels
                     call do_regrid(x2grids,
     &                    sijl(:,:,n,:),tijl(:,:,n,:))
c                     sijl4(:,:,1,:)=sijl(:,:,n,:)
c                     title="toto"
c                  write(unit=777) title,sijl4(:,:,1,:)
                  enddo
                  write(unit=iuout) tijl
               endif
               maxrec=maxrec+1
            enddo
 60      continue
c         close(777)
         endif

         deallocate(sijl4, sijl) 
         deallocate(tijl4, tijl)
         
      else if (cformat .eq. 'lij' .or. cformat .eq. 'LIJ') then
         
         allocate(slij(nlevels,ims,jms,nts),slij4(nlevels,ims,jms,nts))
         allocate(tlij(nlevels,imt,jmt,ntt),tlij4(nlevels,imt,jmt,ntt))         
         
         if (ctitle .eq. 'yes' .or. ctitle .eq. 'YES'
     &        .or. ctitle .eq. 'y' .or. ctitle .eq. 'Y') then
            do 
               read(unit=iuin,end=70) title,slij4
               slij=slij4
               do n=1,nlevels
                  call do_regrid(x2grids,slij(n,:,:,:),tlij(n,:,:,:))
               enddo
               tlij4=tlij
               write(unit=iuout) title,tlij4
               maxrec=maxrec+1
            enddo
 70         continue

         else if (ctitle .eq. 'no' .or. ctitle .eq. 'NO'
     &           .or. ctitle .eq. 'n' .or. ctitle .eq. 'N') then
            open(777,FILE='tc90',FORM='unformatted', STATUS='unknown')
            do 
               if (creal.ne.'real8') then
                  read(unit=iuin,end=80) slij4
                  slij=slij4
                  do n=1,nlevels
                     call do_regrid(x2grids,
     &                    slij(n,:,:,:),tlij(n,:,:,:))
                  enddo
                  tlij4=tlij
                  write(unit=iuout) tlij4
               else
                  read(unit=iuin,end=80) slij
                  do n=1,nlevels
                     call do_regrid(x2grids,
     &                    slij(n,:,:,:),tlij(n,:,:,:))
                     slij4(1,:,:,:)=slij(1,:,:,:)
                     title="toto"
                     write(unit=777) title,slij4(1,:,:,:)
                  enddo
                  write(unit=iuout) tlij
               endif
               maxrec=maxrec+1
            enddo
 80         continue
         endif
         close(777)
         deallocate(slij4, slij) 
         deallocate(tlij4, tlij)
         
      endif

      close(iuin)
      close(iuout)

 
      write(6,*) "remapped file contains:"
      write(6,100) maxrec," records"
      write(6,200) "record format ",trim(cformat)
      write(6,200) trim(cfields)," fields per record"
      write(6,200) trim(clevels)," levels per array"
      write(6,200) "records contains title? ", ctitle
      write(6,200) "file contains main title? ", cmaintitle

 100  format(3X, I3, A)
 200  format(4X, A, A)


      end program ll2cs



