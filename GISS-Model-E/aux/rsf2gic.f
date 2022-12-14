      program extract_gic_from_rsf_file
C**** extract GIC records from an rsf_file
C**** compile with: gmake rsf2gic.o
C**** f90 -o rsf2gic rsf2gic.o ../model/ *.o -O2 -64 -mips4 \
C****                      -OPT:reorg_comm=off -w2 -listing
C**** Note that since it uses modules and routines from the model, it
C**** must be compiled after the model
      USE MODEL_COM, only : irsfic,iowrite
      IMPLICIT NONE
      CHARACTER infile*60, outfile*60
      INTEGER IARGC,n,N1,ioerr,kunit,iaction,it

      IF (IARGC().lt.1) THEN
        PRINT*,"Extract GIC_file from an rsf_file"
        PRINT*,"rsf2gic rsf_file_name"
        STOP
      END IF

      CALL GETARG(1,infile)

      N1=1
      if(index(infile,'/').gt.0) then
        do n=1,60
           if(infile(n:n).eq.'/') n1=n+1
        end do
      end if
      outfile='GIC.'//trim(infile(n1+12:60))//'.'//infile(n1:n1+7)

      kunit=9

      print*,trim(infile)

      iaction=irsfic
      ioerr=-1
      call io_rsf  (trim(infile),it,iaction,ioerr)
      if (ioerr.eq.1) then
        WRITE(6,*) "I/O ERROR IN RSF FILE: KUNIT=",kunit
       stop 'no action'
      end if

      OPEN(kunit,FILE=trim(outfile),
     *     FORM="UNFORMATTED",STATUS="UNKNOWN")

      iaction=iowrite
      ioerr=-1
      call io_ocean  (kunit,iaction,ioerr)
      call io_seaice (kunit,iaction,ioerr)
      call io_earth  (kunit,iaction,ioerr)
      call io_soils  (kunit,iaction,ioerr)
      call io_landice(kunit,iaction,ioerr)
      if (ioerr.eq.1) WRITE(6,*) "I/O ERROR IN GIC FILE: KUNIT="
     *     ,kunit
      close (kunit)

      print*,"New GIC file written out to ",trim(outfile)
      stop
      end
