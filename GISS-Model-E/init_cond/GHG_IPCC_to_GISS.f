      real conc_in(35), conc_out(6)
      character*80 filein

      call getarg(1,filein)
      open(1,file=filein,form='formatted')

      call getarg(2,filein)
      read(filein,*) nhead

      do n=1,nhead
      read(1,*)
      end do

      CFC12_to_CFC11   = .32/.25        ! conversion to CFC-11 equiv.
      HFC134a_to_CFC11 = .16/.25
      do
        read(1,*,end=100) iyr,conc_in
        conc_out(1) = conc_in(3) + .5   ! CO2
        conc_out(2) = conc_in(5)/1000.  ! N2O    ppb->ppm
        conc_out(3) = conc_in(4)/1000.  ! CH4    ppb->ppm
        conc_out(4) = conc_in(20)/1000. ! CFC-11 ppt->ppb
        conc_out(5) = conc_in(21)/1000. ! CFC-12 ppt->ppb
        conc_out(6) = 1e-3*             ! others (as CFC-11 equivalent)
     *   ( conc_in(6)*HFC134a_to_CFC11 +   ! extra Kyoto Protocol gases
     *     (conc_in(7)-conc_in(21))*CFC12_to_CFC11 - conc_in(20) )
                                           ! extra Montreal Prot. gases
        write(2,'(i4,f7.1,f7.4,4f7.3)') iyr,conc_out
      end do

  100 stop
      end
