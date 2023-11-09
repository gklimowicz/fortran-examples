c-----------------------------------------------------------------------------
      real c1,c2,c3,c4,c5,c6,c7,c8,c9,pref,sclkap,refprs
c
ccc      data c1,c2,c3,c4,c5,c6,c7/
c --- coefficients for sigma-0 (based on Brydon & Sun fit: -2:30,30:38)
ccc     . -1.36471E-01, 4.68181E-02, 8.07004E-01,-7.45353E-03,-2.94418E-03,
ccc    .  3.43570E-05, 3.48658E-05/
ccc      data pref/0./
c --- coefficients for sigma-2 (based on Brydon & Sun fit: -2:30,30:38)
css     .  9.77093E+00,-2.26493E-02, 7.89879E-01,-6.43205E-03,-2.62983E-03,
css     .  2.75835E-05, 3.15235E-05/
c     data c1,c2,c3,c4,c5,c6,c7,c8,c9/  ! T:[-2,32],S:[16,38] c1-c9 fit for sig2
c    .  9.887030E+00,-1.571861E-02, 7.827997E-01,-6.575396E-03,
c    . -2.913760E-03, 2.974136E-05, 3.248280E-05, 1.062302E-04,
c    .  3.524143E-06/
c     data c1,c2,c3,c4,c5,c6,c7,c8,c9/  ! T:[-2,30],S:[18,38] c1-c9 fit for sig2
c    .  9.893322E+00,-1.500683E-02, 7.823698E-01,-6.630868E-03,
c    . -2.936529E-03, 3.068372E-05, 3.343679E-05, 1.135806E-04,
c    .  3.620535E-06/

c --- for sigma1: Jackett et al. 2006, J. of atm & oceanic technology
      data c1,c2,c3,c4,c5,c6,c7,c8,c9/	! T:[-2,30],S:[30,38]
     .  5.019935E+00, 1.470290E-02, 7.904674E-01,-6.861022E-03,
     . -3.051459E-03, 3.099320E-05, 3.050299E-05, 1.046036E-04, 
     . 5.287567E-06/
      data pref/1.e7/

c --- for sigma2: Jackett et al. 2006, J. of atm & oceanic technology
c     data c1,c2,c3,c4,c5,c6,c7,c8,c9/	! T:[-2,30],S:[30,38]
c    .  9.929853E+00,-1.990743E-02, 7.813083E-01,-6.343472E-03,
c    . -2.883841E-03, 2.746757E-05, 2.897260E-05, 1.160426E-04, 
c    .  5.007679E-06/
c     data pref/2.e7/
c
c --- sub-coefficients for locally referenced sigma
c --- based on Jackett et al. 2006, J. of atm & oceanic technology.
      real*8, parameter, dimension(9) ::   ! T:[-2,30],S:[30,38], P:[0:4000] in dbar
     .  alphap = (/-6.313707E-03, 5.066326E-02, 7.999489E-01,   
     .             -7.403005E-03,-3.225958E-03, 3.471955E-05,   
     .              3.210467E-05, 9.259891E-05, 5.580351E-06 /)
     . ,betap  = (/ 5.083339E-03,-3.664164E-05,-9.637964E-06,
     .              5.540090E-07, 1.778390E-07,-3.824842E-09,
     .             -1.636328E-09, 1.227641E-08,-2.989838E-10 /)
     . ,gammap = (/-5.762706E-08, 6.781354E-10, 1.588334E-10,
     .             -1.212103E-11,-3.390138E-12, 9.942362E-14,
     .              3.514305E-14,-2.772681E-13, 6.323460E-15 /)
c
c> Revision history:
c>
c> Mar. 1999 - made thermobaric coefficients conform to Brydon & Sun
c> Nov. 2002 - updated thermobaric coefficients
c> May  2003 - further expansion of thermobaric options
c> July 2005 - added coeeficients for in-situ density (alphap,betap,gammap)
c-----------------------------------------------------------------------------
