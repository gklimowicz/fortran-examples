      block data blkprf
c
      include 'kprf_scalars.h'
c
      data dp00  /   3.0  / ! deep    z-level spacing minimum thickness (m)
      data locsig/ .false./ ! locally-referenced pot. density for stability
      data cb    /   3.e-3/ ! coefficient of quadratic bottom friction  
      data tmljmp/   0.2  / ! equivalent temperature jump across mixed-layer (degC)
      data rigr  /   0.25 / ! PWP:     critical gradient richardson number
      data ribc  /   0.65 / ! PWP:     critical bulk     richardson number
      data rinfty/   0.7  / ! KPP:     maximum  gradient richardson number (shear inst.)  
      data ricr  /   0.45 / ! KPP:     critical bulk     richardson number
      data bldmin/  20.0  / ! KPP:     minimum surface boundary layer thickness
      data bldmax/1200.0  / ! KPP:     maximum surface boundary layer thickness
      data cekman/   0.7  / ! KPP/KT:  scale factor for Ekman depth 
      data cmonob/   1.0  / ! KPP:     scale factor for Monin-Obukov depth
      data bblkpp/ .false./ ! KPP:     activate bottom boundary layer
      data shinst/ .true. / ! KPP:     activate shear instability mixing
      data dbdiff/ .true. / ! KPP:     activate double diffusion  mixing
      data nonloc/ .true. / ! KPP:     activate nonlocal b. layer mixing
      data latdiw/ .true. / ! K-PROF:  activate lat.dep. int.wave mixing
      data botdiw/ .true. / ! GISS:    activate bot.enhan.int.wav mixing
      data difsmo/ .false./ ! K-PROF:  activate horiz smooth diff coeffs
      data difm0 /  50.e-4/ ! KPP:     max viscosity   due to shear instability (m**2/s)  
      data difs0 /  50.e-4/ ! KPP:     max diffusivity due to shear instability (m**2/s)  
      data difmiw/   3.e-5/ ! KPP:     background/internal wave viscosity       (m**2/s)  
      data difsiw/   1.e-5/ ! KPP:     background/internal wave diffusivity     (m**2/s)  
      data dsfmax/  10.e-4/ ! KPP:     salt fingering diffusivity factor        (m**2/s)  
      data rrho0 /   1.9  / ! KPP:     salt fingering rp=(alpha*delT)/(beta*delS)         
      data cs    /  98.96 / ! KPP:     value for nonlocal flux term                       
      data cstar /  10.0  / ! KPP:     value for nonlocal flux term                       
      data cv    /   0.0  / ! KPP:     buoyancy frequency ratio (0.0 to use a fn. of N)   
      data c11   /   5.0  / ! KPP:     value for turb velocity scale                      
      data hblflg/   1    / ! KPP:     b. layer interpolation flag (1=lin.,2=quad.)       
      end
