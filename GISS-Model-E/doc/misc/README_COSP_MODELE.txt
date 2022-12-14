###############################################################################
                    COSP (CFMIP Observation Simulator Package)

Home Page: http://cfmip.metoffice.com/COSP.html
Code Base: http://code.google.com/p/cfmip-obs-sim/
COSP Documentation: 
    http://code.google.com/p/cfmip-obs-sim/downloads/detail?name=COSP_user_manual.v1.3.1.pdf

Current Version 1.3.2 (April 2011)

Overview: Bodas-Salcedo, A., and Coauthors, 2011: COSP: Satellite simulation 
    software for model assessment. Bull. Amer. Meteor. Soc., 92, 1023-1043.
    doi: 10.1175/2011BAMS2856.1

    From Bodas-Salcede et al. 2011

        COSP is a flexible software tool that enables the simulation from model
        variables of data from several satelliteborne active and passive sensors.
        It facilitates the use of satellite data to evaluate models in a 
        consistent way. The flexibility of COSP makes it suitable for use in many
        types of numerical models, from high-resolution models (~1-km 
        resolution) to coarse-resolution models, such as the GCMs used in climate
        modeling, and the scales in be- tween used in weather forecast and 
        regional models. The fact that COSP includes several simulators under the
        same interface facilitates the implementation of a range of simulators in
        models. Another advantage of COSP-and in general, the simulator 
        approach--is that it facilitates model intercomparison, not only 
        model-satellite comparison (e.g., comparisons of cloud properties
        simulated by GCMs and CRMs).
        ...
        The current version of COSP includes simulators for datasets produced by
        the following instruments (Table 1): the cloud profiling radar (CPR) on
        board CloudSat (Stephens et al. 2002), the Cloud-Aerosol Lidar with 
        Orthogonal Polarization (CALIOP) lidar onboard Cloud-Aerosol Lidar and 
        Infrared Pathfinder Satellite Observations (CALIPSO; Winker et al. 2010),
        the ISCCP (Rossow and Schiffer 1999), the MISR (Diner et al. 2005), and 
        the Moderate Resolution Imaging Spectroradiometer (MODIS) (King et al. 
        2003).

    * Instrument simulators:
      ISCCP simulator
      Radar simulator (QuickBeam)
      Lidar and PARASOL simulators (ACTSIM)
      MISR simulator
      MODIS simulator
      RTTOV simulator

Created: Nov 2011 Mike Bauer
Updated:
###############################################################################

INSTALL:

    Step: Where to install. 
    
    To minimize complications between the code repositories for modelE and COSP
    it is not advisable to clone the COSP code base directly into that of modelE.
    Rather, checkout COSP to a directory outside of the modelE code base
    and move the needed files into modelE. 
    
    *************************************************************************** 
    *                              IMPORTANT NOTE                             *
    * Please be aware that COSP was already a reserved term in modelE prior   *
    * to the introduction of COSP.                                            *
    *                                                                         *
    *  !@var  COSP, COSV cos of latitude at primary, secondary latitudes      *
    *                                                                         *
    * As a result, CFMIP is used to refer to COSP when referring to COSP from *
    * the context of using it inside of modelE. Conversely, COSP is used when *
    * referring to the stand alone COSP product.                              *
    ***************************************************************************

    Step: Acquire a current version of modelE

        See "Check out, compile and run ModelE" 
        https://modelingguru.nasa.gov/docs/DOC-1616

        > git clone username@simplex.giss.nasa.gov:/giss/gitrepo/modelE.git 

        This should create the $MHOME/modelE in the current directory (referred to
        here as $MHOME).

        Switch to the appropriate branch
        > cd modelE
        > git checkout branch_name

    Step: Acquire a current version of COSP from the modelE code base.
    
        This should be done from $MHOME as well, creating the directory CFMIP.

        > git clone username@simplex.giss.nasa.gov:/giss/gitrepo/modelE.git

        Note: See below for installation/update procedures of COSP itself.
    
    Step: Modified modelE files:
        /model/CLOUDS2_COM.f
        /model/CLOUDS2.f
        /model/CLOUDS2_DRV.f

        COSP has some settings that alter what sort of information it needs
        from modelE (see initializing COSP below). The following modifications
        allow for some choices by the user to be automatically implemented.

        There are a number of COPS options in the rundeck. 
        
        Define CFMIP preprocessor directives, evoke SUBDD (to save output) and
        trigger CFMIP as a model component.

        Preprocessor Options:
        #define CFMIP       - Enable COSP
        #define CFMIP_PFLUX ! Use precipitation flux (otherwise hydrometeor mixing ratio)
        #define CFMIP_USERE ! Use model hydrometeor effective radii (otherwise COSP calculated)
        #define CACHED_SUBDD ! Required for COSP

        CFMIP_PFLUX - Use precipitation flux (otherwise hydrometeor mixing ratio)
        CFMIP_USERE - Use model hydrometeor effective radii (otherwise COSP calculated)
        CACHED_SUBDD - Used to store output

        Object modules:
        SUBDD                               ! call subdd

        Components:
        ../../CFMIP

        Point to the proper namelist files (be sure this are consistent with
        the preprocessor options above!). Again, you may have to use a template 
        and create a special set of namelist files for your run. 

        Data input files:
        CFMIPIN=cosp_input_nl.txt ! CFMIP/COSP options and data
        CFMIPOUT=cosp_output_nl.txt ! CFMIP/COSP options and data
        
        CFMIPIN - Namelist tells COSP about the model and how COSP is to process
            the data:
            + NPOINTS
                The number of horizontal grid points of the model (e.g., jm*im).
                Example:
                    NPOINTS = 12960 2 x 2.5 degree grid [144,90]
            + NPOINTS_IT
                Controls the *maximum* number of grid points to be processed per 
                iteration of the COSP simulator. Lower values reduce COSP memory
                usage but increase processing time.
                Example:
                    NPOINTS_IT = 1000 or im (whole longitude)   
            + NCOLUMNS
                Sets the number of subcolumns in SCOPS (Subgrid Cloud Overlap
                Profile Sampler). Lower values reduce both COSP memory usage 
                and processing time.
                
                The COSP documentation in one place recommends that NCOLUMNS
                be ~ model resolution (in degrees) x 100, but not less than
                ~50. That is, for a 1x1 deg model => NCOLUMN=100.

                However, in another place it says:                 
                    The simulator uses a Monte Carlo method for sampling
                    various columns within each model gridbox. The number
                    of columns is set by the value of NCOLUMNS. The value
                    that you want to set NCOLUMNS to depends on the
                    accuracy you want and amount of averaging you are doing
                    on the outputs.

                    The recommended rule of thumb recommended by the
                    authors is that you should aim for something like 2400
                    samples to keep statistical noise to a reasonable
                    level.

                    For example, if you are doing no averaging, (i.e. you
                    are be calling the simulator once on instantaneous
                    model variables and looking directly at the results),
                    you should set might expect that you need to set NCOLUMNS
                    to something around 2400.

                    If you are looking at daily means, and are calculating
                    this by averaging 8 3-hourly calls to the simulator,
                    NCOLUMNS should be set to 300 (2400/8).

                    If, say, you are looking at monthly means, and are
                    calling the simulator, say, every 15 hours, NCOLUMNS should
                    be set to 50 (2400/(24*30/15)).

                    If you are looking at monthly means, and are
                    calculating this by averaging 8 3-hourly calls per day to
                    the simulator, NCOLUMNS should be set to 10 (2400/(8*30)).
                Example:
                    NCOLUMNS = 200 for 2 x 2.5 degree grid
            + NLEVELS: Sets the number of full model levels.
                Example:
                    NLEVELS =  40 for modelE F40.

            The following settings concerning the way statistical outputs from 
            COSP are handled in the vertical dimension.

            + NLR: Sets the number of levels used when USE_VGRID is .true..
                Example:
                    NLR = 40 for CFMIP-2
            + USE_VGRID: If .true. COSP output is on a fixed evenly spaced 
                altitude grid of size NLR. If .false. NLEVELS (i.e., model 
                levels) is used instead.
                Example:
                    USE_VGRID = .true. for CFMIP-2
            + CSAT_VGRID: If .true. the CloudSat standard grid of 40 evenly 
                spaced levels (480 m) is used for relevant statistical outputs.
                Only used if USE_VGRID is also .true.
                Example:
                    CSAT_VGRID = .true. for CFMIP-2   

            The following settings are related to the simulations. Only the 
            options likely to be altered by the user are listed here. See the
            COSP documentation for any remaining questions.
            
            + USE_MIE_TABLES: Use precomputed mie lookup tables. CFMIP suggests
                that not to use this option.
                Example:
                    USE_MIE_TABLES = 0 for CFMIP-2 

            + USE_REFF: Use effective radius in the radar simulator. This should
                always be turned on. Note that the preprocessor command 
                CFMIP_USERE controls what effective radii are used by COSP. In
                short the options are to either give COSP model derived values or
                let COSP calculate them using a complicated set of assumptions
                and the model data it does have (See HCLASS Table). Experiments 
                with modelE suggest that using the model effective radii works 
                best, with the optimum setting being USE_REFF=.true. and 
                USE_MODEL_RE enabled.
                Example:
                    USE_REFF = .true. for CFMIP-2
            + USE_PRECIPITATION_FLUXES: Pass model precipitation fluxes into 
                COSP where they are then converted into hydrometeor mixing ratios
                and then derived radar reflectivities. This involves a number
                of assumptions for each hydrometeor type (See HCLASS Table). The
                alternative is to provide the model hydrometeor mixing ratios 
                directly and then derive the radar reflectivities. However, even
                in this case COSP uses the assumptions for each hydrometeor 
                type (See HCLASS Table) to determine the discrete drop size 
                distribution. Experiments with modelE suggest that allowing COSP
                to do all the conversions (i.e., use precipitation fluxes) gives
                a better result. This may be because the resulting mixing ratios
                are more consistent with the assumptions that go into creating
                the discrete drop size distribution and radar reflectivities. 
                ***Be sure that USE_PFLUX and USE_PRECIPITATION_FLUXES agree***
                Example:
                    USE_PRECIPITATION_FLUXES = .true. for CFMIP-2
            + OVERLAP: Controls the subgrid vertical cloud overlap distribution 
                assumed by SCOPS (Subgrid Cloud Overlap Profile Sampler).
                    1 = Maximum overlap
                    2 = Random overlap
                    3 = Maximum/random overlap
                Maximum overlap is applied to the convective cloud, and 
                maximum/random is used for large-scale cloud.
                Example:
                    OVERLAP = 3 for CFMIP-2
            + ISCCP_TOPHEIGHT: Controls cloud top calculations in the ISCCP 
                simulator.
                    1 = Adjust top height using both a computed infrared 
                        brightness temperature and the visible optical depth to 
                        adjust cloud top pressure. Note that this calculation is
                        most appropriate to compare to ISCCP data during sunlit
                        hours.
                    2 = Do not adjust top height. That is, cloud top pressure is
                        the actual cloud top pressure in the model.
                    3 = Adjust top height using only the computed infrared 
                        brightness temperature. Note that this calculation is 
                        most appropriate to compare to ISCCP IR only algortihm
                        (i.e. you can compare to nighttime ISCCP data with this
                        option).
                Example:
                    ISCCP_TOPHEIGHT = 1 for CFMIP-2
            + ISCCP_TOPHEIGHT_DIRECTION: Controls the direction for finding 
                atmosphere pressure level with interpolated temperature equal to
                the radiance determined cloud-top temperature in the ISCCP 
                simulator.
                    1 = Find the *lowest* altitude (highest pressure) level with
                        interpolated temperature equal to the radiance 
                        determined cloud-top temperature.
                    2 = Find the *highest* altitude (lowest pressure) level with
                        interpolated temperature equal to the radiance 
                        determined cloud-top temperature. This is the default 
                        value since V4.0 of the ISCCP simulator.
                Only applicable if ISCCP_TOPHEIGHT is 1 or 3
                Example:
                    ISCCP_TOPHEIGHT_DIRECTION = 2 for CFMIP-2

        CFMIPOUT - Namelist that tells COSP the simulations to run and variables
            to save.

            ! Simulator flags
            
            + Lradar_sim use CloudSat
            + Llidar_sim use CALIPSO PARASOL
            + Lisccp_sim use ISCCP
            + Lmisr_sim  use MISR
            + Lmodis_sim use MODIS
            + Lrttov_sim use RTTOV

            ! Output variables
            
            !- ISCCP
            + Lalbisccp         Mean Cloud Albedo
            + Lboxptopisccp     Cloud Top Pressure in Each Column
            + Lboxtauisccp      Optical Depth in Each Column
            + Lpctisccp         Mean Cloud Top Pressure
            + Lclisccp          Cloud Area Fraction
            + Ltauisccp         Mean Optical Depth
            + Lcltisccp         Total Cloud Fraction
            + Lmeantbisccp      Mean all-sky 10.5 micron brightness temperature
            + Lmeantbclrisccp   Mean clear-sky 10.5 micron brightness temperature            

            !- MISR
            + LclMISR           Cloud Fraction

            !- MODIS
            + Lcltmodis         Total Cloud Fraction
            + Lclwmodis         Liquid Cloud Fraction
            + Lclimodis         Ice Cloud Fraction
            + Lclhmodis         High Level Cloud Fraction
            + Lclmodis          Mid Level Cloud Fraction
            + Lcllmodis         Low Level Cloud Fraction
            + Ltautmodis        Total Cloud Optical Thickness
            + Ltauwmodis        Liquid Cloud Optical Thickness
            + Ltauimodis        Ice Cloud Optical Thickness
            + Ltautlogmodis     Total Cloud Optical Thickness (Log10 Mean)
            + Ltauwlogmodis     Liquid Cloud Optical Thickness (Log10 Mean)
            + Ltauilogmodis     Ice Cloud Optical Thickness (Log10 Mean)
            + Lreffclwmodis     Liquid Cloud Particle Size
            + Lreffclimodis     Ice Cloud Particle Size
            + Lpctmodis         Cloud Top Pressure
            + Llwpmodis         Cloud Liquid Water Path
            + Liwpmodis         Cloud Ice Water Path

            !- CALIPSO
            + Latb532           Lidar Attenuated Total Backscatter (532 nm)
            + LcfadLidarsr532   Scattering Ratio (Cloud Frequency Altitude Diagrams)
            + Lclcalipso        Cloud Area Fraction
            + Lclhcalipso       High Level Cloud Fraction
            + Lclmcalipso       Mid Level Cloud Fraction
            + Lcllcalipso       Low Level Cloud Fraction
            + Lcltcalipso       Total Cloud Fraction
            + LparasolRefl      PARASOL Reflectance
            !- Use lidar and radar
            + Lclcalipso2       CALIPSO Cloud Fraction Undetected by CloudSat
            + Lcltlidarradar    Lidar and Radar Total Cloud Fraction

            !- CloudSat
            + Lcfaddbze94       Radar Reflectivity (Cloud Frequency Altitude Diagrams)
            + Ldbze94           Radar Reflectivity            

        There are templates for these namelist files in COSP_GISS_Tools which
        can be copied and modified as needed and placed in /cmrun and denoted
        in the rundeck with the CFMIPIN and CFMIPOUT noted above.

        NOTE: Due to limits with the SUBDD sub-daily diagnostics routine some of
        the above output is not stored in its native shape. That is, SUBDD only
        allows for 2 and 3D arrays, whereas some COSP arrays are 4 or 5D. For these
        we simply collapse the extra dimensions down to 3D, which then requires
        post-processing. Also, file size limitations with netcdf necessitate that
        COSP data is stored on a daily basis, rather than monthly. This depends
        on the model resolution and the time frequency at which COSP is called.
        Also, the primary COSP output variables that will cause these problems
        are LparasolRefl, LcfadLidarsr532 and especially Lcfaddbze94. In cases
        where these variables are not requested and/or the COSP sampling frequency
        is low SUBDD can be safely directed to store monthly output (see 
        write_daily_files in the rundeck).

        Alter the parameters. Note do not run COSP with the old isccp simulator
        enabled (isccp_diags=0). COSP requires a non-zero sub-daily 
        diagnostics NSUBDD but saves at its own Nsubdd_for_cfmip
        
        ! parameters that affect at most diagn. output:
        SUBDD=' '       ! no sub-daily frequency diags
        NSUBDD=24        ! saving sub-daily diags every NSUBDD*DTsrc/3600. hour(s)
        Nsubdd_for_cfmip=24 ! saving sub-daily COSP diags every NSUBDD-th physics time step (1/2 hr)            
        isccp_diags=0   ! use =0 to save cpu time, but you lose some key diagnostics
        write_daily_files=1 ! SUBDD saves daily rather than monthly.

================================================================================

    Possible Alterations to COSP (cosp_constants.F90 and HCLASS Table). 
        COSP makes some assumptions about the nature of the hydrometeor 
        microphysics/size distribution. These can be controls by a number of
        settings in the cosp_constants.F90 file. These value almost exclusively
        impact the radar and lidar simulations.

        Key to Hydrometeor Type below.
        -----------------------------
        LSL = large_scale_cloud_liquid
        LSI = large_scale_cloud_ice
        LSR = large_scale_rain
        LSS = large_scale_snow
        CVL = convective_cloud_liquid
        CVI = convective_cloud_ice
        CVR = convective_cloud_rain
        CVS = convective_cloud_snow
        LSG = large_scale_cloud_graupel

        Step 1: In cosp_constants.F90
        ^^^^^^
        
        Set the distribution type for each hydrometeor.
        
        This is done with the HCLASS_TYPE table.
        
        HCLASS_TYPE = 
            Set to 1 for modified gamma distribution, 
                   2 for exponential distribution, 
                   3 for power law distribution, 
                   4 for monodisperse distribution, 
                   5 for lognormal distribution. 
                  -1 to ignore this hydrometeor.
        
        ****REQUIRE MODELE VALUES****
        Default COSP Values:
                        LSL    LSI    LSR   LSS   CVL   CVI   CVR   CVS   LSG
        HCLASS_TYPE      5,     1,     2,    2,    5,    1,    2,    2,    2
        HCLASS_PHASE     0,     1,     0,    1,    0,    1,    0,    1,    1

        These defaults have cloud liquid governed by a lognormal distribution and
        cloud ice by a modified gamma distribution. All precipitation follows an
        exponential distribution.
        
        * Setting anything to -1 means no radar data at all for that. Even when 
            a non-zero mixing ratio and an effective radius is provided for that
            hydrometeor.
        * Because GISS lacks LSG this has be defaulted to -1.
        * Phase (0 - liquid, 1 - ice).
                   
        Step 2: 
        ^^^^^^
        For hydrometers with power law distributions (i.e., HCLASS_TYPE=3) a
        minimum/maximum drop size limits (um) is required.

        Default COSP Values:
                       LSL    LSI    LSR   LSS   CVL   CVI   CVR   CVS   LSG        
        HCLASS_DMIN     -1,    -1,    -1,   -1,   -1,   -1,   -1,   -1,   -1
        HCLASS_DMAX     -1,    -1,    -1,   -1,   -1,   -1,   -1,   -1,   -1
        
        Step 3: 
        ^^^^^^
        Hydrometeor mass must be specified as either a constant value that applies
        to all particles of a given class or expressed as a function of particle 
        diameter (see Eq 4 in COSP_user_manual).
        
        HCLASS_APM: The alpha_x coefficient in equation 4 in the COSP Manual [kg m^(Beta*m)]
        HCLASS_BPM: The beta_x coefficient in equation 4 in the COSP Manual [unitless]
        HCLASS_RHO: Alternate constant density hydrometeor [kg m^-3]
        
        The Mass-Diameter Relationship:
            D = Partical Diameter
            Partical Mass() = HCLASS_APM*D^(HCLASS_BPM)
            
            For liquid HCLASS_BPM close to 3. For ice, < 3.
        
        To specify a constant density for particles, set HCLASS_APM and HCLASS_BPM
        to -1 and HCLASS_RHO to the constant value desired. Then in dsd.F90
            apm = (pi/6)*rho_c
            bpm = 3.

        To let mass vary as a function of diameter, specify values for HCLASS_APM 
        and HCLASS_BPM and set HCLASS_RHO to -1. 
        
        ** HCLASS_APM and HCLASS_BPM are used for all HCLASS_TYPE distributions and
            are used to determine the discrete drop size distribution (calling dsd 
            from dsd.F90) in radar_simulator.f90 **
        
        Default COSP Values:
                       LSL    LSI    LSR   LSS   CVL   CVI   CVR   CVS   LSG          
        HCLASS_APM/    524, 110.8,   524,   -1,  524,110.8,  524,   -1,   -1
        HCLASS_BPM/     3,   2.91,     3,   -1,    3, 2.91,    3,   -1,   -1
        HCLASS_RHO/    -1,     -1,    -1,  100,   -1,   -1,   -1,  100,  400
        
        The default values assume that liquid particles are spherical drops with 
        a density equal to that for liquid water (based on volume of a sphere 
        4/3*Pi*(D/2)^3). All ice particles use a lesser value. Snow on the 
        other hand, uses a constant density.  

        Step 4: 
        ^^^^^^
        
        COSP requires an effective radius as input for CALIPSO and CloudSat. These
        can either be provided by the model or calculated by COSP. The calculated
        values are are 30 um for the lidar, and the values defined in HCLASS_P1 for
        CloudSat. If the model effective radii are passed as an input parameter 
        the values in HCLASS_P1, HCLASS_P2, HCLASS_P3 are used instead to determine
        the proper discrete drop size distribution. That is, in all cases some
        assumptions are made by COSP when finding the discrete drop size 
        distribution.

        For HCLASS_TYPE = 1
            The modified gamma distribution requires that one of the parameters 
            (HCLASS_P1,HCLASS_P2) be specified; the other should be set to -1. 
            The user must also specify a value for HCLASS_P3.
        
            HCLASS_P1 - Sets the total particle number concentration 
                N_t/rho_a (kg^-1), where rho_a is the density of air in the radar volume.
        
            HCLASS_P2 - Sets the particle mean diameter D (um).
        
            HCLASS_P3 - Sets the distribution width, a_x + 1.
        
        For HCLASS_TYPE = 2
            The exponential distribution requires that one of the parameters 
            (HCLASS_P1, HCLASS_P2, HCLASS_P3) be specified; the remaining two 
            must be set to -1.

            HCLASS_P1 - Sets a constant intercept parameter N_0(m^-4) and the
                the slope parameter lambda is calculated.
        
            HCLASS_P2 - Sets the slope parameter lambda (um^-1) and intercept parameter
                is calculated.
        
            HCLASS_P3 - Set to 2 to indicate that lambda should be evaluated as a 
                function of temperature. Only useful for ice particles.

        For HCLASS_TYPE = 3
            The power law distribution requires that only HCLASS_P1 be specified and
            HCLASS_P2, HCLASS_P3 be set to 0 (NOT -1). It is critical that the user 
            specify reasonable values for HCLASS_DMIN and HCLASS_DMAX.
        
            HCLASS_P1 - Either set this to the value of a constant power law 
                parameter b_r, set to -2 to evaluate b_r according to the method of 
                Ryan (2000) for cirrus type clouds, or set to -3 to evaluate with
                the same method but for frontal type clouds. This method is useful
                only for ice particles.

        For HCLASS_TYPE = 4 
            The monodisperse distribution sets particles to a uniform size and  
            concentration. Only HCLASS_P1 is specified while HCLASS_P2 and 
            HCLASS_P3 are set to 0 (NOT -1).
     
            HCLASS_P1 - Set to a constant diameter D_0 (um) 
        
        For HCLASS_TYPE = 5
            The lognormal distribution, defined in terms of particle radius r 
            rather than diameter for consistency with the CloudSat 2B-LWC algorithm,
            requires that one of the parameters (HCLASS_P1,HCLASS_P2) be specified; 
            the others should be set to -1. The user must also specify a value for
            HCLASS_P3.
     
            HCLASS_P1 - Sets the total particle number concentration 
                N_t/rho_a (kg^-1), where rho_a is the density of air in the radar volume.
        
            HCLASS_P2 - Sets the the geometric mean particle radius r_g (um).
        
            HCLASS_P3 - Sets the natural logarithm of the geometric standard 
                deviation, ln(sigma_g).
     
        Default COSP Values:
                       LSL    LSI    LSR   LSS   CVL   CVI   CVR   CVS   LSG       
        HCLASS_TYPE     5,     1,     2,    2,    5,    1,    2,    2,    2
        HCLASS_P1/     -1,    -1,  8.e6, 3.e6,   -1,   -1, 8.e6, 3.e6, 4.e6
        HCLASS_P2/      6,    40,    -1,   -1,    6,   40,   -1,   -1,   -1
        HCLASS_P3/    0.3,     2,    -1,   -1,  0.3,    2,   -1,   -1,   -1    

        Step 5: 
        ^^^^^^
        
        COSP, if use_precipitation_fluxes, requires these settings to convert the
        flux into mixing ratios (by COSP_PRECIP_MXRATIO in cosp_utils.F90).

        Default COSP Values:
        ! Microphysical settings for the precipitation flux to mixing ratio conversion
                         LSL    LSI       LSR       LSS   CVL    CVI       CVR       CVS      LSG
        data N_ax/       -1.,   -1.,     8.e6,     3.e6,  -1.,   -1.,     8.e6,     3.e6,     4.e6/
        data N_bx/       -1.,   -1.,      0.0,      0.0,  -1.,   -1.,      0.0,      0.0,      0.0/
        data alpha_x/    -1.,   -1.,      0.0,      0.0,  -1.,   -1.,      0.0,      0.0,      0.0/
        data c_x/        -1.,   -1.,    842.0,     4.84,  -1.,   -1.,    842.0,     4.84,     94.5/
        data d_x/        -1.,   -1.,      0.8,     0.25,  -1.,   -1.,      0.8,     0.25,      0.5/
        data g_x/        -1.,   -1.,      0.5,      0.5,  -1.,   -1.,      0.5,      0.5,      0.5/
        data a_x/        -1.,   -1.,    524.0,    52.36,  -1.,   -1.,    524.0,    52.36,   209.44/
        data b_x/        -1.,   -1.,      3.0,      3.0,  -1.,   -1.,      3.0,      3.0,      3.0/
        data gamma_1/    -1.,   -1., 17.83725, 8.284701,  -1.,   -1., 17.83725, 8.284701, 11.63230/
        data gamma_2/    -1.,   -1.,      6.0,      6.0,  -1.,   -1.,      6.0,      6.0,      6.0/
        data gamma_3/    -1.,   -1.,      2.0,      2.0,  -1.,   -1.,      2.0,      2.0,      2.0/
        data gamma_4/    -1.,   -1.,      6.0,      6.0,  -1.,   -1.,      6.0,      6.0,      6.0/
        
        OVERVIEW and GUIDANCE:

        Experiments with modelE suggest that the best option is to use the 
        model's implied precipitation flux and calculated effective radii. 
        
        The CloudSat simulator is most sensitive to how the ice phase is handled.
        For example, turning off all liquid HCLASS_TYPEs has minimal impact
        on the resultant CloudSat data.
        
        Setting all HCLASS_TYPEs to follow a modified gamma distribution, rather
        than just the cloud ice as is the default, for example returns a very 
        different CloudSat picture than does setting all HCLASS_TYPEs to follow
        a lognormal distribution, rather than just the cloud liquid as is the 
        default, and neither of these closely resembles the CloudSat picture 
        using the defaults.

        On the other hand, CloudSat's sensitivity to precipitation means that
        setting all HCLASS_TYPEs to follow an exponential distribution, rather
        than just precipitation as is the default, returns CloudSat picture 
        that closely resembles that from using the defaults.

================================================================================
    How to acquire the COSP-modelE interface and tools.

        Included are scripts for updating the COSP code base, i.e., if there is
        a newer version of COSP available, and other tools for working with COSP
        output.

        These are kept in the directory $MHOME/COSP_GISS_Tools.

         > mkdir COSP_GISS_Tools

        *FIX* Where are these to be stored and made available?

        Acquire a current version of COSP from within $MHOME/COSP_GISS_Tools, 
        creating the directory COSP
        
        svn checkout http://cfmip-obs-sim.googlecode.com/svn/stable/current/trunk COSP
    
        Copy the necessary COSP files into the modelE code base. This is easiest
        by calling the python script COSP_GISS_Tools/move_cosp_to_modele.py 
        which does the following:

            > python COSP_GISS_Tools/move_cosp_to_modele.py COSP CFMIP

            a) Make a CFMIP directory $MHOME (assuming one is not there, otherwise
               the current files are **overwritten**!)
            b) Copy needed COSP files over to CFMIP. 

            Necessary COSP files (as of version 1.3.2):
                cosp_types.F90
                cosp_constants.F90
                cosp_utils.F90
                radar_simulator_types.f90
                load_hydrometeor_classes.f90
                cosp_modis_simulator.F90
                modis_simulator.F90
                scops.f
                prec_scops.f
                cosp_simulator.F90
                cosp_radar.F90
                gases.f90
                zeff.f90
                pf_to_mr.f
                array_lib.f90
                atmos_lib.f90
                radar_simulator.f90
                dsd.f90
                format_input.f90
                math_lib.f90
                optics_lib.f90
                mrgrnk.f90
                cosp_isccp_simulator.F90
                icarus.f
                cosp_stats.F90
                cosp_misr_simulator.F90
                MISR_simulator.f
                cosp.F90
                cosp_lidar.F90
                lidar_simulator.F90
                congvec.f
                llnl_stats.F90
                lmd_ipsl_stats.F90
                cosp_defs.h                       
            * NOTES
                - CFMIP ignores the subdirectory structure native to COSP.
                - COSP uses a mixture of *.f90 and *.F90 extensions, these are 
                  converted to *.F90.
                - COSP_GISS_Tools/Makefile is added to CFMIP. This file allows 
                  modelE to compile COSP as a component (CFMIP).
                - COSP_GISS_Tools/cfmip_drv.F90 is added to CFMIP. This file 
                  contains the COSP-modelE driver/interface. 

            Modify cosp_constants.F90 HCLASS Table. At the very least set the
            LSG (Large-scale grauple) element to -1 as modelE does not have 
            this hydrometeor type. See below for other possible alterations. 
            This is done automatically by move_cosp_to_modele.py.

                           LSL    LSI    LSR   LSS   CVL   CVI   CVR   CVS   LSG
            HCLASS_TYPE      5,     1,     2,    2,    5,    1,    2,    2,    -1
================================================================================
================================================================================
================================================================================
