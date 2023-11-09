#include "rundeck_opts.h"

!@auth M. Kelley
!@ver 1.0
!@sum  subdd_mod and subroutines in this file provide functionality
!@+    facilitating the output of model diagnostics at sub-daily
!@+    frequencies.  The model currently writes the diagnostics
!@+    to disk in "raw" form; the scaleacc utility must be used
!@+    to extract groups of outputs (see "Terminology" below for the
!@+    definition of "group").  Model diagnostics not defined on the
!@+    model's horizontal grid are not supported. Currently, only the
!@+    atmospheric model can use this package.
!@+
!@usage Two interfaces are available, for
!@+       I. outputs whose groupings and other metadata are defined
!@+          during model initialization
!@+      II. outputs whose registration is deferred until the stage
!@+          of model execution at which the first time-slice of the
!@+          output data is saved/accumulated.
!@+
!@+     The type II interface is intended for outputs that do not
!@+     naturally fall into one of the pre-existing categories of
!@+     type I, and for which the creation of a new category
!@+     would yield little or no long-term benefit.  For example,
!@+     special (temporary) diagnostics to investigate parameterization-
!@+     specific factors determining the behavior of a particular model
!@+     component are easily introduced as type II, as the deferred
!@+     registration permits a more compact coding footprint than I.
!@+     For standard fields (e.g. SLP, SAT, precipitation), type I
!@+     is easier to use over time, via its single set of rundeck
!@+     parameters.  Automatic vertical regridding of outputs to
!@+     constant-pressure levels is currently only possible via I
!@+     (will be added to II soon).
!@+
!@+     Requests for type I outputs are made through rundeck
!@+     strings SUBDD, SUBDD1, ..., following the traditional
!@+     subdaily diagnostics framework.  However, the parsing
!@+     of the requests has been extended to allow optional
!@+     specification of the output frequency Nsubdd for individual
!@+     fields using the syntax varname:NN where integer NN is the
!@+     value of Nsubdd to be in effect for varname.  Requests for
!@+     which NN is not specified will default to the global Nsubdd.
!@+     By default, outputs are time-averaged over each time interval;
!@+     instantaneous output can be requested by appending a trailing
!@+     "i" (e.g. varname:i or varname:NNi).  For a time-averaged
!@+     quantity, the number of timesteps per output file (Nday*days_per_file)
!@+     must be divisible by its Nsubdd.
!@+
!@+     To predefine a new output field, two steps are sufficient when
!@+     it is appropriate to place that field in a pre-existing output
!@+     category (see Terminology for the definition of "category"):
!@+       (1)  Specify metadata for that field in the declaration
!@+            procedure for the pre-existing category.  This is
!@+            done by setting the elements of the info_type structure.
!@+       (2)  In the model component in which the output field is
!@+            available, add a section to a CACHED_SUBDD "select case"
!@+            block corresponding to the pre-existing category.
!@+            These blocks pass gridded arrays to the inc_subdd
!@+            procedure, which copies the contents of the arrays
!@+            into a master database of all outputs (incrementing
!@+            accumulated quantities when appropriate).  See the
!@+            next paragraph for an example of such a block.
!@+     To define a new category, it is necessary to
!@+       (1)  Create a new procedure to define the metadata
!@+            for its output fields, and invoke the procedure
!@+            in the appropriate section of parse_subdd (where
!@+            the calls to ijh_defs et al. are made).
!@+       (2)  In the part(s) of the model in which outputs in this
!@+            category are to be collected, add a do-loop and
!@+            select-case block like the following for category
!@+            'xyz' (there are many examples to clone):
!@+                integer ngroups, grpids(subdd_ngroups)
!@+                type(subdd_type), pointer :: subdd
!@+                ! find groups belonging to category 'xyz'
!@+                call find_groups('xyz',grpids,ngroups)
!@+                do igrp=1,ngroups ! loop over groups in this category
!@+                   subdd => subdd_groups(grpids(igrp)) ! for brevity
!@+                   do k=1,subdd%ndiags ! loop over diags in this group
!@+                     select case (subdd%name(k))
!@+                     case ('some_var') ! array already available
!@+                       call inc_subdd(subdd,k,some_array)
!@+                     case ('complicated_var') ! need to compute output
!@+                       temp_array = (....)
!@+                       call inc_subdd(subdd,k,temp_array)
!@+                     end select
!@+                   end do
!@+                end do
!@+
!@+     Type II outputs are declared/saved/accumulated as follows:
!@+       call inc_subdd(
!@+                 'some_var',       ! name of this field
!@+          real*8  some_array,      ! data to be saved/accumulated
!@+          integer Nsubdd,          ! output frequency for this field
!@+          logical instant,         ! T for snapshots, F for averages
!@+         optional units='xxxx',    ! units specification
!@+         optional long_name='xxxx' ! description
!@+                     )
!@+     The metadata and control parameters are only examined during
!@+     the FIRST call to inc_subdd for a given output field.  No coding
!@+     is in place yet to check whether they were changed afterward.
!@+     Arrays or rank 2 and 3 are accepted.  In the latter case, the
!@+     size of this extra dimension can be automatically determined,
!@+     but if it does not correspond to the 3rd index of the array it is
!@+     necessary to specify (via optional argument jdim) which of the
!@+     dimensions is the 2nd horizontal dimension.  Optional arguments
!@+     dim3name and coordvalues can be used to name this extra dimension
!@+     in output files and provide a coordinate axis if appropriate.
!@+       call inc_subdd(
!@+                  'some_var',
!@+                  ......
!@+                  jdim=3,          ! only needed if extra dim is the first
!@+                  dim3name='zzz',  ! name instead of size_some_var
!@+                  coordvalues=zzz  ! 1D array of length size(dim3)
!@+                     )
!@+     The extra dimension need not be "vertical"; in some circumstances
!@+     it may be useful to bundle collections of 2D fields into 3D output
!@+     arrays.
!@+     Currently, each type II output variable constitues a separate
!@+     group (of size 1) and hence will be placed in its own output
!@+     file by scaleacc.  Another consequence of this arrangement
!@+     is that the checkpointing of model state will write a CACHED_SUBDD
!@+     array having the name 'some_var' if the above example were present
!@+     in the code, so a name collision will occur if 'some_var' is already
!@+     the checkpoint-file name of a model array. Please choose your output
!@+     names keeping this in mind (i.e. 't', 'q', 'p', 'u', 'v' will not
!@+     work, sorry) until I resolve the issue.  This is not an issue for
!@+     type I outputs.
!@+     If snapshots are to be written, inc_subdd can still be called
!@+     every timestep; it will simply return immediately during the
!@+     intervening timesteps.  In case you wish to avoid calculating
!@+     some expensive diagnostic every timestep, keep in mind that
!@+     snapshots are taken when mod(itime+1,Nsubdd)==0.
!@+
!@+   Terminology:
!@+   Categories occupy the middle position of a 3-level taxonomy
!@+   based on string identifiers (names).
!@+     (a) shape: a name for the grid on which a field is
!@+         defined.  Currently recognized names:
!@+            i. aijh  horizontal grid only
!@+           ii. aijph horizontal grid + constant-pressure levels
!@+          iii. aijlh horizontal grid + model levels
!@+           iv. gijlh horizontal grid + soil levels
!@+     (b) category: the name given to a collection of outputs
!@+         of a particular shape, e.g. those from a single model
!@+         component.  The identifier must be unique globally, i.e.
!@+         two categories having different shapes must have different
!@+         category names.
!@+     (c) group: an instance of a given category for
!@+         a particular output frequency and sampling type
!@+         (e.g. instantaneous versus time-averaged).  A group
!@+         may contain multiple output fields.  The names of
!@+         groups are generated by combining the names of
!@+         their parent categories with their sampling information.
!@+   The hierarchy is stored in flattened form, i.e. as a list of groups
!@+   having shape and category attributes.  Each group is placed into
!@+   a separate output file by the scaleacc procedure.
!@+
!@+   Other notes:
!@+   -  Arrays passed to inc_subdd must be horizontally dimensioned to
!@+      include the 1-row "halo".
!@+   -  To define a time-averaged type-I field as a ratio of sums,
!@+      set the dname element of the info_type structure to the
!@+      sname of the field serving as the denominator.  For
!@+      example: if there is a cloud fraction output named 'cldf',
!@+      and one wishes to weight a water content diagnostic by
!@+      cloud fraction, one can set dname='cldf' along with the
!@+      other metadata for water content.  The call to inc_subdd()
!@+      for water content would then pass an array equal to in-cloud
!@+      water content multiplied by cloud fraction.  It is the
!@+      responsibility of the user to request that 'cldf' be
!@+      output, and at the same time frequency as water content;
!@+      execution will stop if dname is set but the corresponding
!@+      field does not exist in the numerator''s group.

      module subdd_mod
      use cdl_mod, only : cdl_type
      use mdiag_com, only : sname_strlen,units_strlen,lname_strlen
      use resolution, only : lm
      implicit none

!@dbparam write_daily_files flag indicating whether to output data
!@+   once per day (the default is once per month)
      logical :: write_one_file = .false.
      logical :: write_daily_files = .false.
      logical :: write_monthly_files ! = .not. write_daily_files
!@dbparam days_per_file output files are written every days_per_file days
!@+   if set in the rundeck or if write_daily_files is true, in which
!@+   case days_per_file=1.  Otherwise files are written once per month.
      integer :: days_per_file

!@var vinterp_using_timeavgs flag indicating whether to use time-averaged
!@+       pressure to calculate vertical interpolation weights for
!@+       time-averaged constant-pressure outputs.  Default: no.
      logical :: vinterp_using_timeavgs

      integer, parameter :: subdd_nsched=3,
     &     sched_src=1,sched_rad=2,sched_inst=3

!@param reduc_avg,reduc_min,reduc_max indices for declaring the type
!@+   of reduction operation for time-reduced outputs.  The default
!@+   is to average.  At this time, such declarations are made in
!@+   source code along with the other metadata specifications for
!@+   a particular output, but the parsing of requests in the rundeck
!@+   will be changed to allow runtime choices.
      integer, parameter :: reduc_avg=1,reduc_min=2,reduc_max=3

!@dbparam subdd_npres number of subdd constant-pressure levels
      integer :: subdd_npres=0
!@dbparam subdd_pres subdd constant-pressure levels
!@var subdd_pk subdd_pres to the power kapa
      real*8, allocatable :: subdd_pres(:),subdd_pk(:)
!@var AIJPh_l1,AIJPh_l2 vertical interpolation coefficients for AIJPh
      real*8, dimension(:,:,:), allocatable :: AIJPh_l1,AIJPh_l2
!@dbparam LmaxSUBDD: the max L when writing "ALL" levels
      INTEGER :: LmaxSUBDD = LM

      integer, parameter ::
     &     num_subdd_str = 18, subdd_strlen = 64,
     &     subddt_len = num_subdd_str*subdd_strlen

      character(len=subddt_len) :: subddt

!@type subdd_type a derived type holding control parameters,
!@+    output metadata, and the timeseries for a given group
      type subdd_type
!@var initialized whether a given instance of this type has been initialized
      logical :: initialized=.false.
!@dbparam Nsubdd: DT_save_SUBDD =  Nsubdd*DTsrc sub-daily diag freq.
      INTEGER :: nsubdd = 0
!@var nperiod number of subdd accumulation periods per output file
      integer :: nperiod=0
!@var subdd_period index of the current subdd accumulation period
      integer :: subdd_period
!@var ndiags number of outputs in this group
      integer :: ndiags=0
!@var dtime length of each output period (hr)
!@var timelast time of last accumulation/store in this group and output file
!@+            (hrs since start of run)
      real*8 :: dtime=0.,timelast=0.
!@var is_inst whether this diag group holds snapshots
      logical :: is_inst=.false.
!@var accum_this_step always true for time-reduction output groups,
!@+   and true every nsubdd timesteps for snapshot groups
      logical :: accum_this_step=.false.
!@var nacc number of accumulations for each output period and schedule
      integer, allocatable :: nacc(:,:) ! (nperiod,subdd_nsched)
!@var catshape unique name for the shape (type) of the category to
!@+   which this group belongs
!@var catshape_accum normally the same as catshape, except for
!@+   constant-pressure outputs with deferred vertical regridding
!@+   from model levels
      character(len=sname_strlen) :: catshape,catshape_accum
!@var dsize3_input for 3D output, size of third dimension of arrays
!@+   passed to inc_subdd
      integer :: dsize3_input
!@var category unique name for the category to which this group belongs
      character(len=sname_strlen) :: category
!@var grpname name of this output group
      character(len=sname_strlen) :: grpname
!@var strdimlen, accshape, rsfshape strings used to define acc-arrays in
!@+   output files
      character(len=sname_strlen) :: strdimlen
      character(len=lname_strlen) :: accshape
      character(len=lname_strlen) :: rsfshape=''
!@var split_dim the dimension in v4d or v5d corresponding to output quantity
!@+   (needed for postprocessing)
      integer :: split_dim
!@var tile_dim_out needs to be to the dim after j for the cubed-sphere case
      integer :: tile_dim_out=0
!@var next_soloid a saved ID for optimizing string lookups
      integer :: next_soloid=0
!@var denom index of the denominator for outputs that are ratios
!@var sched the index of the accumulation schedule for each output
!@var reduc indicates the type of time reduction for each output
!@+   (average, min, max).  Has no effect if is_inst=.true.
      integer, allocatable :: denom(:),sched(:),reduc(:)
!@var scale scale factor for each output
      real*8, allocatable :: scale(:)
!@var name short name for each output
      character(len=sname_strlen), allocatable :: name(:)
!@var cdl0 consolidated metadata for this group without time data
!@var cdl  like cdl0 but with time info
      type(cdl_type) :: cdl0,cdl
!@var v4d accumulation array for 2D subdaily diagnostics
      real*8, allocatable :: v4d(:,:,:,:)   !(:,:,nperiod,ndiags)
!@var v5d accumulation array for 3D subdaily diagnostics
      real*8, allocatable :: v5d(:,:,:,:,:) !(:,:,:,nperiod,ndiags)
!@var v5dout is written to the final output file instead of v5d when
!@+   grid transformations are applied to v5d (e.g. when outputs
!@+   are stored on model levels and later remapped to constant-pressure
!@+   levels using time-averaged pressure).
      real*8, allocatable :: v5dout(:,:,:,:,:) !(:,:,:,nperiod,ndiags)
      end type subdd_type

!@var cdl_ijt template for consolidated metadata for files
!@+   containing a time dimension.  Once the monthly-mean
!@+   files are given a time dimension, this template will
!@+   be obsolete.
      type(cdl_type) :: cdl_ijt

!@param subdd_ngroups_max maximum number of output groups per run
!@+     (increase as necessary)
#ifdef COSP_SIM
      integer, parameter :: subdd_ngroups_max=40
#else
      integer, parameter :: subdd_ngroups_max=30
#endif
c SUSA only for MEEEEE
c       integer, parameter :: subdd_ngroups_max=50
c SUSA
!@var subdd_ngroups the number of output groups for this run
      integer :: subdd_ngroups=0
!@var subdd_solovar1 the index starting from which groups are solo variables
!@var lastsolo the last solo variable that was incremented
      integer :: subdd_solovar1=0,lastsolo=0
!@var subdd_groups a handle for the collection of output groups
      !type(subdd_type), target :: subdd_groups(subdd_ngroups_max)
      type(subdd_type), allocatable, target :: subdd_groups(:)

!@var rsf_save saves the checkpoint filename so that deferred-registration
!@+   groups can be initialized after the model has finished restarting
      character(len=80) :: rsf_save='NOTAFILE'

!@param namedd_strlen the maximum length of names of subdaily
!@+     output fields
      integer, parameter :: namedd_strlen=sname_strlen

!@itimei_subdd timestep counter at beginning of execution
!@+           (only matters for write_one_file case)
      integer :: itimei_subdd

!@type info_type a derived type for metadata about output fields
      type info_type
      character(len=sname_strlen) :: sname='no output'
      character(len=units_strlen) :: units='no output'
      character(len=lname_strlen) :: lname='no output'
      character(len=sname_strlen) :: dname='no denom'
      real*8 :: scale = 1d0
      integer :: sched=sched_src
      integer :: reduc=reduc_avg
      end type info_type

      interface inc_subdd
        module procedure inc_subdd_2d
        module procedure inc_subdd_3d
        module procedure inc_subdd_solo_2d
        module procedure inc_subdd_solo_3d
        module procedure inc_subdd_solo_4d
      end interface inc_subdd

      contains

      function info_type_(sname,units,lname,dname,scale,sched,reduc)
c this is a homebrew info_type constructor for older compilers that
c require all elements of a derived type to be specified when invoking
c a structure constructor
      type(info_type) :: info_type_
      character(len=*), optional :: sname
      character(len=*), optional :: units
      character(len=*), optional :: lname
      character(len=*), optional :: dname
      real*8, optional :: scale
      integer, optional :: sched
      integer, optional :: reduc
c
      character(len=sname_strlen) :: sname_
      character(len=units_strlen) :: units_
      character(len=lname_strlen) :: lname_
      character(len=sname_strlen) :: dname_
      real*8 :: scale_
      integer :: sched_
      integer :: reduc_
c
      sname_='no output'
      units_='no output'
      lname_='no output'
      dname_='no denom'
      scale_ = 1d0
      sched_ = sched_src
      reduc_ = reduc_avg
c
      if(present(sname)) sname_ = sname
      if(present(units)) units_ = units
      if(present(lname)) lname_ = lname
      if(present(dname)) dname_ = dname
      if(present(scale)) scale_ = scale
      if(present(sched)) sched_ = sched
      if(present(reduc)) reduc_ = reduc

      if(reduc_ .ne. reduc_avg) sched_ = sched_inst
c
      info_type_%sname = sname_
      info_type_%units = units_
      info_type_%lname = lname_
      info_type_%dname = dname_
      info_type_%scale = scale_
      info_type_%sched = sched_
      info_type_%reduc = reduc_
      end function info_type_

      subroutine find_groups(str,grpids,ngroups)
      implicit none
      character(len=*) :: str
      integer, dimension(subdd_ngroups) :: grpids
      integer :: ngroups
c
      integer :: idgrp
c
      ngroups = 0
      do idgrp=1,subdd_ngroups
        if(.not.subdd_groups(idgrp)%accum_this_step) cycle
        if(trim(subdd_groups(idgrp)%category).ne.trim(str)) cycle
        ngroups = ngroups + 1
        grpids(ngroups) = idgrp
      enddo
      return
      end subroutine find_groups

      subroutine inc_subdd_2d(subdd,k,arr)
      use domain_decomp_atm, only : grid
      use geom, only : imaxj
      implicit none
      type(subdd_type) :: subdd
      integer :: k
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) :: arr
c
      integer :: i,j,ip
c
      ip = subdd%subdd_period
      if(subdd%reduc(k).eq.reduc_avg) then
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,imaxj(j)
          subdd%v4d(i,j,ip,k) = subdd%v4d(i,j,ip,k) + arr(i,j)
        enddo
        enddo
      elseif(subdd%reduc(k).eq.reduc_min) then
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,imaxj(j)
          subdd%v4d(i,j,ip,k) = min(subdd%v4d(i,j,ip,k),arr(i,j))
        enddo
        enddo
      elseif(subdd%reduc(k).eq.reduc_max) then
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,imaxj(j)
          subdd%v4d(i,j,ip,k) = max(subdd%v4d(i,j,ip,k),arr(i,j))
        enddo
        enddo
      endif
      return
      end subroutine inc_subdd_2d

      subroutine inc_subdd_3d(subdd,k,arr,jdim)
      use domain_decomp_atm, only : grid
      implicit none
      type(subdd_type) :: subdd
      integer :: k
      real*8, dimension(:,:,:) :: arr
      integer, optional :: jdim
c
      integer :: jdim_
      integer :: ni,nj,sizes(3),dim
c
      jdim_ = 2
      if(present(jdim)) jdim_ = jdim
      if(jdim_.eq.1) then
        call stop_model('inc_subdd_3d: jdim==1 is invalid',255)
      endif
      ni = 1+grid%i_stop_halo-grid%i_strt_halo
      nj = 1+grid%j_stop_halo-grid%j_strt_halo
      do dim=1,3
        sizes(dim) = size(arr,dim)
      enddo
      if(ni.ne.sizes(jdim_-1) .or. nj.ne.sizes(jdim_)) then
        call stop_model('inc_subdd_3d: horz dimsize mismatch',255)
      endif
      select case (trim(subdd%catshape_accum))
      case ('aijlh', 'gijlh')
        if(jdim_ .eq. 2) then
          if(subdd%dsize3_input.ne.sizes(3)) then
            call stop_model('inc_subdd_3d: xijlh dimsize mismatch',255)
          endif
          call inc_subdd_aijlh(subdd,k,arr)
        else
          if(subdd%dsize3_input.ne.sizes(1)) then
            call stop_model('inc_subdd_3d: xlijh dimsize mismatch',255)
          endif
          call inc_subdd_alijh(subdd,k,arr)
        endif
      case ('aijph')
        if(jdim_ .eq. 2) then
          if(subdd%dsize3_input.eq.sizes(3)) then
            call inc_subdd_cp_ijl(subdd,k,arr)
          elseif(subdd_npres.eq.sizes(3)) then
            call inc_subdd_cp(subdd,k,arr)
          else
            call stop_model('inc_subdd_3d: unrecognized size3',255)
          endif
        else
          if(subdd%dsize3_input.eq.sizes(1)) then
            call inc_subdd_cp_lij(subdd,k,arr)
          elseif(subdd_npres.eq.sizes(1)) then
            !
          else
            call stop_model('inc_subdd_3d: unrecognized size1',255)
          endif
        endif
      case default
        call stop_model('inc_subdd_3d: unrecognized catshape',255)
      end select
      return
      end subroutine inc_subdd_3d

      subroutine inc_subdd_aijlh(subdd,k,arr)
      use domain_decomp_atm, only : grid
      use geom, only : imaxj
      implicit none
      type(subdd_type) :: subdd
      integer :: k
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,
     &                  subdd%dsize3_input) ::
     &     arr
c
      integer :: i,j,l,ip
c
      ip = subdd%subdd_period
      do l=1,size(subdd%v5d,3)
      do j=grid%j_strt,grid%j_stop
      do i=grid%i_strt,imaxj(j)
        subdd%v5d(i,j,l,ip,k) = subdd%v5d(i,j,l,ip,k) + arr(i,j,l)
      enddo
      enddo
      enddo
      return
      end subroutine inc_subdd_aijlh

      subroutine inc_subdd_alijh(subdd,k,arr)
      use domain_decomp_atm, only : grid
      use geom, only : imaxj
      implicit none
      type(subdd_type) :: subdd
      integer :: k
      real*8, dimension(subdd%dsize3_input,
     &                  grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo) ::
     &     arr
c
      integer :: i,j,l,ip
c
      ip = subdd%subdd_period
      do l=1,size(subdd%v5d,3)
      do j=grid%j_strt,grid%j_stop
      do i=grid%i_strt,imaxj(j)
        subdd%v5d(i,j,l,ip,k) = subdd%v5d(i,j,l,ip,k) + arr(l,i,j)
      enddo
      enddo
      enddo
      return
      end subroutine inc_subdd_alijh

      subroutine inc_subdd_cp(subdd,k,arr)
      use domain_decomp_atm, only : grid
      use geom, only : imaxj
      implicit none
      type(subdd_type) :: subdd
      integer :: k
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,subdd_npres)
     &     :: arr
c
      integer :: i,j,l,ip
c
      ip = subdd%subdd_period

      do l=1,subdd_npres
      do j=grid%j_strt,grid%j_stop
      do i=grid%i_strt,imaxj(j)
        subdd%v5d(i,j,l,ip,k) = subdd%v5d(i,j,l,ip,k) + arr(i,j,l)
      enddo
      enddo
      enddo
      return
      end subroutine inc_subdd_cp

      subroutine inc_subdd_cp_ijl(subdd,k,arr)
      use domain_decomp_atm, only : grid
      use geom, only : imaxj
      use resolution, only : lm
      implicit none
      type(subdd_type) :: subdd
      integer :: k
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,lm) ::
     &     arr
c
      integer :: i,j,l,ip,ldn,lup
      real*8 :: wtdn,wtup
c
      ip = subdd%subdd_period

      do l=1,subdd_npres
      do j=grid%j_strt,grid%j_stop
      do i=grid%i_strt,imaxj(j)
        ldn = aijph_l1(i,j,l)
        lup = aijph_l2(i,j,l)
        wtdn = aijph_l1(i,j,l)-ldn
        wtup = aijph_l2(i,j,l)-lup
        subdd%v5d(i,j,l,ip,k) = subdd%v5d(i,j,l,ip,k) +
     &         (wtdn*arr(i,j,ldn) +wtup*arr(i,j,lup))
      enddo
      enddo
      enddo
      return
      end subroutine inc_subdd_cp_ijl

      subroutine inc_subdd_cp_lij(subdd,k,arr)
      use domain_decomp_atm, only : grid
      use geom, only : imaxj
      use resolution, only : lm
      implicit none
      type(subdd_type) :: subdd
      integer :: k
      real*8, dimension(lm,grid%i_strt_halo:grid%i_stop_halo,
     &                     grid%j_strt_halo:grid%j_stop_halo) ::
     &     arr
c
      integer :: i,j,l,ip,ldn,lup
      real*8 :: wtdn,wtup
c
      ip = subdd%subdd_period

      do l=1,subdd_npres
      do j=grid%j_strt,grid%j_stop
      do i=grid%i_strt,imaxj(j)
        ldn = aijph_l1(i,j,l)
        lup = aijph_l2(i,j,l)
        wtdn = aijph_l1(i,j,l)-ldn
        wtup = aijph_l2(i,j,l)-lup
        subdd%v5d(i,j,l,ip,k) = subdd%v5d(i,j,l,ip,k) +
     &         (wtdn*arr(ldn,i,j) +wtup*arr(lup,i,j))
      enddo
      enddo
      enddo
      return
      end subroutine inc_subdd_cp_lij

      subroutine inc_subdd_solo_3d(vname,arr,nsubdd,inst,
     &     jdim,units,long_name,dim3name,coordvalues)
      use domain_decomp_atm, only : grid,get=>getdomainbounds
      use geom, only : imaxj
      use cdl_mod, only : add_var,add_varline,add_dim,add_coord
      implicit none
      character(len=*) :: vname
      real*8, dimension(:,:,:) :: arr
      integer :: nsubdd
      logical :: inst
      integer, optional :: jdim
      character(len=*), optional :: units,long_name,dim3name
      real*8, optional :: coordvalues(:)
c
      type(subdd_type), pointer :: subdd
      integer :: jdim_
      integer :: ni,nj,sizes(3),dim,ip,idgrp,dsize
      integer :: i,j,l
      integer :: i_0h,j_0h,i_1h,j_1h
      character(len=64) :: dimstr,dname
      character(len=8) :: cnperiod,cdsize
c
      call find_solovar(vname,idgrp)

      subdd => subdd_groups(idgrp)

      call get(grid,i_strt_halo=i_0h,i_stop_halo=i_1h,
     &              j_strt_halo=j_0h,j_stop_halo=j_1h)

      jdim_ = 2
      if(present(jdim)) jdim_ = jdim
      if(jdim_.eq.1) then
        call stop_model('inc_subdd_solo: jdim==1 is invalid',255)
      endif

      if(.not.subdd%initialized) then

        call create_solo(subdd,nsubdd,vname,inst)
        cnperiod=''
        write(cnperiod,'(i8)') subdd%nperiod

        if(jdim_.eq.2) then
          dsize = size(arr,3)
        else
          dsize = size(arr,1)
        endif
        cdsize=''
        write(cdsize,'(i8)') dsize
        dname = 'size_'//trim(adjustl(cdsize))
        subdd%accshape =
     &       'dist_im,dist_jm,'//trim(dname)//
     &       ',nperiod_'//trim(adjustl(cnperiod))
        subdd%split_dim = 5
        allocate(subdd%v5d(I_0H:I_1H,J_0H:J_1H,dsize,subdd%nperiod,1))
        subdd%v5d = 0.

        ! read from rsf if necessary
        call read_subdd_rsf1(subdd)

          ! define metadata for this output

        if(present(dim3name)) then
          dname = dim3name
        !else
        !  dname = 'size_'//trim(vname)
        endif

        subdd%cdl0 = cdl_ijt

#ifdef CUBED_SPHERE
        subdd%tile_dim_out = 4
        dimstr='(time,tile,'//trim(dname)//',y,x) ;'
#else
        dimstr='(time,'//trim(dname)//',lat,lon) ;'
#endif


        call add_var(subdd%cdl0,'float '//trim(vname)//trim(dimstr))

        if(subdd%is_inst) then
        ! global attribute
          call add_varline(subdd%cdl0,':sampling = "instantaneous" ;')
        endif

        if(present(units)) then
          call add_varline(subdd%cdl0,trim(vname)//
     &         ':units = "'//trim(units)//'" ;')
        endif
        if(present(long_name)) then
          call add_varline(subdd%cdl0,trim(vname)//
     &         ':long_name = "'//trim(long_name)//'" ;')
        endif
        if(present(coordvalues)) then
          if(size(coordvalues).lt.dsize) then
            call stop_model('inc_subdd_solo: coordvalues too small',255)
          endif
          call add_coord(subdd%cdl0,trim(dname),dsize,
     &         coordvalues=coordvalues(1:dsize))
        else
          call add_dim(subdd%cdl0,trim(dname),dsize)
        endif

        call set_subdd_period1(subdd)

        subdd%initialized = .true.
      endif ! initialized?
      subdd_groups(lastsolo)%next_soloid = idgrp
      lastsolo = idgrp
c
      if(.not.subdd%accum_this_step) return
c
      ni = 1+i_1h-i_0h
      nj = 1+j_1h-j_0h
      do dim=1,3
        sizes(dim) = size(arr,dim)
      enddo
      if(ni.ne.sizes(jdim_-1) .or. nj.ne.sizes(jdim_)) then
        call stop_model('inc_subdd_solo: horz dimsize mismatch',255)
      endif
      ip = subdd%subdd_period

      if(jdim_ .eq. 2) then
        if(size(subdd%v5d,3).ne.sizes(3)) then
          call stop_model('inc_subdd_solo: aijlh dimsize mismatch',255)
        endif
        do l=1,sizes(3)
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,imaxj(j)
          subdd%v5d(i,j,l,ip,1) = subdd%v5d(i,j,l,ip,1) +
     &         arr(1+i-i_0h,1+j-j_0h,l)
        enddo
        enddo
        enddo
      else
        if(size(subdd%v5d,3).ne.sizes(1)) then
          call stop_model('inc_subdd_solo: alijh dimsize mismatch',255)
        endif
        do l=1,sizes(1)
        do j=grid%j_strt,grid%j_stop
        do i=grid%i_strt,imaxj(j)
          subdd%v5d(i,j,l,ip,1) = subdd%v5d(i,j,l,ip,1) +
     &         arr(l,1+i-i_0h,1+j-j_0h)
        enddo
        enddo
        enddo
      endif

      return
      end subroutine inc_subdd_solo_3d

      subroutine inc_subdd_solo_4d(vname,arr,nsubdd,inst,
     &     jdim,units,long_name,
     &     dim3name,dim4name,
     &     coordvalues3,coordvalues4)
      use domain_decomp_atm, only : grid,get=>getdomainbounds
      use geom, only : imaxj
      use cdl_mod, only : add_var,add_varline,add_dim,add_coord
      implicit none
      character(len=*) :: vname
      real*8, dimension(:,:,:,:) :: arr
      integer :: nsubdd
      logical :: inst
      integer, optional :: jdim
      character(len=*), optional :: units,long_name,dim3name,dim4name
      real*8, optional :: coordvalues3(:),coordvalues4(:)
c
      type(subdd_type), pointer :: subdd
      integer :: jdim_
      integer :: ni,nj,sizes(4),dim,ip,idgrp,dsize
      integer :: i,j,k,l,m,kdim,ldim
      integer :: i_0h,j_0h,i_1h,j_1h
      character(len=64) :: dimstr,dname,dname3,dname4,dname43
      character(len=8) :: cnperiod,cdsize
c
      call find_solovar(vname,idgrp)

      subdd => subdd_groups(idgrp)

      call get(grid,i_strt_halo=i_0h,i_stop_halo=i_1h,
     &              j_strt_halo=j_0h,j_stop_halo=j_1h)

      jdim_ = 2
      if(present(jdim)) jdim_ = jdim
      if(jdim_.eq.1) then
        call stop_model('inc_subdd_solo: jdim==1 is invalid',255)
      endif

      do dim=1,4
        sizes(dim) = size(arr,dim)
      enddo

      if(jdim_ .eq. 2) then
        kdim = 3
        ldim = 4
      elseif(jdim_.eq.3) then
        kdim = 1
        ldim = 4
      else
        kdim = 1
        ldim = 2
      endif

      if(.not.subdd%initialized) then

        call create_solo(subdd,nsubdd,vname,inst)
        cnperiod=''
        write(cnperiod,'(i8)') subdd%nperiod

        dsize = sizes(kdim)*sizes(ldim)
        cdsize=''
        write(cdsize,'(i8)') dsize
        dname = 'size_'//trim(adjustl(cdsize))
        subdd%accshape =
     &       'dist_im,dist_jm,'//trim(dname)//
     &       ',nperiod_'//trim(adjustl(cnperiod))
        subdd%split_dim = 5
        allocate(subdd%v5d(I_0H:I_1H,J_0H:J_1H,dsize,subdd%nperiod,1))
        subdd%v5d = 0.

        ! read from rsf if necessary
        call read_subdd_rsf1(subdd)

          ! define metadata for this output

        if(present(dim3name)) then
          dname3 = dim3name
        else
          cdsize=''
          write(cdsize,'(i8)') sizes(kdim)
          dname3 = 'size_'//trim(adjustl(cdsize))
        endif
        if(present(dim4name)) then
          dname4 = dim4name
        else
          cdsize=''
          write(cdsize,'(i8)') sizes(ldim)
          dname4 = 'size_'//trim(adjustl(cdsize))
        endif
        dname43=trim(dname4)//','//trim(dname3)

        subdd%cdl0 = cdl_ijt

#ifdef CUBED_SPHERE
        subdd%tile_dim_out = 4 ! dims 3,4 are fused in acc array
        dimstr='(time,tile,'//trim(dname43)//',y,x) ;'
#else
        dimstr='(time,'//trim(dname43)//',lat,lon) ;'
#endif


        call add_var(subdd%cdl0,'float '//trim(vname)//trim(dimstr))

        if(subdd%is_inst) then
        ! global attribute
          call add_varline(subdd%cdl0,':sampling = "instantaneous" ;')
        endif

        if(present(units)) then
          call add_varline(subdd%cdl0,trim(vname)//
     &         ':units = "'//trim(units)//'" ;')
        endif
        if(present(long_name)) then
          call add_varline(subdd%cdl0,trim(vname)//
     &         ':long_name = "'//trim(long_name)//'" ;')
        endif

        if(present(coordvalues3)) then
          if(size(coordvalues3).lt.sizes(kdim)) then
           call stop_model('inc_subdd_solo: coordvalues3 too small',255)
          endif
          call add_coord(subdd%cdl0,trim(dname3),sizes(kdim),
     &         coordvalues=coordvalues3(1:sizes(kdim)))
        else
          call add_dim(subdd%cdl0,trim(dname3),sizes(kdim))
        endif
        if(present(coordvalues4)) then
          if(size(coordvalues4).lt.sizes(ldim)) then
           call stop_model('inc_subdd_solo: coordvalues4 too small',255)
          endif
          call add_coord(subdd%cdl0,trim(dname4),sizes(ldim),
     &         coordvalues=coordvalues4(1:sizes(ldim)))
        else
          call add_dim(subdd%cdl0,trim(dname4),sizes(ldim))
        endif

        call set_subdd_period1(subdd)

        subdd%initialized = .true.
      endif ! initialized?
      subdd_groups(lastsolo)%next_soloid = idgrp
      lastsolo = idgrp
c
      if(.not.subdd%accum_this_step) return
c
      ni = 1+i_1h-i_0h
      nj = 1+j_1h-j_0h
      if(ni.ne.sizes(jdim_-1) .or. nj.ne.sizes(jdim_)) then
        call stop_model('inc_subdd_solo: horz dimsize mismatch',255)
      endif
      ip = subdd%subdd_period

      if(size(subdd%v5d,3) .ne.sizes(kdim)*sizes(ldim)) then
        call stop_model('inc_subdd_solo_4d: dimsize mismatch',255)
      endif

      m = 0
      do l=1,sizes(ldim)
      do k=1,sizes(kdim)
        m = m + 1
        if(jdim_ .eq. 2) then
          do j=grid%j_strt,grid%j_stop
          do i=grid%i_strt,imaxj(j)
            subdd%v5d(i,j,m,ip,1) = subdd%v5d(i,j,m,ip,1) +
     &           arr(1+i-i_0h,1+j-j_0h,k,l)
          enddo
          enddo
        elseif(jdim_ .eq. 3) then
          do j=grid%j_strt,grid%j_stop
          do i=grid%i_strt,imaxj(j)
            subdd%v5d(i,j,m,ip,1) = subdd%v5d(i,j,m,ip,1) +
     &           arr(k,1+i-i_0h,1+j-j_0h,l)
          enddo
          enddo
        else
          do j=grid%j_strt,grid%j_stop
          do i=grid%i_strt,imaxj(j)
            subdd%v5d(i,j,m,ip,1) = subdd%v5d(i,j,m,ip,1) +
     &           arr(k,l,1+i-i_0h,1+j-j_0h)
          enddo
          enddo
        endif
      enddo
      enddo

      return
      end subroutine inc_subdd_solo_4d

      subroutine inc_subdd_solo_2d(vname,arr,nsubdd,inst,
     &     units,long_name)
      use domain_decomp_atm, only : grid,get=>getdomainbounds
      use geom, only : imaxj
      use cdl_mod, only : add_var,add_varline
      implicit none
      character(len=*) :: vname
      real*8, dimension(:,:) :: arr
      integer :: nsubdd
      logical :: inst
      character(len=*), optional :: units,long_name
c
      type(subdd_type), pointer :: subdd
      integer :: ni,nj,ip,idgrp
      integer :: i,j
      integer :: i_0h,j_0h,i_1h,j_1h
      character(len=32) :: dimstr
      character(len=8) :: cnperiod
c
      call find_solovar(vname,idgrp)

      subdd => subdd_groups(idgrp)

      call get(grid,i_strt_halo=i_0h,i_stop_halo=i_1h,
     &              j_strt_halo=j_0h,j_stop_halo=j_1h)

      if(.not.subdd%initialized) then

        call create_solo(subdd,nsubdd,vname,inst)
        cnperiod=''
        write(cnperiod,'(i8)') subdd%nperiod

        subdd%accshape =
     &       'dist_im,dist_jm,nperiod_'//trim(adjustl(cnperiod))
        subdd%split_dim = 4
        allocate(subdd%v4d(I_0H:I_1H,J_0H:J_1H,subdd%nperiod,1))
        subdd%v4d = 0.

        ! read from rsf if necessary
        call read_subdd_rsf1(subdd)

          ! define metadata for this output
        subdd%cdl0 = cdl_ijt
#ifdef CUBED_SPHERE
        subdd%tile_dim_out = 3
        dimstr='(time,tile,y,x) ;'
#else
        dimstr='(time,lat,lon) ;'
#endif
        call add_var(subdd%cdl0,'float '//trim(vname)//trim(dimstr))
        if(subdd%is_inst) then
        ! global attribute
          call add_varline(subdd%cdl0,':sampling = "instantaneous" ;')
        endif
        if(present(units)) then
          call add_varline(subdd%cdl0,trim(vname)//
     &         ':units = "'//trim(units)//'" ;')
        endif
        if(present(long_name)) then
          call add_varline(subdd%cdl0,trim(vname)//
     &         ':long_name = "'//trim(long_name)//'" ;')
        endif

        call set_subdd_period1(subdd)
        subdd%initialized = .true.
      endif ! initialized?
      subdd_groups(lastsolo)%next_soloid = idgrp
      lastsolo = idgrp
c
      if(.not.subdd%accum_this_step) return
c
      ni = 1+i_1h-i_0h
      nj = 1+j_1h-j_0h
      if(ni.ne.size(arr,1) .or. nj.ne.size(arr,2)) then
        call stop_model('inc_subdd_solo: horz dimsize mismatch',255)
      endif
      ip = subdd%subdd_period
      do j=grid%j_strt,grid%j_stop
      do i=grid%i_strt,imaxj(j)
        subdd%v4d(i,j,ip,1) = subdd%v4d(i,j,ip,1) +
     &       arr(1+i-i_0h,1+j-j_0h)
      enddo
      enddo

      return
      end subroutine inc_subdd_solo_2d

      subroutine find_solovar(varname,ind)
      implicit none
      character(len=*) :: varname
      integer :: ind
      integer :: idgrp
      ind = subdd_ngroups+1
      if(subdd_solovar1==0) then
        subdd_solovar1 = ind
        lastsolo = subdd_solovar1
        subdd_groups(ind)%next_soloid = ind ! just an initial guess
      endif
      do idgrp=subdd_groups(lastsolo)%next_soloid,subdd_ngroups
        if(trim(subdd_groups(idgrp)%grpname).eq.trim(varname)) then
          ind = idgrp
          return
        endif
      enddo
      do idgrp=subdd_solovar1,subdd_ngroups
        if(trim(subdd_groups(idgrp)%grpname).eq.trim(varname)) then
          ind = idgrp
          return
        endif
      enddo
      if(ind.gt.subdd_ngroups_max) then
        call stop_model('too many groups in CACHED_SUBDD',255)
      endif
      subdd_ngroups = ind
      subdd_groups(ind)%next_soloid = ind ! just an initial guess
      return
      end subroutine find_solovar

      subroutine create_group(
     &     subdd,nsubdd,
     &     catshape,category,grpname,
     &     namedd,kdd,is_inst,
     &     diaglist,listlen,ndiags,dsize3_input)
      use model_com, only : nday,dtsrc
      use model_com, only : itimei,itimee
      use resolution, only : lm ! temporary?
      use cdl_mod, only : add_var,add_varline,add_coord
      use domain_decomp_atm, only : grid,get=>getdomainbounds
      implicit none
      type(subdd_type) :: subdd
      integer :: nsubdd
      character(len=sname_strlen) :: catshape,category,grpname
      character(len=namedd_strlen) :: namedd(kdd)
      integer :: kdd,listlen,ndiags,dsize3_input
      logical :: is_inst
      type(info_type) :: diaglist(listlen)
c
      integer :: i,j,k,kk,l,nperiod,dsize3,dsize3out
      integer, parameter :: kmax=1024
      ! overdimensioned temporary instances of metadata to collect declarations
      real*8, dimension(kmax) :: scale_tmp
      integer, dimension(kmax) :: sched_tmp,reduc_tmp,denom_tmp
      character(len=sname_strlen), dimension(kmax) :: name_tmp
      character(len=sname_strlen), dimension(kmax) :: dname_tmp
      character(len=units_strlen), dimension(kmax) :: units_tmp
      character(len=lname_strlen), dimension(kmax) :: lname_tmp
      character(len=32) :: dimstr
      character(len=256) :: errstr
      real*8, allocatable :: lvlarr(:)
      integer :: i_0h,i_1h,j_0h,j_1h

      k = 0
      do i=1,kdd
      do j=1,listlen
        if(trim(diaglist(j)%sname).eq.trim(namedd(i))) then
          k = k + 1
           name_tmp(k) = diaglist(j)%sname
          lname_tmp(k) = diaglist(j)%lname
          dname_tmp(k) = diaglist(j)%dname
          units_tmp(k) = diaglist(j)%units
          scale_tmp(k) = diaglist(j)%scale
          sched_tmp(k) = diaglist(j)%sched
          if(is_inst) sched_tmp(k) = sched_inst
          reduc_tmp(k) = diaglist(j)%reduc
          exit
        endif
      enddo
      enddo
      ndiags = k

      if(ndiags.eq.0) return

      ! find the indices of denominators from their names
      denom_tmp(:) = 0
      do k=1,ndiags
        if(trim(dname_tmp(k)).eq.'no denom') cycle
        do kk=1,ndiags
          if(trim(dname_tmp(k)).eq.trim(name_tmp(kk))) then
            denom_tmp(k) = kk
            exit
          endif
        enddo
        if(denom_tmp(k).eq.0) then
          errstr='field '//trim(dname_tmp(k))//' needed as '//
     &         'the denominator for subdd output '//trim(name_tmp(k))
          call stop_model(trim(errstr),255)
        endif
      enddo

      if(write_one_file) then
        !nperiod = ceiling(real(itimee-itimei)/real(nsubdd))
        nperiod = ceiling(real(itimee-itimei_subdd)/real(nsubdd))
      else
        nperiod = ceiling(real(days_per_file*nday)/real(nsubdd))
      endif

      call get(grid,i_strt_halo=i_0h,i_stop_halo=i_1h,
     &              j_strt_halo=j_0h,j_stop_halo=j_1h)

      subdd%cdl0 = cdl_ijt

      subdd%grpname = grpname
      subdd%strdimlen = 'k'//trim(grpname)

      subdd%nsubdd = nsubdd
      subdd%dtime = DTsrc*subdd%nsubdd/3600d0
      subdd%nperiod = nperiod
      allocate(subdd%nacc(nperiod,subdd_nsched))
      subdd%is_inst = is_inst
      subdd%catshape = catshape
      subdd%catshape_accum = catshape
      subdd%category = category

      dsize3out = 0

      if(is_inst) then
        ! global attribute
        call add_varline(subdd%cdl0,':sampling = "instantaneous" ;')
      endif

      select case (catshape)
      case ('aijh')
        dsize3 = 0
        subdd%accshape = 'dist_im,dist_jm,nperiod_'//trim(grpname)
#ifdef CUBED_SPHERE
        subdd%tile_dim_out = 3
        dimstr='(time,tile,y,x) ;'
#else
        dimstr='(time,lat,lon) ;'
#endif

      case('aijlh','gijlh')
        if(catshape.eq.'aijlh') then
          dsize3 = lmaxsubdd
          subdd%accshape =
     &       'dist_im,dist_jm,lmaxsubdd,nperiod_'//trim(grpname)
        else
          dsize3 = dsize3_input
          subdd%accshape =
     &         'dist_im,dist_jm,lm_'//trim(grpname)//
     &         ',nperiod_'//trim(grpname)
        endif
#ifdef CUBED_SPHERE
        subdd%tile_dim_out = 4
        dimstr='(time,tile,level,y,x) ;'
#else
        dimstr='(time,level,lat,lon) ;'
#endif
        allocate(lvlarr(dsize3))
        do l=1,dsize3
          lvlarr(l) = l
        enddo
        call add_coord(subdd%cdl0,'level',size(lvlarr),
     &       coordvalues=lvlarr)
        deallocate(lvlarr)

      case('aijph')
        if( (.not.is_inst) .and. vinterp_using_timeavgs) then
          ! deferred vertical regridding of model-level accumulations
          dsize3 = lm
          subdd%rsfshape =
     &         'dist_im,dist_jm,lm,nperiod_'//trim(grpname)
          dsize3out = subdd_npres
          subdd%catshape_accum = 'aijlh'
        else
          dsize3 = subdd_npres
        endif
        subdd%accshape =
     &       'dist_im,dist_jm,subdd_npres,nperiod_'//trim(grpname)
#ifdef CUBED_SPHERE
        subdd%tile_dim_out = 4
        dimstr='(time,tile,p,y,x) ;'
#else
        dimstr='(time,p,lat,lon) ;'
#endif
        call add_coord(subdd%cdl0,'p',subdd_npres,units='mb',
     &         coordvalues=subdd_pres)
        ! extra position for a denominator containing the count of the
        ! number of times each pressure level was accumulated
        ndiags = ndiags + 1
        k = ndiags
        name_tmp(k) = 'no output'
        lname_tmp(k) = 'no output'
        units_tmp(k) = 'no output'
        scale_tmp(k) = 1.
        reduc_tmp(k) = reduc_avg
        sched_tmp(k) = sched_tmp(k-1)
        denom_tmp(1:ndiags-1) = ndiags
        denom_tmp(ndiags)=0

      case default
        call stop_model('unrecognized catshape in create_group',255)

      end select

      subdd%dsize3_input = dsize3_input
      subdd%ndiags = ndiags
      allocate(subdd%scale(ndiags))
      allocate(subdd%sched(ndiags))
      allocate(subdd%reduc(ndiags))
      allocate(subdd%name(ndiags))
      allocate(subdd%denom(ndiags))
      subdd%scale(:) = scale_tmp(1:ndiags)
      subdd%sched(:) = sched_tmp(1:ndiags)
      subdd%name(:) = name_tmp(1:ndiags)
      subdd%reduc(:) = reduc_tmp(1:ndiags)
      subdd%denom(:) = denom_tmp(1:ndiags)

      if(dsize3.le.0) then
        subdd%split_dim = 4
        allocate(
     &       subdd%v4d(I_0H:I_1H,J_0H:J_1H,nperiod,ndiags))
      else
        subdd%split_dim = 5
        allocate(
     &       subdd%v5d(I_0H:I_1H,J_0H:J_1H,dsize3,nperiod,ndiags))
        if(dsize3out.gt.0) then
          allocate(
     &       subdd%v5dout(I_0H:I_1H,J_0H:J_1H,dsize3out,nperiod,ndiags))
          subdd%v5dout = 0.
        endif
      endif

      do k=1,ndiags
        if(trim(name_tmp(k)).eq.'no output') cycle
        call add_var(subdd%cdl0,
     &       'float '//trim(name_tmp(k))//trim(dimstr),
     &       units=trim(units_tmp(k)),
     &       long_name=trim(lname_tmp(k)))
        if(denom_tmp(k).ne.0) then
          call add_varline(subdd%cdl0,trim(name_tmp(k))//
     &         ':missing_value = -1.e30f ;')
        endif
      enddo

      if(subdd%catshape_accum.eq.'aijph') deallocate(subdd%sched) ! not needed

      subdd%initialized = .true.

      return
      end subroutine create_group

      subroutine create_solo(subdd,nsubdd,vname,is_inst)
      use model_com, only : nday,dtsrc
      use model_com, only : itimei,itimee
      implicit none
      type(subdd_type) :: subdd
      integer :: nsubdd
      character(len=*) :: vname
      logical :: is_inst
c
      integer :: ndiags,nperiod
      character(len=32) :: dname
      integer :: i_0h,i_1h,j_0h,j_1h

      ndiags = 1

      if(write_one_file) then
        !nperiod = ceiling(real(itimee-itimei)/real(nsubdd))
        nperiod = ceiling(real(itimee-itimei_subdd)/real(nsubdd))
      else
        nperiod = ceiling(real(days_per_file*nday)/real(nsubdd))
      endif

      subdd%grpname = vname
      subdd%strdimlen = 'k'//trim(vname)

      subdd%nsubdd = nsubdd
      subdd%dtime = DTsrc*subdd%nsubdd/3600d0
      subdd%nperiod = nperiod

      subdd%catshape = 'solo'
      subdd%category = vname

      allocate(subdd%nacc(nperiod,subdd_nsched))
      subdd%nacc = 0
      subdd%nacc(:,sched_inst) = 1

      subdd%ndiags = ndiags
      allocate(subdd%scale(ndiags))
      allocate(subdd%sched(ndiags))
      allocate(subdd%reduc(ndiags))
      allocate(subdd%name(ndiags))
      allocate(subdd%denom(ndiags))
      subdd%scale(:) = 1.

      subdd%name(1) = vname
      subdd%reduc(:) = reduc_avg
      subdd%denom(:) = 0
      subdd%reduc = reduc_avg

      subdd%is_inst = is_inst

      if(is_inst) then
        subdd%sched(:) = sched_inst
      else
        subdd%sched(:) = sched_src
      endif

      return
      end subroutine create_solo

      end module subdd_mod

      subroutine parse_subdd
!@sum parse_subdd parse sub daily diag requests
!@+   and declare subdaily diag metadata and allocate space
!@+   for requested outputs
!@auth M. Kelley
      use model_com, only : dtsrc,nday,itime
      use resolution, only : lm
      use constant, only : kapa
      use diag_com, only : cdl_ij_template
      use cdl_mod, only : add_var,add_coord,add_varline,add_unlimdim
      use domain_decomp_atm, only : grid,get=>getdomainbounds
      use dictionary_mod
      use subdd_mod, only : write_daily_files,days_per_file,
     &     vinterp_using_timeavgs,write_monthly_files,write_one_file,
     &     create_group,subdd_ngroups,subdd_ngroups_max,subdd_groups,
     &     cdl_ijt,info_type,namedd_strlen,sname_strlen,itimei_subdd,
     &     lmaxsubdd,subdd_npres,subdd_pres,subdd_pk,aijph_l1,aijph_l2,
     &     subddt
      use ghy_com, only: ngm
      implicit none
      integer :: i,j,k,l,kk,listlen,idcat
      integer :: avg_and_inst
      character(len=32) :: dimstr,freqstr
      character(len=32) :: subdd_timeunitstr
      character(len=3) :: c3
      character(len=4) :: c4
      integer :: kdd_thisfreq,nvars_found,ifreq,ngroups
!@dbparam Nsubdd: DT_save_SUBDD =  Nsubdd*DTsrc sub-daily diag freq.
      INTEGER :: Nsubdd = 0
!@var kddmax maximum number of sub-daily diags outputs
      INTEGER, PARAMETER :: kddmax = 110
!@var namedd array of names of sub-daily diags
      character(len=namedd_strlen), DIMENSION(kddmax) :: namedd
!@var kdd total number of sub-daily diags
      INTEGER :: kdd
      character(len=namedd_strlen), DIMENSION(kddmax) :: namedd_thisfreq
      integer, dimension(kddmax) :: ddfreq
      logical :: is_inst(kddmax),instcat,found
      integer, dimension(8), parameter :: allowed_hrfreqs_timeavg=
     &     (/ 1, 2, 3, 4, 6, 8, 12, 24 /)
      integer, dimension(11) :: allowed_freqs_timeavg
      integer, parameter :: nmax_possible=1024,ncats_max=15
      character(len=sname_strlen) :: catshape,grpname
      character(len=sname_strlen), dimension(ncats_max) ::
     &     catshapes,categories
      integer, dimension(ncats_max) :: diaglens
     &     ,input_sizes3 ! size of 3rd dim of 3d arrays passed to inc_subdd
      type(info_type), dimension(nmax_possible,ncats_max) ::
     &     diaglists
      integer :: write_daily_files_int,vinterp_using_timeavgs_int
      integer :: write_one_file_int

      integer :: i_0h,i_1h,j_0h,j_1h
      character(len=80) :: errmsg

      call sync_param( "Nsubdd",Nsubdd)
      if(Nsubdd.le.0) return

      itimei_subdd = itime

      allocate(subdd_groups(subdd_ngroups_max))

      if(is_set_param('write_one_file')) then
        call get_param( "write_one_file" ,write_one_file_int)
        write_one_file = write_one_file_int == 1
        if(write_one_file) then
          days_per_file = 1 ! not really used
        endif
      elseif(is_set_param('days_per_file')) then
        call get_param( "days_per_file" ,days_per_file)
        write_daily_files = .true.
      else
        write_daily_files_int = 0
        call sync_param( "write_daily_files" ,write_daily_files_int)
        write_daily_files = write_daily_files_int == 1
        if(write_daily_files) then
          days_per_file = 1
        else
          days_per_file = 31
        endif
      endif
      write_monthly_files = .not.
     &     (write_daily_files .or. write_one_file)

      vinterp_using_timeavgs_int = 0
      call sync_param("vinterp_using_timeavgs" ,
     &     vinterp_using_timeavgs_int)
      vinterp_using_timeavgs = vinterp_using_timeavgs_int == 1

      allowed_freqs_timeavg(1) = 1
      allowed_freqs_timeavg(2:9) = allowed_hrfreqs_timeavg(:)*(nday/24)
      allowed_freqs_timeavg(10) = days_per_file
      allowed_freqs_timeavg(11) = days_per_file*nday

      call sync_param( "LmaxSUBDD",LmaxSUBDD)

      call get_subdd_strings('database',-9999,subddt)

c
c count/parse names
c
      kdd = 0
      listlen = len_trim(subddt)
      do k=1,listlen-1
        if(subddt(k:k).ne.' ' .and. subddt(k+1:k+1).eq.' ')
     &       kdd = kdd + 1
      enddo
      if(listlen.gt.0) kdd = kdd + 1
      if (kdd.gt.kddmax) call stop_model
     *     ("Increase kddmax: No. of sub-daily diags too big",255)
      read(subddt,*) namedd(1:kdd)
      do k=1,kdd
        j = index(namedd(k),':')
        i = len_trim(namedd(k))
        is_inst(k) = .false.
        ddfreq(k) = nsubdd
        if(j.gt.0 .and. j.lt.i) then
          if(j.gt.1 .and. j.lt.i) then
            freqstr = namedd(k)(j+1:i)
            l = len_trim(freqstr)
            if(freqstr(l:l).eq.'i') then
              is_inst(k) = .true.
              freqstr = freqstr(1:l-1)
              if(l.gt.1) read(freqstr,*) ddfreq(k)
            else
              read(freqstr,*) ddfreq(k)
            endif
            namedd(k) = namedd(k)(1:j-1)
          else
            call stop_model ("sub-daily parse error",255)
          endif
        endif
      enddo

      do k=1,kdd
        if(ddfreq(k).lt.1 .or. ddfreq(k).gt.days_per_file*nday) then
          errmsg ='parse_subdd: invalid frequency for field '//namedd(k)
          call stop_model (errmsg,255)
        endif
        if(.not.is_inst(k)) then
          if(.not.any(allowed_freqs_timeavg.eq.ddfreq(k))) then
            errmsg = 'parse_subdd: invalid time-averaging frequency '//
     &           'for field '//namedd(k)
            call stop_model (errmsg,255)
          endif
        endif
      enddo

      cdl_ijt = cdl_ij_template

      call add_unlimdim(cdl_ijt,'time')

      call get_subdd_timeunitstr(subdd_timeunitstr)
      call add_var(cdl_ijt,'double time(time);',
     &     units=subdd_timeunitstr)
      call add_varline(cdl_ijt,'time:calendar = "noleap" ;')

      if(is_set_param('subdd_npres')) then
        ! for now, all ijph using the same interp coeffs
        call get_param( 'subdd_npres' ,subdd_npres)
        allocate(subdd_pres(subdd_npres),subdd_pk(subdd_npres))
        call get_param( 'subdd_pres' ,subdd_pres, subdd_npres)
        subdd_pk(:) = subdd_pres(:)**kapa
        call get(grid,i_strt_halo=i_0h,i_stop_halo=i_1h,
     &                j_strt_halo=j_0h,j_stop_halo=j_1h)
        allocate(AIJPh_l1(I_0H:I_1H,J_0H:J_1H,subdd_npres),
     &           AIJPh_l2(I_0H:I_1H,J_0H:J_1H,subdd_npres))
      endif

c
c Compose the lists of known outputs and their shapes/categories.
c
c Some category names are identical to the names of their
c corresponding shapes, but this is purely coincidental.
c If/when separate categories are created for various components,
c as has already been done for tracers, most category
c names will differ from their corresponding shape names.
c To create multiple categories for outputs of a given shape,
c add (calls to) the analogs of ijh_defs et al.

      catshapes(:) = ''
      k = 0

      k = k + 1
      catshapes(k) = 'aijh'; categories(k) = 'aijh'
      input_sizes3(k) = 0
      call ijh_defs (diaglists(1,k),nmax_possible,diaglens(k))

      k = k + 1
      catshapes(k) = 'aijph'; categories(k) = 'aijph'
      input_sizes3(k) = lm
      call ijph_defs(diaglists(1,k),nmax_possible,diaglens(k))

      k = k + 1
      catshapes(k) = 'aijlh'; categories(k) = 'aijlh'
      input_sizes3(k) = lm
      call ijlh_defs(diaglists(1,k),nmax_possible,diaglens(k))

      k = k + 1
      catshapes(k) = 'aijh'; categories(k) = 'rijh'
      input_sizes3(k) = 0
      call rijh_defs (diaglists(1,k),nmax_possible,diaglens(k))

      k = k + 1
      catshapes(k) = 'gijlh'; categories(k) = 'gijlh'
      input_sizes3(k) = ngm
      call gijlh_defs(diaglists(1,k),nmax_possible,diaglens(k))

      k = k + 1
      catshapes(k) = 'aijlh'; categories(k) = 'rijlh'
      input_sizes3(k) = lm
      call rijlh_defs(diaglists(1,k),nmax_possible,diaglens(k))

      k = k + 1
      catshapes(k) = 'aijlh'; categories(k) = 'sijlh'
      input_sizes3(k) = lm
      call sijlh_defs(diaglists(1,k),nmax_possible,diaglens(k))

      k = k + 1
      catshapes(k) = 'aijh'; categories(k) = 'cijh'
      input_sizes3(k) = 0
      call cijh_defs(diaglists(1,k),nmax_possible,diaglens(k))

      k = k + 1
      catshapes(k) = 'aijlh'; categories(k) = 'cijlh'
      input_sizes3(k) = lm
      call cijlh_defs(diaglists(1,k),nmax_possible,diaglens(k))

#ifdef SCM
      k = k + 1
      catshapes(k) = 'aijlh'; categories(k) = 'fijlh'
      input_sizes3(k) = lm
      call fijlh_defs(diaglists(1,k),nmax_possible,diaglens(k))
#endif

#ifdef TRACERS_ON
      k = k + 1
      catshapes(k) = 'aijh'; categories(k) = 'taijh'
      input_sizes3(k) = 0
      call tijh_defs (diaglists(1,k),nmax_possible,diaglens(k))

      k = k + 1
      catshapes(k) = 'aijlh'; categories(k) = 'taijlh'
      input_sizes3(k) = lm
      call tijlh_defs(diaglists(1,k),nmax_possible,diaglens(k))

      k = k + 1
      catshapes(k) = 'aijph'; categories(k) = 'taijph'
      input_sizes3(k) = lm
      call tijph_defs(diaglists(1,k),nmax_possible,diaglens(k))
#endif

c
c check whether each requested diagnostic is in the list
c of declared possible outputs
c
      do k=1,kdd
        found = .false.
        name_check_loop: do idcat=1,ncats_max
        do kk=1,diaglens(idcat)
          if(trim(namedd(k)).eq.trim(diaglists(kk,idcat)%sname)) then
            found = .true.
            exit name_check_loop
          endif
        enddo
        enddo name_check_loop
        if(.not. found) call stop_model('subdd request '//
     &       trim(namedd(k))//' not found',255)
      enddo


c
c compose the list of actual output groups
c
      ngroups = 0
      do ifreq=1,nday*days_per_file
      do avg_and_inst=1,2
        instcat = avg_and_inst == 2
        kdd_thisfreq = 0
        do k=1,kdd
          if(ddfreq(k).ne.ifreq) cycle
          if(instcat) then
            if(.not.is_inst(k)) cycle
          else
            if(     is_inst(k)) cycle
          endif
          kdd_thisfreq = kdd_thisfreq + 1
          namedd_thisfreq(kdd_thisfreq) = namedd(k)
        enddo
        if(kdd_thisfreq.eq.0) cycle
c
        write(c3,'(i3)') ifreq
        c4 = adjustl(c3)
        if(instcat) c4=trim(c4)//'i'
c
        do idcat=1,ncats_max
          catshape = catshapes(idcat)
          if(len_trim(catshape).eq.0) exit
          grpname = trim(categories(idcat))//trim(c4)
          call create_group(
     &         subdd_groups(ngroups+1),ifreq,
     &         catshape,categories(idcat),grpname,
     &         namedd_thisfreq,kdd_thisfreq,instcat,
     &         diaglists(1,idcat),diaglens(idcat),
     &         nvars_found,input_sizes3(idcat))
          if(nvars_found.gt.0) ngroups = ngroups + 1
        enddo
      enddo ! avg_and_inst
      enddo ! ifreq
      subdd_ngroups = ngroups

c      if(grid%gid.eq.0) then
c        write(6,*) 'ngroups ',ngroups
c        do k=1,ngroups
c          write(6,*) 'group '
c     &         ,subdd_groups(k)%grpname
c     &         ,subdd_groups(k)%nsubdd
c     &         ,subdd_groups(k)%is_inst
c     &         ,subdd_groups(k)%name(1:subdd_groups(k)%ndiags)
c        enddo
c      endif
c      call stop_model('got here',255)

      return
      end subroutine parse_subdd

      subroutine get_subdd_strings(from,fid,subddt)
      use subdd_mod, only : subddt_len,num_subdd_str,subdd_strlen
      use dictionary_mod, only : sync_param
      use domain_decomp_atm, only : grid
      use pario, only : read_attr
      implicit none
      character(len=*) :: from
      integer :: fid
!@var subddt = subdd + subdd1,2,3 = all variables for sub-daily diags
      character(len=subddt_len) :: subddt
!
!@dbparam subdd string contains variables to save for sub-daily diags
!@+
!@dbparam subdd1-subdd9, subd10-subd17 additional strings of requested output

      character(len=6) :: strnames(num_subdd_str)
      character(len=subdd_strlen) :: str
      integer :: n,idum

      strnames = (/
     &     'subdd ','subdd1','subdd2','subdd3','subdd4','subdd5',
     &     'subdd6','subdd7','subdd8','subdd9','subd10','subd11',
     &     'subd12','subd13','subd14','subd15','subd16','subd17'
     &     /)

      subddt = ''
      do n=1,num_subdd_str
        str = ''
        if(from.eq.'database') then
          call sync_param( trim(strnames(n)), str)
        elseif(from.eq.'rsf') then
          call read_attr(grid,fid,'cparam',trim(strnames(n)),idum,str)
        endif
        subddt = trim(subddt)//' '//trim(str)
      enddo
      subddt = adjustl(subddt)

      end subroutine get_subdd_strings

      subroutine ijh_defs(arr,nmax,decl_count)
c
c 2D outputs
c
      use model_com, only : dtsrc,nday
      use constant, only : rhow,bygrav
      use TimeConstants_mod, only: SECONDS_PER_DAY
      use subdd_mod, only : info_type,sched_rad,reduc_min,reduc_max
! info_type_ is a homemade structure constructor for older compilers
      use subdd_mod, only : info_type_
      implicit none
      integer :: nmax,decl_count
      type(info_type) :: arr(nmax)
c
c note: next() is a locally declared function to increment decl_count
c

      decl_count = 0

      arr(next()) = info_type_(
     &  sname = 'tsavg',
     &  lname = 'SURFACE AIR TEMPERATURE',
     &  units = 'C'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'tsmin',
     &  lname = 'Minimum Daily Surface Temperature',
     &  units = 'C',
     &  reduc = reduc_min
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'tsmax',
     &  lname = 'Maximum Daily Surface Temperature',
     &  units = 'C',
     &  reduc = reduc_max
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'fireCount',
     &  lname = 'FIRE COUNT FOR DYN BIOBURN',
     &  units = 'm-2 s-1'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'FLAMM',
     &  lname = 'VEGETATION FLAMMABILITY',
     &  units = 'none'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'FLAMM_prec',
     &  lname = 'prec for FLAMMABILITY',
     &  units = 'mm/day'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'FLAMM_rh',
     &  lname = 'rh for FLAMMABILITY',
     &  units = 'none'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'FVDEN',
     &  lname = 'FIRE MODEL VEGETATION DENSITY',
     &  units = 'none'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'f_ignCG',
     &  lname = 'FRAC OF BB DUE TO CG LIGT IGN. ONLY',
     &  units = 'm-2 s-1'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'wsavg',
     &  lname = 'SURFACE WIND SPEED',
     &  units = 'm/s'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'wsmax',
     &  lname = 'Maximum Daily Surface Wind Speed',
     &  units = 'm/s',
     &  reduc = reduc_max
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'prec',
     &  lname = 'PRECIPITATION',
     &  units = 'mm/day',
     &  scale = real(nday,kind=8)
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'evap',
     &  lname = 'EVAPORATION',
     &  units = 'mm/day',
     &  scale = real(nday,kind=8)
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'snowfall',
     &  lname = 'SNOW FALL (H2O EQUIV)',
     &  units = 'mm/day',
     &  scale = real(nday,kind=8)
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'slp',
     &  lname = 'SEA LEVEL PRESSURE',
     &  units = 'mb'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'pw',
     &  lname = 'PRECIPITABLE WATER',
     &  units = 'kg/m^2'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'us',
     &  lname = 'SURFACE ZONAL WIND',
     &  units = 'm/s'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'vs',
     &  lname = 'SURFACE MERIDIONAL WIND',
     &  units = 'm/s'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'qs',
     &  lname = 'SURFACE SPECIFIC HUMIDITY',
     &  units = 'kg/kg'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'rs',
     &  lname = 'Surface Relative Humidity',
     &  units = '%',
     &  scale = 1d2
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'rsmin',
     &  lname = 'Minimum Daily Surface Relative Humidity',
     &  units = '%',
     &  scale = 1d2,
     &  reduc = reduc_min
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'rsmax',
     &  lname = 'Maximum Daily Surface Relative Humidity',
     &  units = '%',
     &  scale = 1d2,
     &  reduc = reduc_max
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'mcp',
     &  lname = 'MOIST CONVECTIVE PRECIPITATION',
     &  units = 'mm/day',
     &  scale = real(nday,kind=8)
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'ssp',
     &  lname = 'STRATIFORM PRECIPITATION',
     &  units = 'mm/day',
     &  scale = real(nday,kind=8)
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'sst',
     &  lname = 'SEA SURFACE TEMPERATURE',
     &  units = 'C'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'SIT',
     &  lname = 'Surface Sea/Lake Ice Temperature',
     &  units = 'C'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'LIT',
     &  lname = 'Surface Land Ice Temperature',
     &  units = 'C'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'GT1',
     &  lname = 'Level 1 Ground Temperature, Land',
     &  units = 'C'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'FOICE',
     &  lname = 'Ocean Ice Cover',
     &  units = 'fraction'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'FOOPN',
     &  lname = 'Ice-Free Ocean Cover',
     &  units = 'fraction'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'FLKICE',
     &  lname = 'Lake Ice Cover',
     &  units = 'fraction'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'FLKOPN',
     &  lname = 'Ice-Free Lake Cover',
     &  units = 'fraction'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'FLICE',
     &  lname = 'Land Ice Cover',
     &  units = 'fraction'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'FLOPN',
     &  lname = 'Ice-Free Land Cover',
     &  units = 'fraction'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'p_surf',
     &  lname = 'SURFACE PRESSURE',
     &  units = 'hPa (mb)'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'z_surf',
     &  lname = 'surface height',
     &  units = 'm',
     &  scale = bygrav
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'gtempr',
     &  lname = 'SKIN RADIATIVE TEMPERATURE',
     &  units = 'K'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'gtemp',
     &  lname = 'SKIN TEMPERATURE',
     &  units = 'K'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'ustar',
     &  lname = 'FRICTION VELOCITY',
     &  units = 'm/s'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'pblht',
     &  lname = 'planetary boundary layer height',
     &  units = 'm'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'shflx',
     &  lname = 'SENSIBLE HEAT FLUX',
     &  units = 'W/m^2'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'lhflx',
     &  lname = 'LATENT HEAT FLUX',
     &  units = 'W/m^2'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'pwv',
     &  lname = 'PRECIPITABLE WATER VAPOR',
     &  units = 'kg/m^2'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'puq',
     &  lname = 'East-west humidity flux (vert sum)',
     &  units = 'kg m^-1 s^-1'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'pvq',
     &  lname = 'North-south humidity flux (vert sum)',
     &  units = 'kg m^-1 s^-1'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'lwp',
     &  lname = 'LIQUID WATER PATH',
     &  units = 'kg/m^2'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'iwp',
     &  lname = 'ICE WATER PATH',
     &  units = 'kg/m^2'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'column_fmse',
     &  lname = 'column of frozen moist static energy',
     &  units = 'J/m^2'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'snowdp',
     &  lname = 'SNOW DEPTH',
     &  units = 'mm'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'c_iwp',
     &  lname = 'CLOUD ICE WATER PATH',
     &  units = 'kg/m^2'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'c_lwp',
     &  lname = 'CLOUD LIQUID WATER PATH',
     &  units = 'kg/m^2'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'mc_lwp',
     &  lname = 'CONVECTIVE CLOUD LIQUID WATER PATH',
     &  units = 'kg/m^2'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'dcnvf',
     &  lname = 'DEEP CONV CLOUD FREQUENCY',
     &  units = '%',
     &  scale = 1d2
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'scnvf',
     &  lname = 'SHALLOW CONV CLOUD FREQUENCY',
     &  units = '%',
     &  scale = 1d2
     &     )
c
c NOT IMPLEMENTED YET
c      arr(next()) = info_type_(
c     &  sname = 'olrsrf',
c     &  lname = 'Outgoing Longwave Radiation at TOA '//
c     &         'accounting for tsurf variability not seen by RADIA',
c     &  units = 'W/m^2'
c     &     )
c
      arr(next()) = info_type_(
     &  sname = 'snowd',
     &  lname = 'Snow depth',
     &  units = 'mm water equiv.',
     &  scale = 1d3/rhow
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'snowc',
     &  lname = 'Snow cover',
     &  units = 'fraction of grid area'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'icef',
     &  lname = 'Ice Fraction Over Open Water',
     &  units = '%',
     &  scale = 1d2
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'rnft',
     &  lname = 'Total runoff',
     &  units = 'mm/day',
     &  scale = SECONDS_PER_DAY/dtsrc
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'smst',
     &  lname = 'Near Surface Soil Moisture',
     &  units = 'kg/m^3'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'cdnc_ls',
     &  lname = 'CDNC stratif. warm clouds',
     &  units = '#/cm^3'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'cdnc_mc',
     &  lname = 'CDNC convect. warm clouds',
     &  units = '#/cm^3'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'cdnc_RB',
     &  lname = 'CDNC large scale screened after Bennartz',
     &  units = '#/cm^3'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'dzwm',
     &  lname = 'height warm conv clouds',
     &  units = 'm'
     &     )
      arr(next()) = info_type_(
     &  sname = 'dzim',
     &  lname = 'height ice conv clouds',
     &  units = 'm'
     &     )
      arr(next()) = info_type_(
     &  sname = 'dzws',
     &  lname = 'height warm large scale clouds',
     &  units = 'm'
     &     )
      arr(next()) = info_type_(
     &  sname = 'dzis',
     &  lname = 'height ice large scale clouds',
     &  units = 'm'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'ctp_mc',
     &  lname = 'Convective Cloud top pressure',
     &  units = 'Pa'
     &     )
c
#ifdef CFMIP3_SUBDD
      arr(next()) = info_type_(
     &  sname = 'tauus',
     &  lname = 'U COMPON OF MOMENTUM SRF DRAG',
     &  units = 'g/m*s^2',
     &  scale = 1d3
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'tauvs',
     &  lname = 'V COMPON OF MOMENTUM SRF DRAG',
     &  units = 'g/m*s^2',
     &  scale = 1d3
     &     )
#endif
c
      arr(next()) = info_type_(
     &  sname = 'r_w_mc',
     &  lname = 'Warm-Cloud effective Radius convective',
     &  units = 'um'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'r_i_mc',
     &  lname = 'Ice-Cloud effective Radius convective',
     &  units = 'um'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'r_w_ls',
     &  lname = 'Warm-Cloud effective Radius Large Scale',
     &  units = 'um'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'r_i_ls',
     &  lname = 'Ice-Cloud effective Radius Large scale',
     &  units = 'um'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'pn',
     &  lname = 'Number Concentration of dg > 0.1 um',
     &  units = '#/m^2'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'apn',
     &  lname = 'Activated Particles Number Concentration',
     &  units = '#/m^2'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'ptrop',
     &  lname = 'Tropopause pressure',
     &  units = 'mb'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'ttrop',
     &  lname = 'Tropopause temperature',
     &  units = 'K'
     &     )
      return
      contains
      integer function next()
      decl_count = decl_count + 1
      next = decl_count
      end function next
      end subroutine ijh_defs

      subroutine ijph_defs(arr,nmax,decl_count)
c
c 3D constant-pressure outputs
c
      use constant, only : bygrav
      use subdd_mod, only : info_type
! info_type_ is a homemade structure constructor for older compilers
      use subdd_mod, only : info_type_
      implicit none
      integer :: nmax,decl_count
      type(info_type) :: arr(nmax)
c
c note: next() is a locally declared function to increment decl_count
c

      decl_count = 0

c
      arr(next()) = info_type_(
     &  sname = 'ucp',
     &  lname = 'east-west velocity',
     &  units = 'm/s'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'vcp',
     &  lname = 'north-south velocity',
     &  units = 'm/s'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'tcp',
     &  lname = 'temperature',
     &  units = 'K'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'qcp',
     &  lname = 'specific humidity',
     &  units = 'kg/kg'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'zcp',
     &  lname = 'height',
     &  units = 'm',
     &  scale = bygrav
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'vortcp',
     &  lname = 'relative vorticity',
     &  units = '1/s'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'wcp',
     &  lname = 'Vertical Velocity',
     &  units = 'm/s'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'rhcp',
     &  lname = 'relative humidity',
     &  units = 'percent'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'cfraccp',
     &  lname = 'Convective Cloud Fraction',
     &  units = 'percent',
     &  scale = 1d2
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'stfraccp',
     &  lname = 'Stratiform Cloud Fraction',
     &  units = 'percent',
     &  scale = 1d2
     &     )

#ifdef TES_SIM
      arr(next()) = info_type_(
     &   sname = 'HDORaw',
     &   lname = 'Raw model HDO',
     &   units = 'kg/kg'
     &     )

      arr(next()) = info_type_(
     &   sname = 'H2ORaw',
     &   lname = 'Raw model H2O',
     &   units = 'kg/kg'
     &     )

      arr(next()) = info_type_(
     &   sname = 'nTES',
     &   lname = 'Number of TES measurements',
     &   units = 'Count'
     &     )

      arr(next()) = info_type_(
     &   sname = 'nTESGoodC',
     &   lname = 'Number of good categorical points',
     &   units = 'Count'
     &     )

      arr(next()) = info_type_(
     &  sname = 'HDOC',
     &  lname = 'HDO from categorical TES operator',
     &  dname = 'nTESGoodC',
     &  units = 'kg/kg'
     &     )

      arr(next()) = info_type_(
     &  sname = 'H2OC',
     &  lname = 'H2O from categorica TES operator',
     &  dname = 'nTESGoodC',
     &  units = 'kg/kg'
     &     )

      arr(next()) = info_type_(
     & sname = 'nTESGoodR',
     & lname = 'Number of good retrieval points',
     & units = 'Count'
     &     )

      arr(next()) = info_type_(
     & sname = 'HDOR',
     & lname = 'HDO from retrieval-based TES operator',
     & dname = 'nTESGoodR',
     & units = 'kg/kg'
     &     )

      arr(next()) = info_type_(
     & sname = 'H2OR',
     & lname = 'H2O from retrieval-based TES operator',
     & dname = 'nTESGoodR',
     & units = 'kg/kg'
     &    )
#endif

      return
      contains
      integer function next()
      decl_count = decl_count + 1
      next = decl_count
      end function next
      end subroutine ijph_defs

      subroutine ijlh_defs(arr,nmax,decl_count)
c
c 3D model-level outputs
c
      use subdd_mod, only : info_type,sched_rad
! info_type_ is a homemade structure constructor for older compilers
      use subdd_mod, only : info_type_
#ifdef CFMIP3_SUBDD
      use model_com, only : dtsrc
#endif
      use constant, only : bygrav,kapa
      implicit none
      integer :: nmax,decl_count
      type(info_type) :: arr(nmax)
c
c note: next() is a locally declared function to increment decl_count
c

      decl_count = 0
c
      arr(next()) = info_type_(
     &  sname = 't',
     &  lname = 'TEMPERATURE',
     &  units = 'K'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'th',
     &  lname = 'potential temperature',
     &  units = 'K',
     &  scale = 1000.**kapa
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'q',
     &  lname = 'SPECIFIC HUMIDITY',
     &  units = 'kg/kg'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'rh',
     &  lname = 'RELATIVE HUMIDITY',
     &  units = '%'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'qcl',
     &  lname = 'cloud water',
     &  units = 'kg/kg'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'qci',
     &  lname = 'cloud ice',
     &  units = 'kg/kg'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'cldss',
     &  lname = 'cloud fraction, stratiform',
     &  units = '-'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'cldmc',
     &  lname = 'cloud fraction, convective',
     &  units = '-'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'z',
     &  lname = 'HEIGHT',
     &  units = 'm',
     &  scale = bygrav
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'u',
     &  lname = 'east-west velocity',
     &  units = 'm/s'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'v',
     &  lname = 'north-south velocity',
     &  units = 'm/s'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'ub',
     &  lname = 'east-west velocity (b-grid)',
     &  units = 'm/s'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'vb',
     &  lname = 'north-south velocity (b-grid)',
     &  units = 'm/s'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'w',
     &  lname = 'vertical velocity',
     &  units = 'm/s'
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'swhr',
     &  lname = 'Shortwave Radiative Heating Rate',
     &  units = 'K/day',
     &  sched = sched_rad
     &     )
c
      arr(next()) = info_type_(
     &  sname = 'lwhr',
     &  lname = 'Longwave Radiative Heating Rate',
     &  units = 'K/day',
     &  sched = sched_rad
     &     )
c
c  indirect effect diagnostics
c
c
      arr(next()) = info_type_(
     &  sname = 'cdncls3d',
     &  lname = 'CDNC warm large scale clouds',
     &  units = 'cm-3'
     &     )

      arr(next()) = info_type_(
     &  sname = 'cdncmc3d',
     &  lname = 'CDNC warm convective clouds',
     &  units = 'cm-3'
     &     )

      arr(next()) = info_type_(
     &  sname = 'ijl_cfwm',
     &  lname = 'cloud fraction warm conv. clouds',
     &  units = '-'
     &     )

      arr(next()) = info_type_(
     &  sname = 'ijl_cfws',
     &  lname = 'cloud fraction warm large scale clouds',
     &  units = '-'
     &     )

      arr(next()) = info_type_(
     &  sname = 'p_3d',
     &  lname = 'pressure on model levels',
     &  units = 'mb'
     &     )
      arr(next()) = info_type_(
     &  sname = 'rh_3d',
     &  lname = 'Rel. Humidity on model levels',
     &  units = '%'
     &     )
      arr(next()) = info_type_(
     &  sname = 'prec_3d',
     &  lname = 'Precip on model levels',
     &  units = 'kg m-2 s-1'
     &     )
      arr(next()) = info_type_(
     &  sname = 'lwcls_3d',
     &  lname = 'LWP large scale cl. on model levels',
     &  units = 'g/m-3'
     &     )
      arr(next()) = info_type_(
     &  sname = 'lwcmc_3d',
     &  lname = 'LWP conv clds on model levels',
     &  units = 'g/m-3'
     &     )
      arr(next()) = info_type_(
     &  sname = 'iwc_3d',
     &  lname = 'IWP on model levels',
     &  units = 'g/m-3'
     &     )
      arr(next()) = info_type_(
     &  sname = 'r_wls_3d',
     &  lname = 'R eff warm ls on model levels',
     &  units = 'um'
     &     )
      arr(next()) = info_type_(
     &  sname = 'r_ils_3d',
     &  lname = 'R eff ice ls on model levels',
     &  units = 'um'
     &     )
      arr(next()) = info_type_(
     &  sname = 'r_wmc_3d',
     &  lname = 'R eff warm conv cld on model levels',
     &  units = 'um'
     &     )
      arr(next()) = info_type_(
     &  sname = 'r_imc_3d',
     &  lname = 'R eff ice conc cld on model levels',
     &  units = 'um'
     &     )
      arr(next()) = info_type_(
     &  sname = 'cod_w_3d',
     &  lname = 'COD warm cld on model levels',
     &  units = '-'
     &     )
      arr(next()) = info_type_(
     &  sname = 'cod_i_3d',
     &  lname = 'COD ice cld on model levels',
     &  units = '-'
     &     )
      arr(next()) = info_type_(
     &  sname = 'ccn01_3d',
     &  lname = 'CCN 0.1% on model levels',
     &  units = 'cm-3'
     &     )
      arr(next()) = info_type_(
     &  sname = 'ccn02_3d',
     &  lname = 'CCN 0.2% on model levels',
     &  units = 'cm-3'
     &     )
      arr(next()) = info_type_(
     &  sname = 'ccn05_3d',
     &  lname = 'CCN 0.1% on model levels',
     &  units = 'cm-3'
     &     )
      arr(next()) = info_type_(
     &  sname = 'pn_3d',
     &  lname = 'NC Dg>100nm on model levels',
     &  units = 'cm-3'
     &     )
      arr(next()) = info_type_(
     &  sname = 'apn_3d',
     &  lname = 'APN on model levels',
     &  units = 'cm-3'
     &     )
#ifdef CFMIP3_SUBDD
       arr(next()) = info_type_(
     &  sname = 'mcamfx',
     &  lname = 'MC Air Mass Flux',
     &  units = 'kg/s',
     &  scale = 100.*bygrav/dtsrc
     &     )
#endif
      return
      contains
      integer function next()
      decl_count = decl_count + 1
      next = decl_count
      end function next
      end subroutine ijlh_defs


      subroutine get_subdd_vinterp_coeffs
      use geom, only : imaxj
      use resolution, only : lm
      use atm_com, only : pmid
      use domain_decomp_atm, only : grid,get=>getdomainbounds
      use subdd_mod, only : subdd_ngroups,subdd_groups,subdd_type,
     &     AIJPh_l1,AIJPh_l2,subdd_npres,inc_subdd
      implicit none
      integer j_0, j_1, i_0, i_1
      integer :: i,j,l, k,ip, idgrp
      type(subdd_type), pointer :: subdd
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo,
     &                  grid%j_strt_halo:grid%j_stop_halo,subdd_npres)
     &     :: onezero
      logical :: any_aijph
      any_aijph = .false.
      do idgrp=1,subdd_ngroups
        subdd => subdd_groups(idgrp)
        if(subdd%catshape.eq.'aijph') then
          any_aijph = .true.
          exit
        endif
      enddo
      if(.not. any_aijph) return
      call get(grid, i_strt=i_0,i_stop=i_1, j_strt=j_0,j_stop=j_1)

      call get_subdd_vinterp_coeffs1(pmid,AIJPh_l1,AIJPh_l2,onezero)

      ! increment counters
      do idgrp=1,subdd_ngroups
        subdd => subdd_groups(idgrp)
        if(.not.subdd%accum_this_step) cycle
        if(subdd%catshape.ne.'aijph') cycle
        k = subdd%ndiags
        if(subdd%catshape_accum.eq.'aijlh') then
          call inc_subdd(subdd,k,pmid,jdim=3)
        else
          call inc_subdd(subdd,k,onezero)
        endif
      enddo
      return
      end subroutine get_subdd_vinterp_coeffs

      subroutine get_subdd_vinterp_coeffs1(pmid,AIJPh_l1,AIJPh_l2,
     &     onezero)
      use geom, only : imaxj
      use resolution, only : lm
      use domain_decomp_atm, only : grid,get=>getdomainbounds
      use domain_decomp_atm, only : hasSouthPole,hasNorthPole
      use subdd_mod, only : subdd_pres,subdd_npres
      implicit none
      real*8, dimension(lm,grid%i_strt_halo:grid%i_stop_halo
     &                    ,grid%j_strt_halo:grid%j_stop_halo) :: pmid
      real*8, dimension(grid%i_strt_halo:grid%i_stop_halo
     &                 ,grid%j_strt_halo:grid%j_stop_halo,subdd_npres)
     &     :: AIJPh_l1,AIJPh_l2,onezero
c
      real*8, dimension(subdd_npres) :: lwdn,lwup,logp_subdd
      real*8, dimension(lm) :: logp
      real*8 :: wtdn,wtup
      integer j_0, j_1, i_0, i_1, i,j,l, ldn,lup
      call get(grid, i_strt=i_0,i_stop=i_1, j_strt=j_0,j_stop=j_1)
      logp_subdd = log(subdd_pres)
      do j=j_0,j_1
      do i=i_0,imaxj(j)
        logp(:) = log(pmid(:,i,j))
        !call xintrp(pmid(1,i,j),lm,subdd_pres,subdd_npres,lwdn,lwup)
        call xintrp(logp,lm,logp_subdd,subdd_npres,lwdn,lwup)
        AIJPh_l1(i,j,:) = lwdn(:)
        AIJPh_l2(i,j,:) = lwup(:)
      enddo
      enddo
      if(hasSouthPole(grid)) then
        do l=1,subdd_npres
          AIJPh_l1(2:i_1,j_0,l) = AIJPh_l1(1,j_0,l)
          AIJPh_l2(2:i_1,j_0,l) = AIJPh_l2(1,j_0,l)
        enddo
      endif
      if(hasNorthPole(grid)) then
        do l=1,subdd_npres
          AIJPh_l1(2:i_1,j_1,l) = AIJPh_l1(1,j_1,l)
          AIJPh_l2(2:i_1,j_1,l) = AIJPh_l2(1,j_1,l)
        enddo
      endif
      do l=1,subdd_npres
      do j=j_0,j_1
      do i=i_0,i_1
        ldn = aijph_l1(i,j,l)    ;  lup = aijph_l2(i,j,l)
       wtdn = aijph_l1(i,j,l)-ldn; wtup = aijph_l2(i,j,l)-lup
        if(wtdn+wtup.gt.0.) then
          onezero(i,j,l) = 1.
        else
          onezero(i,j,l) = 0.
        endif
      enddo
      enddo
      enddo
      contains
      subroutine xintrp(xsrc,nsrc,xdst,ndst,lw1,lw2)
      implicit none
      integer :: nsrc,ndst
      real*8, dimension(nsrc) :: xsrc
      real*8, dimension(ndst) :: xdst
      real*8, dimension(ndst) :: lw1,lw2
      real*8 :: wt1,wt2,dx,s
      integer :: ldst,lsrc,lout
      s = sign(1d0,xsrc(2)-xsrc(1))
      lout = 1
      do while(s*xdst(lout).lt.s*xsrc(1))
        lw1(lout) = 1
        lw2(lout) = 1
        lout = lout+1
      enddo
      lsrc = 1
      do ldst=lout,ndst
        if(s*xdst(ldst).gt.s*xsrc(nsrc)) exit
        do while(s*xsrc(lsrc+1).lt.s*xdst(ldst))
          lsrc = lsrc+1
        enddo
        if(xsrc(lsrc+1).eq.xdst(ldst)) then
          lw1(ldst) = real(lsrc+1,kind=8)+.5
          lw2(ldst) = real(lsrc+1,kind=8)+.5
        else
          dx = xsrc(lsrc+1)-xsrc(lsrc)
          wt2 = (xdst(ldst)-xsrc(lsrc))/dx
          wt1 = 1.-wt2
          lw1(ldst) = real(lsrc  ,kind=8)+wt1
          lw2(ldst) = real(lsrc+1,kind=8)+wt2
        endif
      enddo
      do lout=ldst,ndst
        lw1(lout) = nsrc
        lw2(lout) = nsrc
      enddo
      return
      end subroutine xintrp
      end subroutine get_subdd_vinterp_coeffs1

      subroutine set_subdd_period
!@sum set_subdd_period get index of the current subdd accumulation period
      use subdd_mod, only : subdd_groups,subdd_ngroups
      implicit none
      integer :: k
      do k=1,subdd_ngroups
        call set_subdd_period1(subdd_groups(k))
      enddo
      return
      end subroutine set_subdd_period

      subroutine set_subdd_period1(subdd)
!@sum set_subdd_period get index of the current subdd accumulation period
      use model_com, only : itime,itimei,nday,dtsrc,modelEclock
      use subdd_mod, only : sched_src,subdd_type
     &     ,write_one_file,write_monthly_files,days_per_file
     &     ,itimei_subdd
      implicit none
      integer :: istep,subdd_period,jdate
      type(subdd_type) :: subdd
      real*8 :: timenow
      timenow = real(1+itime-itimei,kind=8)*DTsrc/3600d0
      if(subdd%is_inst) then
        subdd%accum_this_step = mod(1+(itime-itimei),subdd%nsubdd)==0
        if(subdd%accum_this_step) then
          subdd%timelast = timenow
          do subdd_period=1,subdd%nperiod
            if(subdd%nacc(subdd_period,sched_src) == 0) then
              subdd%nacc(subdd_period,sched_src) = 1
              exit
            endif
          enddo
        endif
        subdd_period = max(1,sum(subdd%nacc(:,sched_src)))
      else
        if(write_one_file) then
          !istep = itime-itimei
          istep = itime-itimei_subdd
        elseif(write_monthly_files) then
          call modelEclock%get(date=jdate)
          istep = (jdate-1)*nday + mod(itime,nday)
        else
          istep = mod(itime,days_per_file*nday)
        endif
        subdd_period = 1+istep/subdd%nsubdd
        subdd%accum_this_step = .true.
        subdd%timelast = timenow - .5d0*subdd%dtime
        subdd%nacc(subdd_period,sched_src) =
     &       subdd%nacc(subdd_period,sched_src) + 1
      endif
      subdd%subdd_period = subdd_period
      return
      end subroutine set_subdd_period1

      subroutine get_subdd_timeunitstr(timeunitstr)
      use model_com, only : itimei,iyear1,nday
      implicit none
      character(len=32) :: timeunitstr
      integer :: year1,mon1,day1,jdate1,hour1
      character(len=4) :: amon1,ystr
      character(len=2) :: mstr,dstr,hstr
      call getdte(
     &     itimei,nday,iyear1,year1,mon1,day1,jdate1,hour1,amon1)
      write(ystr,'(i4.4)') year1
      write(mstr,'(i2.2)') mon1
      write(dstr,'(i2.2)') jdate1
      write(hstr,'(i2.2)') hour1
        ! note: assuming simulations all start at HH:00 UTC
      timeunitstr = 'hours since '//
     &     ystr//'-'//mstr//'-'//dstr//' '//hstr//':00 UTC'
      return
      end subroutine get_subdd_timeunitstr

      SUBROUTINE reset_cached_subdd
!@sum reset_cached_subdd resets cached_subdd accumulations
      use subdd_mod, only : sched_inst,subdd_groups,subdd_ngroups
     &     ,reduc_min,reduc_max
      IMPLICIT NONE
      INTEGER k,l
      do k=1,subdd_ngroups
        subdd_groups(k)%nacc = 0
        subdd_groups(k)%nacc(:,sched_inst) = 1
        if(allocated(subdd_groups(k)%v4d)) then
          subdd_groups(k)%v4d=0.
          do l=1,size(subdd_groups(k)%v4d,4)
            if(subdd_groups(k)%reduc(l).eq.reduc_min) then
              subdd_groups(k)%v4d(:,:,:,l) = +1d30
            elseif(subdd_groups(k)%reduc(l).eq.reduc_max) then
              subdd_groups(k)%v4d(:,:,:,l) = -1d30
            endif
          enddo
        elseif(allocated(subdd_groups(k)%v5d)) then
          subdd_groups(k)%v5d=0.
        endif
      enddo
      return
      end subroutine reset_cached_subdd

      subroutine def_rsf_subdd_acc(fid,r4_on_disk)
      use subdd_mod, only : subdd_groups,subdd_ngroups
      use domain_decomp_atm, only : grid
      use pario, only : defvar
      use mdiag_com, only : sname_strlen,lname_strlen
      implicit none
      integer :: fid
      logical :: r4_on_disk
      character(len=sname_strlen) :: grpname
      character(len=lname_strlen) :: dimstr,arrshape
      integer :: k

      do k=1,subdd_ngroups
        grpname = subdd_groups(k)%grpname
        if(.not.r4_on_disk) then
          dimstr = 'nacc_'//trim(grpname)//'('//
     &         'nperiod_'//trim(grpname)//',subdd_nsched)'
          call defvar(grid,fid,subdd_groups(k)%nacc,trim(dimstr))
        endif
        if(.not.r4_on_disk .and.
     &       len_trim(subdd_groups(k)%rsfshape).gt.0) then
          arrshape = subdd_groups(k)%rsfshape
        else
          arrshape = subdd_groups(k)%accshape
        endif
        dimstr =
     &       trim(grpname)//'('//
     &       trim(arrshape)//','//
     &       trim(subdd_groups(k)%strdimlen)//')'
        if(allocated(subdd_groups(k)%v4d)) then
          call defvar(grid,fid,subdd_groups(k)%v4d,trim(dimstr),
     &         r4_on_disk=r4_on_disk)
        elseif(allocated(subdd_groups(k)%v5d)) then
          if(r4_on_disk .and. allocated(subdd_groups(k)%v5dout)) then
            call defvar(grid,fid,subdd_groups(k)%v5dout,trim(dimstr),
     &           r4_on_disk=r4_on_disk)
          else
            call defvar(grid,fid,subdd_groups(k)%v5d,trim(dimstr),
     &           r4_on_disk=r4_on_disk)
          endif
        endif
      enddo

      return
      end subroutine def_rsf_subdd_acc

      subroutine write_subdd_accdata(fid,iaction)
      use resolution, only : im,jm
      use model_com, only : ioread,iowrite,iowrite_single
      use subdd_mod, only : subdd_groups,subdd_ngroups
      use domain_decomp_atm, only : grid,hasSouthPole,hasNorthPole
      use pario, only : write_dist_data,read_dist_data
      use pario, only : write_data,read_data
      use mdiag_com, only : sname_strlen
      implicit none
      integer fid   !@var fid unit number of read/write
      integer iaction !@var iaction flag for reading or writing to file
      integer :: k
      character(len=sname_strlen) :: grpname
      do k=1,subdd_ngroups
        grpname = subdd_groups(k)%grpname
        if(iaction.eq.iowrite) then
          call write_data(grid,fid,'nacc_'//trim(grpname),
     &         subdd_groups(k)%nacc)
        endif
        if(allocated(subdd_groups(k)%v4d)) then
          call write_dist_data(grid,fid,trim(grpname),
     &         subdd_groups(k)%v4d)
        elseif(allocated(subdd_groups(k)%v5d)) then
          if(iaction.eq.iowrite_single .and.
     &         allocated(subdd_groups(k)%v5dout)) then
            call write_dist_data(grid,fid,trim(grpname),
     &           subdd_groups(k)%v5dout)
          else
            call write_dist_data(grid,fid,trim(grpname),
     &           subdd_groups(k)%v5d)
          endif
        endif
      enddo
      return
      end subroutine write_subdd_accdata

      subroutine prep_subdd_acc
      use resolution, only : im,jm,lm
      use subdd_mod, only : subdd_npres,subdd_groups,subdd_ngroups
      use domain_decomp_atm, only : grid,hasSouthPole,hasNorthPole,
     &     get=>getdomainbounds
      use mdiag_com, only : sname_strlen
      implicit none
      integer :: kdn,kup
      real*8 :: wtdn,wtup
      integer :: i,j,k,l,m,n
      real*8, dimension(:,:,:), allocatable ::
     &     pmid,AIJPh_k1,AIJPh_k2,onezero
      integer j_0, j_1, i_0, i_1
      integer j_0h, j_1h, i_0h, i_1h
      call get(grid, i_strt=i_0,i_stop=i_1, j_strt=j_0,j_stop=j_1)
      call get(grid, i_strt_halo=i_0h,i_stop_halo=i_1h,
     &               j_strt_halo=j_0h,j_stop_halo=j_1h)
      do n=1,subdd_ngroups
        if(allocated(subdd_groups(n)%v4d)) then
          do l=1,size(subdd_groups(n)%v4d,4) ! pole fill
          do k=1,size(subdd_groups(n)%v4d,3)
            if(hasSouthPole(grid))
     &           subdd_groups(n)%v4d(2:im, 1,k,l) =
     &           subdd_groups(n)%v4d(   1, 1,k,l)
            if(hasNorthPole(grid))
     &           subdd_groups(n)%v4d(2:im,jm,k,l) =
     &           subdd_groups(n)%v4d(   1,jm,k,l)
          enddo
          enddo
          if(allocated(subdd_groups(n)%sched)) then ! scale
            do l=1,size(subdd_groups(n)%v4d,4)
            do k=1,size(subdd_groups(n)%v4d,3)
              subdd_groups(n)%v4d(:,:,k,l) =
     &        subdd_groups(n)%v4d(:,:,k,l) / max(1,
     &             subdd_groups(n)%nacc(k,subdd_groups(n)%sched(l)) )
            enddo
            enddo
          endif
        elseif(allocated(subdd_groups(n)%v5d)) then
          do m=1,size(subdd_groups(n)%v5d,5) ! pole fill
          do l=1,size(subdd_groups(n)%v5d,4)
          do k=1,size(subdd_groups(n)%v5d,3)
            if(hasSouthPole(grid))
     &           subdd_groups(n)%v5d(2:im, 1,k,l,m) =
     &           subdd_groups(n)%v5d(   1, 1,k,l,m)
            if(hasNorthPole(grid))
     &           subdd_groups(n)%v5d(2:im,jm,k,l,m) =
     &           subdd_groups(n)%v5d(   1,jm,k,l,m)
          enddo
          enddo
          enddo
          if(allocated(subdd_groups(n)%sched)) then
            do m=1,size(subdd_groups(n)%v5d,5) ! scale
            do l=1,size(subdd_groups(n)%v5d,4)
            do k=1,size(subdd_groups(n)%v5d,3)
              subdd_groups(n)%v5d(:,:,k,l,m) =
     &        subdd_groups(n)%v5d(:,:,k,l,m) / max(1,
     &             subdd_groups(n)%nacc(l,subdd_groups(n)%sched(m)) )
            enddo
            enddo
            enddo
          endif
          if(allocated(subdd_groups(n)%v5dout)) then
            ! Currently, the only case for which this is true is
            ! deferred vertical regridding to CP from model levels
            ! using time-averaged pressure.
            allocate(
     &           pmid(lm,i_0h:i_1h,j_0h:j_1h)
     &          ,AIJPh_k1(i_0h:i_1h,j_0h:j_1h,subdd_npres)
     &          ,AIJPh_k2(i_0h:i_1h,j_0h:j_1h,subdd_npres)
     &          ,onezero(i_0h:i_1h,j_0h:j_1h,subdd_npres)
     &           )
            do l=1,size(subdd_groups(n)%v5dout,4) ! loop over times
              m=size(subdd_groups(n)%v5dout,5)
              do j=j_0,j_1
              do i=i_0,i_1
                pmid(:,i,j) = subdd_groups(n)%v5d(i,j,:,l,m)
              enddo
              enddo
              call get_subdd_vinterp_coeffs1(pmid,AIJPh_k1,AIJPh_k2,
     &             onezero)
              subdd_groups(n)%v5dout(:,:,:,l,m) = onezero(:,:,:)
              do m=1,size(subdd_groups(n)%v5dout,5)-1
              do k=1,subdd_npres
              do j=j_0,j_1
              do i=i_0,i_1
                kdn = aijph_k1(i,j,k)    ;  kup = aijph_k2(i,j,k)
               wtdn = aijph_k1(i,j,k)-kdn; wtup = aijph_k2(i,j,k)-kup
               subdd_groups(n)%v5dout(i,j,k,l,m) =
     &              wtdn*subdd_groups(n)%v5d(i,j,kdn,l,m)
     &             +wtup*subdd_groups(n)%v5d(i,j,kup,l,m)
              enddo
              enddo
              enddo
              enddo
            enddo
            deallocate(pmid,AIJPh_k1,AIJPh_k2,onezero)
          endif
        endif
      enddo
      return
      end subroutine prep_subdd_acc

      subroutine defmeta_subdd(fid,subdd)
      use model_com, only : jdendofm,idacc
      use mdiag_com, only : monacc,sname_strlen
      use subdd_mod, only : subdd_type,write_monthly_files,sched_src
      use cdl_mod, only : defvar_cdl,add_vardata
      use pario, only : defvar,write_attr
      use domain_decomp_atm, only : grid
      implicit none
      integer :: fid         !@var fid file id
      type(subdd_type) :: subdd
c
      integer :: k,l,int_dummy
      integer :: nperiod_thismo
      real*8, allocatable :: time_axis(:)
      character(len=sname_strlen) :: grpname,strdimlen

      if(.not. allocated(subdd%name)) return

      grpname = subdd%grpname
      strdimlen = subdd%strdimlen

      if(subdd%is_inst) then
        nperiod_thismo = sum(subdd%nacc(:,sched_src))
      else
        nperiod_thismo = subdd%nperiod
        if(write_monthly_files) then
          do k=1,12
            if(monacc(k).eq.1) then
              nperiod_thismo =
     &             subdd%nperiod*(jdendofm(k)-jdendofm(k-1))/31
              exit
            endif
          enddo
        endif
      endif
      allocate(time_axis(nperiod_thismo))
      if(idacc(1).gt.1) then  ! normal end-of-month write
        time_axis(nperiod_thismo) = subdd%timelast
        do k=nperiod_thismo-1,1,-1
          time_axis(k) = time_axis(k+1) - subdd%dtime
        enddo
      else                    ! IC write
        time_axis(1) = subdd%timelast
        do k=2,nperiod_thismo
          time_axis(k) = time_axis(k-1) + subdd%dtime
        enddo
      endif

      subdd%cdl = subdd%cdl0
      call add_vardata(subdd%cdl,'time',time_axis,
     &     fmtstr='f13.3')  ! appropriate for hours
      deallocate(time_axis)

      call defvar(grid,fid,subdd%scale,
     &     'scale_'//trim(grpname)//'('//trim(strdimlen)//')')
      call defvar(grid,fid,subdd%denom,
     &     'denom_'//trim(grpname)//'('//trim(strdimlen)//')')
      call defvar(grid,fid,subdd%name,
     &     'sname_'//trim(grpname)//'(sname_strlen,'//
     &     trim(strdimlen)//')')
      call defvar_cdl(grid,fid,subdd%cdl,
     &     'cdl_'//trim(grpname)//'(cdl_strlen,kcdl_'//
     &     trim(grpname)//')')
      call defvar(grid,fid,int_dummy,'ntime_'//trim(grpname))
      call write_attr(grid,fid,trim(grpname),
     &     'split_dim',subdd%split_dim)
      if(subdd%tile_dim_out.gt.0) then
        call write_attr(grid,fid,trim(grpname),
     &       'tile_dim_out',subdd%tile_dim_out)
      endif

      return
      end subroutine defmeta_subdd

      subroutine writemeta_subdd(fid,subdd)
      use subdd_mod, only : subdd_type
      use pario, only : write_data
      use cdl_mod, only : write_cdl
      use domain_decomp_atm, only : grid
      use mdiag_com, only : sname_strlen
      implicit none
      integer :: fid         !@var fid file id
      type(subdd_type) :: subdd
c
      integer :: ntime
      character(len=sname_strlen) :: grpname,strdimlen

      if(.not. allocated(subdd%name)) return

      grpname = subdd%grpname
      strdimlen = subdd%strdimlen

      ntime = 1

      call write_data(grid,fid,'ntime_'//trim(grpname),
     &     ntime)
      call write_data(grid,fid,'scale_'//trim(grpname),
     &     subdd%scale)
      call write_data(grid,fid,'denom_'//trim(grpname),
     &     subdd%denom)
      call write_data(grid,fid,'sname_'//trim(grpname),
     &     subdd%name)
      call write_cdl(grid,fid,'cdl_'//trim(grpname),
     &     subdd%cdl)

      return
      end subroutine writemeta_subdd

      subroutine write_subdd_accfile (fname)
      use model_com, only : iowrite_single,xlabel
      use domain_decomp_atm, only : grid
      use pario, only : par_open,par_enddef,par_close,defvar
      use pario, only : write_attr,write_dist_data
      use subdd_mod, only : subdd_groups,subdd_ngroups
      use geom, only : axyp
#ifdef CUBED_SPHERE
      use geom, only : lon2d_dg,lat2d_dg,lonbds,latbds
#endif
      implicit none
!@var fname base name of file to be read or written
      character(len=*) :: fname
      character(len=200) :: tmpname
      logical :: r4_on_disk=.true.
      integer :: fid
      integer :: k
c
      tmpname = trim(fname)//'.nc'
      fid = par_open(grid,trim(tmpname),'create')
c
      call write_attr(grid,fid,'global','xlabel',xlabel)
      call defvar(grid,fid,k,'idacc') ! not really used
      call def_rsf_subdd_acc(fid,r4_on_disk)
      do k=1,subdd_ngroups
        call defmeta_subdd(fid,subdd_groups(k))
      enddo
#ifdef CUBED_SPHERE
      call defvar(grid,fid,lon2d_dg,'lon(dist_im,dist_jm)')
      call defvar(grid,fid,lat2d_dg,'lat(dist_im,dist_jm)')
      call defvar(grid,fid,lonbds,'lonbds(four,dist_im,dist_jm)')
      call defvar(grid,fid,latbds,'latbds(four,dist_im,dist_jm)')
#endif
      call defvar(grid,fid,axyp,'axyp(dist_im,dist_jm)')
c
      call par_enddef(grid,fid)
c
      do k=1,subdd_ngroups
        call writemeta_subdd(fid,subdd_groups(k))
      enddo
      call prep_subdd_acc
      call write_subdd_accdata(fid,iowrite_single)
#ifdef CUBED_SPHERE
      call write_dist_data(grid,fid,'lon',lon2d_dg)
      call write_dist_data(grid,fid,'lat',lat2d_dg)
      call write_dist_data(grid,fid,'lonbds',lonbds,jdim=3)
      call write_dist_data(grid,fid,'latbds',latbds,jdim=3)
#endif
      call write_dist_data(grid,fid,'axyp',axyp)
c
      call par_close(grid,fid)
c
c Reset accumulations after every acc write
c
      call reset_cached_subdd
c
      return
      end subroutine write_subdd_accfile

      subroutine read_subdd_rsf(fname)
      use subdd_mod, only : subdd_groups,subdd_ngroups,rsf_save
      use subdd_mod, only : subddt,subddt_len
      use domain_decomp_atm, only : grid,am_i_root
      use pario, only : par_open,par_close
      use pario, only : read_dist_data,read_data,read_attr
      use dictionary_mod, only : get_param,is_set_param
      use mdiag_com, only : sname_strlen
      implicit none
      character(len=*) :: fname
      integer fid   !@var fid unit number of read/write
      integer :: n
      character(len=sname_strlen) :: grpname
      character(len=subddt_len) :: subddt_rsf
      logical :: do_read(4)
      if(subdd_ngroups.le.0) return
      fid = par_open(grid,trim(fname),'read')
      ! todo: add additional difference checks
      call get_subdd_strings('rsf',fid,subddt_rsf)
      do_read = (/
     &     trim(subddt) == trim(subddt_rsf),
     &     rsf_equals_db_int('nsubdd'),
     &     rsf_equals_db_int('subdd_npres'),
     &     rsf_equals_db_int('lmaxsubdd')
     &     /)
      if(all(do_read)) then
        do n=1,subdd_ngroups
          grpname = subdd_groups(n)%grpname
          call read_data(grid,fid,'nacc_'//trim(grpname),
     &         subdd_groups(n)%nacc, bcast_all=.true.)
          if(allocated(subdd_groups(n)%v4d)) then
            call read_dist_data(grid,fid,trim(grpname),
     &           subdd_groups(n)%v4d)
          elseif(allocated(subdd_groups(n)%v5d)) then
            call read_dist_data(grid,fid,trim(grpname),
     &           subdd_groups(n)%v5d)
          endif
        enddo
      else
        if(am_i_root()) then
          write(6,*) 'WARNING: skipping rsf read of subdd info '//
     &         'because request list has changed'
          write(6,*) do_read
        endif
      endif
      call par_close(grid,fid)
      rsf_save = fname
      return
      contains
      logical function rsf_equals_db_int(parname)
      character(len=*) :: parname
      !
      integer :: val_rsf,val_db

      val_rsf = 0; val_db = 0
      call read_attr(grid,fid,'iparam',trim(parname),n,val_rsf)
      if(is_set_param(trim(parname))) then
        call get_param(trim(parname) ,val_db)
      endif
      rsf_equals_db_int = val_rsf == val_db
      end function rsf_equals_db_int
      end subroutine read_subdd_rsf

      subroutine read_subdd_rsf1(subdd)
      use subdd_mod, only : subdd_type,rsf_save
      use domain_decomp_atm, only : grid
      use pario, only : par_open,par_close
      use pario, only : read_dist_data,read_data
      use mdiag_com, only : sname_strlen
      implicit none
      type(subdd_type) :: subdd
      integer fid   !@var fid unit number of read/write
      character(len=sname_strlen) :: grpname
      if(trim(rsf_save).eq.'NOTAFILE') return
      fid = par_open(grid,trim(rsf_save),'read')
      grpname = subdd%grpname
      call read_data(grid,fid,'nacc_'//trim(grpname),
     &     subdd%nacc, bcast_all=.true.)
      if(allocated(subdd%v4d)) then
        call read_dist_data(grid,fid,trim(grpname),subdd%v4d)
      elseif(allocated(subdd%v5d)) then
        call read_dist_data(grid,fid,trim(grpname),subdd%v5d)
      endif
      call par_close(grid,fid)
      return
      end subroutine read_subdd_rsf1
