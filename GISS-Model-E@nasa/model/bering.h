c-----------------------------------------------------------------------------
c --- opening the bering strait requires information exchange across a
c --- 'u' face represented in 2 different locations in the tri-pole grid.
c --- the 2 locations are the northern tip of the bering 'inlet' (truncated
c --- bering channel) on the pacific side and the southern (i.e., upper)
c --- tip of the bering inlet in the panam pole patch
c
c --- ipacn,jpac:  grid point north of inlet head on pacific side
c --- ipacs,jpac:  grid point south of inlet head on pacific side
c --- iatln,jatl:  grid point north of inlet head on arctic ocean side
c --- iatls,jatl:  grid point south of inlet head on arctic ocean side
c
c --- thus, the pairs [(ipacn,jpac),(iatln,jatl)],[(ipacs,jpac),(iatls,jatl)]
c --- refer to identical grid cells in physical space.

#ifdef HYCOM2deg
c --- 2deg hycom
      integer,parameter :: ipacn=67,ipacs=68,jpac= 95
      integer,parameter :: iatln= 2,iatls= 1,jatl=156
#endif
#ifdef HYCOM1deg
c --- 1deg hycom refined (389x360) or unrefined (359x360)
      integer,parameter :: ipacn=137,ipacs=138,jpac=189
      integer,parameter :: iatln= 2,iatls= 1,jatl=312
#endif

      logical,parameter :: beropn=.true.    !  true if bering strait open
c-----------------------------------------------------------------------------
