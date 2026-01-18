#!/bin/sh
#
# Copyright (C) 2008-2015 Michele Martone
# 
# This file is part of librsb.
# 
# librsb is free software; you can redistribute it and/or modify it
# under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
# 
# librsb is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
# License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with librsb; see the file COPYING.
# If not, see <http://www.gnu.org/licenses/>.

CMD=" ./rsbench  -oa -Ob --dense 100 --compare-competitors --verbose -R -qH --no-want-ancillary-execs -n 2"
#rm -vfRI scorep-*
rm -vfR ./scorep-*_*_*

export SCOREP_METRIC_PAPI=PAPI_L2_TCM,PAPI_L1_TCM
$CMD

export SCOREP_TOTAL_MEMORY=163840000 # memory may not suffice
export SCOREP_FILTERING_FILE=scorep_filter.filt
cat > $SCOREP_FILTERING_FILE << EOF
SCOREP_REGION_NAMES_BEGIN
  EXCLUDE *
  INCLUDE rsb__do_spmv_uaua rsb__BCSR_spmv_uaua_double_H__tN_r1_c1_uu_sU_dE_uG rsb__BCSR_spmv_uaua_double_C__tN_r1_c1_uu_sU_dE_uG main rsb__mkl_csr_spmv rsb__mkl_coo_spmv
SCOREP_REGION_NAMES_END
EOF
export SCOREP_ENABLE_TRACING=true
$CMD

#export SCOREP_ENABLE_PROFILING=true
#$CMD
