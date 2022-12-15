#!/bin/bash
#
# Copyright (C) 2008-2019 Michele Martone
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

# postprocess benchmark results.
#SBATCH -J librsb-job
#SBATCH -D .
#SBATCH -o librsb-job.%J.%N.out  
#SBATCH -e librsb-job.%J.%N.err
#SBATCH --export=NONE
#SBATCH --no-requeue
#SBATCH --get-user-env
#SBATCH --mail-type=end
#SBATCH --time=00:30:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --partition=test
#SBATCH --export=NONE
. /etc/profile
. /etc/profile.d/modules.sh
module load slurm_setup
module load gnuplot
module list
set -ex
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export KMP_AFFINITY="verbose,granularity=fine,compact,1,0"
set -vx
printenv | grep -i KMP | sort
printenv | grep -i SLURM | sort
export PATH=`pwd`:$PATH
# scripts/benchmark.sh cs.mtx
which rsbench
TRF=`mktemp librsb.reference.XXX.tmp`
touch ${TRF}
rsbench -oa -Ob --bench cs.mtx \
	--performance-record-name-prepend librsb-job.${SLURM_JOBID}.
	--verbose --verbose
BRF=`find -name '*.rpr' -newer ${TRF}`
rm ${TRF}
rsbench --read-performance-record ${BRF} > ${BRF/rpr/txt}
ls -ltr ${BRF/rpr/txt}
RSB_PR_WLTC=2 RSB_PR_SR=0 \
 rsbench --read-performance-record ${BRF} > ${BRF/rpr/tex}
exit
