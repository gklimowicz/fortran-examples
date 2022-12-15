#!/bin/bash
#
# Copyright (C) 2008-2021 Michele Martone
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

# This script is intended for the librsb developer usage.

# This script gets a bunch of useful test matrices.

# our default matrix set is :
#bayer02.mtx  coater2.mtx  crystk03.mtx  ex11.mtx  lhr10.mtx  memplus.mtx  orani678.mtx  raefsky4.mtx  wang4.mtx

# Sources for other matrices:
# https://sparse.tamu.edu/matrices/list_by_nnz.html

#https://sparse.tamu.edu/MM/Pothen/commanche_dual.tar.gz
#https://sparse.tamu.edu/MM/Simon/venkat01.tar.gz
#https://sparse.tamu.edu/MM/Simon/venkat25.tar.gz
#https://sparse.tamu.edu/MM/Simon/venkat50.tar.gz

set -e
set -x

MATRICES="\
https://sparse.tamu.edu/MM/Grund/bayer02.tar.gz		\
https://sparse.tamu.edu/MM/Brethour/coater2.tar.gz	\
https://sparse.tamu.edu/MM/Boeing/crystk03.tar.gz	\
https://sparse.tamu.edu/MM/FIDAP/ex11.tar.gz		\
https://sparse.tamu.edu/MM/Mallya/lhr10.tar.gz		\
https://sparse.tamu.edu/MM/Hamm/memplus.tar.gz		\
https://sparse.tamu.edu/MM/HB/orani678.tar.gz		\
https://sparse.tamu.edu/MM/Simon/raefsky4.tar.gz	\
https://sparse.tamu.edu/MM/Simon/raefsky3.tar.gz	\
https://sparse.tamu.edu/MM/Wang/wang4.tar.gz"

# matrices used in Buluc et al,
# "Reduced-Bandwidth Multithreaded Algorithms for Sparse Matrix-Vector Multiplication"
#  https://doi.org/10.1109/IPDPS.2011.73
#  https://dblp.uni-trier.de/rec/bibtex/conf/ipps/BulucWOD11
MATRICES_CSB="\
https://sparse.tamu.edu/MM/Sandia/ASIC_320k.tar.gz	\
https://sparse.tamu.edu/MM/FEMLAB/sme3Dc.tar.gz		\
https://sparse.tamu.edu/MM/Wissgott/parabolic_fem.tar.gz \
https://sparse.tamu.edu/MM/Mittelmann/cont11_l.tar.gz	\
https://sparse.tamu.edu/MM/Rucci/Rucci1.tar.gz		\
https://sparse.tamu.edu/MM/Norris/torso1.tar.gz		\
https://sparse.tamu.edu/MM/Zaoui/kkt_power.tar.gz	\
https://sparse.tamu.edu/MM/Rajat/rajat31.tar.gz		\
https://sparse.tamu.edu/MM/GHS_psdef/ldoor.tar.gz	\
https://sparse.tamu.edu/MM/Oberwolfach/bone010.tar.gz"

MATRICES_HUGE_BIN_SYM="\
https://sparse.tamu.edu/MM/DIMACS10/delaunay_n24.tar.gz	\
https://sparse.tamu.edu/MM/DIMACS10/rgg_n_2_23_s0.tar.gz	\
https://sparse.tamu.edu/MM/DIMACS10/kron_g500-logn21.tar.gz	\
https://sparse.tamu.edu/MM/DIMACS10/rgg_n_2_24_s0.tar.gz"

MATRICES_HUGE_REAL_SYM="\
https://sparse.tamu.edu/MM/Janna/Flan_1565.tar.gz	\
https://sparse.tamu.edu/MM/Janna/Cube_Coup_dt0.tar.gz	\
https://sparse.tamu.edu/MM/Janna/Cube_Coup_dt6.tar.gz	\
https://sparse.tamu.edu/MM/Janna/Bump_2911.tar.gz	\
https://sparse.tamu.edu/MM/Schenk/nlpkkt160.tar.gz	\
https://sparse.tamu.edu/MM/Fluorem/HV15R.tar.gz	\
https://sparse.tamu.edu/MM/Janna/Queen_4147.tar.gz	\
https://sparse.tamu.edu/MM/Schenk/nlpkkt200.tar.gz	\
https://sparse.tamu.edu/MM/Schenk/nlpkkt240.tar.gz"

MATRICES_HUGE_LAW="\
https://sparse.tamu.edu/MM/LAW/indochina-2004.tar.gz	\
https://sparse.tamu.edu/MM/LAW/uk-2002.tar.gz		\
https://sparse.tamu.edu/MM/LAW/arabic-2005.tar.gz	\
https://sparse.tamu.edu/MM/LAW/uk-2005.tar.gz		\
https://sparse.tamu.edu/MM/LAW/webbase-2001.tar.gz	\
https://sparse.tamu.edu/MM/LAW/it-2004.tar.gz		\
https://sparse.tamu.edu/MM/LAW/sk-2005.tar.gz"

# matrices used in Martone et al,
# "Assembling Recursively Stored Sparse Matrices"
#  https://dblp.uni-trier.de/rec/bibtex/conf/imcsit/MartoneFPT10
MATRICES_CANA_ASSEMBLING_2010_PAPER="\
https://sparse.tamu.edu/MM/Buss/12month1.tar.gz		\
https://sparse.tamu.edu/MM/Schenk_AFE/af_shell10.tar.gz	\
https://sparse.tamu.edu/MM/vanHeukelum/cage15.tar.gz		\
https://sparse.tamu.edu/MM/Mittelmann/cont11_l.tar.gz		\
https://sparse.tamu.edu/MM/DNVS/fcondp2.tar.gz			\
https://sparse.tamu.edu/MM/JGD_GL7d/GL7d19.tar.gz		\
https://sparse.tamu.edu/MM/GHS_psdef/ldoor.tar.gz		\
https://sparse.tamu.edu/MM/Mittelmann/neos.tar.gz		\
https://sparse.tamu.edu/MM/Pajek/patents.tar.gz		\
https://sparse.tamu.edu/MM/Mittelmann/rail2586.tar.gz		\
https://sparse.tamu.edu/MM/JGD_Relat/relat9.tar.gz		\
https://sparse.tamu.edu/MM/FEMLAB/sme3Dc.tar.gz		\
https://sparse.tamu.edu/MM/Gleich/wb-edu.tar.gz"

# matrices used in Martone et al,
# "Use of Hybrid Recursive {CSR/COO} Data Structures in Sparse Matrices-Vector Multiplication"
#  https://dblp.uni-trier.de/rec/bibtex/conf/imcsit/MartoneFGPT10
MATRICES_CANA_RSBCOO_2010_PAPER_SQUARE="\
https://sparse.tamu.edu/MM/Bourchtein/atmosmodl.tar.gz		\
https://sparse.tamu.edu/MM/Vavasis/av41092.tar.gz		\
https://sparse.tamu.edu/MM/vanHeukelum/cage15.tar.gz		\
https://sparse.tamu.edu/MM/Mallya/lhr71.tar.gz			\
https://sparse.tamu.edu/MM/Pajek/patents.tar.gz		\
https://sparse.tamu.edu/MM/Simon/raefsky3.tar.gz		\
https://sparse.tamu.edu/MM/Rajat/rajat31.tar.gz		\
https://sparse.tamu.edu/MM/Bova/rma10.tar.gz			\
https://sparse.tamu.edu/MM/FEMLAB/sme3Dc.tar.gz		\
https://sparse.tamu.edu/MM/Norris/torso1.tar.gz		\
https://sparse.tamu.edu/MM/Simon/venkat01.tar.gz		\
https://sparse.tamu.edu/MM/Gleich/wb-edu.tar.gz"
# and:
MATRICES_CANA_RSBCOO_2010_PAPER_NON_SQUARE="\
https://sparse.tamu.edu/MM/Buss/12month1.tar.gz		\
https://sparse.tamu.edu/MM/JGD_Groebner/c8_mat11_I.tar.gz	\
https://sparse.tamu.edu/MM/JGD_GL7d/GL7d19.tar.gz		\
https://sparse.tamu.edu/MM/Mittelmann/neos.tar.gz		\
https://sparse.tamu.edu/MM/Mittelmann/rail2586.tar.gz		\
https://sparse.tamu.edu/MM/JGD_Relat/rel9.tar.gz		\
https://sparse.tamu.edu/MM/JGD_Relat/relat9.tar.gz		\
https://sparse.tamu.edu/MM/Rucci/Rucci1.tar.gz			\
https://sparse.tamu.edu/MM/Mittelmann/spal_004.tar.gz		
https://sparse.tamu.edu/MM/Meszaros/tp-6.tar.gz"

# matrices used in Martone et al,
# "On BLAS Operations with Recursively Stored Sparse Matrices"
#  https://doi.org/10.1109/SYNASC.2010.72 => https://ieeexplore.ieee.org/document/5715268
#  https://dblp.uni-trier.de/rec/bibtex/conf/synasc/MartoneFPT10
MATRICES_SYNASC_2010_BLASOPS_PAPER="\
https://sparse.tamu.edu/MM/Rucci/Rucci1.tar.gz			\
https://sparse.tamu.edu/MM/FEMLAB/sme3Dc.tar.gz		\
https://sparse.tamu.edu/MM/BenElechi/BenElechi1.tar.gz		\
https://sparse.tamu.edu/MM/PARSEC/Ga41As41H72.tar.gz		\
https://sparse.tamu.edu/MM/Schenk_AFE/af_shell10.tar.gz	\
https://sparse.tamu.edu/MM/Oberwolfach/bone010.tar.gz		\
https://sparse.tamu.edu/MM/GHS_psdef/ldoor.tar.gz		\
https://sparse.tamu.edu/MM/Botonakis/FEM_3D_thermal2.tar.gz	\
https://sparse.tamu.edu/MM/Botonakis/FEM_3D_thermal1.tar.gz	\
https://sparse.tamu.edu/MM/Hollinger/g7jac180.tar.gz		\
https://sparse.tamu.edu/MM/Schenk_ISEI/ohne2.tar.gz		\
https://sparse.tamu.edu/MM/Norris/torso1.tar.gz		\
"

# matrices used in Martone et al,
# "On the Usage of 16 Bit Indices in Recursively Stored Sparse Matrices"
#  https://doi.org/10.1109/SYNASC.2010.77 => https://ieeexplore.ieee.org/document/5715269
MATRICES_SYNASC_2010_16BIT_PAPER_SYMMETRIC="\
https://sparse.tamu.edu/MM/Schenk_AFE/af_shell10.tar.gz	\
https://sparse.tamu.edu/MM/BenElechi/BenElechi1.tar.gz		\
https://sparse.tamu.edu/MM/Oberwolfach/bone010.tar.gz		\
https://sparse.tamu.edu/MM/GHS_psdef/crankseg_1.tar.gz		\
https://sparse.tamu.edu/MM/Boeing/ct20stif.tar.gz		\
https://sparse.tamu.edu/MM/Koutsovasilis/F1.tar.gz		\
https://sparse.tamu.edu/MM/DNVS/fcondp2.tar.gz			\
https://sparse.tamu.edu/MM/Zaoui/kkt_power.tar.gz		\
https://sparse.tamu.edu/MM/GHS_psdef/ldoor.tar.gz		\
https://sparse.tamu.edu/MM/Andrianov/mip1.tar.gz		\
https://sparse.tamu.edu/MM/ND/nd24k.tar.gz			\
https://sparse.tamu.edu/MM/GHS_psdef/s3dkq4m2.tar.gz		\
"
# and:
MATRICES_SYNASC_2010_16BIT_PAPER_GENERAL="\
https://sparse.tamu.edu/MM/Buss/12month1.tar.gz		\
https://sparse.tamu.edu/MM/Bourchtein/atmosmodl.tar.gz		\
https://sparse.tamu.edu/MM/Vavasis/av41092.tar.gz		\
https://sparse.tamu.edu/MM/vanHeukelum/cage15.tar.gz		\
https://sparse.tamu.edu/MM/Mittelmann/cont11_l.tar.gz		\
https://sparse.tamu.edu/MM/JGD_GL7d/GL7d19.tar.gz		\
https://sparse.tamu.edu/MM/Mallya/lhr71.tar.gz			\
https://sparse.tamu.edu/MM/Mittelmann/neos.tar.gz		\
https://sparse.tamu.edu/MM/Pajek/patents.tar.gz		\
https://sparse.tamu.edu/MM/Simon/raefsky3.tar.gz		\
https://sparse.tamu.edu/MM/Mittelmann/rail2586.tar.gz		\
https://sparse.tamu.edu/MM/Rajat/rajat31.tar.gz		\
https://sparse.tamu.edu/MM/JGD_Relat/rel9.tar.gz		\
https://sparse.tamu.edu/MM/JGD_Relat/relat9.tar.gz		\
https://sparse.tamu.edu/MM/Bova/rma10.tar.gz			\
https://sparse.tamu.edu/MM/Rucci/Rucci1.tar.gz			\
https://sparse.tamu.edu/MM/FEMLAB/sme3Dc.tar.gz		\
https://sparse.tamu.edu/MM/Mittelmann/spal_004.tar.gz		\
https://sparse.tamu.edu/MM/Norris/torso1.tar.gz		\
https://sparse.tamu.edu/MM/Meszaros/tp-6.tar.gz		\
https://sparse.tamu.edu/MM/Simon/venkat01.tar.gz		\
https://sparse.tamu.edu/MM/Gleich/wb-edu.tar.gz		\
"

# matrices used in Martone et al,
# "An Improved Sparse Matrix-Vector Multiply Based on Recursive Sparse Blocks Layout"
#  https://doi.org/10.1007/978-3-642-29843-1_69 =>
#  https://link.springer.com/chapter/10.1007%2F978-3-642-29843-1_69 
#  https://dblp.uni-trier.de/rec/bibtex/conf/lssc/MartonePF11
MATRICES_SCICOM11_2011_PAPER="\
https://sparse.tamu.edu/MM/vanHeukelum/cage15.tar.gz		\
https://sparse.tamu.edu/MM/Freescale/circuit5M_dc.tar.gz	\
https://sparse.tamu.edu/MM/Lee/fem_hifreq_circuit.tar.gz	\
https://sparse.tamu.edu/MM/JGD_GL7d/GL7d19.tar.gz		\
https://sparse.tamu.edu/MM/Pajek/patents.tar.gz		\
https://sparse.tamu.edu/MM/Fluorem/RM07R.tar.gz		\
https://sparse.tamu.edu/MM/TSOPF/TSOPF_RS_b2383.tar.gz		\
https://sparse.tamu.edu/MM/Gleich/wikipedia-20070206.tar.gz"

# matrices used in Martone,
# "Efficient multithreaded untransposed, transposed or symmetric sparse matrix-vector multiplication with the Recursive Sparse Blocks format"
#  https://doi.org/10.1016/j.parco.2014.03.008 => https://www.sciencedirect.com/science/article/pii/S0167819114000386?via%3Dihub
#  ( preprint at https://pure.mpg.de/pubman/faces/ViewItemOverviewPage.jsp?itemId=item_2053189 )
MATRICES_PARALLEL_COMPUTING_2014_PAPER_SYMMETRIC="\
https://sparse.tamu.edu/MM/GHS_psdef/audikw_1.tar.gz			\
https://sparse.tamu.edu/MM/Oberwolfach/bone010.tar.gz			\
https://sparse.tamu.edu/MM/DIMACS10/channel-500x100x100-b050.tar.gz	\
https://sparse.tamu.edu/MM/Janna/Cube_Coup_dt6.tar.gz			\
https://sparse.tamu.edu/MM/DIMACS10/delaunay_n24.tar.gz		\
https://sparse.tamu.edu/MM/Dziekonski/dielFilterV3real.tar.gz		\
https://sparse.tamu.edu/MM/DIMACS10/europe_osm.tar.gz			\
https://sparse.tamu.edu/MM/Janna/Flan_1565.tar.gz			\
https://sparse.tamu.edu/MM/Janna/Geo_1438.tar.gz			\
https://sparse.tamu.edu/MM/Dziekonski/gsm_106857.tar.gz		\
https://sparse.tamu.edu/MM/LAW/hollywood-2009.tar.gz			\
https://sparse.tamu.edu/MM/Janna/Hook_1498.tar.gz			\
https://sparse.tamu.edu/MM/DIMACS10/kron_g500-logn21.tar.gz		\
https://sparse.tamu.edu/MM/Janna/Long_Coup_dt6.tar.gz			\
https://sparse.tamu.edu/MM/Schenk/nlpkkt160.tar.gz			\
https://sparse.tamu.edu/MM/Schenk/nlpkkt200.tar.gz			\
https://sparse.tamu.edu/MM/Schenk/nlpkkt240.tar.gz			\
https://sparse.tamu.edu/MM/DIMACS10/rgg_n_2_23_s0.tar.gz		\
https://sparse.tamu.edu/MM/DIMACS10/rgg_n_2_24_s0.tar.gz		\
https://sparse.tamu.edu/MM/DIMACS10/road_usa.tar.gz			\
https://sparse.tamu.edu/MM/Janna/Serena.tar.gz"

# and:
MATRICES_PARALLEL_COMPUTING_2014_PAPER_UNSYMMETRIC="\
https://sparse.tamu.edu/MM/LAW/arabic-2005.tar.gz			\
https://sparse.tamu.edu/MM/JGD_GL7d/GL7d19.tar.gz			\
https://sparse.tamu.edu/MM/Fluorem/HV15R.tar.gz			\
https://sparse.tamu.edu/MM/LAW/indochina-2004.tar.gz			\
https://sparse.tamu.edu/MM/JGD_Relat/relat9.tar.gz			\
https://sparse.tamu.edu/MM/Fluorem/RM07R.tar.gz			\
https://sparse.tamu.edu/MM/LAW/uk-2002.tar.gz"

# matrices used in Martone and Bacchio,
# "PyRSB: Portable Performance on Multithreaded Sparse BLAS Operations"
#  https://conference.scipy.org/proceedings/scipy2021/bib/martone_bacchio_pyrsb.bib => https://conference.scipy.org/proceedings/scipy2021/martone_bacchio_pyrsb.html
MATRICES_SCIPY_PYRSB_2021_PAPER_SYMMETRIC="$MATRICES_PARALLEL_COMPUTING_2014_PAPER_SYMMETRIC"
MATRICES_SCIPY_PYRSB_2021_PAPER_UNSYMMETRIC="$MATRICES_PARALLEL_COMPUTING_2014_PAPER_UNSYMMETRIC"

[[ -d "$1" ]] && { cd "$1" || exit -1 ; }

for m in $MATRICES
do
	mbn=`basename $m`
	mn=${mbn//.tar.gz/}
	mfn=$mn.mtx
	mtn=$mn.mtx.gz
	mtu=${m//.tar.gz/}
	mtu=${mtu//MM/matrices}.html # matrix descriptive www page
	echo "Getting matrix $mn ; for more info, $mtu "

#	file based
#	[ -f $mbn ] || wget $m
#	tar xzf $mbn $mn/$mfn -O > $mfn
	
	# pipe based
	# [ -f $mfn ] || wget $m -O - | tar xzf - $mn/$mfn -O > $mfn || exit -1
	[ -f $mtn ] || { wget $m -O - | tar xzf - $mn/$mfn -O | gzip - > $mtn.tmp && test "`zcat $mtn.tmp | wc -l`" -gt 3 && mv $mtn.tmp $mtn || { rm -f $mtn.tmp; false; } } || exit -1
done
