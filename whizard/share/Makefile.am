## Makefile.am -- Makefile for WHIZARD data files
##
## Process this file with automake to produce Makefile.in
#
# Copyright (C) 1999-2022 by 
#     Wolfgang Kilian <kilian@physik.uni-siegen.de>
#     Thorsten Ohl <ohl@physik.uni-wuerzburg.de>
#     Juergen Reuter <juergen.reuter@desy.de>
#     with contributions from
#     cf. main AUTHORS file
#
# WHIZARD is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by 
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# WHIZARD is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
########################################################################

## Subdirectories to configure
SUBDIRS = \
  doc \
  susy \
  models \
  cuts \
  beam-sim \
  tests \
  examples \
  muli \
  SM_tt_threshold_data \
  gui

## Data files needed for running WHIZARD
## to be installed with the 'models' prefix kept
modelsdir = $(pkgdatadir)/models
dist_models_DATA = \
  models/THDM.mdl \
  models/THDM_CKM.mdl \
  models/AltH.mdl \
  models/GravTest.mdl \
  models/HSExt.mdl \
  models/Littlest_Eta.mdl \
  models/Littlest_Tpar.mdl \
  models/Littlest.mdl \
  models/MSSM.mdl \
  models/MSSM_CKM.mdl \
  models/MSSM_Grav.mdl \
  models/MSSM_Hgg.mdl \
  models/NMSSM.mdl \
  models/NMSSM_CKM.mdl \
  models/NMSSM_Hgg.mdl \
  models/NoH_rx.mdl \
  models/PSSSM.mdl \
  models/QCD.mdl \
  models/QED.mdl \
  models/Simplest.mdl \
  models/Simplest_univ.mdl \
  models/SM_ac_CKM.mdl \
  models/SM_ac.mdl \
  models/SM_dim6.mdl \
  models/SM_CKM.mdl \
  models/SM_rx.mdl \
  models/SM_ul.mdl \
  models/SM_top.mdl \
  models/SM_top_anom.mdl \
  models/SM_tt_threshold.mdl \
  models/SM.mdl \
  models/SM_Higgs.mdl \
  models/SM_Higgs_CKM.mdl \
  models/SM_hadrons.mdl \
  models/SSC.mdl \
  models/SSC_2.mdl \
  models/SSC_AltT.mdl \
  models/Template.mdl \
  models/UED.mdl \
  models/WZW.mdl \
  models/Xdim.mdl \
  models/Zprime.mdl \
  models/Threeshl.mdl \
  models/Threeshl_nohf.mdl \
  models/Test.mdl \
  models/Test_schemes.mdl

## SLHA parameter input files
## to be installed with the 'susy' prefix kept
susydir = $(pkgdatadir)/susy
dist_susy_DATA = \
  susy/sps1a.slha \
  susy/sps1ap_decays.slha \
  susy/nmssm.slha 

## Files containing predefined cut sets
## to be installed with the 'cuts' prefix kept
cutsdir = $(pkgdatadir)/cuts
dist_cuts_DATA = \
  cuts/default_cuts.sin

## Files containing predefined beam simulation files
## to be installed with the 'beam-sim' prefix kept
beamsimdir = $(pkgdatadir)/beam-sim
dist_beamsim_DATA = \
  beam-sim/uniform_spread_2.5%.dat

## Files containing complete examples
examplesdir = $(pkgdatadir)/examples
dist_examples_DATA = \
  examples/Z-lineshape.sin \
  examples/W-endpoint.sin \
  examples/casc_dec.sin \
  examples/Zprime.sin \
  examples/LEP_higgs.sin \
  examples/LEP_cc10.sin \
  examples/eeww_polarized.sin \
  examples/circe1.sin \
  examples/HERA_DIS.sin \
  examples/NLO_eettbar_GoSam.sin \
  examples/NLO_eettbar_OpenLoops.sin \
  examples/NLO_NLL_matched.sin \
  examples/LHC_VBS_likesign.sin \
  examples/manual_api_example.c \
  examples/manual_api_example.cc \
  examples/manual_api_example.f90 \
  examples/manual_api_example_1.py \
  examples/manual_api_example_2.py \
  examples/DrellYanMatchingP.sin \
  examples/DrellYanMatchingW.sin \
  examples/DrellYanNoMatchingP.sin \
  examples/DrellYanNoMatchingW.sin \
  examples/EEMatching2P.sin \
  examples/EEMatching2W.sin \
  examples/EEMatching3P.sin \
  examples/EEMatching3W.sin \
  examples/EEMatching4P.sin \
  examples/EEMatching4W.sin \
  examples/EEMatching5P.sin \
  examples/EEMatching5W.sin \
  examples/EENoMatchingP.sin \
  examples/EENoMatchingW.sin

## The data files for the included PDF sets
pdfdir = $(pkgdatadir)/pdf_builtin
dist_pdf_DATA = \
	pdf_builtin/cteq6l.tbl \
	pdf_builtin/cteq6l1.tbl \
	pdf_builtin/cteq6m.tbl \
	pdf_builtin/cteq6d.tbl \
	pdf_builtin/qed6-10gridp.dat \
	pdf_builtin/qed6-10gridn.dat \
        pdf_builtin/mstw2008lo.00.dat \
	pdf_builtin/mstw2008nlo.00.dat \
	pdf_builtin/mstw2008nnlo.00.dat \
	pdf_builtin/ct10.00.pds \
	pdf_builtin/CJ12_max_00.tbl \
	pdf_builtin/CJ12_mid_00.tbl \
	pdf_builtin/CJ12_min_00.tbl \
	pdf_builtin/CJ15LO_00.tbl \
	pdf_builtin/CJ15NLO_00.tbl \
	pdf_builtin/mmht2014lo.00.dat \
	pdf_builtin/mmht2014nlo.00.dat \
	pdf_builtin/mmht2014nnlo.00.dat \
	pdf_builtin/CT14llo.pds \
	pdf_builtin/CT14lo.pds \
	pdf_builtin/CT14n.00.pds \
	pdf_builtin/CT14nn.00.pds \
	pdf_builtin/CT18n.00.pds \
	pdf_builtin/CT18nn.00.pds

## Files containing precompiled cross sections for muli
## to be installed with the 'muli' prefix kept
mulidir = $(pkgdatadir)/muli
dist_muli_DATA = \
  muli/dsigma_cteq6ll.LHpdf.xml \
  muli/pdf_norm_cteq6ll.LHpdf.xml

thresholddir = $(pkgdatadir)/SM_tt_threshold_data
dist_threshold_DATA = \
  SM_tt_threshold_data/download_data.sh \
  SM_tt_threshold_data/threshold_virtual.f90 \
  SM_tt_threshold_data/threshold.f90
