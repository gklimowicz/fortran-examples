#!/bin/bash
# compile regridding tools
gmake ll2cs
gmake ncll2cs
gmake input

# CDN
./remap.pl -par regridc1x1.par -in AL30RL360X180.ext -out CD_CS90

# SOILCARB
./remap.pl -par regride2x2.5.par -in soilcarb_top30cm_nmaps_2x2.5bin.dat -out soilcarb_C90

#     CROPS uses a land mask:  
#     CROPS_288X180N.ext using Reto's data & programs in
#     dirac:/archive/u/rruedy/GISS/crops.tar
#     original data defined on a  1/2 x 1/2 degree grid.
#     from Ramankutti and Foley and covers the period 1700-1992
#     http://www.earthsystematlas.org/explanations/crop_data/crop_data_technical.html
./remap.pl -par regrid1x1Q.par -in CROPS_288X180N.ext -out CROPS_CS90

#     top_index_360x180.ij.ext has been extended using tools from
#     /discover/nobackup/dgueyffi/create_top_index_inputfile
./remap.pl -par regridc1x1.par -in top_index_360x180.ij.ext -out top_index_CS90 


#     regrid all other quantities which need pre/post processing
./inputll2cs
