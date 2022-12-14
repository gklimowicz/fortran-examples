from rundeck_section import *
from rundeck_section_defs import *


Title="E002ia_foo.R GISS Model E Tom Clune 02/28/2011"

Description="""E002ia_foo:  This is a major update to the old E1M20 rundeck 
to correspond to AR5 defaults."""

# global definitions
res="M20"

GlobalCPPOptions=[
    "#define NEW_IO"
]

Sections=[
    Grid_M20, # resolution-dependent source files
    Main,
    Atm(res),
    TrAdv_QUS3D,
    Surface,
    OceanPrescribed,
    VegetationEnt
    ]
  


#Label and Namelist:  (next 2 lines)
#E002ia_foo (ModelE1 4x5, 20 lyrs, 1850 atm/ocn)

Namelist="""
 &INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1949,MONTHE=12,DATEE=2,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
/"""


prt_rundeck(Title,Description,GlobalCPPOptions,Sections,Namelist,res)
