from rundeck_section import *
from rundeck_section_defs import *


Title="E005ia_foo.R GISS Model E  1850 ocn/atm          larissa        04/15/2010"

Description="""!! E4F40 is for NIsurf=1 (U00a=0.72; U00b=1.68)

!! delete lines starting with '!!' unless E4F40 prepares a q-flux ocean run
!! E4qsF40.R GISS Model E  1850 atm, ocn: q-flux 65m             rar 07/15/2009

!! E4qsF40 = E4F40 with 65m q-flux ocean
E005ia_foo = modelE as frozen in April 2010:
modelE1 (3.0) 2x2.5 hor. grid with 40 lyrs, top at .1 mb (+ 3 rad.lyrs)
atmospheric composition from year 1850
ocean data: prescribed, 1876-1885 climatology
uses turbulence scheme (no dry conv), grav.wave drag
time steps: dynamics 3.75 min leap frog; physics 30 min.; radiation 2.5 hrs
filters: U,V in E-W and N-S direction (after every physics time step)
         U,V in E-W direction near poles (after every dynamics time step)
         sea level pressure (after every physics time step)
"""

# global definitions
res="F40"

GlobalCPPOptions=[
    "#define NEW_IO"
]

Sections=[
    Grid_F40, # resolution-dependent source files
    Main,
    Atm(res),
    TrAdv_TQUS,
    StratDyn,
    Surface,
    OceanPrescribed,
    VegetationEnt,
    RundeckSection(
        name("local settings"),
        parameters_text("""
U00a=0.54      ! above 850mb w/o MC region; tune this first to get 30-35% high clouds
U00b=1.00      ! below 850mb and MC regions; then tune this to get rad.balance
""")
        )
    ]
  


#Label and Namelist:  (next 2 lines)
#E002ia_foo (ModelE1 4x5, 20 lyrs, 1850 atm/ocn)

Namelist="""
 &INPUTZ
 YEARI=1949,MONTHI=12,DATEI=1,HOURI=0, ! pick IYEAR1=YEARI (default) or < YEARI
 YEARE=1949,MONTHE=12,DATEE=2,HOURE=0,     KDIAG=12*0,9,
 ISTART=2,IRANDI=0, YEARE=1949,MONTHE=12,DATEE=1,HOURE=1,
!! suggested settings for E4qsF40:
!! YEARI=1901,MONTHI=1,DATEI=1,HOURI=0,
!! YEARE=1931,MONTHE=1,DATEE=1,HOURE=0,     KDIAG=12*0,9,
!! ISTART=8,IRANDI=0, YEARE=1901,MONTHE=1,DATEE=1,HOURE=1,
/"""


prt_rundeck(Title,Description,GlobalCPPOptions,Sections,Namelist,res)
