# IRIX64 - specific options

CPP = /lib/cpp -P
CPPFLAGS = -DMACHINE_SGI

# if COMPILER not defined use default IRIX64
COMPILER ?= IRIX64

# hack to deal with missing ranlib
RANLIB = echo

