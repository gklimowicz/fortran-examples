# Darwin - specific options

CPP = cpp -P -traditional
CPPFLAGS = -DMACHINE_MAC

# this is a hack to work around Xcode/MacPorts bug with Xcode 11.x.x
# (the compiler can't find the system include files)
### For now commenting it out - use explicit flags in ~/.modelErc
# XCODE_VERSION = $(word 2,$(shell xcodebuild -version))
#XCODE_VERSION_MAJOR = $(word 1,$(subst ., ,$(XCODE_VERSION)))
# ifeq ($(XCODE_VERSION_MAJOR),11)
#   CPATH_HACK=CPATH=/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/include 
# endif

