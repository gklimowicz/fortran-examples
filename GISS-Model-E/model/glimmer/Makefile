# modelE/ top-level directory
GISSCLIM_DIR=$(shell pwd)/../..
GLIMMER_DIR=$(GISSCLIM_DIR)/../glimmer

default_target :
	cd $(GLIMMER_DIR); make install
	cd $(GISSCLIM_DIR)/model; mkdir -p mod
	cd $(GISSCLIM_DIR)/model/mod; ln -s ../../../glimmer/include/*.mod .
	cd $(GISSCLIM_DIR)/model/glimmer; ln -s ../../../glimmer/lib/*.mod .


depend :
	cd $(GLIMMER_DIR); make install
	cd $(GLIMMER_DIR)/include; ls -1 *.mod | sed 's/$$/:/' > $(GISSCLIM_DIR)/model/glimmer/.depend 

