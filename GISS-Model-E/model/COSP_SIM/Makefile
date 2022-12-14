# makefile for the component COSP/CFMIP
.SUFFIXES:

SUPPORTED_OPTIONS = COSP_PFLUX COSP_USERE COSP_DEBUG

FSRCS = MISR_simulator.f icarus.f
FSRCS += pf_to_mr.f prec_scops.f scops.f

F90SRCS = COSP_drv.F90
F90SRCS += array_lib.F90 atmos_lib.F90 calc_Re.F90 cosp.F90 cosp_constants.F90
F90SRCS += cosp_isccp_simulator.F90 cosp_lidar.F90 cosp_misr_simulator.F90
F90SRCS += cosp_modis_simulator.F90 cosp_radar.F90 cosp_simulator.F90 cosp_stats.F90
F90SRCS += cosp_types.F90 cosp_utils.F90 dsd.F90 format_input.F90 gases.F90
F90SRCS += lidar_simulator.F90 llnl_stats.F90 lmd_ipsl_stats.F90 math_lib.F90
F90SRCS += modis_simulator.F90 mrgrnk.F90 optics_lib.F90 radar_simulator.F90
F90SRCS += radar_simulator_init.F90 radar_simulator_types.F90 scale_LUTs_io.F90 zeff.F90

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
include $(CONFIG_DIR)/rules.mk
