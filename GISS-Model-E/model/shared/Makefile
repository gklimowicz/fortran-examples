# makefile for the component shared
.SUFFIXES:

SUPPORTED_OPTIONS = ESMF

FSRCS = orbpar.f

F90SRCS  = AbstractAttribute.F90 \
	ArrayBundle_mod.F90 \
	AttributeDictionary.F90 \
	AttributeHashMap.F90 \
	AttributeReference.F90 \
	Attributes.F90 \
	Constants_mod.F90 \
	CubicEquation_mod.F90 \
	Dictionary_mod.F90 \
	FileManager_mod.F90 \
	GaussianQuadrature.F90 \
	GenericType_mod.F90 \
	Geometry_mod.F90 \
	GetTime_mod.F90 \
	KeyValuePair_mod.F90 \
	Parser_mod.F90 \
	MathematicalConstants.F90 \
	PlanetaryParams.F90\
	PlanetParams_mod.F90 \
	PolynomialInterpolator.F90 \
	Precision_mod.F90 \
	Random_mod.F90 \
	RootFinding_mod.F90 \
	RunTimeControls_mod.F90 \
	SpecialFunctions.F90 \
	stop_model.F90 \
	StringUtilities_mod.F90 \
	System.F90 \
	SystemTimers_mod.F90 \
	TimeConstants.F90 \
	Time.F90 \
	Utilities.F90 \
	SystemTools.F90 \
	dast.F90

CSRCS = system_tools.c

# Time Management
F90SRCS += KindParameters.F90 \
	Rational.F90 \
	BaseTime.F90 \
	TimeInterval.F90 \
	AbstractTimeStamp.F90 \
	AnniversaryDate.F90 \
	CalendarDate.F90 \
	CalendarMonth.F90 \
	AbstractCalendar.F90 \
	FixedCalendar.F90 \
	JulianCalendar.F90 \
	PlanetaryCalendar.F90 \
	Time.F90 \
	OrbitUtilities.F90 \
	AbstractOrbit.F90 \
	FixedOrbit.F90 \
	Earth365DayOrbit.F90 \
	ParameterizedEarthOrbit.F90 \
	PlanetaryOrbit.F90 \
	ModelClock.F90

# automatically generated files
CLEAN_LIST = AttributeDictionary.F90 Attributes.F90 RunTimeControls_mod.F90

ifndef GISSCLIM_DIR
       GISSCLIM_DIR := $(shell pwd)/../..
endif
ifndef CONFIG_DIR
       CONFIG_DIR := $(GISSCLIM_DIR)/config
endif

include $(CONFIG_DIR)/base.mk
include $(CONFIG_DIR)/rules.mk



