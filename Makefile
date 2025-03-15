# Define directories
SRCDIR=src
BUILDDIR=build

# Define default target directories and run directories
ifeq ($(MAKECMDGOALS),)
  OBJDIR=$(BUILDDIR)/obj/run
  BINDIR=$(BUILDDIR)/bin/run
else
  OBJDIR=$(BUILDDIR)/obj
  BINDIR=$(BUILDDIR)/bin
endif

# Define source files and objects
SOURCES=$(wildcard $(SRCDIR)/*.cc)
OBJECTS=$(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)

MARCH:=$(shell root-config --arch)
CXX:=  $(shell root-config --cxx)
ifneq ($(findstring g,$(CXX)),g)
 ifneq ("$(findstring oneAPI,$(INTELDIR))","")
   ifndef NOICX
     CXX:=icpx
   endif
 endif
endif

# Set default optimization flags
OPTFLAGS?=-O3
CERNLIB?=

#----
VERSION6     := $(shell $(ROOTSYS)/bin/root-config --version | cut -b1-1)
#----
VERSION      := $(shell $(ROOTSYS)/bin/root-config --version | cut -b1-4)
ifeq ($(VERSION),5.27)
  VERSIONP=
else
  #VERSIONP= $(VERSION)_7
  VERSIONP= $(VERSION)_$(USLC)
endif

ifneq ("$(wildcard /etc/redhat-release)","")
  RH_Version := $(shell cat /etc/redhat-release|awk -F'[^0-9]+' '{ print $$2 }')
  VERSIONP = $(VERSION)_$(RH_Version)
endif

OBJS=$(SRC:%.cc=%.o)

############
CXXFLAGS=$(OPTFLAGS) -qopenmp -std=c++11 -Wabi-tag
ifeq ($(findstring icpx,$(CXX)),icpx)
  CXXFLAGS=$(OPTFLAGS) -qopenmp -std=c++11 -axSSE4.2,SSE4.1,SSSE3,AVX,CORE-AVX2,CORE-AVX512 -fp-model precise -Wabi
else ifeq ($(findstring g,$(CXX)),g)
  CXXFLAGS=$(OPTFLAGS) -std=c++11 -Wall -Wabi-tag
endif 
CXXFLAGS+= -D_GLIBCXX_USE_CXX11_ABI=0

############
INCLUDEDIRS:=-Iinclude -isystem $(AMSWD)/include -isystem $(ROOTSYS)/include -isystem $(AMSWD)/GenFit/include

############
LIBAUXS:=-L/cvmfs/ams.cern.ch/Offline/CERN/NagLib -lnag64 -L/cvmfs/ams.cern.ch/Offline/CERN/2005/lib -lmathlib
LIBAUXS+= -lifcore -lifcoremt -lsvml -lintlc -lirng -static-libgcc
ifeq ($(findstring g,$(CXX)),g)
  LIBAUXS:=-L/cvmfs/ams.cern.ch/Offline/CERN/NagLib -lnaggcc64 -L/cvmfs/ams.cern.ch/Offline/CERN/2005.gcc64.44/lib -lmathlib
  LIBAUXS+= -lgfortran -lm -lgomp
endif
ifeq ($(VERSION6),6)
  LIBAUXS:=-lCling -lz $(LIBAUXS)
endif

############
DEFINES:=-D_PGTRACK_ -DUSEGENFIT -D__ROOTSHAREDLIBRARY__
ifeq ($(VERSION6),6)
  DEFINES+= -DVERSION6
endif
DEFINES+= -D_IONL1PRESCALE_ -D_USELAYERRES_
############
APP0=
APP1=
APP2=
APP3=
APP4=
APP5=
APP6=
APP7=
APP8=
ifdef USEPASS7
  DEFINES+= -D_USEPASS7_
  APP0:=P8
endif
ifdef USEEVENTORDER
  DEFINES+= -D_USEEVENTORDER_
  APP1:=ORD
endif
#-------
ifdef USEPRHEION
  DEFINES+= -D_PRHEPRESCALE3_
  APP5:=$(APP5)H
endif
ifdef USEPRL1
  DEFINES+= -D_PRL1PRESCALE_
  APP5:=$(APP5)P
endif
ifdef USEHEL1
  DEFINES+= -D_HEL1PRESCALE_
  APP5:=$(APP5)A
endif
ifdef USEHEINNER
  DEFINES+= -D_HEINNERPRESCALE_
  APP5:=$(APP5)B
endif
ifdef USEIONINNER
  DEFINES+= -D_IONINNERPRESCALE_
  APP5:=$(APP5)I
endif
ifdef USESAVENEG
  DEFINES+= -D_SAVENEGSCALE_
  APP5:=$(APP5)N
endif
#-------
ifdef USENEWL1L9G
  DEFINES+= -D_USENEWL1L9G_
  APP2:=NEWG
endif
ifdef USEMCTKRAW
  DEFINES+= -D_USEMCTKRAW_
  APP3:=MCTKRAW
endif
ifdef USEONEEV
  DEFINES+= -D_USEONEEV_
  APP4:=USEONEEV
endif
ifdef USEADDTKHIT
  DEFINES+= -D_USEADDTKHIT_
  APP6:=MTKHIT
endif
ifdef USECALIB
  DEFINES+= -D_USECALIB_
  APP7:=CALIB
endif
ifdef USENOLINEARCOR
  DEFINES+= -D_USENOLINEARCOR_
  APP7:=NOLINEAR
endif
ifdef USEKALMANFIT
  DEFINES+= -D_USEKALMANFIT_
  APP8:=KM
endif

# Define directories and targets
EXE=$(BINDIR)/betaFitter_ROOT$(VERSION6)SLC6$(APP0)$(APP1)$(APP5)$(APP2)$(APP3)$(APP4)$(APP6)$(APP7)$(APP8)

# Define static library path and lib directory
NTUPLE_LIBDIR=$(AMSWD)/lib/$(MARCH)$(VERSIONP)
NTUPLE_LIB=$(NTUPLE_LIBDIR)/libntuple_slc6_PG.a

.PHONY: all clean init debug

# Default target should be 'all' instead of 'main'
all: init $(EXE)

# Initialize build directories
init:
	@mkdir -p $(OBJDIR)
	@mkdir -p $(BINDIR)

# Compilation rules
$(OBJDIR)/%.o: $(SRCDIR)/%.cc
	$(CXX) $(CXXFLAGS) $(DEFINES) $(INCLUDEDIRS) -c $< -o $@

$(EXE): $(OBJECTS) $(NTUPLE_LIB)
	$(CXX) -o $@ $(OBJECTS) $(CXXFLAGS) $(DEFINES) -L$(NTUPLE_LIBDIR) -lntuple_slc6_PG$(LIB_SUFFIX) $(shell root-config --libs) -lMinuit -lTMVA -lNetx -lGeom -lMathMore -lEG -lTreePlayer -lMLP -lXMLIO $(LIBAUXS) $(CERNLIB)

clean:
	rm -rf $(BUILDDIR)/*

# Handle debug target
debug: OPTFLAGS=-O0 -g
debug: CERNLIB=-L/cvmfs/ams.cern.ch/Offline/CERN/2005/lib -lpacklib -lmathlib -lkernlib
debug: LIB_SUFFIX=_debug
debug: all
