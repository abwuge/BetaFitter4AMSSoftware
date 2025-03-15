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
SOURCES = $(wildcard $(SRCDIR)/*.cc) $(wildcard $(SRCDIR)/*.C)
OBJECTS = $(SOURCES:$(SRCDIR)/%.cc=$(OBJDIR)/%.o)
OBJECTS := $(OBJECTS:$(SRCDIR)/%.C=$(OBJDIR)/%.o)

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
# Define directories and targets
EXE=$(BINDIR)/betaFitter

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

$(OBJDIR)/%.o: $(SRCDIR)/%.C
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
