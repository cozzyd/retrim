CXX=g++
LDFLAGS=`root-config --libs`
CXXFLAGS= `root-config --cflags` -fPIC -O3 -g 

LIBDIR=lib
BUILDDIR=build
INCLUDEDIR=include
BINDIR=bin


OBJS := $(addprefix $(BUILDDIR)/, retrimCint.o\
																	CollisionReader.o TrackMaker.o ReadSav3D.o \
																	HistogramReader.o PeriodicTable.o  TableReader.o\
                                  IonizationModel.o EnergyLossModel.o)

INCLUDES := $(addprefix $(INCLUDEDIR)/, $(shell ls $(INCLUDEDIR)))


all: shared

shared: $(LIBDIR)/libretrim.so

$(LIBDIR)/libretrim.so: $(OBJS) | $(LIBDIR) 
	@echo Building shared library $@
	@$(CXX) $(LDFLAGS) $(OBJS) -g -shared -o $@

$(OBJS): | $(BUILDDIR)

$(LIBDIR): 
	mkdir -p $(LIBDIR)

$(BUILDDIR): 
	mkdir -p $(BUILDDIR)

$(BUILDDIR)/retrimCint.cc: $(INCLUDES) LinkDef.h | $(BUILDDIR)
	@echo Running rootcint
	@rootcint -f $(BUILDDIR)/retrimCint.cc -c $(INCLUDES) LinkDef.h


$(BUILDDIR)/%.o: src/%.cc $(INCLUDES) Makefile | $(BUILDDIR) 
	@echo Compiling  $< 
	@$(CXX) $(CXXFLAGS) -I./include  -o $@ -c $< 

$(BUILDDIR)/%.o: build/%.cc $(INCLUDES) Makefile | $(BUILDDIR) 
	@echo Compiling  $< 
	@$(CXX) $(CXXFLAGS) -I. -o $@ -c $< 



clean:
	rm -rf build
	rm -rf lib



