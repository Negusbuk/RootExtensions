#
# TARGET is the basename of the file containing main()
# This should be a separate file 

TARGET	= RootExtensions

#
# MODULES are all .cc/.hh file combinations containing your
# own classes except the ones which have to put in a
# shared library
# The ROOT linkage has to be specified in 'LinkDef.hh'

MODULES	= TMPalette TMGraph TMGraphErrors TMGraph2D TMGraph2DErrors \
	  TArcArrow

SUBDIRS = CSV2Root
# 
# Starting from here no changes should be necessary
# 

ARCHITECTURE := $(shell uname)

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

# c++-compiler
CXX           = g++
CXXFLAGS      = -Wall -fPIC
LD            = g++
LDFLAGS       = -O2
SOFLAGS       = -fPIC -O2 -shared

ifeq ($(ARCHITECTURE),Darwin)
OSXVERSION    = $(shell sw_vers -productVersion | cut -d . -f 2)
OSXTARGET     = 10.$(OSXVERSION)
CXX           = MACOSX_DEPLOYMENT_TARGET=$(OSXTARGET) g++
LD            = MACOSX_DEPLOYMENT_TARGET=$(OSXTARGET) g++
LDFLAGS       = -O2
SOFLAGS       = -dynamiclib -single_module -undefined dynamic_lookup
endif

CXXFLAGS     += $(ROOTCFLAGS)

RLIBMAP	 := rlibmap
RLIBDEP	 := libGraf.so

ALLDEPEND = $(addsuffix .d,$(MODULES))
EXISTDEPEND = $(shell find . -name \*.d -type f -print)

all:  depend lib$(TARGET).so lib$(TARGET).rootmap sub

sub:
	@for dir in $(SUBDIRS); do (cd $$dir; make); done

lib$(TARGET).so: $(addsuffix .o,$(MODULES)) $(TARGET)Dict.o
	@echo "Linking shared library $@"
	@$(CXX) -fPIC -O2 -shared -dynamiclib -undefined dynamic_lookup $^ -o $@

lib$(TARGET).rootmap: lib$(TARGET).so LinkDef.hh
	@echo "Creating library rootmap $@"
	@$(RLIBMAP) -f -o $@ -d lib$(TARGET).so $(RLIBDEP) -c LinkDef.hh

$(TARGET)Dict.cc: $(addsuffix .hh,$(MODULES)) LinkDef.hh
	@echo "Generating dictionary $@"
	rootcint -f $(TARGET)Dict.cc -c $(CPPFLAGS) -p $(addsuffix .hh,$(MODULES)) LinkDef.hh

%.d: %.cc
	@echo Making dependency for file $< ...
	@set -e;\
	$(CXX) -M $(CPPFLAGS) $(CXXFLAGS)  $< |\
	sed 's!$*\.o!& $@!' >$@;\
	[ -s $@ ] || rm -f $@

%.o: %.cc
	@echo "Compiling $<"
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

install: all
	@for dir in $(SUBDIRS); do (cd $$dir; make install); done

depend: $(ALLDEPEND)

doc:
	doxygen documentation/Doxyfile
	cd html
	git push origin gh-pages
	cd ..

clean:
	@rm -f $(addsuffix .o,$(MODULES))
	@rm -f *.so
	@rm -f *Dict.*
	@rm -f *.rootmap
	@rm -f *.d
	@rm -f *~

ifneq ($(EXISTDEPEND),)
-include $(EXISTDEPEND)
endif
