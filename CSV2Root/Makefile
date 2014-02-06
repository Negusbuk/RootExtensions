#
# TARGET is the basename of the file containing main()
# This should be a separate file 

TARGET	= 	CSV2Root

#
# MODULES are all .cc/.hh file combinations containing your
# own classes except the ones which have to put in a
# shared library

MODULES	=	CSV2Root

# 
# Starting from here no changes should be necessary
# 

ARCHITECTURE := $(shell uname)

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

# c-compiler
CC            = gcc
CFLAGS        = -Wall

# c++-compiler
CXX           = g++
CXXFLAGS      = -Wall -fPIC -std=c++0x
ifeq ($(JUSTROOT),1)
CFLAGS       += -DJUSTROOT=1
CXXFLAGS     += -DJUSTROOT=1
endif

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

ALLDEPEND = $(addsuffix .d,$(MODULES))
EXISTDEPEND = $(shell find . -name \*.d -type f -print)

all:  depend $(TARGET)

$(TARGET): $(addsuffix .o,$(MODULES))
	@echo "Linking binary $@"
	@$(LD) $^ -o $@ $(ROOTLIBS)

install: all
	cp $(TARGET) ~/bin

%.d: %.cc
	@echo Making dependency for file $< ...
	@set -e;\
	$(CXX) -M $(CPPFLAGS) $(CXXFLAGS) $< |\
	sed 's!$*\.o!& $@!' >$@;\
	[ -s $@ ] || rm -f $@

%.o: %.cc
	@echo "Compiling $<"
	@$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

depend: $(ALLDEPEND)

clean:
	@rm -f $(addsuffix .o,$(MODULES))
	@rm -f $(TARGET)
	@rm -f *.d
	@rm -f *\~

ifneq ($(EXISTDEPEND),)
-include $(EXISTDEPEND)
endif
