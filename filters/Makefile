#############################################################################

ifeq ($(strip $(BOOST_ROOT)),)
	BOOST_ROOT = /usr/local/include
endif

CC:=g++
SL:=g++ -shared

CFLAGS:= -O2 -fPIC -Wl,--no-as-needed

# ROOT includes and libs
ROOTINC:=`root-config --libs` -lMinuit -lCore 
ROOTLIB:=`root-config --cflags`

# AraROOT includes and libs
ARAROOTINC:=-DARA_UTIL_EXISTS -I$(ARA_UTIL_INSTALL_DIR)/include
ARAROOTLIB:=-L$(ARA_UTIL_INSTALL_DIR)/lib -lAraEvent -lRootFftwWrapper

# Total includes and libs
INCLUDES:=$(ROOTINC) $(ARAROOTINC)
LIBS:=$(ROOTLIB) $(ARAROOTLIB)

# Find all .cxx files (FILTERS) and create variable with object suffix (.o, OBJECTS) for these files
FILTERS=$(basename $(wildcard *.cxx))
OBJECTS=$(addsuffix .o,$(FILTERS))

# End result is shared library libASFilters.so
all : libASFilters.so

%.o : %.cxx
	@echo ""
	@echo -e "\tCompiling: "$<
	@echo -e "\tCommand:\n"
	$(CC) $(CFLAGS) $(LIBS) $(INCLUDES) -c $< -o $@
	@echo ""

libASFilters.so : $(OBJECTS)
	@echo ""
	@echo -e "\tMaking Shared Library: "$@
	@echo -e "\tCommand:\n"
	$(SL) $(CFLAGS) $(LIBS) $(INCLUDES) $(OBJECTS) -o $@
	@echo "Done."

clean:
	@rm -f libASFilters.so
	@rm -f $(OBJECTS)

#############################################################################
