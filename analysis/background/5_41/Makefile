.DELETE_ON_ERROR:

ROOTCONFIG  := root-config
ROOTCFLAGS  := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTGLIBS   := $(shell $(ROOTCONFIG) --glibs)

CXX         := g++

INCLUDES    := -I${EICDIRECTORY}/include
DEFINES     := -DWITHRINT
CXXFLAGS    := -O2 -Wall -fPIC $(ROOTCFLAGS) $(INCLUDES) $(DEFINES)

LD          := g++
LDFLAGS     := -O2 $(ROOTLDFLAGS)

LIBS        := $(ROOTGLIBS)                               \
	       -L${EICDIRECTORY}/lib -leicsmear
SRC_FILES   := $(wildcard *.cxx)
EXE_FILES   := $(SRC_FILES:%.cxx=%) 


.PHONY: all clean

all: $(EXE_FILES)


%: %.o
	$(LD) $(LDFLAGS) $^ $(LIBS) -o $@

%.o: %.cxx
	$(CXX) $(CXXFLAGS) -c $< -o $@


clean:
	rm -f $(EXE_FILES:%=%.o) $(EXE_FILES)

