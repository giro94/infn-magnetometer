# Grab the targets and sources as two batches
SOURCES = $(wildcard src/*.cxx)
TARGETS = $(patsubst src/%.cxx,bin/%,$(SOURCES))

# Debugging/profiling flags
#DEBUG = -g -pg

# Figure out the architecture
UNAME_S = $(shell uname -s)


# Clang compiler
ifeq ($(UNAME_S), Darwin)
        CXX = clang++
        CC  = clang
        FLAGS = -std=c++0x -O3 $(DEBUG)
        CLHEP_INCLUDE = -I/usr/local/Cellar/clhep/2.4.7.1/include/
        CLHEP_LIB  = -L/usr/local/Cellar/clhep/2.4.7.1/lib/ -lCLHEP-Random-2.4.7.1
endif



FLAGS += $(shell root-config --cflags)
LIBS  += $(shell root-config --libs)

# More third party libs
#LIBS  += -lboost_filesystem
INCLUDE += -I ../../../
INCLUDE += $(CLHEP_INCLUDE)
LIBS += $(CLHEP_LIB)
LIBS += -lgsl -lgslcblas

all: $(TARGETS)

bin/%: src/%.cxx
	$(CXX) $(FLAGS) $(INCLUDE) -o $@ $< $(LIBS)

clean:
	rm -f $(OBJECTS) $(TARGETS)
