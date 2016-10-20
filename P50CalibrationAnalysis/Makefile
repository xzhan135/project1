# Makefile for analysis template.

.SUFFIXES:
.SUFFIXES: .c .cc .cpp .o

ifndef PG4_CODE
$(error PG4_CODE is not set)
endif

#CXX = g++
CXXFLAGS = --std=c++11 -fPIC `root-config --cflags` -Wall -I. -I/include \
	-I/$(PG4_CODE)/include/Output/ -I/$(PG4_CODE)/include/EventGen/ \
	-I${G4LEDATA}/../../../../include/Geant4/
LDFLAGS = `root-config --libs` -L./ -L${PG4_BUILD}/lib/ -lEventLib

ODIR = obj/
LIBDIR = lib/
LIBEXP = $(LIBDIR)/libMCANA.a

OBJS = $(addprefix $(ODIR)/, SetupDetResponseTree.o)

$(LIBEXP): $(OBJS)
	ar rs $@ $(OBJS)

$(ODIR)/%.o : %.cc %.hh
	$(CXX) $(CXXFLAGS) $< -c -o $@

% : %.cc
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) -o $@ $<

.PHONY: clean

clean:
	rm -f $(ODIR)/*.o *~ $(OBJS)
