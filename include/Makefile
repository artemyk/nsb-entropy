CPPOPTSF = cppopts
INCDIRSF = incdirs
LDIRSF = ldirs
LOCINC = ./
CXXFNAME = cxxname

CXX  = `cat $(CXXFNAME)`
CXXFLAGS = `cat $(INCDIRSF)` -I $(LOCINC)  
CXXFLAGS += -c -pedantic -Wall -Wpointer-arith
CXXFLAGS += -Wwrite-strings -Wconversion -Woverloaded-virtual
CXXFLAGS += `cat $(CPPOPTSF)`

LD   = $(CXX)
LDFLAGS = `cat $(LDIRSF)` 
LIBS = -lgmp -lstdc++ -lgsl -lgslcblas -lm

TEST_EXEC = test

HEADS = utility.h utility_2.h utility_3.h compiler.h Makefile $(CPPOPTSF) $(INCDIRSF) $(LDIRSF) machine_select

test: tmp.o random.o dierfc.o specfun.o
	$(CXX) $(LDFLAGS) $(LIBS) dierfc.o random.o tmp.o 
	mv a.out $(TEST_EXEC)

tmp.o: tmp.cc $(HEADS)
	$(CXX) $(CXXFLAGS)  tmp.cc

random.o: random.cpp $(HEADS) random.h specfunctions/specfunctions.h
	$(CXX) $(CXXFLAGS)  random.cpp

dierfc.o: specfunctions/dierfc.c specfunctions/specfunctions.h $(HEADS)
	$(CXX) $(CXXFLAGS)  specfunctions/dierfc.cc

specfun.o: specfun.cc specfun.h $(HEADS)
	$(CXX) $(CXXFLAGS)  specfun.cc

$(CPPOPTSF) $(INCDIRSF) $(LDIRSF): Makefile machine_select
	chmod 755 machine_select
	./machine_select

.PHONY: clean

clean:
	rm $(CPPOPTSF) $(INCDIRSF) $(LDIRSF)
	rm *o
	rm $(TEST_EXEC)
