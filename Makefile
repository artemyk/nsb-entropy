BPATH  = ./build
LOCINC = ./include
CPPOPTSF = build/cppopts
INCDIRSF = build/incdirs
LDIRSF = build/ldirs
CXXFNAME = build/cxxname

CXX  = `cat $(CXXFNAME)`
CXXFLAGS = `cat $(INCDIRSF)`  -I $(LOCINC) 
#CXXFLAGS += -g
CXXFLAGS += -Wall -Wpointer-arith -pedantic
CXXFLAGS += -Wwrite-strings -Wconversion -Woverloaded-virtual
CXXFLAGS += `cat $(CPPOPTSF)`


LD   = $(CXX)
LDFLAGS = `cat $(LDIRSF)` 
LIBS = -lgmp -lstdc++ -lgsl -lgslcblas 


EDHEAD = EntrData.h EntrData_except.h Makefile $(CPPOPTSF) $(INCDIRSF) $(LDIRSF) machine_select
INCHEAD = $(LOCINC)/utility.h $(LOCINC)/utility_2.h $(LOCINC)/utility_3.h
INCHEAD += $(LOCINC)/compiler.h $(LOCINC)/random.h
INCHEAD += $(LOCINC)/vect_math.h 
NODESRC = node_access.cc node_constr.cc node_count.cc

OBJS  = $(BPATH)/mpi_map.o $(BPATH)/node_access.o \
	$(BPATH)/node_constr.o $(BPATH)/node_count.o \
	$(BPATH)/Dict_m2o.o $(BPATH)/counts.o $(BPATH)/parsers.o\
	$(BPATH)/recorder.o $(BPATH)/simple_recorder.o \
	$(BPATH)/parallel_recorder.o \
	$(BPATH)/EntrData.o $(BPATH)/nsb.o $(BPATH)/nsb_calc.o \
	$(BPATH)/BxiK_interp.o $(BPATH)/parallel_dict.o\
	$(BPATH)/nsb_stat.o  
INCOBJS = $(LOCINC)/dierfc.o $(LOCINC)/random.o $(LOCINC)/specfun.o
EXEC  = $(BPATH)/nsb-entropy


$(BPATH)/$(EXEC): $(OBJS) $(LOCINC)/test
	$(LD) $(LDFLAGS) $(LIBS) $(OBJS) $(INCOBJS)
	mv a.out $(EXEC)

$(LOCINC)/test: $(LOCINC)
	cd $(LOCINC); make
	cd ..

$(BPATH)/parallel_recorder.o: $(EDHEAD) $(INCHEAD) \
	parallel_recorder.h recorder.h parallel_recorder.cc
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/parallel_recorder.o parallel_recorder.cc

$(BPATH)/parallel_dict.o: $(EDHEAD) $(INCHEAD) \
	parallel_dict.cc simple_dicts.h
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/parallel_dict.o parallel_dict.cc

$(BPATH)/BxiK_interp.o: $(EDHEAD)  $(INCHEAD) BxiK_interp.h $(LOCINC)/interp.h BxiK_interp.cc
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/BxiK_interp.o BxiK_interp.cc

$(BPATH)/nsb.o: $(EDHEAD)  $(INCHEAD) nsb.h BxiK_interp.h nsb.cc \
	$(LOCINC)/specfunctions/specfunctions.h \
	$(LOCINC)/integration.h $(LOCINC)/specfun.h 
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/nsb.o nsb.cc 

$(BPATH)/nsb_calc.o: $(EDHEAD)  $(INCHEAD) nsb.h BxiK_interp.h \
	nsb_calc.cc $(LOCINC)/integration.h \
	$(LOCINC)/specfunctions/specfunctions.h  $(LOCINC)/specfun.h 
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/nsb_calc.o nsb_calc.cc

$(BPATH)/nsb_stat.o: $(EDHEAD) node.h mpi_map.h parsers.h Entropy.h nsb.h nsb_stat.cc
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/nsb_stat.o nsb_stat.cc

$(BPATH)/counts.o: $(EDHEAD) counts.h counts.cc 
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/counts.o counts.cc

$(BPATH)/EntrData.o: $(EDHEAD) EntrData.cc 
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/EntrData.o EntrData.cc

$(BPATH)/recorder.o: $(EDHEAD) recorder.h recorder.cc 
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/recorder.o recorder.cc

$(BPATH)/simple_recorder.o: $(EDHEAD) recorder.h simple_recorder.h simple_recorder.cc 
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/simple_recorder.o simple_recorder.cc

$(BPATH)/parsers.o: $(EDHEAD) parsers.h nsb.h parsers.cc 
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/parsers.o parsers.cc

$(BPATH)/Dict_m2o.o: $(EDHEAD) simple_dicts.h Dict_m2o.cc 
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/Dict_m2o.o Dict_m2o.cc

$(BPATH)/mpi_map.o: $(EDHEAD)  mpi_map.h mpi_map.cc
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/mpi_map.o mpi_map.cc

$(BPATH)/node_access.o: $(EDHEAD) node.h node_access.cc 
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/node_access.o node_access.cc

$(BPATH)/node_constr.o: $(EDHEAD) node.h node_constr.cc 
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/node_constr.o node_constr.cc

$(BPATH)/node_count.o: $(EDHEAD) node.h node_count.cc
	$(CXX) $(CXXFLAGS) -c -o$(BPATH)/node_count.o node_count.cc

$(CPPOPTSF) $(INCDIRSF) $(LDIRSF): Makefile machine_select
	chmod 755 machine_select
	rm -rf $(BPATH)
	mkdir $(BPATH)
	./machine_select

.PHONY: clean

clean:
	rm $(OBJS)
	rm $(EXEC)
	rm $(CPPOPTSF) $(INCDIRSF) $(LDIRSF)
	make -C $(LOCINC) -f Makefile clean

