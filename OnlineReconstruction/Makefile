ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)

CXX = nvcc
CXXFLAGS = -g
COMPILE = $(CXX) $(CXXFLAGS) -std=c++11 -m64 -I/home/hjiang/root/include -c
LD = $(CXX)
LDFLAGS = 
LIBS = -L/home/hjiang/root/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lm -ldl

PROGRAM = online_reconstruction

OBJS = online_reconstruction.o LoadInput.o LoadInputDict.o

$(PROGRAM) : $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

online_reconstruction.o : online_reconstruction.cu 
	$(COMPILE) online_reconstruction.cu -o online_reconstruction.o

LoadInput.o : LoadInput.cxx 
	$(COMPILE) LoadInput.cxx -o LoadInput.o

LoadInputDict.o : LoadInputDict.cxx 
	$(COMPILE) LoadInputDict.cxx -o LoadInputDict.o

clean:
	rm -f *.o *~ $(PROGRAM)
