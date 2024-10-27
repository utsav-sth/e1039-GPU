ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)

CXX = nvcc
CXXFLAGS = -g --cudart=static -ccbin g++
DEFS = -DUSE_CUBLAS
LD = $(CXX)
LDFLAGS = 
LIBS = -L/home/cayuso/Products/test_build_6.14/lib -L/home/cayuso/Documents/Gpufit-build/Gpufit -lGpufit -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lm -ldl -lcublas -lcurand

PROGRAM = online_reconstruction

OBJS = online_reconstruction.o LoadInput.o LoadInputDict.o

COMPILE = $(CXX) $(CXXFLAGS) -std=c++14 -m64 $(DEFS) \
	-I/home/cayuso/Products/test_build_6.14/include -I/home/cayuso/Documents/Gpufit/Gpufit -c

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
