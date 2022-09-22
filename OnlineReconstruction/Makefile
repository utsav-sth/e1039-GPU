ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)

CXX = nvcc
CXXFLAGS = -g --cudart=static -ccbin g++
DEFS = -DUSE_CUBLAS
LD = $(CXX)
LDFLAGS = 
LIBS = -L/home/fermipdoc/Products/test_build_6.18/lib -L/home/fermipdoc/work/Gpufit-build/Gpufit -lGpufit -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lm -ldl -lcublas -lcurand

PROGRAM = online_reconstruction

OBJS = online_reconstruction.o LoadInput.o LoadInputDict.o

COMPILE = $(CXX) $(CXXFLAGS) -std=c++14 -m64 $(DEFS) \
	-I/home/fermipdoc/Products/test_build_6.18/include -I/home/fermipdoc/Documents/Gpufit/Gpufit -I/home/fermipdoc/work/e1039-GPU/OnlineReconstruction -c


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
