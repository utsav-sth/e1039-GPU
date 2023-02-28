ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)

CXX = nvcc
CXXFLAGS = -g --cudart=static -ccbin g++ -rdc=true
#CXXFLAGS = --cudart=static -ccbin g++
DEFS = -DUSE_CUBLAS
LD = $(CXX)
LDFLAGS = 
LIBS = -L${ROOTSYS}/lib -L${GPUFIT}/Gpufit-build/Gpufit -lGpufit -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lm -ldl -lcublas -lcurand
# -linterface_main 
#-lphool -lfun4all -lg4detectors -lktracker

PROGRAM = online_reconstruction

OBJS = online_reconstruction.o LoadInput.o LoadInputDict.o OROutput.o OROutputDict.o

COMPILE = $(CXX) $(CXXFLAGS) -std=c++14 -m64 $(DEFS) \
	-I${ROOTSYS}/include -I${GPUFIT}/Gpufit/Gpufit -I${PWD} -c
#\
#	-I/home/programs/e1039-core/interface_main -c
#	-I/home/programs/e1039-core/simulation/g4detectors \
#	-I/home/programs/e1039-core/packages/reco/ktracker -c
#	-I/home/programs/e1039-core/framework/fun4all \
#	-I/home/programs/e1039-core/framework/phool \

$(PROGRAM) : $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

online_reconstruction.o : online_reconstruction.cu 
	$(COMPILE) online_reconstruction.cu -o online_reconstruction.o

LoadInput.o : LoadInput.cxx 
	$(COMPILE) LoadInput.cxx -o LoadInput.o

LoadInputDict.o : LoadInputDict.cxx 
	$(COMPILE) LoadInputDict.cxx -o LoadInputDict.o

OROutput.o : OROutput.cxx 
	$(COMPILE) OROutput.cxx -o OROutput.o

OROutputDict.o : OROutputDict.cxx 
	$(COMPILE) OROutputDict.cxx -o OROutputDict.o

clean:
	rm -f *.o *~ $(PROGRAM)
