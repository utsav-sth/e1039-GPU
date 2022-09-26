ROOTCFLAGS = $(shell root-config --cflags)
ROOTLIBS   = $(shell root-config --libs)

CXX = nvcc
CXXFLAGS = -g --cudart=static -ccbin g++
DEFS = -DUSE_CUBLAS
LD = $(CXX)
LDFLAGS = 
#LIBS = -L/home/fermipdoc/Products/test_build_6.18/lib -L/home/fermipdoc/work/Gpufit-build/Gpufit -lGpufit -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lm -ldl -lcublas -lcurand
LIBS = -L/data2/e1039/share/root/lib -L/home/programs/Gpufit-build/Gpufit -L/home/programs/core-inst/lib -lGpufit -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -lm -ldl -lcublas -lcurand -lktracker_interface -lkfitter -lgenfit2 -lsqgenfit -lphfield -lphgeom -linterface_main -lfun4all -lgeom_svc -lphg4hit -lphool

PROGRAM = online_reconstruction

OBJS = online_reconstruction.o #LoadInput.o LoadInputDict.o

#COMPILE = $(CXX) $(CXXFLAGS) -std=c++14 -m64 $(DEFS) \
	-I/home/fermipdoc/Products/test_build_6.18/include -I/home/fermipdoc/Documents/Gpufit/Gpufit -c
COMPILE = $(CXX) $(CXXFLAGS) -std=c++14 -m64 $(DEFS) -I/home/programs/e1039-GPU/OnlineReconstruction \
        -I/data2/e1039/share/root/include -I/home/programs/Gpufit/Gpufit \
        -I/home/programs/e1039-core/packages/geom_svc \
        -I/home/programs/e1039-core/packages/reco/SQGenFit \
        -I/home/programs/e1039-core/packages/reco/interface \
        -I/home/programs/e1039-core/packages/reco/kfitter \
        -I/home/programs/e1039-core/packages/reco/ktracker \
        -I/home/programs/e1039-core/packages/reco/patterndb -c

$(PROGRAM) : $(OBJS)
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $(PROGRAM)

online_reconstruction.o : online_reconstruction.cu 
	$(COMPILE) online_reconstruction.cu -o online_reconstruction.o

# LoadInput.o : LoadInput.cxx 
# 	$(COMPILE) LoadInput.cxx -o LoadInput.o

# LoadInputDict.o : LoadInputDict.cxx 
# 	$(COMPILE) LoadInputDict.cxx -o LoadInputDict.o

clean:
	rm -f *.o *~ $(PROGRAM)
