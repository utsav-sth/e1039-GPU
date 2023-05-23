make clean
rm *Dict.cxx *.pcm
rootcint -f LoadInputDict.cxx LoadInput.h LoadInputLinkDef.h
rootcint -f OROutputDict.cxx OROutput.h OROutputLinkDef.h
make && ./online_reconstruction DST.root basic_geometry.txt output_mc.root 0
