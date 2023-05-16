make clean
rm *Dict.cxx *.pcm
rootcint -f LoadInputDict.cxx LoadInput.h LoadInputLinkDef.h
rootcint -f OROutputDict.cxx OROutput.h OROutputLinkDef.h
make && ./online_reconstruction digit_012525_R007.root basic_geometry.txt output.root 1
