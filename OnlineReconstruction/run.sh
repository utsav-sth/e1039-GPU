make clean
rm LoadInputDict.cxx *.pcm
rootcint -f LoadInputDict.cxx LoadInput.h LoadInputLinkDef.h
rootcint -f OROutputDict.cxx OROutput.h OROutputLinkDef.h
#make && ./online_reconstruction digit_028705_009.root basic_geometry.txt output.root 1
make && ./online_reconstruction DST.root basic_geometry.txt output.root
