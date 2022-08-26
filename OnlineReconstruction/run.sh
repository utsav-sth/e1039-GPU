make clean
rm LoadInputDict.cxx *.pcm
rootcint -f LoadInputDict.cxx LoadInput.h LoadInputLinkDef.h
make && ./online_reconstruction digit_028705_009.root output_ex.txt
