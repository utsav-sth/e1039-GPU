void BuildORGPUGeometryDB(const char* surveyfile, const char* alignfile)
{
  
  
  fstream _align_mille;
  _align_mille.open(alignfile.c_str(), ios::in);

  float deltaZ, rotZ, deltaW, resolution;
  
  for(int i = 1; i <= nChamberPlanes; i++){
    _align_mille.getline(buf, 100);
    istringstream stringBuf(buf);
    
    stringBuf >> deltaZ >> rotZ >> deltaW >> resolution;
    
    // planes[i].deltaX = planes[i].deltaW*planes[i].costheta;
    // planes[i].deltaY = planes[i].deltaW*planes[i].sintheta;
    // planes[i].update();
    //if(planes[i].resolution < RESOLUTION_DC) planes[i].resolution = RESOLUTION_DC;
  }
  
  
  
}
