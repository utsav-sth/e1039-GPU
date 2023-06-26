#define nChamberPlanes 30
#define nHodoPlanes 16
#define nPropPlanes 8
#define nDetectors (nChamberPlanes+nHodoPlanes+nPropPlanes+1)
#define RESOLUTION_DC 1.6

void BuildORGPUGeometryDB(const char* surveyfile, const char* aligncham, const char* alignhodo, const char* alignprop, const char* outgeomfile)
{
  //"Member" Arrays declaration:
  // as defined at l.64-119 of GeomSvc.h
  
  //Detector identifier
  int detectorID[nDetectors];
  std::string detectorName[nDetectors];
  int planeType[nDetectors];   //X = 1, angleFromVert > 0 = 2, angleFromVert < 0 = 3, Y = 4

  memset(detectorID, 0, sizeof(int)*nDetectors);
  memset(planeType, 0, sizeof(int)*nDetectors);
  
  //Ideal properties
  int nElements[nDetectors];
  double spacing[nDetectors];
  double cellWidth[nDetectors];
  double xoffset[nDetectors];
  double overlap[nDetectors];
  double angleFromVert[nDetectors];
  double sintheta[nDetectors];
  double costheta[nDetectors];
  double tantheta[nDetectors];
  
  memset(nElements, 0, sizeof(int)*nDetectors);
  memset(spacing, 0, sizeof(double)*nDetectors);
  memset(cellWidth, 0, sizeof(double)*nDetectors);
  memset(xoffset, 0, sizeof(double)*nDetectors);
  memset(overlap, 0, sizeof(double)*nDetectors);
  memset(angleFromVert, 0, sizeof(double)*nDetectors);
  memset(sintheta, 0, sizeof(double)*nDetectors);
  memset(costheta, 0, sizeof(double)*nDetectors);
  memset(tantheta, 0, sizeof(double)*nDetectors);
  
  //Survey info
  double x0[nDetectors];     //x0, y0, z0 define the center of detector
  double y0[nDetectors];
  double z0[nDetectors];
  double x1[nDetectors];     //x1, y1 define the lower/left edge of detector
  double y1[nDetectors];
  double z1[nDetectors];     //z1 is the upstream z position of the detector
  double x2[nDetectors];     //x2, y2 define the upper/right edge of detector
  double y2[nDetectors];
  double z2[nDetectors];     //z2 is the downstream z position of the detector
  double thetaX[nDetectors];
  double thetaY[nDetectors];
  double thetaZ[nDetectors];

  memset(x0, 0, sizeof(double)*nDetectors);
  memset(y0, 0, sizeof(double)*nDetectors);
  memset(z0, 0, sizeof(double)*nDetectors);
  memset(x1, 0, sizeof(double)*nDetectors);
  memset(y1, 0, sizeof(double)*nDetectors);
  memset(z1, 0, sizeof(double)*nDetectors);
  memset(x2, 0, sizeof(double)*nDetectors);
  memset(y2, 0, sizeof(double)*nDetectors);
  memset(z2, 0, sizeof(double)*nDetectors);
  memset(thetaX, 0, sizeof(double)*nDetectors);
  memset(thetaY, 0, sizeof(double)*nDetectors);
  memset(thetaZ, 0, sizeof(double)*nDetectors);

  //Alignment info
  double deltaX[nDetectors];
  double deltaY[nDetectors];
  double deltaZ[nDetectors];
  double deltaW[nDetectors];             //for chambers and hodos
  double deltaW_module[nDetectors][9];   //for prop. tubes only
  double rotX[nDetectors];
  double rotY[nDetectors];
  double rotZ[nDetectors];
  double resolution[nDetectors];

  memset(deltaX, 0, sizeof(double)*nDetectors);
  memset(deltaY, 0, sizeof(double)*nDetectors);
  memset(deltaZ, 0, sizeof(double)*nDetectors);
  memset(deltaW, 0, sizeof(double)*nDetectors);
  memset(deltaW_module, 0, sizeof(double)*nDetectors*9);
  memset(rotX, 0, sizeof(double)*nDetectors);
  memset(rotY, 0, sizeof(double)*nDetectors);
  memset(rotZ, 0, sizeof(double)*nDetectors);
  memset(resolution, 0, sizeof(double)*nDetectors);

  //Final position/rotation
  double xc[nDetectors];
  double yc[nDetectors];
  double zc[nDetectors];
  double wc[nDetectors];
  double rX[nDetectors];
  double rY[nDetectors];
  double rZ[nDetectors];

  memset(xc, 0, sizeof(double)*nDetectors);
  memset(yc, 0, sizeof(double)*nDetectors);
  memset(zc, 0, sizeof(double)*nDetectors);
  memset(wc, 0, sizeof(double)*nDetectors);
  memset(rX, 0, sizeof(double)*nDetectors);
  memset(rY, 0, sizeof(double)*nDetectors);
  memset(rZ, 0, sizeof(double)*nDetectors);
  
  double planeWidth[nDetectors];
  double planeHeight[nDetectors];
  
  memset(planeWidth, 0, sizeof(double)*nDetectors);
  memset(planeHeight, 0, sizeof(double)*nDetectors);
  
  // //Geometric setup
  // TVectorD nVec[nDetectors];// = TVectorD(3);             //Perpendicular to plane
  // TVectorD uVec[nDetectors];// = TVectorD(3);             //measuring direction
  // TVectorD vVec[nDetectors];// = TVectorD(3);             //wire direction
  // TVectorD xVec[nDetectors];// = TVectorD(3);             //X direction
  // TVectorD yVec[nDetectors];// = TVectorD(3);             //Y direction
  // TMatrixD rotM[nDetectors];// = TMatrixD(3, 3);          //rotation matrix

  // for(int k = 1; k<nDetectors; k++){
  //   nVec[k] = TVectorD(3); 
  //   uVec[k] = TVectorD(3); 
  //   vVec[k] = TVectorD(3); 
  //   xVec[k] = TVectorD(3); 
  //   yVec[k] = TVectorD(3); 
  //   rotM[k] = TMatrixD(3, 3);
  // }

  // opening and reading geom scheme
  // corresponds to GeomSvc::initPlaneDbSvc()
  // excpet with a text file for the moment.
  
  fstream _geom_schema;
  _geom_schema.open(surveyfile);
  
  int detectorID_;
  int stationID_;
  int componentID_;
  string detectorName_;
  int myPrimeIs_;
  string geantName_;
  double spacing_;
  double cellWidth_;
  double overlap_;
  int nElements_;
  int lowElementID_;
  double angleFromVert_;
  double xoffset_;
  double planeWidth_;
  double planeHeight_;
  double x0_;
  double y0_;
  double z0_;
  double thetaX_;
  double thetaY_;
  double thetaZ_;
  double resolution_;
  
  char buf[500];
  for(int k = 0; k <= nChamberPlanes+nHodoPlanes; k++){
    _geom_schema.getline(buf, 500);
    if (buf[0] == '#') continue;
    istringstream stringBuf(buf);
    
    stringBuf >> detectorID_ >> stationID_ >> componentID_ >> detectorName_ >> myPrimeIs_
	      >> geantName_ >> spacing_ >> cellWidth_ >> overlap_
	      >> nElements_ >> lowElementID_ >> angleFromVert_ >> xoffset_
	      >> planeWidth_ >> planeHeight_ >> x0_ >> y0_ >> z0_
	      >> thetaX_ >> thetaY_ >> thetaZ_ >> resolution_;

    detectorID[detectorID_] = detectorID_;
    //stationID[detectorID_] = stationID_;
    //componentID[detectorID_] = componentID_;
    detectorName[detectorID_] = detectorName_;
    //myPrimeIs[detectorID_] = myPrimeIs_;
    //geantName[detectorID_] = geantName_;
    spacing[detectorID_] = spacing_;
    cellWidth[detectorID_] = cellWidth_;
    overlap[detectorID_] = overlap_;
    nElements[detectorID_] = nElements_;
    //lowElementID[detectorID_] = lowElementID_;
    angleFromVert[detectorID_] = angleFromVert_;
    xoffset[detectorID_] = xoffset_;
    planeWidth[detectorID_] = planeWidth_;
    planeHeight[detectorID_] = planeHeight_;
    x0[detectorID_] = x0_;
    y0[detectorID_] = y0_;
    z0[detectorID_] = z0_;
    thetaX[detectorID_] = thetaX_;
    thetaY[detectorID_] = thetaY_;
    thetaZ[detectorID_] = thetaZ_;
    resolution[detectorID_] = resolution_;

    if(detectorID_<=30)
      cout << detectorName_ << "\t" << nElements_ << "\t" << spacing_ << "\t" << cellWidth_ << "\t"
	   << angleFromVert_ << "\t" << xoffset_ << "\t" << planeWidth_ << "\t" << planeHeight_ << "\t"
	   << x0_ << "\t" << y0_ << "\t" << z0_ << "\t" << thetaX_ << "\t" << thetaY_ << "\t" << thetaZ_ << endl;
    
    x1[detectorID_] = x0[detectorID_]-planeWidth[detectorID_]*0.5;
    x2[detectorID_] = x0[detectorID_]+planeWidth[detectorID_]*0.5;
    y1[detectorID_] = y0[detectorID_]-planeHeight[detectorID_]*0.5;
    y2[detectorID_] = y0[detectorID_]+planeHeight[detectorID_]*0.5;
    
    //update

    xc[detectorID_] = x0[detectorID_]+deltaX[detectorID_];
    yc[detectorID_] = y0[detectorID_]+deltaY[detectorID_];
    zc[detectorID_] = z0[detectorID_]+deltaZ[detectorID_];

    rX[detectorID_] = thetaX[detectorID_] + rotX[detectorID_];
    rY[detectorID_] = thetaY[detectorID_] + rotY[detectorID_];
    rZ[detectorID_] = thetaZ[detectorID_] + rotZ[detectorID_];

    sintheta[detectorID_] = sin(angleFromVert[detectorID_] + rZ[detectorID_]);
    costheta[detectorID_] = cos(angleFromVert[detectorID_] + rZ[detectorID_]);
    tantheta[detectorID_] = tan(angleFromVert[detectorID_] + rZ[detectorID_]);

    wc[detectorID_] = xc[detectorID_]*costheta[detectorID_] + yc[detectorID_]*sintheta[detectorID_];

    bool isprime = false;
    if(detectorName[detectorID_].find("X") != string::npos || detectorName[detectorID_].find("T") != string::npos || detectorName[detectorID_].find("B") != string::npos)
      {
	planeType[detectorID_] = 1;
      }
    else if((detectorName[detectorID_].find("U") != string::npos || detectorName[detectorID_].find("V") != string::npos) && angleFromVert[detectorID_] > 0)
      {
	planeType[detectorID_] = 2;
      }
    else if((detectorName[detectorID_].find("U") != string::npos || detectorName[detectorID_].find("V") != string::npos) && angleFromVert[detectorID_] < 0)
      {
	planeType[detectorID_] = 3;
      }
    else if(detectorName[detectorID_].find("Y") != string::npos || detectorName[detectorID_].find("L") != string::npos || detectorName[detectorID_].find("R") != string::npos)
      {
	planeType[detectorID_] = 4;
      }
    if(detectorName[detectorID_].find("p") != string::npos)isprime = true;
    
    //cout << "detector " << detectorID[detectorID_] << " type " << planeType[detectorID_] << " isprime " << isprime << endl;
  }

  cout << "loaded survey parameters" << endl;
  
  fstream _align_mille;
  _align_mille.open(aligncham, ios::in);
  
  for(int k = 1; k <= nChamberPlanes; k++){
    _align_mille.getline(buf, 100);
    istringstream stringBuf(buf);
    
    stringBuf >> deltaZ[k] >> rotZ[k] >> deltaW[k] >> resolution[k];
    resolution[k] = 0.04;
    deltaX[k] = deltaW[k]*costheta[k];
    deltaY[k] = deltaW[k]*sintheta[k];

    xc[k] = x0[k]+deltaX[k];
    yc[k] = y0[k]+deltaY[k];
    zc[k] = z0[k]+deltaZ[k];

    rX[k] = thetaX[k] + rotX[k];
    rY[k] = thetaY[k] + rotY[k];
    rZ[k] = thetaZ[k] + rotZ[k];

    sintheta[k] = sin(angleFromVert[k] + rZ[k]);
    costheta[k] = cos(angleFromVert[k] + rZ[k]);
    tantheta[k] = tan(angleFromVert[k] + rZ[k]);

    wc[k] = xc[k]*costheta[k] + yc[k]*sintheta[k];
    
    //if(k<10)cout << x0[k] << " "<< y0[k] << " "<< z0[k] << " " << xc[k] << " "<< yc[k] << " "<< zc[k] << " " << endl;
    /*
    rotM[k][0][0] = cos(rZ[k])*cos(rY[k]);
    rotM[k][0][1] = cos(rZ[k])*sin(rX[k])*sin(rY[k]) - cos(rX[k])*sin(rZ[k]);
    rotM[k][0][2] = cos(rX[k])*cos(rZ[k])*sin(rY[k]) + sin(rX[k])*sin(rZ[k]);
    rotM[k][1][0] = sin(rZ[k])*cos(rY[k]);
    rotM[k][1][1] = sin(rZ[k])*sin(rX[k])*sin(rY[k]) + cos(rZ[k])*cos(rX[k]);
    rotM[k][1][2] = sin(rZ[k])*sin(rY[k])*cos(rX[k]) - cos(rZ[k])*sin(rX[k]);
    rotM[k][2][0] = -sin(rY[k]);
    rotM[k][2][1] = cos(rY[k])*sin(rX[k]);
    rotM[k][2][2] = cos(rY[k])*cos(rX[k]);
    
    uVec[k][0] = cos(angleFromVert[k]);
    uVec[k][1] = sin(angleFromVert[k]);
    uVec[k][2] = 0.;
    uVec[k] = rotM[k]*uVec[k];

    vVec[k][0] = -sin(angleFromVert[k]);
    vVec[k][1] = cos(angleFromVert[k]);
    vVec[k][2] = 0.;
    vVec[k] = rotM[k]*vVec[k];

    xVec[k][0] = 1.;
    xVec[k][1] = 0.;
    xVec[k][2] = 0.;
    xVec[k] = rotM[k]*xVec[k];

    yVec[k][0] = 0.;
    yVec[k][1] = 1.;
    yVec[k][2] = 0.;
    yVec[k] = rotM[k]*yVec[k];

    nVec[k][0] = uVec[k][1]*vVec[k][2] - vVec[k][1]*uVec[k][2];
    nVec[k][1] = uVec[k][2]*vVec[k][0] - vVec[k][2]*uVec[k][0];
    nVec[k][2] = uVec[k][0]*vVec[k][1] - vVec[k][0]*uVec[k][1];
    */
  }

  for(int k = 1; k <= nChamberPlanes; k+=2){
    resolution[k] = RESOLUTION_DC*0.5*(resolution[k]+resolution[k+1]);
    resolution[k+1] = resolution[k];
  }

  fstream _align_hodo;
  _align_hodo.open(alignhodo, ios::in);
  for(int k = nChamberPlanes+1; k <= nChamberPlanes+nHodoPlanes; k++){
    _align_mille.getline(buf, 100);
    istringstream stringBuf(buf);
    
    stringBuf >> deltaW[k];
    
    deltaX[k] = deltaW[k]*costheta[k];
    deltaY[k] = deltaW[k]*sintheta[k];

    xc[k] = x0[k]+deltaX[k];
    yc[k] = y0[k]+deltaY[k];
    zc[k] = z0[k]+deltaZ[k];

    rX[k] = thetaX[k] + rotX[k];
    rY[k] = thetaY[k] + rotY[k];
    rZ[k] = thetaZ[k] + rotZ[k];

    sintheta[k] = sin(angleFromVert[k] + rZ[k]);
    costheta[k] = cos(angleFromVert[k] + rZ[k]);
    tantheta[k] = tan(angleFromVert[k] + rZ[k]);

    wc[k] = xc[k]*costheta[k] + yc[k]*sintheta[k];
 }
  
  
  fstream _align_prop;
  _align_prop.open(alignprop, ios::in);

  
  
  cout << "loaded alignment" << endl;

  float p1x[nDetectors];
  float p1y[nDetectors];
  float p1z[nDetectors];
  float dp1x[nDetectors];
  float dp1y[nDetectors];
  float dp1z[nDetectors];
  float deltapx[nDetectors];
  float deltapy[nDetectors];
  float deltapz[nDetectors];

  //calculate
  for(int k = 1; k<nDetectors; k++){
    TVectorD detCenter(3);
    detCenter[0] = xc[k];
    detCenter[1] = yc[k];
    detCenter[2] = zc[k];

    TVectorD nVec = TVectorD(3);             //Perpendicular to plane
    TVectorD uVec = TVectorD(3);             //measuring direction
    TVectorD vVec = TVectorD(3);             //wire direction
    TVectorD xVec = TVectorD(3);             //X direction
    TVectorD yVec = TVectorD(3);             //Y direction
    TMatrixD rotM = TMatrixD(3, 3);          //rotation matrix
    
    rotM[0][0] = cos(rZ[k])*cos(rY[k]);
    rotM[0][1] = cos(rZ[k])*sin(rX[k])*sin(rY[k]) - cos(rX[k])*sin(rZ[k]);
    rotM[0][2] = cos(rX[k])*cos(rZ[k])*sin(rY[k]) + sin(rX[k])*sin(rZ[k]);
    rotM[1][0] = sin(rZ[k])*cos(rY[k]);
    rotM[1][1] = sin(rZ[k])*sin(rX[k])*sin(rY[k]) + cos(rZ[k])*cos(rX[k]);
    rotM[1][2] = sin(rZ[k])*sin(rY[k])*cos(rX[k]) - cos(rZ[k])*sin(rX[k]);
    rotM[2][0] = -sin(rY[k]);
    rotM[2][1] = cos(rY[k])*sin(rX[k]);
    rotM[2][2] = cos(rY[k])*cos(rX[k]);
    
    uVec[0] = cos(angleFromVert[k]);
    uVec[1] = sin(angleFromVert[k]);
    uVec[2] = 0.;
    uVec = rotM*uVec;

    vVec[0] = -sin(angleFromVert[k]);
    vVec[1] = cos(angleFromVert[k]);
    vVec[2] = 0.;
    vVec = rotM*vVec;

    xVec[0] = 1.;
    xVec[1] = 0.;
    xVec[2] = 0.;
    xVec = rotM*xVec;

    yVec[0] = 0.;
    yVec[1] = 1.;
    yVec[2] = 0.;
    yVec = rotM*yVec;

    nVec[0] = uVec[1]*vVec[2] - vVec[1]*uVec[2];
    nVec[1] = uVec[2]*vVec[0] - vVec[2]*uVec[0];
    nVec[2] = uVec[0]*vVec[1] - vVec[0]*uVec[1];

    // cout << nElements[k] << endl;
    // for(int l = 1; l<=nElements[k]; l++){
    //   TVectorD ep1(3), ep2(3);
    //   if(planeType[k] != 4) //special treatment for Y planes
    //  	{
    //  	  double cellLength = fabs(y2[k] - y1[k])/cos(angleFromVert[k]);
    // 	  double hspacing = spacing[k]/cos(angleFromVert[k]);
    // 	  double hoffset = xoffset[k]/cos(angleFromVert[k]);
    //  	  short sign = -1;
    //  	  ep1 = detCenter + ((l - (nElements[k]+1.)/2.)*hspacing + hoffset)*xVec + 0.5*sign*cellLength*vVec;
    //  	  sign = +1;
    //  	  ep2 = detCenter + ((l - (nElements[k]+1.)/2.)*hspacing + hoffset)*xVec + 0.5*sign*cellLength*vVec;
    //  	}
    // }

    int i_element = 1;
    short sign = -1;
    double cellLength = fabs(y2[k] - y1[k])/cos(angleFromVert[k]);
    double hspacing = spacing[k]/cos(angleFromVert[k]);
    double hoffset = xoffset[k]/cos(angleFromVert[k]);

    TVectorD ep1_w1(3), ep2_w1(3);
    sign = -1;
    ep1_w1 = detCenter + ((i_element - (nElements[k]+1.)/2.)*hspacing + hoffset)*xVec + 0.5*sign*cellLength*vVec;
    sign = +1;
    ep2_w1 = detCenter + ((i_element - (nElements[k]+1.)/2.)*hspacing + hoffset)*xVec + 0.5*sign*cellLength*vVec;
    
    TVectorD ep1_wN(3), ep2_wN(3);
    i_element = nElements[k];
    sign = -1;
    ep1_wN = detCenter + ((i_element - (nElements[k]+1.)/2.)*hspacing + hoffset)*xVec + 0.5*sign*cellLength*vVec;
    sign = +1;
    ep2_wN = detCenter + ((i_element - (nElements[k]+1.)/2.)*hspacing + hoffset)*xVec + 0.5*sign*cellLength*vVec;

    //if(detectorID[k]<3)
      //cout << xc[k] << " "<< yc[k] << " "<< zc[k] << " " << cellLength << " " << hspacing << " " << hoffset << endl;
      
    if(planeType[k]==4){
      cellLength = fabs(x2[k] - x1[k])/sin(angleFromVert[k]);
      hspacing = spacing[k]/sin(angleFromVert[k]);
      hoffset = xoffset[k]/sin(angleFromVert[k]);

      sign = -1;
      ep1_w1 = detCenter + ((i_element - (nElements[k]+1.)/2.)*hspacing + hoffset)*yVec + 0.5*sign*cellLength*vVec;
      sign = +1;
      ep2_w1 = detCenter + ((i_element - (nElements[k]+1.)/2.)*hspacing + hoffset)*yVec + 0.5*sign*cellLength*vVec;
      
      i_element = nElements[k];
      sign = -1;
      ep1_wN = detCenter + ((i_element - (nElements[k]+1.)/2.)*hspacing + hoffset)*yVec + 0.5*sign*cellLength*vVec;
      sign = +1;
      ep2_wN = detCenter + ((i_element - (nElements[k]+1.)/2.)*hspacing + hoffset)*yVec + 0.5*sign*cellLength*vVec;
    }
    p1x[k] = ep1_w1[0];
    p1y[k] = ep1_w1[1];
    p1z[k] = ep1_w1[2];
    
    dp1x[k] = (ep1_wN[0]-ep1_w1[0])/(nElements[k]-1);
    dp1y[k] = (ep1_wN[1]-ep1_w1[1])/(nElements[k]-1);
    dp1z[k] = (ep1_wN[2]-ep1_w1[2])/(nElements[k]-1);

    if(k<=30)cout <<std::setprecision(10) << ep1_w1[0] << " " << ep1_w1[1] << " " << ep1_w1[2] << "; " << ep1_wN[0] << " " << ep1_wN[1] << " " << ep1_wN[2] << "; " << dp1x[k] << " "<< dp1y[k] << " "<< dp1z[k] << endl;
    
    // deltapx[k] = ep2_w1[0]-ep1_w1[0];
    // deltapy[k] = ep2_w1[1]-ep1_w1[1];
    // deltapz[k] = ep2_w1[2]-ep1_w1[2];

    //if(k<=6)cout << ep2_w1[0]-ep1_w1[0] << " " << ep2_w1[1]-ep1_w1[1] << " " << ep2_w1[2]-ep1_w1[2] << " " << sqrt( (ep2_w1[0]-ep1_w1[0])*(ep2_w1[0]-ep1_w1[0]) + (ep2_w1[1]-ep1_w1[1])*(ep2_w1[1]-ep1_w1[1]) + (ep2_w1[2]-ep1_w1[2])*(ep2_w1[2]-ep1_w1[2]) ) << endl;
    
    deltapx[k] = 0;
    deltapy[k] = 0;
    deltapz[k] = 0;

    for(i_element = 1; i_element<=nElements[k]; i_element++){
      TVectorD ep1(3), ep2(3);
      if(planeType[k] != 4){ //special treatment for Y planes
      	
	sign = -1;
	ep1 = detCenter + ((i_element - (nElements[k]+1.)/2.)*hspacing + hoffset)*xVec + 0.5*sign*cellLength*vVec;
	sign = +1;
	ep2 = detCenter + ((i_element - (nElements[k]+1.)/2.)*hspacing + hoffset)*xVec + 0.5*sign*cellLength*vVec;
      }else{
	sign = -1;
	ep1_w1 = detCenter + ((i_element - (nElements[k]+1.)/2.)*hspacing + hoffset)*yVec + 0.5*sign*cellLength*vVec;
	sign = +1;
	ep2_w1 = detCenter + ((i_element - (nElements[k]+1.)/2.)*hspacing + hoffset)*yVec + 0.5*sign*cellLength*vVec;
      }
      deltapx[k]+= ep2[0]-ep1[0];
      deltapy[k]+= ep2[1]-ep1[1];
      deltapz[k]+= ep2[2]-ep1[2];

    }
    
    deltapx[k]/= nElements[k];
    deltapy[k]/= nElements[k];
    deltapz[k]/= nElements[k];

    //if(k<=6)cout << deltapx[k] << " " << deltapy[k] << " " << deltapz[k] << " " << sqrt( deltapx[k]*deltapx[k] + deltapy[k]*deltapy[k] + deltapz[k]*deltapz[k] ) << endl; 
  }

  ofstream out(outgeomfile);
  
  out << " ##  z_plane   N_elem  CellWidth   Spacing xoffset scalex  x_0   x1       x2        costheta        scaley      y_0        y1     y2          sintheta     resolution    p1x     p1y     p1z      deltapx deltapy deltapz  dp1x      dp1y        dp1z           deltaW " << endl;
  out << "## chambers" << endl;
  for(int k = 1; k<=nChamberPlanes; k++){
    out << std::setprecision(8) << k << "\t" << zc[k] << "\t" << nElements[k] << "\t"
	 << cellWidth[k] << "\t" << spacing[k]  << "\t"
	 << xoffset[k] << "\t" << planeWidth[k] << "\t"
	 << x0[k] << "\t" << x1[k] << "\t" << x2[k] << "\t"
	 << costheta[k] << "\t" << planeHeight[k] << "\t"
	 << y0[k] << "\t" << y1[k] << "\t" << y2[k] << "\t"
	 << sintheta[k] << "\t" << resolution[k] << "\t"
	 << p1x[k] << "\t" << p1y[k] << "\t" << p1z[k] << "\t"  
	 << deltapx[k] << "\t" << deltapy[k] << "\t" << deltapz[k] << "\t"  
	 << dp1x[k] << "\t" << dp1y[k] << "\t" << dp1z[k] << "\t" << deltaW[k] << endl;
  }
  out << "## hodoscopes" << endl;
  for(int k = nChamberPlanes+1; k<=nChamberPlanes+nHodoPlanes; k++){
    out << std::setprecision(8) << k << "\t" << zc[k] << "\t" << nElements[k] << "\t"
	<< cellWidth[k] << "\t" << spacing[k]  << "\t"
	<< xoffset[k] << "\t" << planeWidth[k] << "\t"
	<< x0[k] << "\t" << x1[k] << "\t" << x2[k] << "\t"
	<< costheta[k] << "\t" << planeHeight[k] << "\t"
	<< y0[k] << "\t" << y1[k] << "\t" << y2[k] << "\t"
	<< sintheta[k] << "\t" << resolution[k] << "\t"
	<< p1x[k] << "\t" << p1y[k] << "\t" << p1z[k] << "\t"  
	<< deltapx[k] << "\t" << deltapy[k] << "\t" << deltapz[k] << "\t"  
	<< dp1x[k] << "\t" << dp1y[k] << "\t" << dp1z[k] << "\t" << deltaW[k] << endl;
  }  
  // out << "## prop tubes" << endl;
  // for(int k = nChamberPlanes+nHodoPlanes+1; k<=nChamberPlanes+nHodoPlanes+nPropPlanes; k++){
  //   out << std::setprecision(8) << k << "\t" << zc[k] << "\t" << nElements[k] << "\t"
  // 	 << cellWidth[k] << "\t" << spacing[k]  << "\t"
  // 	 << xoffset[k] << "\t" << planeWidth[k] << "\t"
  // 	 << x0[k] << "\t" << x1[k] << "\t" << x2[k] << "\t"
  // 	 << costheta[k] << "\t" << planeHeight[k] << "\t"
  // 	 << y0[k] << "\t" << y1[k] << "\t" << y2[k] << "\t"
  // 	 << sintheta[k] << "\t" << resolution[k] << "\t"
  // 	 << p1x[k] << "\t" << p1y[k] << "\t" << p1z[k] << "\t"  
  // 	 << deltapx[k] << "\t" << deltapy[k] << "\t" << deltapz[k] << "\t"  
  // 	 << dp1x[k] << "\t" << dp1y[k] << "\t" << dp1z[k] << "\t";
  //   for(int l = 0; l<9; l++)out << " " << deltaW_module[k][l];
  //   out << endl;
  // }  
  
}
