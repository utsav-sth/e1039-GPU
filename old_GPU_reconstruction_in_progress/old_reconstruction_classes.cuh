class gFitParams {
public:
      int max_iterations;
      short nparam;
      
      float parameter_limits_min[5];
      float parameter_limits_max[5];
      
      float tolerance;
};


class gFitArrays {
public:
      int npoints;
      float drift_dist[nChamberPlanes]; // hit drift distance
      float resolution[nChamberPlanes]; // detector resolution
      
      float p1x[nChamberPlanes];// x bottom end point of the wire hit 
      float p1y[nChamberPlanes];// y bottom end point of the wire hit 
      float p1z[nChamberPlanes];// z bottom end point of the wire hit 
      
      float deltapx[nChamberPlanes];// x distance between bottom and top end points of the wire hit 
      float deltapy[nChamberPlanes];// y distance between bottom and top end points of the wire hit 
      float deltapz[nChamberPlanes];// z distance between bottom and top end points of the wire hit 
      
      float prev_parameters[5];
      float output_parameters[5];
      float output_parameters_errors[5];
      float chi2prev;
      float chi2;

      float values[nChamberPlanes];
      float derivatives[5*nChamberPlanes];
      float gradients[5];
      float hessians[25];
      
      float scaling_vector[5];
      float deltas[5];
      
      float lambda;//scale factor
      
      //float calc_matrix[25];
      //float abs_row[25];
      //int abs_row_index[25];

      int n_iter;
      short iter_failed;
      short finished;
      short state;
      short skip;
      short singular;
      
      //float x_array[nChamberPlanes];// x position arrays
      //float y_array[nChamberPlanes];// y position arrays
      //float z_array[nChamberPlanes];// z position arrays
      //float dx_array[nChamberPlanes];// x position uncertainty
      //float dy_array[nChamberPlanes];// x position uncertainty

      float A[25];// matrix: max size 5x5, but we can use this unique array for all possible sizes
      float Ainv[25];// matrix
      float B[5];// input vector

      //float output_parameters_steps[5];
      //float doublederivatives[5];
};


//Output class
class gSW {
public:
	gSW(){
		EventID = -1;
		nAH = 0;
		nTracklets = 0;
	}
	thrust::pair<int, int> hitpairs_x[100];
	thrust::pair<int, int> hitpairs_u[100];
	thrust::pair<int, int> hitpairs_v[100];
	int hitidx1[100];
	int hitidx2[100];
	short hitflag1[100];
	short hitflag2[100];
	int EventID;
	int nAH;
	int nTracklets;
	gTracklet AllTracklets[TrackletSizeMax];
	short nTKL_stID[7];//0: D0; 1: D1; 2: D2; 3: D3p; 4: D3m; 5: back partial; 6: global
};

