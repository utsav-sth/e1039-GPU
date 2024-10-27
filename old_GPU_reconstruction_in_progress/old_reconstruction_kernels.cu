// tracklet in station builder: 
__global__ void gkernel_TrackletinStation(gEvent* ic, gSW* oc, gFitArrays* fitarrays, int stID, const gPlane* planes, gFitParams* fitparams) {
	//int blockId = blockIdx.x + blockIdx.y * gridDim.x + gridDim.x * gridDim.y * blockIdx.z;
	//int index = threadIdx.x + blockId * blockDim.x;
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	
	stID--;
	oc[index].EventID = ic[index].EventID;
	oc[index].nAH = ic[index].nAH;
	
	//if(10000<ic[index].EventID && ic[index].EventID<10050){
	//	for(int m = 0; m<30; m++){
	//		if(planes[m].u_win!=0)printf("index= %d, m = %d, u_win = %1.6f, costheta = %1.6f\n", index, m, planes[m].u_win, planes[m].costheta);
	//	}
	//}
	// loop on hits
	//if( (ic[index].EventID)>10000 && (ic[index].EventID)<10100 ){//just look at a subset with something in it
	//	printf("core idx %d, evt %d: reduced AllHits value : %d\n", (index), ic[index].EventID, (Nhits));
	//}
	// answer is yes, we still have the info from the previous function i.e. running after eR we still benefit from hit reduction;
	// was worth checking, just in case...

	//we don't need pairs of *HITS* necessarily, we just need pairs of indices...
	//thrust::pair<int, int> hitpairs_x[100];
	//thrust::pair<int, int> hitpairs_u[100];
	//thrust::pair<int, int> hitpairs_v[100];

	int nx = make_hitpairs_in_station(ic, oc[index].hitpairs_x, oc[index].hitidx1, oc[index].hitidx2, oc[index].hitflag1, oc[index].hitflag2, stID, 0, planes);
	int nu = make_hitpairs_in_station(ic, oc[index].hitpairs_u, oc[index].hitidx1, oc[index].hitidx2, oc[index].hitflag1, oc[index].hitflag2, stID, 1, planes);
	int nv = make_hitpairs_in_station(ic, oc[index].hitpairs_v, oc[index].hitidx1, oc[index].hitidx2, oc[index].hitflag1, oc[index].hitflag2, stID, 2, planes);
	
	short uidx = stID==0? stID*6+1 : stID*6+5;
	short vidx = stID==0? stID*6+5 : stID*6+1;
	
	//bool print = false;
	//if(0<=ic[index].EventID && ic[index].EventID<20){
	//	print = true;
	//printf("evt %d, nx = %d, nu = %d, nv = %d, ucostheta(plane %d) = %1.6f, uwin(plane %d) = %1.6f\n", ic[index].EventID, nx, nu, nv, uidx, planes[uidx].costheta, uidx, planes[uidx].u_win);
	//	printf("evt %d, nx = %d, nu = %d, nv = %d\n", ic[index].EventID, nx, nu, nv);
	//}
	//one has to have at least one hit in x, u, v
	if(nx==0 || nu==0 || nv==0)return;
	int n_tkl = oc[index].nTracklets;

	int nhits_tkl = 0;
	int npts = 0;

	REAL fixedpoint[3] = {0, 0, 0};
		
	//X-U combinations first
	for(int i = 0; i< nx; i++){
		double xpos = oc[index].hitpairs_x[i].second>=0 ? 0.5*(ic[index].AllHits[ oc[index].hitpairs_x[i].first ].pos+ic[index].AllHits[ oc[index].hitpairs_x[i].second ].pos): ic[index].AllHits[ oc[index].hitpairs_x[i].first ].pos;
		double umin = xpos*planes[uidx].costheta-planes[uidx].u_win;
		double umax = umin+2*planes[uidx].u_win;
		//if(print){
		//	printf("evt %d, xpos = %1.6f, umin = %1.6f, umax = %1.6f\n", ic[index].EventID, xpos, umin, umax);
		//	printf("evt %d, x1 pos = %1.6f, x2 pos =%1.6f\n", ic[index].EventID, ic[index].AllHits[ oc[index].hitpairs_x[i].first ].pos, 
		//		oc[index].hitpairs_x[i].second >=0 ? ic[index].AllHits[ oc[index].hitpairs_x[i].second ].pos : -1000000);
		//}
		for(int j = 0; j< nu; j++){
			double upos = oc[index].hitpairs_u[j].second>=0 ? 0.5*(ic[index].AllHits[ oc[index].hitpairs_u[j].first ].pos+ic[index].AllHits[ oc[index].hitpairs_u[j].second ].pos): ic[index].AllHits[ oc[index].hitpairs_u[j].first ].pos;

			if(upos<umin || upos>umax)continue;
			//if(print)printf("evt %d, %1.6f <? upos = %1.6f <? %1.6f \n", ic[index].EventID, umin, upos, umax);
			//we chose by convention to start the numbering of the "planes" object arrays to 0 instead of 1, this is why we have to subtract 1 to detectorID to get the correct information
			double z_x = oc[index].hitpairs_x[i].second>=0 ? planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID ].z_mean : planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID ].z;
			double z_u = oc[index].hitpairs_u[j].second>=0 ? planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID ].z_mean : planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID ].z;
			double z_v = planes[vidx].z_mean;
			//if(ic[index].EventID==0)printf("detid x = %d, detid u = %d, z_x = %1.6f, z_u = %1.6f, z_v = %1.6f\n", ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID, ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID, z_x, z_u, z_v);
			double v_win1 = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID ].v_win_fac1;
			double v_win2 = fabs(z_u+z_v-2*z_x)*planes[ vidx ].v_win_fac2;
			double v_win3 = fabs(z_v-z_u)*planes[ vidx ].v_win_fac3;
			double v_win = v_win1+v_win2+v_win3+2*planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID ].spacing;

			double vmin = 2*xpos*planes[uidx].costheta-upos-v_win;
			double vmax = vmin+2*v_win;
			//if(ic[index].EventID==0)printf("vmin = %1.6f, vmax = %1.6f, vwin = %1.6f, vwin1 = %1.6f, vwin2 = %1.6f, vwin3 = %1.6f\n", vmin, vmax, v_win, v_win1, v_win2, v_win3);
			for(int k = 0; k< nv; k++){
				double vpos = oc[index].hitpairs_v[k].second>=0 ? 0.5*(ic[index].AllHits[ oc[index].hitpairs_v[k].first ].pos+ic[index].AllHits[ oc[index].hitpairs_v[k].second ].pos): ic[index].AllHits[ oc[index].hitpairs_v[k].first ].pos;
				//if(ic[index].EventID<20)printf("evt %d: vmin = %1.6f <? vpos = %1.6f <? vmax = %1.6f\n", ic[index].EventID, vmin, vpos, vmax);
				if(vpos<vmin || vpos>vmax)continue;

				oc[index].AllTracklets[n_tkl].stationID = stID;
				nhits_tkl = 0;
				npts = 0;

				if(oc[index].hitpairs_x[i].first>=0){
					if(ic[index].EventID==0)printf("x first hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID, ic[index].AllHits[ oc[index].hitpairs_x[i].first ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_x[i].first ];
					oc[index].AllTracklets[n_tkl].nXHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;
					
					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_x[i].first ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_x[i].first ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_x[i].first ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_x[i].first ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].first ].detectorID ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_x[i].second>=0){
					if(ic[index].EventID==0)printf("x second hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID, ic[index].AllHits[ oc[index].hitpairs_x[i].second ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_x[i].second ];
					oc[index].AllTracklets[n_tkl].nXHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;

					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_x[i].second ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_x[i].second ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_x[i].second ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_x[i].second ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_x[i].second ].detectorID ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_u[j].first>=0){
					if(ic[index].EventID==0)printf("u first hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID, ic[index].AllHits[ oc[index].hitpairs_u[j].first ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_u[j].first ];
					oc[index].AllTracklets[n_tkl].nUHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;
					
					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_u[j].first ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_u[j].first ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[i].first ].detectorID ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[i].first ].detectorID ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_u[j].first ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_u[j].first ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].first ].detectorID ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_u[j].second>=0){
					if(ic[index].EventID==0)printf("u second hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID, ic[index].AllHits[ oc[index].hitpairs_u[j].second ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_u[j].second ];
					oc[index].AllTracklets[n_tkl].nUHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;

					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_u[j].second ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_u[j].second ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[i].second ].detectorID ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[i].second ].detectorID ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_u[j].second ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_u[j].second ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_u[j].second ].detectorID ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_v[k].first>=0){
					if(ic[index].EventID==0)printf("v first hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID, ic[index].AllHits[ oc[index].hitpairs_v[k].first ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_v[k].first ];
					oc[index].AllTracklets[n_tkl].nVHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;
					
					fitarrays[index].drift_dist[npts] = 0.;//ic[index].AllHits[ oc[index].hitpairs_v[k].first ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_v[k].first ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[i].first ].detectorID ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[i].first ].detectorID ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_v[k].first ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_v[k].first ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].first ].detectorID ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(oc[index].hitpairs_v[k].second>=0){
					if(ic[index].EventID==0)printf("v second hit plane = %d, elem %d, z = %1.4f\n", ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID, ic[index].AllHits[ oc[index].hitpairs_v[k].second ].elementID, planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID ].p1z_w1);
					oc[index].AllTracklets[n_tkl].hits[nhits_tkl]=ic[index].AllHits[ oc[index].hitpairs_v[k].second ];
					oc[index].AllTracklets[n_tkl].nVHits++;
					oc[index].AllTracklets[n_tkl].hitsign[nhits_tkl]=0;

					fitarrays[index].drift_dist[npts] = 0;//ic[index].AllHits[ oc[index].hitpairs_v[k].second ].driftDistance;
					fitarrays[index].resolution[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID ].spacing/3.4641f;//.resolution;
					fitarrays[index].p1x[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID ].p1x_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID ].dp1x * (ic[index].AllHits[ oc[index].hitpairs_v[k].second ].elementID-1) ;
					fitarrays[index].p1y[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[i].second ].detectorID ].p1y_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[i].second ].detectorID ].dp1y * (ic[index].AllHits[ oc[index].hitpairs_v[k].second ].elementID-1) ;
					fitarrays[index].p1z[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID ].p1z_w1 + planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID ].dp1z * (ic[index].AllHits[ oc[index].hitpairs_v[k].second ].elementID-1) ;
					fitarrays[index].deltapx[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID ].deltapx;
					fitarrays[index].deltapy[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID ].deltapy;
					fitarrays[index].deltapz[npts] = planes[ ic[index].AllHits[ oc[index].hitpairs_v[k].second ].detectorID ].deltapz;
					npts++;
					nhits_tkl++;
				}
				if(npts<4)continue;
				if(oc[index].AllTracklets[n_tkl].nXHits<1 || oc[index].AllTracklets[n_tkl].nUHits<1 || oc[index].AllTracklets[n_tkl].nVHits<1)continue;
				

				fitarrays[index].output_parameters[0] = 0.;
				fitarrays[index].output_parameters[1] = 0.;
				fitarrays[index].output_parameters[2] = 0.;
				fitarrays[index].output_parameters[3] = 0.;
				/*
				get_straighttrack_fixedpoint(npts, fitarrays[index].drift_dist, fitarrays[index].resolution, 
								fitarrays[index].p1x, fitarrays[index].p1y, fitarrays[index].p1z, 
								fitarrays[index].deltapx, fitarrays[index].deltapy, fitarrays[index].deltapz,
								fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B,
								fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, 
								fixedpoint, fitarrays[index].chi2);

				fitarrays[index].lambda = 0.01f;
				fitarrays[index].chi2prev = 1.e20;
				
				gpufit_algorithm_fitter(fitparams[0].max_iterations,
										fitarrays[index].n_iter, fitarrays[index].iter_failed, fitarrays[index].finished,
										fitarrays[index].state, fitarrays[index].skip, fitarrays[index].singular,
										fitparams[0].tolerance, fitparams[0].nparam,
										fitparams[0].parameter_limits_min, fitparams[0].parameter_limits_max,  
										fitarrays[index].output_parameters, fitarrays[index].prev_parameters, 
										fitarrays[index].chi2, fitarrays[index].chi2prev,
										fitarrays[index].values, fitarrays[index].derivatives, fitarrays[index].gradients,
										fitarrays[index].hessians, fitarrays[index].scaling_vector, 
										fitarrays[index].deltas, fitarrays[index].lambda,
										fitarrays[index].Ainv,
										//fitarrays[index].calc_matrix, fitarrays[index].abs_row, fitarrays[index].abs_row_index,
										npts, fitarrays[index].drift_dist, fitarrays[index].resolution, 
										fitarrays[index].p1x, fitarrays[index].p1y, fitarrays[index].p1z, 
										fitarrays[index].deltapx, fitarrays[index].deltapy, fitarrays[index].deltapz);					
				*/
				oc[index].AllTracklets[n_tkl].x0 = fitarrays[index].output_parameters[0];
				oc[index].AllTracklets[n_tkl].y0 = fitarrays[index].output_parameters[1];
				oc[index].AllTracklets[n_tkl].tx = fitarrays[index].output_parameters[2];
				oc[index].AllTracklets[n_tkl].ty = fitarrays[index].output_parameters[3];
				oc[index].AllTracklets[n_tkl].err_x0 = fitarrays[index].output_parameters_errors[0];
				oc[index].AllTracklets[n_tkl].err_y0 = fitarrays[index].output_parameters_errors[1];
				oc[index].AllTracklets[n_tkl].err_tx = fitarrays[index].output_parameters_errors[2];
				oc[index].AllTracklets[n_tkl].err_ty = fitarrays[index].output_parameters_errors[3];
				oc[index].AllTracklets[n_tkl].chisq = fitarrays[index].chi2;
				
				/*
				if(ic[index].EventID==0)printf("track: x0 = %1.6f +- %1.6f, y0 = %1.6f +- %1.6f, tx = %1.6f +- %1.6f, ty = %1.6f +- %1.6f; chi2 = %1.6f\n", 
					oc[index].AllTracklets[n_tkl].x0, oc[index].AllTracklets[n_tkl].err_x0, 
					oc[index].AllTracklets[n_tkl].y0, oc[index].AllTracklets[n_tkl].err_y0, 
					oc[index].AllTracklets[n_tkl].tx, oc[index].AllTracklets[n_tkl].err_tx, 
					oc[index].AllTracklets[n_tkl].ty, oc[index].AllTracklets[n_tkl].err_ty,
					oc[index].AllTracklets[n_tkl].chisq);
				*/

				if(!match_tracklet_to_hodo(oc[index].AllTracklets[n_tkl], stID, ic, planes))continue;

				
				if(n_tkl>TrackletSizeMax){
					continue;
					//printf("evt %d: n_tkl = %d > %d\n", oc[index].EventID, n_tkl, TrackletSizeMax);
				}
				n_tkl++;
				oc[index].nTKL_stID[stID]++;

			}
		}
	}
	
	//if(print)printf("evt %d number of tracklets %d\n", oc[index].EventID, n_tkl);
	//printf("%d ", n_tkl);
	oc[index].nTracklets = n_tkl;
}




// kernel to combine tracklets into back partial tracks
__global__ void gkernel_BackPartialTracks(gEvent* ic, gSW* oc, gFitArrays* fitarrays, const gPlane* planes, gFitParams* fitparams){
	int index = threadIdx.x + blockIdx.x * blockDim.x;
	short stID = 5;
	//for the time being, declare in function;
	REAL x_pos[4];
	REAL z_pos[4];
	REAL a, b, sum, det, sx, sy, sxx, syy, sxy;
	short nprop, iprop;
	REAL xExp;
	short nhits_X2, nhits_X3;
	int n_tkl = oc[index].nTracklets;
	int nhits_2 = 0, nhits_3 = 0;
	
	REAL fixedpoint[3] = {0, 0, 0};

	//tracklets in station 3  (D3p, stID 4-1, D3m, stID 5-1)
	for(int i = oc[index].nTKL_stID[2]; i<oc[index].nTKL_stID[2]+oc[index].nTKL_stID[3]+oc[index].nTKL_stID[4]; i++){
		nhits_2 = oc[index].AllTracklets[i].nXHits+oc[index].AllTracklets[i].nUHits+oc[index].AllTracklets[i].nVHits;
		nhits_X2 = 0;
		for(int k = 0; k<nhits_2; k++){
			if(oc[index].AllTracklets[i].hits[k].detectorID==geometry::dets_x[0] || oc[index].AllTracklets[i].hits[k].detectorID==geometry::dets_x[1])
			{
				x_pos[nhits_X2] = oc[index].AllTracklets[i].hits[k].pos;
				z_pos[nhits_X2] = planes[oc[index].AllTracklets[i].hits[k].detectorID].z;
				nhits_X2++;
			}
		}
		
		//tracklets in station 2 (D2, stID 3-1)
		for(int j = 0; j<oc[index].nTKL_stID[2]; j++){
			if(fabs(oc[index].AllTracklets[i].tx - oc[index].AllTracklets[j].tx) > TX_MAX || fabs(oc[index].AllTracklets[i].ty - oc[index].AllTracklets[i].ty) > TY_MAX)continue;
			
			nhits_3 = oc[index].AllTracklets[j].nXHits+oc[index].AllTracklets[j].nUHits+oc[index].AllTracklets[j].nVHits;
			nhits_X3 = nhits_X2;
			for(int k = 0; k<nhits_3; k++){
				if(oc[index].AllTracklets[j].hits[k].detectorID==geometry::dets_x[2] || oc[index].AllTracklets[j].hits[k].detectorID==geometry::dets_x[3] ||
			   	oc[index].AllTracklets[j].hits[k].detectorID==geometry::dets_x[4] || oc[index].AllTracklets[j].hits[k].detectorID==geometry::dets_x[5])
				{
					x_pos[nhits_X3] = oc[index].AllTracklets[i].hits[k].pos;
					z_pos[nhits_X3] = planes[oc[index].AllTracklets[i].hits[k].detectorID].z;
					nhits_X3++;
				}
			}
			
			chi2_simplefit(nhits_X2+nhits_X3, z_pos, x_pos, a, b, sum, det, sx, sy, sxx, syy, sxy);
			if(fabs(a)>X0_MAX || fabs(b)>TX_MAX)continue;
			//prop matching
			nprop = 0;
			for(short ip = 0; ip<4; ip++){
				iprop = 48+ip;
				xExp = a*planes[iprop].z+b;
				//loop on hits to find prop
				for(int k = 0; k<ic[index].nAH; k++){
					if(ic[index].AllHits[k].detectorID==iprop){
						if(fabs(ic[index].AllHits[k].pos-xExp)<5.08f){
							nprop++;
							break;
						}
					}
				}
				if(nprop>0)break;
			}
			if(nprop==0)continue;
				
			//Add tracklets: first, get number of hits in each
			oc[index].AllTracklets[n_tkl].stationID = stID;
			oc[index].AllTracklets[n_tkl].nXHits = oc[index].AllTracklets[i].nXHits+oc[index].AllTracklets[j].nXHits;
			oc[index].AllTracklets[n_tkl].nUHits = oc[index].AllTracklets[i].nUHits+oc[index].AllTracklets[j].nUHits;
			oc[index].AllTracklets[n_tkl].nVHits = oc[index].AllTracklets[i].nVHits+oc[index].AllTracklets[j].nVHits;
			
			//add in the new tracklet the hits from the two existing tracklets;
			//also fill the fit arrays
			for(int k = 0; k<nhits_2; k++){
				oc[index].AllTracklets[n_tkl].hits[k]=oc[index].AllTracklets[i].hits[k];
				oc[index].AllTracklets[n_tkl].hitsign[k]=oc[index].AllTracklets[i].hitsign[k];
				
				fitarrays[index].drift_dist[k] = 0.;//oc[index].AllTracklets[n_tkl].hits[k].driftDistance;
				fitarrays[index].resolution[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].spacing/3.4641f;//.resolution;
				fitarrays[index].p1x[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].p1x_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].dp1x * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1);
				fitarrays[index].p1y[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].p1y_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].dp1y * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1);
				fitarrays[index].p1z[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].p1z_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].dp1z * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1);
				fitarrays[index].deltapx[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].deltapx;
				fitarrays[index].deltapy[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].deltapy;
				fitarrays[index].deltapz[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].deltapz;
			}
				
			for(int k = nhits_2; k<nhits_2+nhits_3; k++){
				oc[index].AllTracklets[n_tkl].hits[k]=oc[index].AllTracklets[j].hits[k-nhits_2];
				oc[index].AllTracklets[n_tkl].hitsign[k]=oc[index].AllTracklets[j].hitsign[k-nhits_2];
					fitarrays[index].drift_dist[k] = 0.;//oc[index].AllTracklets[n_tkl].hits[k].driftDistance;
				fitarrays[index].resolution[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].spacing/3.4641f;//.resolution;
				fitarrays[index].p1x[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].p1x_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].dp1x * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1) ;
				fitarrays[index].p1y[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].p1y_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].dp1y * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1) ;
				fitarrays[index].p1z[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].p1z_w1 + planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].dp1z * (oc[index].AllTracklets[n_tkl].hits[k].elementID-1) ;
				fitarrays[index].deltapx[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].deltapx;
				fitarrays[index].deltapy[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].deltapy;
				fitarrays[index].deltapz[k] = planes[ oc[index].AllTracklets[n_tkl].hits[k].detectorID ].deltapz;
			}
			
			//then fit:
				
			get_straighttrack_fixedpoint(nhits_2+nhits_3, fitarrays[index].drift_dist, fitarrays[index].resolution, 
							fitarrays[index].p1x, fitarrays[index].p1y, fitarrays[index].p1z, 
							fitarrays[index].deltapx, fitarrays[index].deltapy, fitarrays[index].deltapz,
							fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B,
							fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, 
							fixedpoint, fitarrays[index].chi2);
			
			fitarrays[index].output_parameters[0] = b;
			fitarrays[index].output_parameters[2] = a;
			
			// do a rough estimation of the starting fit
			//fitarrays[index].output_parameters[0] = (oc[index].AllTracklets[i].x0+oc[index].AllTracklets[j].x0)*0.5f;
			//fitarrays[index].output_parameters[1] = (oc[index].AllTracklets[i].y0+oc[index].AllTracklets[j].y0)*0.5f;
			//fitarrays[index].output_parameters[2] = (oc[index].AllTracklets[i].tx+oc[index].AllTracklets[j].tx)*0.5f;
			//fitarrays[index].output_parameters[3] = (oc[index].AllTracklets[i].ty+oc[index].AllTracklets[j].ty)*0.5f;
			fitarrays[index].lambda = 0.01f;
			fitarrays[index].chi2prev = 1.e20;
			
			//include fit here:
			/*
			gpufit_algorithm_fitter(fitparams[0].max_iterations,
								fitarrays[index].n_iter, fitarrays[index].iter_failed, fitarrays[index].finished,
								fitarrays[index].state, fitarrays[index].skip, fitarrays[index].singular,
								fitparams[0].tolerance, fitparams[0].nparam,
								fitparams[0].parameter_limits_min, fitparams[0].parameter_limits_max,  
								fitarrays[index].output_parameters, fitarrays[index].prev_parameters, 
								fitarrays[index].chi2, fitarrays[index].chi2prev,
								fitarrays[index].values, fitarrays[index].derivatives, fitarrays[index].gradients,
								fitarrays[index].hessians, fitarrays[index].scaling_vector, 
								fitarrays[index].deltas, fitarrays[index].lambda,
								fitarrays[index].Ainv, 
								nhits_2+nhits_3, fitarrays[index].drift_dist, fitarrays[index].resolution, 
								fitarrays[index].p1x, fitarrays[index].p1y, fitarrays[index].p1z, 
								fitarrays[index].deltapx, fitarrays[index].deltapy, fitarrays[index].deltapz);					
			
			*/
			
			if(fitarrays[index].chi2>9000)continue;
			
			oc[index].AllTracklets[n_tkl].x0 = fitarrays[index].output_parameters[0];
			oc[index].AllTracklets[n_tkl].y0 = fitarrays[index].output_parameters[1];
			oc[index].AllTracklets[n_tkl].tx = fitarrays[index].output_parameters[2];
			oc[index].AllTracklets[n_tkl].ty = fitarrays[index].output_parameters[3];
			oc[index].AllTracklets[n_tkl].err_x0 = fitarrays[index].output_parameters_errors[0];
			oc[index].AllTracklets[n_tkl].err_y0 = fitarrays[index].output_parameters_errors[1];
			oc[index].AllTracklets[n_tkl].err_tx = fitarrays[index].output_parameters_errors[2];
			oc[index].AllTracklets[n_tkl].err_ty = fitarrays[index].output_parameters_errors[3];
			oc[index].AllTracklets[n_tkl].chisq = fitarrays[index].chi2;
			
			n_tkl++;
			oc[index].nTKL_stID[stID]++;
		}
	}
	
	oc[index].nTracklets = n_tkl;
	
}
