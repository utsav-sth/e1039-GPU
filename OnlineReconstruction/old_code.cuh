#ifdef OLDCODE
				if(straighttrackbuilder[index].hitpairs_u2[bin2+nbins_st2*i_u2].first>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[bin2+nbins_st2*i_u2].first], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_u2[bin2+nbins_st2*i_u2].first;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}

				if(straighttrackbuilder[index].hitpairs_u2[bin2+nbins_st2*i_u2].second>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u2[bin2+nbins_st2*i_u2].second], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_u2[bin2+nbins_st2*i_u2].second;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}
				
				
				if(straighttrackbuilder[index].hitpairs_v2[bin2+nbins_st2*i_v2].first>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[bin2+nbins_st2*i_v2].first], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_v2[bin2+nbins_st2*i_v2].first;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}
				if(straighttrackbuilder[index].hitpairs_v2[bin2+nbins_st2*i_v2].second>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v2[bin2+nbins_st2*i_v2].second], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_v2[bin2+nbins_st2*i_v2].second;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}

				
				if(straighttrackbuilder[index].hitpairs_u3[bin3+nbins_st3*i_u3].first>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3[bin3+nbins_st3*i_u3].first], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_u3[bin3+nbins_st3*i_u3].first;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}
				if(straighttrackbuilder[index].hitpairs_u3[bin3+nbins_st3*i_u3].second>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_u3[bin3+nbins_st3*i_u3].second], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_u3[bin3+nbins_st3*i_u3].second;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}

				nhits_u3 = nhits_uv-nhits_u2-nhits_v2;
				if(nhits_u3==0) continue;

				if(straighttrackbuilder[index].hitpairs_v3[bin3+nbins_st3*i_v3].first>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3[bin3+nbins_st3*i_v3].first], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_v3[bin3+nbins_st3*i_v3].first;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}
				if(straighttrackbuilder[index].hitpairs_v3[bin3+nbins_st3*i_v3].second>=0){
					if(calculate_y_uvhit(y, err_y, ic[index].AllHits[straighttrackbuilder[index].hitpairs_v3[bin3+nbins_st3*i_v3].second], 0, straighttrackbuilder[index].trackXZ, planes)){
					straighttrackbuilder[index].trackYZ.hitlist[nhits_uv] = straighttrackbuilder[index].hitpairs_v3[bin3+nbins_st3*i_v3].second;
					FillFitArrays_UV(nhits_uv, ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[nhits_uv]], fitarrays[index], planes, y, err_y);
					nhits_uv++;
					}
				}
				nhits_v3 = nhits_uv-nhits_u2-nhits_v2-nhits_u3;
				if(nhits_v3==0) continue;
				
				fit_2D_track(nhits_uv, fitarrays[index].y_array, fitarrays[index].z_array, fitarrays[index].dy_array, fitarrays[index].A, fitarrays[index].Ainv, fitarrays[index].B, fitarrays[index].output_parameters, fitarrays[index].output_parameters_errors, fitarrays[index].chi2_2d);

				straighttrackbuilder[index].trackYZ.x_0 = fitarrays[index].output_parameters[0];
				straighttrackbuilder[index].trackYZ.err_x_0 = fitarrays[index].output_parameters_errors[0];
				straighttrackbuilder[index].trackYZ.tx_ = fitarrays[index].output_parameters[1];
				straighttrackbuilder[index].trackYZ.err_tx_ = fitarrays[index].output_parameters_errors[1];
				straighttrackbuilder[index].trackYZ.nhits = nhits_uv;
				
				// now evaluate the back track candidate.
				oc[index].AllTracklets[ntkl].y0 = straighttrackbuilder[index].trackYZ.x_0;
				oc[index].AllTracklets[ntkl].err_y0 = straighttrackbuilder[index].trackYZ.err_x_0;
				oc[index].AllTracklets[ntkl].ty = straighttrackbuilder[index].trackYZ.tx_;
				oc[index].AllTracklets[ntkl].err_ty = straighttrackbuilder[index].trackYZ.err_tx_;
				
				// filter with hodoscope matching before evaluating chi2.
				if(!match_tracklet_to_hodo(oc[index].AllTracklets[ntkl], 2, ic, planes))continue;
				if(!match_tracklet_to_hodo(oc[index].AllTracklets[ntkl], 3, ic, planes) &&
				 !match_tracklet_to_hodo(oc[index].AllTracklets[ntkl], 4, ic, planes))continue;
				
				for(i_hit = nhits_x; i_hit<nhits_x+nhits_uv; i_hit++){
					oc[index].AllTracklets[ntkl].hits[i_hit] = ic[index].AllHits[straighttrackbuilder[index].trackYZ.hitlist[i_hit-nhits_x]];
				}
				oc[index].AllTracklets[ntkl].nUHits = nhits_u2+nhits_u3;
				oc[index].AllTracklets[ntkl].nVHits = nhits_v2+nhits_v3;
				
				resolve_leftright(oc[index].AllTracklets[ntkl], planes, 40.);
				resolve_leftright(oc[index].AllTracklets[ntkl], planes, 150.);
				resolve_single_leftright(oc[index].AllTracklets[oc[index].nTracklets], planes);
			
				refit_backpartialtrack_with_drift(oc[index].AllTracklets[ntkl], fitarrays[index], planes);
				for(i_hit = 0; i_hit<nhits_x; i_hit++){
					straighttrackbuilder[index].trackXZ.hitsign[i_hit] = oc[index].AllTracklets[ntkl].hitsign[i_hit];
				}
				for(i_hit = nhits_x; i_hit<nhits_x+nhits_uv; i_hit++){
					straighttrackbuilder[index].trackYZ.hitsign[i_hit-nhits_x] = oc[index].AllTracklets[ntkl].hitsign[i_hit];
				}
				
				if(best_candyz_only){
					if(fitarrays[index].chi2<chi2min){
						chi2min = fitarrays[index].chi2;
						straighttrackbuilder[index].besttrackYZ.x_0 = oc[index].AllTracklets[ntkl].y0;
						straighttrackbuilder[index].besttrackYZ.err_x_0 = oc[index].AllTracklets[ntkl].err_y0;
						straighttrackbuilder[index].besttrackYZ.tx_ = oc[index].AllTracklets[ntkl].ty;
						straighttrackbuilder[index].besttrackYZ.err_tx_ = oc[index].AllTracklets[ntkl].err_ty;
						straighttrackbuilder[index].besttrackYZ.nhits = nhits_uv;
						nhits_v = nhits_v2+nhits_v3;
						for(i_hit = 0; i_hit<nhits_uv; i_hit++){
							straighttrackbuilder[index].besttrackYZ.hitlist[i_hit] = straighttrackbuilder[index].trackYZ.hitlist[i_hit];
							straighttrackbuilder[index].besttrackYZ.hitsign[i_hit] = straighttrackbuilder[index].trackYZ.hitsign[i_hit];
						}
					}
				}else{
					oc[index].AllTracklets[ntkl] = tkl;
					ntkl++;
				}
#endif
