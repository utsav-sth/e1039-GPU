#ifdef KTRACKER_REC	
	// copy result of event reconstruction from device_gEvent to device_input_TKL
	// this input_tkl should be the information that the device uses to reconstruct the tracklets
	// gpuErrchk( cudaMemcpy(device_input_TKL, device_gEvent, NBytesAllEvent, cudaMemcpyDeviceToDevice));

	// shouldn't this function actually be called? should it be the function that puts together tracklets? and then call the fitting???
	// gkernel_TKL<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_input_TKL, device_output_TKL);

	//for(int m = 1; m<=30; m++){
	//	if(plane[m].u_win!=0)printf("plane, m = %d, u_win = %1.6f, costheta = %1.6f\n", m, plane[m].u_win, plane[m].costheta);
	//	if(device_gPlane[m].u_win!=0)printf("device_gplane, m = %d, u_win = %1.6f, costheta = %1.6f\n", m, device_gPlane[m].u_win, device_gPlane[m].costheta);
	//}
	
	auto start_tkl2 = std::chrono::system_clock::now();
	// I first want to see if indeed we can reuse the "gEvent" pointer
	int stID = 3;// to make explicit that we are requiring station 3
	
	gkernel_TrackletinStation<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output_TKL, device_gFitArrays, stID, device_gPlane, device_gFitParams);
	auto end_tkl2 = std::chrono::system_clock::now();

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	auto start_tkl3 = std::chrono::system_clock::now();
	stID = 4;
	gkernel_TrackletinStation<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output_TKL, device_gFitArrays, stID, device_gPlane, device_gFitParams);

	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );

	stID = 5;
	gkernel_TrackletinStation<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output_TKL, device_gFitArrays, stID, device_gPlane, device_gFitParams);
	//cout << endl;
	// check status of device and synchronize again;
	auto end_tkl3 = std::chrono::system_clock::now();
	
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
	
	gkernel_BackPartialTracks<<<BLOCKS_NUM,THREADS_PER_BLOCK>>>(device_gEvent, device_output_TKL, device_gFitArrays, device_gPlane, device_gFitParams);
	
	gpuErrchk( cudaPeekAtLastError() );
	gpuErrchk( cudaDeviceSynchronize() );
#endif
