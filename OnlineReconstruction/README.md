System requirements:
- At least one of the CUDA-enabled GPUs (list of CUDA-enabled GPUs: https://developer.nvidia.com/cuda-gpus).
- CUDA 10.0
- ROOT 6.14

Files:
Main: online_reconstruction.cu
Auxillary classes: classes for handling root input and output.
	  - LoadInput.h, LoadInput.cxx;
	  - OROutput.h, OROutput.cxx: 
headers: - reconstruction constants.h: carries the constants that will be used throughout the program;
	 - operations.h: definitions for matrix calculations
cuda headers: - reconstruction_classes.cuh: defines the structures used by the program;
     	      - reconstruction_kernels.cuh: defines the main core/global functions (kernels);
	      - reconstruction_helper.cuh: defines all auxilliary functions (device);
	      - trackanalyticalminimizer.cuh: defines functions that will be called for track fitting
	      - trackfittingtools.cuh: defines elementary functions for track fitting (namely - but not restricted to - matrix operations);


Input files:
- Decoded rootfiles with data saved as "SRawEvent" (example: /data3/analysis/production/02/87/digit_028705_009.root).

