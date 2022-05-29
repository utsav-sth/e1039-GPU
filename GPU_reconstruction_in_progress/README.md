System requirements:
- At least one of the CUDA-enabled GPUs (list of CUDA-enabled GPUs: https://developer.nvidia.com/cuda-gpus).
- CUDA 11.0
- ROOT 6.14
- GPUfit library (https://github.com/gpufit/Gpufit)
- cuBLASLt API (https://docs.nvidia.com/cuda/cublas/index.html)
- Boost 1.58 or later (To build the tests)
- Python for building the Python bindings (Python version 2.x or 3.x)

Input files:
- Rootfiles in E906/E1039 format
- Switch data format type to process in "LoadInput.cxx" file
.sh

To run main source code, simply execute file titles "run.sh" which calls the MakeFile first and then defines the data/MC file needed to be processed.

