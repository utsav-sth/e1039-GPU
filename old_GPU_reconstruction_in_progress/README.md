To clone this repository:
- git clone git@github.com:efuchey/e1039-GPU

System requirements:
- At least one of the CUDA-enabled GPUs (list of CUDA-enabled GPUs: https://developer.nvidia.com/cuda-gpus).
- CUDA Toolkit 11.0
- ROOT 6.14
- CMake 3.11 or later
- GPUfit library (https://github.com/gpufit/Gpufit)
- cuBLASLt API (https://docs.nvidia.com/cuda/cublas/index.html)
- Boost 1.58 or later (To build the tests)
- Python for building the Python bindings (Python version 2.x or 3.x)

Input files:
- Rootfiles in E-906/E-1039 format
- Switch data format type to process in "LoadInput.cxx" file

Latest/most advanced version of the GPU analysis code is "online_reconstruction.cu"

To run main executable, simply run file titled "run.sh" which calls the MakeFile first, runs executable and defines the data/MC file needed to be processed



