System requirements:
- At least one of the CUDA-enabled GPUs (list of CUDA-enabled GPUs: https://developer.nvidia.com/cuda-gpus).
- CUDA 10.0
- ROOT 6.14

Available parameters:
- EstnEvtMax: estimated maximum number of events.
- THREADS_PER_BLOCK: threads per block.
- BLOCKS_NUM: the number of blocks (EstnEvtMax/THREADS_PER_BLOCK).
- EstnAHMax: estimated maximum number of all hits per event.
- EstnTHMax: estimated maximum number of trigger hits per event.
- ClusterSizeMax: maximum number of hits allowed in a cluster.

Input files:
- Decoded rootfiles with data saved as "SRawEvent" (example: /data3/analysis/production/02/87/digit_028705_009.root).

