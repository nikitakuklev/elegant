#ifndef GPU_LRWAKE_H
#define GPU_LRWAKE_H

void gpu_determine_bucket_assignments(long np, long idSlotsPerBunch, 
       double P0, double **d_time, long **npBucket, long *nBuckets, long lastNBuckets);

#endif


