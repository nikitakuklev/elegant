#ifndef GPU_LRWAKE_H
#define GPU_LRWAKE_H

void gpu_index_bunch_assignments(long np, long idSlotsPerBunch, 
       double P0, double **d_time, long **npBunch, long *nBunches, long lastNBunches);

#endif


