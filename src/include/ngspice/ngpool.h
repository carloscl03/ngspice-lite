/* ngpool.h - persistent worker pool for the parallel device load loops.
 * See src/maths/ni/ngpool.c. */
#ifndef ngspice_NGPOOL_H
#define ngspice_NGPOOL_H

#include "ngspice/cktdefs.h"

typedef int (*ngpool_item_fn)(int idx, void *arg);

/* job shared by the device loaders: the instance array and the circuit */
struct ngpool_devjob {
    void **inst;
    CKTcircuit *ckt;
};

/* start (or resize) the pool with nthreads threads including the caller;
 * returns the number in use, 0 when the pool is off */
int ngpool_configure(int nthreads);
int ngpool_threads(void);
void ngpool_stop(void);

/* run fn(i, arg) for i in [0, n) on the pool; returns the last non-zero
 * value returned by fn, 0 if all returned 0 */
int ngpool_run(int n, ngpool_item_fn fn, void *arg);
/* same, with an explicit chunk size (1 = one item per grab) */
int ngpool_run_chunk(int n, ngpool_item_fn fn, void *arg, int chunk);

#endif
