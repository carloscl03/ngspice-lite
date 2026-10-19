/* ngpool.c - persistent worker pool for the parallel device load loops.
 *
 * The OpenMP loops in the device loaders open a fork/join region on every
 * Newton iteration and leave the idle threads to libgomp: with its default
 * spin they keep the package at the all-core frequency while the master
 * runs the serial phases (stamping, factorisation, step control), and with
 * passive waiting every region pays a wake-up. Measured on a 59-MOSFET
 * sky130 transient (Ryzen 7 3700U, WSL2): the serial gap between two
 * loads grew from 38 us with one thread to 56 us with four OpenMP
 * threads.
 *
 * This pool keeps the worker threads alive across regions, hands the
 * instances out dynamically in small chunks, and after a short bounded
 * spin puts the workers to sleep on a condition variable, so the serial
 * phases run on an otherwise idle package. Same circuit, four threads:
 * region 29 us (OpenMP: 41), gap 40 us (OpenMP: 56).
 *
 * Thread count: `set num_threads=N` (CKTsetup). Tuning, read once:
 *   NGSPICE_POOL_SPIN=K   pause iterations before sleeping (5000)
 *   NGSPICE_POOL_WAIT=yield|cond   policy after the spin (cond)
 *   NGSPICE_POOL_CHUNK=C  instances per grab (2)
 *   NGSPICE_POOL_PIN=1    pin thread t to CPU 2t (physical cores on SMT boxes)
 *   NGSPICE_NOPOOL=1      keep the OpenMP loops (for comparisons)
 */

#include "ngspice/ngspice.h"
#include "ngspice/ngpool.h"

#ifdef USE_OMP

#include <pthread.h>
#include <stdatomic.h>
#include <sched.h>
#include <stdint.h>
#ifdef __x86_64__
#include <immintrin.h>
#define NGPOOL_PAUSE() _mm_pause()
#else
#define NGPOOL_PAUSE() do { } while (0)
#endif

#define NGPOOL_MAXTHREADS 64

static struct {
    int nthreads;              /* including the caller, which is thread 0 */
    long spin;
    int wait;                  /* 0 = sched_yield loop, 1 = condvar sleep */
    int chunk;
    int pin;
    int disabled;
    int env_read;
    pthread_t th[NGPOOL_MAXTHREADS];
    _Atomic int generation;
    _Atomic int next;
    _Atomic int done;
    _Atomic int error;
    _Atomic int quit;
    ngpool_item_fn fn;
    void *arg;
    int n;
    pthread_mutex_t mu;
    pthread_cond_t cv;
} pool;


static int env_int(const char *name, int dflt)
{
    const char *e = getenv(name);
    return e ? atoi(e) : dflt;
}


static void read_env(void)
{
    const char *w;
    if (pool.env_read)
        return;
    pool.env_read = 1;
    pool.spin = env_int("NGSPICE_POOL_SPIN", 5000);
    pool.chunk = env_int("NGSPICE_POOL_CHUNK", 2);
    if (pool.chunk < 1)
        pool.chunk = 1;
    pool.pin = env_int("NGSPICE_POOL_PIN", 0);
    w = getenv("NGSPICE_POOL_WAIT");
    pool.wait = (w && strcmp(w, "yield") == 0) ? 0 : 1;
    pool.disabled = env_int("NGSPICE_NOPOOL", 0);
    pthread_mutex_init(&pool.mu, NULL);
    pthread_cond_init(&pool.cv, NULL);
}


static void do_work(void)
{
    int ch = pool.chunk;
    for (;;) {
        int i = atomic_fetch_add_explicit(&pool.next, ch, memory_order_relaxed);
        int e;
        if (i >= pool.n)
            break;
        e = i + ch;
        if (e > pool.n)
            e = pool.n;
        for (; i < e; i++) {
            int err = pool.fn(i, pool.arg);
            if (err)
                atomic_store(&pool.error, err);
        }
    }
}


static void pin_self(int tid)
{
    cpu_set_t cs;
    CPU_ZERO(&cs);
    CPU_SET(2 * tid, &cs);
    pthread_setaffinity_np(pthread_self(), sizeof cs, &cs);
}


static void *worker(void *v)
{
    int tid = (int)(intptr_t)v;
    int seen = 0;
    if (pool.pin)
        pin_self(tid);
    for (;;) {
        long i = 0;
        while (atomic_load_explicit(&pool.generation, memory_order_acquire) == seen) {
            if (i < pool.spin) {
                i++;
                NGPOOL_PAUSE();
            } else if (pool.wait == 0) {
                sched_yield();
            } else {
                pthread_mutex_lock(&pool.mu);
                while (atomic_load(&pool.generation) == seen)
                    pthread_cond_wait(&pool.cv, &pool.mu);
                pthread_mutex_unlock(&pool.mu);
            }
        }
        seen = atomic_load(&pool.generation);
        if (atomic_load(&pool.quit))
            break;
        do_work();
        atomic_fetch_add_explicit(&pool.done, 1, memory_order_release);
    }
    return NULL;
}


static void wake_all(void)
{
    if (pool.wait == 1) {
        pthread_mutex_lock(&pool.mu);
        atomic_fetch_add_explicit(&pool.generation, 1, memory_order_release);
        pthread_cond_broadcast(&pool.cv);
        pthread_mutex_unlock(&pool.mu);
    } else {
        atomic_fetch_add_explicit(&pool.generation, 1, memory_order_release);
    }
}


void ngpool_stop(void)
{
    int t;
    if (pool.nthreads < 2)
        return;
    atomic_store(&pool.quit, 1);
    wake_all();
    for (t = 1; t < pool.nthreads; t++)
        pthread_join(pool.th[t], NULL);
    atomic_store(&pool.quit, 0);
    pool.nthreads = 0;
}


int ngpool_configure(int nthreads)
{
    int t;
    read_env();
    if (pool.disabled)
        nthreads = 0;
    if (nthreads > NGPOOL_MAXTHREADS)
        nthreads = NGPOOL_MAXTHREADS;
    if (nthreads < 2)
        nthreads = 0;
    if (nthreads == pool.nthreads)
        return pool.nthreads;
    ngpool_stop();
    if (nthreads == 0)
        return 0;
    pool.nthreads = nthreads;
    if (pool.pin)
        pin_self(0);
    for (t = 1; t < nthreads; t++)
        pthread_create(&pool.th[t], NULL, worker, (void *)(intptr_t)t);
    return nthreads;
}


int ngpool_threads(void)
{
    return pool.nthreads;
}


int ngpool_run(int n, ngpool_item_fn fn, void *arg)
{
    pool.fn = fn;
    pool.arg = arg;
    pool.n = n;
    atomic_store(&pool.next, 0);
    atomic_store(&pool.done, 0);
    atomic_store(&pool.error, 0);
    wake_all();
    do_work();
    while (atomic_load_explicit(&pool.done, memory_order_acquire) < pool.nthreads - 1)
        NGPOOL_PAUSE();
    return atomic_load(&pool.error);
}

#else /* !USE_OMP */

void ngpool_stop(void) { }
int ngpool_configure(int nthreads) { NG_IGNORE(nthreads); return 0; }
int ngpool_threads(void) { return 0; }
int ngpool_run(int n, ngpool_item_fn fn, void *arg)
{
    int i, err = 0;
    for (i = 0; i < n; i++) {
        int e = fn(i, arg);
        if (e)
            err = e;
    }
    return err;
}

#endif
