/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: September 2003 Paolo Nenzi
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "capdefs.h"
#include "ngspice/trandefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/ngpool.h"
#include "ngspice/hash.h"
#include <time.h>
#include <string.h>

/* fork-lite: the capacitor loop is a walk over every instance that reads two
   node voltages and the charge history and writes four matrix entries and two
   RHS entries. With tens of thousands of parasitic capacitors it is one of the
   serial stretches between the parallel BSIM4 regions (SPI with parasitics,
   77 k capacitors: 1.7 ms of evaluation and 4.2 ms of stamping per call, the
   stamping being 460 k scattered adds).

   Same architecture as the fused BSIM4 load: each pool thread evaluates and
   stamps its own static range of instances; a matrix entry or RHS node that
   is touched by more than one range is redirected to a per-thread slot and
   the slots are reduced afterwards (in parallel, by ranges of slots). Entries
   touched by one range only are written directly. Only the shared entries can
   differ from the serial loop, and only in the rounding of the sum order.
   NGSPICE_NOPCAP=1 restores the serial loop; NGSPICE_LOADTRACE=1 prints the
   table sizes and the time per phase at exit. */

static int cl_on = -1, cl_trace = 0;
static long cl_calls = 0;
static double cl_eval = 0, cl_reduce = 0;

static struct {
    int valid, nthreads, ninst;
    int *lo;              /* nthreads + 1 range bounds */
    double **red;         /* 4*ninst: where each matrix entry is added (entry or slot) */
    int *rhsred;          /* 2*ninst: node index, or -(slot+1) */
    int nshared, nrshared;
    double **targets;     /* nshared: matrix entries behind the slots */
    int *rtargets;        /* nrshared: RHS nodes behind the RHS slots */
    double *slots;        /* nthreads x nshared */
    double *rslots;       /* nthreads x nrshared */
} cst;

void CAPStampTablesInvalidate(void) { cst.valid = 0; }

static double cl_now(void)
{
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return ts.tv_sec + 1e-9 * ts.tv_nsec;
}

static void cl_report(void)
{
    if (cl_calls)
        fprintf(stderr, "LOADTRACE cap: calls=%ld threads=%d inst=%d shared entries=%d shared rhs=%d eval+stamp=%.3fs (%.1fus/call) reduce=%.3fs (%.1fus/call)\n",
                cl_calls, cst.nthreads, cst.ninst, cst.nshared, cst.nrshared,
                cl_eval, 1e6 * cl_eval / cl_calls, cl_reduce, 1e6 * cl_reduce / cl_calls);
}

static inline double *cap_ptr(CAPinstance *h, int k)
{
    switch (k) {
    case 0: return h->CAPposPosPtr;
    case 1: return h->CAPnegNegPtr;
    case 2: return h->CAPposNegPtr;
    default: return h->CAPnegPosPtr;
    }
}
static inline int cap_node(CAPinstance *h, int k) { return k ? h->CAPnegNode : h->CAPposNode; }

/* ownership pass over a set of keys: ins(key) -> owner thread or -1 when shared */
static void
cap_tables_build(CAPmodel *model, CKTcircuit *ckt, int T)
{
    CAPinstance **inst = model->CAPinstArray;
    int n = model->CAPinstCount, i, k, t, rhs_n;
    NGHASHPTR owner, slot;
    int *rowner, *rslot;

    free(cst.lo); free(cst.red); free(cst.rhsred); free(cst.targets); free(cst.rtargets); free(cst.slots); free(cst.rslots);
    memset(&cst, 0, sizeof cst);
    cst.valid = 1; cst.nthreads = T; cst.ninst = n;
    if (T < 2 || n == 0)
        return;

    cst.lo = TMALLOC(int, T + 1);
    for (t = 0; t <= T; t++)
        cst.lo[t] = (int)((long) n * t / T);
    cst.red = TMALLOC(double *, 4 * (size_t) n);
    cst.rhsred = TMALLOC(int, 2 * (size_t) n);

    /* matrix entries: owner thread of each pointer, -1 when shared */
    owner = nghash_init_pointer(4 * n);
    for (t = 0; t < T; t++)
        for (i = cst.lo[t]; i < cst.lo[t + 1]; i++)
            for (k = 0; k < 4; k++) {
                double *p = cap_ptr(inst[i], k);
                void *v = nghash_find(owner, p);
                if (!v)
                    nghash_insert(owner, p, (void *)(intptr_t)(t + 1));
                else if ((intptr_t) v != t + 1 && (intptr_t) v != -1) {
                    nghash_delete(owner, p);
                    nghash_insert(owner, p, (void *)(intptr_t) -1);
                }
            }
    slot = nghash_init_pointer(4 * n);
    {
        int cap = 16;
        cst.targets = TMALLOC(double *, cap);
        for (i = 0; i < n; i++)
            for (k = 0; k < 4; k++) {
                double *p = cap_ptr(inst[i], k);
                void *v;
                if ((intptr_t) nghash_find(owner, p) != -1) {
                    cst.red[4 * i + k] = p;
                    continue;
                }
                v = nghash_find(slot, p);
                if (!v) {
                    if (cst.nshared == cap) { cap *= 2; cst.targets = TREALLOC(double *, cst.targets, cap); }
                    cst.targets[cst.nshared] = p;
                    v = (void *)(intptr_t)(cst.nshared + 1);
                    nghash_insert(slot, p, v);
                    cst.nshared++;
                }
                cst.red[4 * i + k] = (double *)(intptr_t) v;   /* slot code, resolved below */
            }
    }
    cst.slots = TMALLOC(double, (size_t) T * (cst.nshared ? cst.nshared : 1));
    for (t = 0; t < T; t++)
        for (i = cst.lo[t]; i < cst.lo[t + 1]; i++)
            for (k = 0; k < 4; k++)
                if ((intptr_t) nghash_find(owner, cap_ptr(inst[i], k)) == -1) {
                    intptr_t code = (intptr_t) cst.red[4 * i + k] - 1;
                    cst.red[4 * i + k] = &cst.slots[(size_t) t * cst.nshared + code];
                }
    nghash_free(owner, NULL, NULL);
    nghash_free(slot, NULL, NULL);

    /* RHS nodes: same thing with plain arrays indexed by node */
    rhs_n = SMPmatSize(ckt->CKTmatrix) + 1;
    rowner = TMALLOC(int, rhs_n);
    rslot = TMALLOC(int, rhs_n);
    for (i = 0; i < rhs_n; i++) { rowner[i] = 0; rslot[i] = 0; }
    for (t = 0; t < T; t++)
        for (i = cst.lo[t]; i < cst.lo[t + 1]; i++)
            for (k = 0; k < 2; k++) {
                int nd = cap_node(inst[i], k);
                if (nd <= 0 || nd >= rhs_n) continue;          /* ground: never stamped */
                if (rowner[nd] == 0) rowner[nd] = t + 1;
                else if (rowner[nd] != t + 1) rowner[nd] = -1;
            }
    {
        int cap = 16;
        cst.rtargets = TMALLOC(int, cap);
        for (i = 0; i < n; i++)
            for (k = 0; k < 2; k++) {
                int nd = cap_node(inst[i], k);
                if (nd <= 0 || nd >= rhs_n || rowner[nd] != -1) {
                    cst.rhsred[2 * i + k] = nd;
                    continue;
                }
                if (rslot[nd] == 0) {
                    if (cst.nrshared == cap) { cap *= 2; cst.rtargets = TREALLOC(int, cst.rtargets, cap); }
                    cst.rtargets[cst.nrshared] = nd;
                    rslot[nd] = ++cst.nrshared;
                }
                cst.rhsred[2 * i + k] = -rslot[nd];
            }
    }
    cst.rslots = TMALLOC(double, (size_t) T * (cst.nrshared ? cst.nrshared : 1));
    free(rowner); free(rslot);
    if (cl_trace)
        fprintf(stderr, "LOADTRACE cap stamp tables: inst=%d threads=%d shared entries=%d shared rhs nodes=%d\n",
                n, T, cst.nshared, cst.nrshared);
}

/* evaluate one instance: charge state and integration -> m*geq, m*ceq */
static inline int
CAPevalOne(CAPinstance *here, CKTcircuit *ckt, int cond1, double *geqM, double *ceqM)
{
    double vcap, geq, ceq;
    int error;

    if(cond1) {
        vcap = here->CAPinitCond;
    } else {
        vcap = *(ckt->CKTrhsOld+here->CAPposNode) -
               *(ckt->CKTrhsOld+here->CAPnegNode) ;
    }
    if(ckt->CKTmode & (MODETRAN | MODEAC)) {
#ifndef PREDICTOR
        if(ckt->CKTmode & MODEINITPRED) {
            *(ckt->CKTstate0+here->CAPqcap) =
                *(ckt->CKTstate1+here->CAPqcap);
        } else { /* only const caps - no poly's */
#endif /* PREDICTOR */
            *(ckt->CKTstate0+here->CAPqcap) = here->CAPcapac * vcap;
            if((ckt->CKTmode & MODEINITTRAN)) {
                *(ckt->CKTstate1+here->CAPqcap) =
                    *(ckt->CKTstate0+here->CAPqcap);
            }
#ifndef PREDICTOR
        }
#endif /* PREDICTOR */
        error = NIintegrate(ckt,&geq,&ceq,here->CAPcapac,
                            here->CAPqcap);
        if(error) return(error);
        if(ckt->CKTmode & MODEINITTRAN) {
            *(ckt->CKTstate1+here->CAPccap) =
                *(ckt->CKTstate0+here->CAPccap);
        }
        *geqM = here->CAPm * geq;
        *ceqM = here->CAPm * ceq;
    } else {
        *(ckt->CKTstate0+here->CAPqcap) = here->CAPcapac * vcap;
        *geqM = 0.0;
        *ceqM = 0.0;
    }
    return(OK);
}

struct CAPfuseJob { CAPinstance **inst; CKTcircuit *ckt; int cond1, T; int err[64]; };

static int
CAPfuseItem(int t, void *arg)
{
    struct CAPfuseJob *j = (struct CAPfuseJob *) arg;
    CKTcircuit *ckt = j->ckt;
    double *rhs = ckt->CKTrhs;
    double *sl = cst.slots + (size_t) t * cst.nshared;
    double *rsl = cst.rslots + (size_t) t * cst.nrshared;
    int i, e = 0;

    memset(sl, 0, sizeof(double) * (size_t) cst.nshared);
    memset(rsl, 0, sizeof(double) * (size_t) cst.nrshared);
    for (i = cst.lo[t]; i < cst.lo[t + 1]; i++) {
        CAPinstance *here = j->inst[i];
        double **red = cst.red + 4 * (size_t) i;
        int *rr = cst.rhsred + 2 * (size_t) i;
        double geqM, ceqM;
        e = CAPevalOne(here, ckt, j->cond1, &geqM, &ceqM);
        if (e) break;
        *red[0] += geqM;
        *red[1] += geqM;
        *red[2] -= geqM;
        *red[3] -= geqM;
        if (rr[0] >= 0) rhs[rr[0]] -= ceqM; else rsl[-rr[0] - 1] -= ceqM;
        if (rr[1] >= 0) rhs[rr[1]] += ceqM; else rsl[-rr[1] - 1] += ceqM;
    }
    j->err[t] = e;
    return 0;
}

struct CAPreduceJob { CKTcircuit *ckt; int T; };

static int
CAPreduceItem(int r, void *arg)
{
    struct CAPreduceJob *j = (struct CAPreduceJob *) arg;
    int T = j->T, t, k;
    int k0 = (int)((long) cst.nshared * r / T), k1 = (int)((long) cst.nshared * (r + 1) / T);
    int i0 = (int)((long) cst.nrshared * r / T), i1 = (int)((long) cst.nrshared * (r + 1) / T);
    double *rhs = j->ckt->CKTrhs;
    for (t = 0; t < T; t++) {
        const double *sl = cst.slots + (size_t) t * cst.nshared;
        const double *rsl = cst.rslots + (size_t) t * cst.nrshared;
        for (k = k0; k < k1; k++)
            *cst.targets[k] += sl[k];
        for (k = i0; k < i1; k++)
            rhs[cst.rtargets[k]] += rsl[k];
    }
    return 0;
}

static inline void
CAPstampOne(CAPinstance *here, CKTcircuit *ckt, double geqM, double ceqM)
{
    *(here->CAPposPosPtr) += geqM;
    *(here->CAPnegNegPtr) += geqM;
    *(here->CAPposNegPtr) -= geqM;
    *(here->CAPnegPosPtr) -= geqM;
    *(ckt->CKTrhs+here->CAPposNode) -= ceqM;
    *(ckt->CKTrhs+here->CAPnegNode) += ceqM;
}

int
CAPload(GENmodel *inModel, CKTcircuit *ckt)
/* actually load the current capacitance value into the
 * sparse matrix previously provided
 */
{
    CAPmodel *model = (CAPmodel*)inModel;
    CAPinstance *here;
    int cond1;
    int error;

    /* check if capacitors are in the circuit or are open circuited */
    if(ckt->CKTmode & (MODETRAN|MODEAC|MODETRANOP) ) {
        /* evaluate device independent analysis conditions */
        cond1=
            ( ( (ckt->CKTmode & MODEDC) &&
                (ckt->CKTmode & MODEINITJCT) )
              || ( ( ckt->CKTmode & MODEUIC) &&
                   ( ckt->CKTmode & MODEINITTRAN) ) ) ;

        if (cl_on < 0) {
            cl_on = getenv("NGSPICE_NOPCAP") == NULL;
            cl_trace = getenv("NGSPICE_LOADTRACE") != NULL;
            if (cl_trace) atexit(cl_report);
        }
        {
            int T = ngpool_threads();
            if (cl_on && T > 1 && T <= 64 && model->CAPinstArray &&
                model->CAPinstCount >= ckt->CKTompMinInst) {
                struct CAPfuseJob j;
                struct CAPreduceJob rj;
                double t0 = cl_trace ? cl_now() : 0, t1;
                int t;
                if (!cst.valid || cst.nthreads != T || cst.ninst != model->CAPinstCount)
                    cap_tables_build(model, ckt, T);
                j.inst = model->CAPinstArray; j.ckt = ckt; j.cond1 = cond1; j.T = T;
                ngpool_run_chunk(T, CAPfuseItem, &j, 1);
                for (t = 0; t < T; t++)
                    if (j.err[t]) return(j.err[t]);
                if (cl_trace) { t1 = cl_now(); cl_eval += t1 - t0; t0 = t1; }
                rj.ckt = ckt; rj.T = T;
                ngpool_run_chunk(T, CAPreduceItem, &rj, 1);
                if (cl_trace) { cl_calls++; cl_reduce += cl_now() - t0; }
                return(OK);
            }
        }

        /*  loop through all the capacitor models */
        for( ; model != NULL; model = CAPnextModel(model)) {

            /* loop through all the instances of the model */
            for (here = CAPinstances(model); here != NULL ;
                    here=CAPnextInstance(here)) {
                double geqM, ceqM;
                error = CAPevalOne(here, ckt, cond1, &geqM, &ceqM);
                if(error) return(error);
                CAPstampOne(here, ckt, geqM, ceqM);
            }
        }
    }
    return(OK);
}
