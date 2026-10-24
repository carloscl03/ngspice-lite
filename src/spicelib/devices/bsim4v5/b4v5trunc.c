/**** BSIM4.5.0 Released by Xuemei (Jane) Xi 07/29/2005 ****/

/**********
 * Copyright 2004 Regents of the University of California. All rights reserved.
 * File: b4trunc.c of BSIM4.5.0.
 * Author: 2000 Weidong Liu
 * Authors: 2001- Xuemei Xi, Mohan Dunga, Ali Niknejad, Chenming Hu.
 * Project Director: Prof. Chenming Hu.
 **********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "bsim4v5def.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"


#include "ngspice/ngpool.h"

/* fork-lite: the truncation-error pass reads the charge states of every
   instance (another walk over the 4.5 kB structs). On the pool each thread
   takes a range and keeps its own minimum; the minimum is order
   independent, so the result is exactly the serial one. */
struct BSIM4v5TruncJob { BSIM4v5instance **inst; CKTcircuit *ckt; int n, T; double ts[64]; };

static void
BSIM4v5truncOne(BSIM4v5instance *here, CKTcircuit *ckt, double *timeStep)
{
    CKTterr(here->BSIM4v5qb,ckt,timeStep);
    CKTterr(here->BSIM4v5qg,ckt,timeStep);
    CKTterr(here->BSIM4v5qd,ckt,timeStep);
    if (here->BSIM4v5trnqsMod)
        CKTterr(here->BSIM4v5qcdump,ckt,timeStep);
    if (here->BSIM4v5rbodyMod)
    {   CKTterr(here->BSIM4v5qbs,ckt,timeStep);
        CKTterr(here->BSIM4v5qbd,ckt,timeStep);
    }
    if (here->BSIM4v5rgateMod == 3)
        CKTterr(here->BSIM4v5qgmid,ckt,timeStep);
}

static int
BSIM4v5TruncItem(int t, void *arg)
{
    struct BSIM4v5TruncJob *j = (struct BSIM4v5TruncJob *) arg;
    int lo = (int)((long) j->n * t / j->T), hi = (int)((long) j->n * (t + 1) / j->T), i;
    double ts = j->ts[t];
    for (i = lo; i < hi; i++)
        BSIM4v5truncOne(j->inst[i], j->ckt, &ts);
    j->ts[t] = ts;
    return 0;
}

int
BSIM4v5trunc(
GENmodel *inModel,
CKTcircuit *ckt,
double *timeStep)
{
BSIM4v5model *model = (BSIM4v5model*)inModel;
BSIM4v5instance *here;
#ifdef USE_OMP
    {
        static int ptrunc = -1;
        int T = ngpool_threads();
        if (ptrunc < 0) ptrunc = getenv("NGSPICE_NOPTRUNC") == NULL;
        if (ptrunc && T > 1 && T <= 64 && model->BSIM4v5InstanceArray && model->BSIM4v5InstCount >= ckt->CKTompMinInst) {
            struct BSIM4v5TruncJob j;
            int t;
            j.inst = model->BSIM4v5InstanceArray; j.ckt = ckt; j.n = model->BSIM4v5InstCount; j.T = T;
            for (t = 0; t < T; t++) j.ts[t] = *timeStep;
            ngpool_run_chunk(T, BSIM4v5TruncItem, &j, 1);
            for (t = 0; t < T; t++)
                if (j.ts[t] < *timeStep) *timeStep = j.ts[t];
            return(OK);
        }
    }
#endif

#ifdef STEPDEBUG
    double debugtemp;
#endif /* STEPDEBUG */

    for (; model != NULL; model = BSIM4v5nextModel(model))
    {    for (here = BSIM4v5instances(model); here != NULL;
	      here = BSIM4v5nextInstance(here))
	 {

#ifdef STEPDEBUG
            debugtemp = *timeStep;
#endif /* STEPDEBUG */
            CKTterr(here->BSIM4v5qb,ckt,timeStep);
            CKTterr(here->BSIM4v5qg,ckt,timeStep);
            CKTterr(here->BSIM4v5qd,ckt,timeStep);
            if (here->BSIM4v5trnqsMod)
                CKTterr(here->BSIM4v5qcdump,ckt,timeStep);
            if (here->BSIM4v5rbodyMod)
            {   CKTterr(here->BSIM4v5qbs,ckt,timeStep);
                CKTterr(here->BSIM4v5qbd,ckt,timeStep);
	    }
	    if (here->BSIM4v5rgateMod == 3)
		CKTterr(here->BSIM4v5qgmid,ckt,timeStep);
#ifdef STEPDEBUG
            if(debugtemp != *timeStep)
	    {  printf("device %s reduces step from %g to %g\n",
                       here->BSIM4v5name,debugtemp,*timeStep);
            }
#endif /* STEPDEBUG */
        }
    }
    return(OK);
}
