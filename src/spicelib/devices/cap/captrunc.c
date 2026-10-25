/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: September 2003 Paolo Nenzi
**********/
/*
 */

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "capdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/suffix.h"
#include "ngspice/ngpool.h"

/* fork-lite: per-thread minimum over a range of instances; the minimum does
   not depend on the order, so the result equals the serial loop exactly */
struct CAPtruncJob { CAPinstance **inst; CKTcircuit *ckt; int n, T; double ts[64]; };

static int
CAPtruncItem(int t, void *arg)
{
    struct CAPtruncJob *j = (struct CAPtruncJob *) arg;
    int lo = (int)((long) j->n * t / j->T), hi = (int)((long) j->n * (t + 1) / j->T), i;
    double ts = j->ts[t];
    for (i = lo; i < hi; i++)
        CKTterr(j->inst[i]->CAPqcap, j->ckt, &ts);
    j->ts[t] = ts;
    return 0;
}

int
CAPtrunc(GENmodel *inModel, CKTcircuit *ckt, double *timeStep)
{
    CAPmodel *model = (CAPmodel*)inModel;
    CAPinstance *here;
    static int ptrunc = -1;
    int T = ngpool_threads();

    if (ptrunc < 0) ptrunc = getenv("NGSPICE_NOPCAP") == NULL;
    if (ptrunc && T > 1 && T <= 64 && model->CAPinstArray &&
        model->CAPinstCount >= ckt->CKTompMinInst) {
        struct CAPtruncJob j;
        int t;
        j.inst = model->CAPinstArray; j.ckt = ckt; j.n = model->CAPinstCount; j.T = T;
        for (t = 0; t < T; t++) j.ts[t] = *timeStep;
        ngpool_run_chunk(T, CAPtruncItem, &j, 1);
        for (t = 0; t < T; t++)
            if (j.ts[t] < *timeStep) *timeStep = j.ts[t];
        return(OK);
    }

    for( ; model!= NULL; model = CAPnextModel(model)) {
        for(here = CAPinstances(model); here != NULL ;
                here = CAPnextInstance(here)) {

            CKTterr(here->CAPqcap,ckt,timeStep);
        }
    }
    return(OK);
}
