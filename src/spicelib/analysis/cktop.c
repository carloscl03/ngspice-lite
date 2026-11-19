/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2000  AlansFixes
Modified: 2005 Paolo Nenzi - Restructured
**********/

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "ngspice/devdefs.h"
#include "vsrc/vsrcdefs.h"
#include "dev.h"
#include "ngspice/sperror.h"
#include "ngspice/cpextern.h"
#include "ngspice/fteext.h"

#ifdef XSPICE
#include "ngspice/enh.h"
#endif

static int dynamic_gmin(CKTcircuit *, long int, long int, int);
static int spice3_gmin(CKTcircuit *, long int, long int, int);
static int new_gmin(CKTcircuit*, long int, long int, int);
static int gillespie_src(CKTcircuit *, long int, long int, int);
static int spice3_src(CKTcircuit *, long int, long int, int);


int
CKTop (CKTcircuit *ckt, long int firstmode, long int continuemode,
       int iterlim)
{
    int converged;

#ifdef HAS_PROGREP
    /* If this is called from dc simulation, don't set "op" */
    if (ckt->CKTcurJob->JOBtype != 2 && (!ckt->CKTnoOpIter || ckt->CKTnumGminSteps >= 1 || ckt->CKTnumSrcSteps >= 1))
        SetAnalyse("op", 0);
#endif

    ckt->CKTmode = firstmode;

    if (!ckt->CKTnoOpIter) {
#ifdef XSPICE
        /* gtri - wbk - add convergence problem reporting flags */
        ckt->enh->conv_debug.last_NIiter_call =
            (ckt->CKTnumGminSteps <= 0) && (ckt->CKTnumSrcSteps <= 0);
#endif
        converged = NIiter (ckt, iterlim);
        if (converged == 0)
            return converged;   /* successfull */
    } else {
        converged = 1;          /* the 'go directly to gmin stepping' option */
    }


    /* no convergence on the first try, so we do something else */
    /* first, check if we should try gmin stepping */

    if (ckt->CKTnumGminSteps >= 1) {
        if (ckt->CKTnumGminSteps == 1) {
            /* only the old gmin */
            if (cp_getvar("dyngmin", CP_BOOL, NULL, 0)) {
                converged = dynamic_gmin(ckt, firstmode, continuemode, iterlim);
            }
            /* first the old gmin, then the new gmin */
            else {
                converged = dynamic_gmin(ckt, firstmode, continuemode, iterlim);
                if(converged != 0) {
                    converged = new_gmin(ckt, firstmode, continuemode, iterlim);
                }

            }
        }
        else {
            converged = spice3_gmin(ckt, firstmode, continuemode, iterlim);
        }
        if (converged == 0) /* If gmin-stepping worked... move out */
            return converged;
    }

    /* ... otherwise try stepping sources ...
     * now, we'll try source stepping - we scale the sources
     * to 0, converge, then start stepping them up until they
     * are at their normal values
     */

    if (ckt->CKTnumSrcSteps >= 1) {
        if (ckt->CKTnumSrcSteps == 1)
            converged = gillespie_src(ckt, firstmode, continuemode, iterlim);
        else
            converged = spice3_src(ckt, firstmode, continuemode, iterlim);
        if (converged == 0) /* If gmin-stepping worked... move out */
            return converged;
    }

    /* If command 'optran' is not given, the function
       returns immediately with the previous 'converged' */
    int prevconverged = converged;
    converged = OPtran(ckt, converged);
    if (converged == 106)
        fprintf(cp_err, "Error: Transient op failed, timestep too small\n\n");
    else if (converged != 0 && converged != prevconverged)
        fprintf(cp_err, "Error: Transient op failed, cause unrecorded\n\n");
    else if (converged == 0)
            return converged;

#ifdef XSPICE
    /* gtri - wbk - add convergence problem reporting flags */
    ckt->enh->conv_debug.last_NIiter_call = MIF_FALSE;
#endif

    fprintf(cp_err, "\nError: The operating point could not be simulated successfully.\n");
    if (ft_stricterror)
        controlled_exit(1);
    fprintf(cp_err, "    Any of the following steps may fail.!\n\n");

    return converged;
}


/* CKTconvTest(ckt)
 *    this is a driver program to iterate through all the various
 *    convTest functions provided for the circuit elements in the
 *    given circuit
 */

int
CKTconvTest (CKTcircuit *ckt)
{
    int i;

    for (i = 0; i < DEVmaxnum; i++) {

        if (DEVices[i] && DEVices[i]->DEVconvTest && ckt->CKThead[i]) {
            int error = DEVices[i]->DEVconvTest (ckt->CKThead[i], ckt);
            if (error)
                return error;
        }

        if (ckt->CKTnoncon) {
            /* printf("convTest: device %s failed\n",
             * DEVices[i]->DEVpublic.name); */
            return OK;
        }
    }

    return OK;
}


/* Dynamic gmin stepping
 * Algorithm by Alan Gillespie
 * Modified 2005 - Paolo Nenzi (extracted from CKTop.c code)
 *
 * return value:
 * 0 -> method converged
 * 1 -> method failed
 *
 * Note that no path out of this code allows ckt->CKTdiagGmin to be
 * anything but CKTgshunt.
 */

static int
dynamic_gmin (CKTcircuit *ckt, long int firstmode,
              long int continuemode, int iterlim)
{
    double OldGmin, gtarget, factor;
    int converged;

    int NumNodes, iters, i;
    double *OldRhsOld, *OldCKTstate0;
    CKTnode *n;

    ckt->CKTmode = firstmode;
    SPfrontEnd->IFerrorf (ERR_INFO, "Starting dynamic gmin stepping");

    NumNodes = 0;
    for (n = ckt->CKTnodes; n; n = n->next)
        NumNodes++;

    OldRhsOld = TMALLOC(double, NumNodes + 1);
    OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1);

    for (n = ckt->CKTnodes; n; n = n->next)
        ckt->CKTrhsOld [n->number] = 0;

    for (i = 0; i < ckt->CKTnumStates; i++)
        ckt->CKTstate0 [i] = 0;

    factor = ckt->CKTgminFactor;
    OldGmin = 1e-2;
    ckt->CKTdiagGmin = OldGmin / factor;
    gtarget = MAX (ckt->CKTgmin, ckt->CKTgshunt);

    for (;;) {
        if (ft_ngdebug)
            fprintf (stderr, "Trying gmin = %12.4E ", ckt->CKTdiagGmin);

        ckt->CKTnoncon = 1;
        iters = ckt->CKTstat->STATnumIter;
        converged = NIiter (ckt, ckt->CKTdcTrcvMaxIter);
        iters = ckt->CKTstat->STATnumIter - iters;

        if (converged == 0) {
            ckt->CKTmode = continuemode;
            if (ft_ngdebug)
                SPfrontEnd->IFerrorf (ERR_INFO, "One successful gmin step");

            if (ckt->CKTdiagGmin <= gtarget)
                break;          /* successfull */

            for (i = 0, n = ckt->CKTnodes; n; n = n->next)
                OldRhsOld[i++] = ckt->CKTrhsOld[n->number];

            memcpy(OldCKTstate0, ckt->CKTstate0,
                   (size_t) ckt->CKTnumStates * sizeof(double));

            if (iters <= (ckt->CKTdcTrcvMaxIter / 4)) {
                factor *= sqrt (factor);
                if (factor > ckt->CKTgminFactor)
                    factor = ckt->CKTgminFactor;
            }

            if (iters > (3 * ckt->CKTdcTrcvMaxIter / 4))
                factor = MAX(sqrt (factor), 1.00005);

            OldGmin = ckt->CKTdiagGmin;

            if (ckt->CKTdiagGmin < factor * gtarget) {
                factor = ckt->CKTdiagGmin / gtarget;
                ckt->CKTdiagGmin = gtarget;
            } else {
                ckt->CKTdiagGmin /= factor;
            }
        } else {
            if (factor < 1.00005) {
                if (ft_ngdebug)
                    SPfrontEnd->IFerrorf (ERR_WARNING, "Last gmin step failed");
                break;          /* failed */
            }
            if (ft_ngdebug)
                SPfrontEnd->IFerrorf (ERR_WARNING, "Further gmin increment");
            factor = sqrt (sqrt (factor));
            ckt->CKTdiagGmin = OldGmin / factor;

            for (i = 0, n = ckt->CKTnodes; n; n = n->next)
                ckt->CKTrhsOld[n->number] = OldRhsOld[i++];

            memcpy(ckt->CKTstate0, OldCKTstate0,
                   (size_t) ckt->CKTnumStates * sizeof(double));
        }
    }

    ckt->CKTdiagGmin = ckt->CKTgshunt;
    FREE (OldRhsOld);
    FREE (OldCKTstate0);

#ifdef XSPICE
    /* gtri - wbk - add convergence problem reporting flags */
    ckt->enh->conv_debug.last_NIiter_call = (ckt->CKTnumSrcSteps <= 0);
#endif

    converged = NIiter (ckt, iterlim);

    if (converged != 0) {
        SPfrontEnd->IFerrorf (ERR_WARNING, "Dynamic gmin stepping failed");
    } else {
        SPfrontEnd->IFerrorf (ERR_INFO, "Dynamic gmin stepping completed");
#ifdef XSPICE
        /* gtri - wbk - add convergence problem reporting flags */
        ckt->enh->conv_debug.last_NIiter_call = MIF_FALSE;
#endif
    }

    return converged;
}


/* Spice3 gmin stepping
 * Modified 2000 - Alan Gillespie (added gshunt)
 * Modified 2005 - Paolo Nenzi (extracted from CKTop.c code)
 *
 * return value:
 * 0 -> method converged
 * 1 -> method failed
 *
 * Note that no path out of this code allows ckt->CKTdiagGmin to be
 * anything but CKTgshunt.
 */

static int
spice3_gmin (CKTcircuit *ckt, long int firstmode,
             long int continuemode, int iterlim)
{
    int converged, i;

    ckt->CKTmode = firstmode;
    SPfrontEnd->IFerrorf (ERR_INFO, "Starting spice3 gmin stepping");

    ckt->CKTdiagGmin =
        (ckt->CKTgshunt == 0) ? ckt->CKTgmin :  ckt->CKTgshunt;

    for (i = 0; i < ckt->CKTnumGminSteps; i++)
        ckt->CKTdiagGmin *= ckt->CKTgminFactor;

    for (i = 0; i <= ckt->CKTnumGminSteps; i++) {
        if (ft_ngdebug)
            fprintf (stderr, "Trying gmin = %12.4E ", ckt->CKTdiagGmin);

        ckt->CKTnoncon = 1;
        converged = NIiter (ckt, ckt->CKTdcTrcvMaxIter);

        if (converged != 0) {
            ckt->CKTdiagGmin = ckt->CKTgshunt;
            if (ft_ngdebug)
                SPfrontEnd->IFerrorf (ERR_WARNING, "gmin step failed");
            break;
        }

        ckt->CKTdiagGmin /= ckt->CKTgminFactor;
        ckt->CKTmode = continuemode;
        if (ft_ngdebug)
            SPfrontEnd->IFerrorf (ERR_INFO, "One successful gmin step");
    }

    ckt->CKTdiagGmin = ckt->CKTgshunt;

#ifdef XSPICE
    /* gtri - wbk - add convergence problem reporting flags */
    ckt->enh->conv_debug.last_NIiter_call = (ckt->CKTnumSrcSteps <= 0);
#endif

    converged = NIiter (ckt, iterlim);

    if (converged == 0) {
        SPfrontEnd->IFerrorf (ERR_INFO, "spice3 gmin stepping completed");

#ifdef XSPICE
        /* gtri - wbk - add convergence problem reporting flags */
        ckt->enh->conv_debug.last_NIiter_call = MIF_FALSE;
#endif

    } else {
        SPfrontEnd->IFerrorf (ERR_WARNING, "spice3 gmin stepping failed");
    }

    return converged;
}

/* just step the real gmin found in every device model */
static int
new_gmin(CKTcircuit* ckt, long int firstmode,
    long int continuemode, int iterlim)
{
    double OldGmin, gtarget, factor, startgmin;
    int converged;

    int NumNodes, iters, i;
    double* OldRhsOld, * OldCKTstate0;
    CKTnode* n;

    ckt->CKTmode = firstmode;
    SPfrontEnd->IFerrorf(ERR_INFO, "Starting true gmin stepping");

    NumNodes = 0;
    for (n = ckt->CKTnodes; n; n = n->next)
        NumNodes++;

    OldRhsOld = TMALLOC(double, NumNodes + 1);
    OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1);

    for (n = ckt->CKTnodes; n; n = n->next)
        ckt->CKTrhsOld[n->number] = 0;

    for (i = 0; i < ckt->CKTnumStates; i++)
        ckt->CKTstate0[i] = 0;

    startgmin = ckt->CKTgmin;
    factor = ckt->CKTgminFactor;
    OldGmin = 1e-2;
    /*ckt->CKTdiagGmin = */ckt->CKTgmin = OldGmin / factor;
    gtarget = MAX(startgmin, ckt->CKTgshunt);

    for (;;) {
        if (ft_ngdebug)
            fprintf(stderr, "Trying gmin = %12.4E ", ckt->CKTgmin);

        ckt->CKTnoncon = 1;
        iters = ckt->CKTstat->STATnumIter;
        converged = NIiter(ckt, ckt->CKTdcTrcvMaxIter);
        iters = ckt->CKTstat->STATnumIter - iters;

        if (converged == 0) {
            ckt->CKTmode = continuemode;
            if (ft_ngdebug)
                SPfrontEnd->IFerrorf(ERR_INFO, "One successful gmin step");

            if (ckt->CKTgmin <= gtarget)
                break;          /* successfull */

            for (i = 0, n = ckt->CKTnodes; n; n = n->next)
                OldRhsOld[i++] = ckt->CKTrhsOld[n->number];

            memcpy(OldCKTstate0, ckt->CKTstate0,
                (size_t)ckt->CKTnumStates * sizeof(double));

            if (iters <= (ckt->CKTdcTrcvMaxIter / 4)) {
                factor *= sqrt(factor);
                if (factor > ckt->CKTgminFactor)
                    factor = ckt->CKTgminFactor;
            }

            if (iters > (3 * ckt->CKTdcTrcvMaxIter / 4))
                factor = MAX(sqrt(factor), 3);

            OldGmin = ckt->CKTgmin;

            if (ckt->CKTgmin < factor * gtarget) {
                factor = ckt->CKTgmin / gtarget;
                /*ckt->CKTdiagGmin = */ckt->CKTgmin = gtarget;
            }
            else {
                /*ckt->CKTdiagGmin = */ckt->CKTgmin /= factor;
            }
        }
        else {
            if (factor < 1.00005) {
                if (ft_ngdebug)
                    SPfrontEnd->IFerrorf(ERR_WARNING, "Last gmin step failed");
                break;          /* failed */
            }
            if (ft_ngdebug)
                SPfrontEnd->IFerrorf(ERR_WARNING, "Further gmin increment");
            factor = sqrt(sqrt(factor));
            /*ckt->CKTdiagGmin = */ckt->CKTgmin = OldGmin / factor;

            for (i = 0, n = ckt->CKTnodes; n; n = n->next)
                ckt->CKTrhsOld[n->number] = OldRhsOld[i++];

            memcpy(ckt->CKTstate0, OldCKTstate0,
                (size_t)ckt->CKTnumStates * sizeof(double));
        }
    }

    /*ckt->CKTdiagGmin = */ckt->CKTgmin = MAX(startgmin, ckt->CKTgshunt);
    FREE(OldRhsOld);
    FREE(OldCKTstate0);

#ifdef XSPICE
    /* gtri - wbk - add convergence problem reporting flags */
    ckt->enh->conv_debug.last_NIiter_call = (ckt->CKTnumSrcSteps <= 0);
#endif

    converged = NIiter(ckt, iterlim);

    if (converged != 0) {
        SPfrontEnd->IFerrorf(ERR_WARNING, "True gmin stepping failed");
    }
    else {
        SPfrontEnd->IFerrorf(ERR_INFO, "True gmin stepping completed");
#ifdef XSPICE
        /* gtri - wbk - add convergence problem reporting flags */
        ckt->enh->conv_debug.last_NIiter_call = MIF_FALSE;
#endif
    }

    return converged;
}


/* Gillespie's Source stepping
 * Modified 2005 - Paolo Nenzi (extracted from CKTop.c code)
 *
 * return value:
 * 0 -> method converged
 * 1 -> method failed
 *
 * Note that no path out of this code allows ckt->CKTsrcFact to be
 * anything but 1.00000.
 */

static int
gillespie_src (CKTcircuit *ckt, long int firstmode,
               long int continuemode, int iterlim)
{
    int converged, i, iters;
    double ConvFact;
    CKTnode *n;
    double gminstart = ckt->CKTgmin;

    NG_IGNORE(iterlim);

    ckt->CKTmode = firstmode;
    SPfrontEnd->IFerrorf (ERR_INFO, "Starting source stepping");

    ckt->CKTsrcFact = 0;
    ConvFact = 0;

    for (n = ckt->CKTnodes; n; n = n->next)
        ckt->CKTrhsOld[n->number] = 0;

    for (i = 0; i < ckt->CKTnumStates; i++)
        ckt->CKTstate0[i] = 0;

    /*  First, try a straight solution with all sources at zero */
    if (ft_ngdebug)
        fprintf (stderr, "Supplies reduced to %8.4f%% ", ckt->CKTsrcFact * 100);
    converged = NIiter (ckt, ckt->CKTdcTrcvMaxIter);

    /*  If this doesn't work, try gmin stepping as well for the first solution */

    if (converged != 0) {
        if (ft_ngdebug)
            fprintf (stderr, "\n");

        ckt->CKTdiagGmin =
            (ckt->CKTgshunt <= 0) ? ckt->CKTgmin : ckt->CKTgshunt;

        for (i = 0; i < 10; i++)
            ckt->CKTdiagGmin *= 10;

        for (i = 0; i <= 10; i++) {
            if (ft_ngdebug)
                fprintf (stderr, "Trying gmin = %12.4E ", ckt->CKTdiagGmin);

#ifdef XSPICE
            /* gtri - wbk - add convergence problem reporting flags */
            ckt->enh->conv_debug.last_NIiter_call = MIF_TRUE;
#endif

            ckt->CKTnoncon = 1;
            converged = NIiter (ckt, ckt->CKTdcTrcvMaxIter);

            if (converged != 0) {
                ckt->CKTdiagGmin = ckt->CKTgshunt;
                if (ft_ngdebug)
                    SPfrontEnd->IFerrorf (ERR_WARNING, "gmin step failed");
#ifdef XSPICE
                /* gtri - wbk - add convergence problem reporting flags */
                ckt->enh->conv_debug.last_NIiter_call = MIF_FALSE;
#endif
                break;
            }

            ckt->CKTdiagGmin /= 10;
            ckt->CKTmode = continuemode;
            if (ft_ngdebug)
                SPfrontEnd->IFerrorf (ERR_INFO, "One successful gmin step");
        }
        ckt->CKTdiagGmin = ckt->CKTgshunt;
    }

    /*  If we've got convergence, then try stepping up the sources  */

    if (converged == 0) {

        double *OldRhsOld, *OldCKTstate0, raise = 0.001;

        int NumNodes = 0;
        for (n = ckt->CKTnodes; n; n = n->next)
            NumNodes++;

        OldRhsOld = TMALLOC(double, NumNodes + 1);
        OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1);

        for (i = 0, n = ckt->CKTnodes; n; n = n->next)
            OldRhsOld[i++] = ckt->CKTrhsOld[n->number];

        memcpy(OldCKTstate0, ckt->CKTstate0,
               (size_t) ckt->CKTnumStates * sizeof(double));

        if (ft_ngdebug)
            SPfrontEnd->IFerrorf (ERR_INFO, "One successful source step");
        ckt->CKTsrcFact = ConvFact + raise;

        do {
            if (ft_ngdebug)
                fprintf (stderr,
                         "Supplies reduced to %8.4f%% ", ckt->CKTsrcFact * 100);

#ifdef XSPICE
            /* gtri - wbk - add convergence problem reporting flags */
            ckt->enh->conv_debug.last_NIiter_call = MIF_TRUE;
#endif

            iters = ckt->CKTstat->STATnumIter;
            converged = NIiter (ckt, ckt->CKTdcTrcvMaxIter);
            iters = ckt->CKTstat->STATnumIter - iters;

            ckt->CKTmode = continuemode;

            if (converged == 0) {
                ConvFact = ckt->CKTsrcFact;

                for (i = 0, n = ckt->CKTnodes; n; n = n->next)
                    OldRhsOld[i++] = ckt->CKTrhsOld[n->number];

                memcpy(OldCKTstate0, ckt->CKTstate0,
                       (size_t) ckt->CKTnumStates * sizeof(double));
                if (ft_ngdebug)
                    SPfrontEnd->IFerrorf (ERR_INFO, "One successful source step");

                ckt->CKTsrcFact = ConvFact + raise;

                if (iters <= (ckt->CKTdcTrcvMaxIter / 4))
                    raise *= 1.5;

                if (iters > (3 * ckt->CKTdcTrcvMaxIter / 4))
                    raise *= 0.5;

                /* if (raise > 0.01)
                 *     raise = 0.01;
                 */

            }
/*            else if (ckt->CKTgmin < 1e-3){
                ckt->CKTdiagGmin = ckt->CKTgmin *= 10;
                fprintf(stderr,
                    "gmin raised to %8.4e\n", ckt->CKTgmin);
            }*/
            else {

                if (ckt->CKTsrcFact - ConvFact < 1e-8)
                    break;

                raise /= 10;

                if (raise > 0.01)
                    raise = 0.01;

                ckt->CKTsrcFact = ConvFact;

                for (i = 0, n = ckt->CKTnodes; n; n = n->next)
                    ckt->CKTrhsOld[n->number] = OldRhsOld[i++];

                memcpy(ckt->CKTstate0, OldCKTstate0,
                       (size_t) ckt->CKTnumStates * sizeof(double));
            }

            if (ckt->CKTsrcFact > 1)
                ckt->CKTsrcFact = 1;

        } while ((raise >= 1e-7) && (ConvFact < 1));

        /* Restore gmin, but leave no extra conductance on the diagonal:
         * like the other stepping routines, return with CKTdiagGmin equal
         * to CKTgshunt.  Leaving it at gmin put 1/gmin from every node to
         * ground into every later analysis of the same circuit (a later
         * op or tran that converged directly), silently changing its
         * solution until some other stepping routine reset it.
         */
        ckt->CKTgmin = gminstart;
        ckt->CKTdiagGmin = ckt->CKTgshunt;
        FREE (OldRhsOld);
        FREE (OldCKTstate0);
    }

    ckt->CKTsrcFact = 1;

    if (ConvFact != 1) {
        ckt->CKTcurrentAnalysis = DOING_TRAN;
        SPfrontEnd->IFerrorf (ERR_WARNING, "source stepping failed");
        return E_ITERLIM;
    } else {
        SPfrontEnd->IFerrorf (ERR_INFO, "Source stepping completed");
        return 0;
    }
}


/* Spice3 Source stepping
 * Modified 2005 - Paolo Nenzi (extracted from CKTop.c code)
 *
 * return value:
 * 0 -> method converged
 * 1 -> method failed
 *
 * Note that no path out of this code allows ckt->CKTsrcFact to be
 * anything but 1.00000.
 */

static int
spice3_src (CKTcircuit *ckt, long int firstmode,
            long int continuemode, int iterlim)
{
    int converged, i;

    NG_IGNORE(iterlim);

    ckt->CKTmode = firstmode;
    SPfrontEnd->IFerrorf (ERR_INFO, "Starting source stepping");

    for (i = 0; i <= ckt->CKTnumSrcSteps; i++) {
        ckt->CKTsrcFact = ((double) i) / ((double) ckt->CKTnumSrcSteps);
#ifdef XSPICE
        /* gtri - wbk - add convergence problem reporting flags */
        ckt->enh->conv_debug.last_NIiter_call = MIF_TRUE;
#endif
        converged = NIiter (ckt, ckt->CKTdcTrcvMaxIter);
        ckt->CKTmode = continuemode;
        if (converged != 0) {
            ckt->CKTsrcFact = 1;
            ckt->CKTcurrentAnalysis = DOING_TRAN;
            SPfrontEnd->IFerrorf (ERR_WARNING, "source stepping failed");
#ifdef XSPICE
            /* gtri - wbk - add convergence problem reporting flags */
            ckt->enh->conv_debug.last_NIiter_call = MIF_FALSE;
#endif
            return converged;
        }
        if (ft_ngdebug)
            SPfrontEnd->IFerrorf (ERR_INFO, "One successful source step");
    }

    if (ft_ngdebug)
        SPfrontEnd->IFerrorf (ERR_INFO, "Source stepping completed");
    ckt->CKTsrcFact = 1;

#ifdef XSPICE
    /* gtri - wbk - add convergence problem reporting flags */
    ckt->enh->conv_debug.last_NIiter_call = MIF_FALSE;
#endif

    return 0;
}



/* Largest level of an independent voltage source: the DC value and, for a
 * source with a waveform, the levels of the waveform (V1/V2 of PULSE and
 * EXP, offset + amplitude of SINE and SFFM, the values of PWL).  Times,
 * frequencies and the like are not levels and are skipped. */
static double
vsrc_level(VSRCinstance *here)
{
    double *c = here->VSRCcoeffs, v = fabs(here->VSRCdcValue);
    int n = here->VSRCfunctionOrder, i;

    if (!c || !here->VSRCcoeffsGiven)
        return v;
    switch (here->VSRCfunctionType) {
    case PULSE:
    case EXP:
        if (n > 0) v = MAX(v, fabs(c[0]));
        if (n > 1) v = MAX(v, fabs(c[1]));
        break;
    case SINE:
    case SFFM:
        if (n > 1) v = MAX(v, fabs(c[0]) + fabs(c[1]));
        else if (n > 0) v = MAX(v, fabs(c[0]));
        break;
    case PWL:
        for (i = 1; i < n; i += 2)
            v = MAX(v, fabs(c[i]));
        break;
    default:
        break;
    }
    return v;
}

static int
uf_find(int *parent, int x)
{
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

/* CKTsourceWindows --
 *
 * Damping window of every node for NIiter(): W = k * Vs + 1 V, where Vs is
 * the largest level of the independent voltage sources in the node's part
 * of the circuit.  The parts are the groups of nodes joined by the
 * terminals of any device instance; ground does not join them, so that,
 * e.g., a source that only carries a parameter to the .control section
 * does not widen the window of the rest of the circuit.  A part without a
 * voltage source of its own gets the window of the largest source of the
 * whole circuit.  The groups are computed once per circuit; the source
 * levels are read at every call (alter may change them).
 *
 * Returns an array indexed by equation number (valid up to CKTmaxEqNum),
 * owned by this routine.
 */
double *
CKTsourceWindows(CKTcircuit *ckt, double k)
{
    static CKTcircuit *c_ckt = NULL;
    static void *c_matrix = NULL;
    static int c_n = -1, *parent = NULL;
    static double *level = NULL, *win = NULL;
    SPICEdev **devs = devices();
    int n = ckt->CKTmaxEqNum, t, i, vtype = CKTtypelook("Vsource");
    double gmax = 0.0;
    GENmodel *m;
    GENinstance *inst;

    if (ckt != c_ckt || n != c_n || (void *) ckt->CKTmatrix != c_matrix) {
        parent = TREALLOC(int, parent, n + 1);
        level = TREALLOC(double, level, n + 1);
        win = TREALLOC(double, win, n + 1);
        for (i = 0; i <= n; i++)
            parent[i] = i;
        for (t = 0; t < DEVmaxnum; t++) {
            int nterm;
            if (!ckt->CKThead[t] || !devs[t] || !devs[t]->DEVpublic.terms)
                continue;
            nterm = *(devs[t]->DEVpublic.terms);
            for (m = ckt->CKThead[t]; m; m = m->GENnextModel)
                for (inst = m->GENinstances; inst; inst = inst->GENnextInstance) {
                    int *node = GENnode(inst), first = 0, j;
                    for (j = 0; j < nterm; j++) {
                        int a = node[j];
                        if (a <= 0 || a > n)
                            continue;
                        if (!first)
                            first = a;
                        else
                            parent[uf_find(parent, a)] = uf_find(parent, first);
                    }
                }
        }
        c_ckt = ckt;
        c_n = n;
        c_matrix = (void *) ckt->CKTmatrix;
    }

    for (i = 0; i <= n; i++)
        level[i] = -1.0;           /* no source in this part (yet) */
    if (vtype >= 0)
        for (m = ckt->CKThead[vtype]; m; m = m->GENnextModel) {
            VSRCinstance *here;
            for (here = VSRCinstances((VSRCmodel *) m); here;
                 here = VSRCnextInstance(here)) {
                double v = vsrc_level(here);
                int r;
                gmax = MAX(gmax, v);
                if (here->VSRCposNode > 0) {
                    r = uf_find(parent, here->VSRCposNode);
                    level[r] = MAX(level[r], v);
                }
                if (here->VSRCnegNode > 0) {
                    r = uf_find(parent, here->VSRCnegNode);
                    level[r] = MAX(level[r], v);
                }
            }
        }
    for (i = 0; i <= n; i++) {
        double l = level[uf_find(parent, i)];
        win[i] = k * (l < 0 ? gmax : l) + 1.0;
    }
    return win;
}
