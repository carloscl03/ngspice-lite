/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

/*
 * NIconvTest(ckt)
 *  perform the convergence test - returns 1 if any of the 
 *  values in the old and new arrays have changed by more 
 *  than absTol + relTol*(max(old,new)), otherwise returns 0
 */

#include "ngspice/ngspice.h"
#include "ngspice/cktdefs.h"
#include "ngspice/smpdefs.h"

extern bool ft_ngdebug;


static double ct_tot_g; static long ct_calls_g;
static void ct_report(void)
{
    if (ct_calls_g)
        fprintf(stderr, "LOADTRACE convtest: calls=%ld total=%.3fs (%.1fus/call)\n", ct_calls_g, ct_tot_g, 1e6 * ct_tot_g / ct_calls_g);
}

int
NIconvTest(CKTcircuit *ckt)
{
    int i; /* generic loop variable */
    int size;  /* size of the matrix */
    CKTnode *node; /* current matrix entry */
    double old;
    double new;
    double tol;
    static int nancount = 0;

    node = ckt->CKTnodes;
    size = SMPmatSize(ckt->CKTmatrix);
#ifdef STEPDEBUG
    for (i=1;i<=size;i++) {
        new =  ckt->CKTrhs [i] ;
        old =  ckt->CKTrhsOld [i] ;
        printf("chk for convergence:   %s    new: %g    old: %g\n",CKTnodName(ckt,i),new,old);
    }
#endif /* STEPDEBUG */
    for (i=1;i<=size;i++) {
        node = node->next;
        new =  ckt->CKTrhs [i] ;
        old =  ckt->CKTrhsOld [i] ;
        if (isnan(new)) {
            if (ft_ngdebug && nancount < 10) {
                fprintf(stderr, "Warning: non-convergence, node %s is nan\n", CKTnodName(ckt, i));
                nancount++;
            }
            else if (ft_ngdebug && nancount == 10) {
                fprintf(stderr, "    non-convergence warnings (nan) limited to 10 node %s\n", CKTnodName(ckt, i));
                nancount++;
            }
            return 1;
        }
        if(node->type == SP_VOLTAGE) {
            tol =  ckt->CKTreltol * (MAX(fabs(old),fabs(new))) +
                    ckt->CKTvoltTol;
            if (fabs(new-old) >tol ) {
                /* fork-lite diagnostic: NGSPICE_VTRACE=t0 tabulates, from time t0 on,
                   the first node failing the voltage test at each iteration */
                {
                    static double vt0 = -1; static long vshown = 0;
                    if (vt0 < 0) { const char *e = getenv("NGSPICE_VTRACE"); vt0 = e ? atof(e) : 1e300; }
                    if (ckt->CKTtime >= vt0 && vshown < 400) {
                        vshown++;
                        fprintf(stderr, "ST VNONCONV t=%.6e h=%.3e it=%d node=%s new=%.9e old=%.9e diff=%.3e tol=%.3e\n",
                                ckt->CKTtime, ckt->CKTdelta, ckt->CKTstat->STATnumIter, CKTnodName(ckt, i), new, old, new-old, tol);
                    }
                }
#ifdef STEPDEBUG
                printf(" non-convergence at node (type=3) %s (fabs(new-old)>tol --> fabs(%g-%g)>%g)\n",CKTnodName(ckt,i),new,old,tol);
                printf("    reltol: %g    voltTol: %g   (tol=reltol*(MAX(fabs(old),fabs(new))) + voltTol)\n",ckt->CKTreltol,ckt->CKTvoltTol);
#endif /* STEPDEBUG */
                ckt->CKTtroubleNode = i;
                ckt->CKTtroubleElt = NULL;
                return(1);
            }
        } else {
            tol =  ckt->CKTreltol * (MAX(fabs(old),fabs(new))) +
                    ckt->CKTabstol;
            if (fabs(new-old) >tol ) {
                {
                    static double it0 = -1; static long ishown = 0;
                    if (it0 < 0) { const char *e = getenv("NGSPICE_VTRACE"); it0 = e ? atof(e) : 1e300; }
                    if (ckt->CKTtime >= it0 && ishown < 300) {
                        ishown++;
                        fprintf(stderr, "ST INONCONV t=%.6e h=%.3e it=%d branch=%s new=%.9e old=%.9e diff=%.3e tol=%.3e\n",
                                ckt->CKTtime, ckt->CKTdelta, ckt->CKTstat->STATnumIter, CKTnodName(ckt, i), new, old, new-old, tol);
                    }
                }
#ifdef STEPDEBUG
                printf(" non-convergence at node (type=%d) %s (fabs(new-old)>tol --> fabs(%g-%g)>%g)\n",node->type,CKTnodName(ckt,i),new,old,tol);
                printf("    reltol: %g    abstol: %g   (tol=reltol*(MAX(fabs(old),fabs(new))) + abstol)\n",ckt->CKTreltol,ckt->CKTabstol);
#endif /* STEPDEBUG */
                ckt->CKTtroubleNode = i;
                ckt->CKTtroubleElt = NULL;
                return(1);
            }
        }
    }

#ifdef NEWCONV
    /* fork-lite diagnostic: time spent in the device convergence tests */
    {
        static int ct_on = -1; static double ct_total; static long ct_calls;
        double ct0 = 0;
        if (ct_on < 0) { ct_on = getenv("NGSPICE_LOADTRACE") != NULL; if (ct_on) atexit(ct_report); }
        if (ct_on) ct0 = SPfrontEnd->IFseconds();
        i = CKTconvTest(ckt);
        if (ct_on) { ct_total += SPfrontEnd->IFseconds() - ct0; ct_calls++; ct_tot_g = ct_total; ct_calls_g = ct_calls; }
    }
    /* The individual testers, called by CKTconvTest, set
     * ckt->CKTnoncon and ckt->CKTtroubleElt appropriately. */
    if (ckt->CKTnoncon != 0) {
        /* fork-lite: the node voltages converged, only device current tests
         * failed. curconv=0 ignores them (ngspice <= 46), 1 insists (47),
         * 2 insists until they have failed curconvmax consecutive times, then accepts until they pass again (see CKTsetup) */
        int give_up = ckt->CKTcurConv == 0 ||
            (ckt->CKTcurConv == 2 && ++ckt->CKTcurConvFails >= ckt->CKTcurConvMax);
        static int shown = 0;
        if (shown < 40 && getenv("NGSPICE_STEPTRACE") && !getenv("NGSPICE_VTRACE")) {
            shown++;
            fprintf(stderr, "ST CURNONCONV t=%.6e elt=%s fails=%d%s\n", ckt->CKTtime,
                    ckt->CKTtroubleElt ? ckt->CKTtroubleElt->GENname : "?",
                    ckt->CKTcurConvFails, give_up ? " (accepted)" : "");
        }
        if (!give_up) {
            ckt->CKTtroubleNode = 0;
            return(1);
        }
        ckt->CKTnoncon = 0;
    } else {
        /* currents converged: insist on them again from now on */
        ckt->CKTcurConvFails = 0;
    }
    /* CKTconvTest early-returns nonzero 'i' on the first error
     * in evaluating convergence (such as parameter out of range) so
     * there may be untested devices that have not yet converged */
    if (i)
	ckt->CKTtroubleNode = 0;
    return(i);
#else /* NEWCONV */
    return(0);
#endif /* NEWCONV */
}
