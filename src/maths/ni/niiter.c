/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
Modified: 2001 AlansFixes
**********/

/*
 * NIiter(ckt,maxIter)
 *
 *  This subroutine performs the actual numerical iteration.
 *  It uses the sparse matrix stored in the circuit struct
 *  along with the matrix loading program, the load data, the
 *  convergence test function, and the convergence parameters
 */

#include "ngspice/ngspice.h"
#include "ngspice/trandefs.h"
#include "ngspice/cktdefs.h"
#include "ngspice/smpdefs.h"
#include "ngspice/sperror.h"
#include "ngspice/fteext.h"

extern double *CKTsourceWindows(CKTcircuit *ckt, double k);

/* Limit the number of 'singular matrix' warnings */
static int msgcount = 0;

/* NIiter() - return value is non-zero for convergence failure */

int
NIiter(CKTcircuit *ckt, int maxIter)
{
    double startTime, *OldCKTstate0 = NULL;
    int error, i, j;

    int iterno = 0;
    int ipass = 0;

    /* some convergence issues that get resolved by increasing max iter */
    if (maxIter < 100)
        maxIter = 100;

    if ((ckt->CKTmode & MODETRANOP) && (ckt->CKTmode & MODEUIC)) {
        SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);
        error = CKTload(ckt);
        if (error)
            return(error);
        return(OK);
    }

#ifdef WANT_SENSE2
    if (ckt->CKTsenInfo) {
        error = NIsenReinit(ckt);
        if (error)
            return(error);
    }
#endif

    if (ckt->CKTniState & NIUNINITIALIZED) {
        error = NIreinit(ckt); /* always returns 0 */
        if (error) {
#ifdef STEPDEBUG
            printf("re-init returned error \n");
#endif
            return(error);
        }
    }

    /* OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1); */

    /* fork-lite experiment: NGSPICE_DAMPWIN=k selects the window damping
     * (see the node damping block below), window = k * Vsrc_max + 1 V. */
    static double dampwin = -1;
    static int dampafter = -1;
    if (dampwin < 0) {
        const char *e = getenv("NGSPICE_DAMPWIN");
        dampwin = e ? atof(e) : 0.0;
    }
    /* fork-lite experiment: NGSPICE_DAMPAFTER=n adds the step damping
     * (absdv/reldv) once n iterations have passed without convergence,
     * e.g. when Newton is caught in a limit cycle. */
    if (dampafter < 0) {
        const char *e = getenv("NGSPICE_DAMPAFTER");
        dampafter = e ? atoi(e) : 0;
    }

    {
        static int ot = -1;
        if (ot < 0) ot = getenv("NGSPICE_OPTRACE") != NULL;
        if (ot) fprintf(stderr, "ST NIITER begin mode=0x%lx gmin=%.2e srcfact=%.3f time=%.4e maxiter=%d\n",
                        (unsigned long) ckt->CKTmode, ckt->CKTdiagGmin, ckt->CKTsrcFact, ckt->CKTtime, maxIter);
    }
    for (;;) {

        ckt->CKTnoncon = 0;

#ifdef NEWPRED
        if (!(ckt->CKTmode & MODEINITPRED))
#endif
        {

            error = CKTload(ckt);
            /* printf("loaded, noncon is %d\n", ckt->CKTnoncon); */
            /* fflush(stdout); */
            iterno++;
            if (error) {
                ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                printf("load returned error \n");
#endif
                FREE(OldCKTstate0);
                return (error);
            }

            /* printf("after loading, before solving\n"); */
            /* CKTdump(ckt); */

            if (!(ckt->CKTniState & NIDIDPREORDER)) {
                error = SMPpreOrder(ckt->CKTmatrix);
                if (error) {
                    ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                    printf("pre-order returned error \n");
#endif
                    FREE(OldCKTstate0);
                    return(error); /* badly formed matrix */
                }
                ckt->CKTniState |= NIDIDPREORDER;
            }

            if ((ckt->CKTmode & MODEINITJCT) ||
                ((ckt->CKTmode & MODEINITTRAN) && (iterno == 1)))
            {
                ckt->CKTniState |= NISHOULDREORDER;
            }

            if (ckt->CKTniState & NISHOULDREORDER) {
                startTime = SPfrontEnd->IFseconds();

#ifdef KLU
                if (ckt->CKTkluMODE) {
                    ckt->CKTmatrix->SMPkluMatrix->KLUloadDiagGmin = 1 ;
                }
#endif

                error = SMPreorder(ckt->CKTmatrix, ckt->CKTpivotAbsTol,
                                   ckt->CKTpivotRelTol, ckt->CKTdiagGmin);
                ckt->CKTstat->STATreorderTime +=
                    SPfrontEnd->IFseconds() - startTime;
                if (error) {
                    /* new feature - we can now find out something about what is
                     * wrong - so we ask for the troublesome entry
                     * Limit the number of messages to 6, if not 'set ngdebug'.
                     */
                    if (ft_ngdebug || msgcount < 6) {
                        SMPgetError(ckt->CKTmatrix, &i, &j);
                        if(eq(NODENAME(ckt, i), NODENAME(ckt, j)))
                            SPfrontEnd->IFerrorf(ERR_WARNING, "singular matrix:  check node %s\n", NODENAME(ckt, i));
                        else
                            SPfrontEnd->IFerrorf(ERR_WARNING, "singular matrix:  check nodes %s and %s\n", NODENAME(ckt, i), NODENAME(ckt, j));
                        msgcount += 1;
                    }
                    ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                    printf("reorder returned error \n");
#endif
                    FREE(OldCKTstate0);
                    return(error); /* can't handle these errors - pass up! */
                }
                ckt->CKTniState &= ~NISHOULDREORDER;
            } else {
                startTime = SPfrontEnd->IFseconds();

#ifdef KLU
                if (ckt->CKTkluMODE) {
                    ckt->CKTmatrix->SMPkluMatrix->KLUloadDiagGmin = 1 ;
                }
#endif

                error = SMPluFac(ckt->CKTmatrix, ckt->CKTpivotAbsTol,
                                 ckt->CKTdiagGmin);
                ckt->CKTstat->STATdecompTime +=
                    SPfrontEnd->IFseconds() - startTime;

#ifdef KLU
                if ((ckt->CKTkluMODE) && (error == E_SINGULAR)) {

                    /* Francesco Lannutti - 25 Aug 2020
                     * If the matrix is numerically singular during ReFactorization, take the same matrix and factor it from scratch in the same iteration.
                     * This is my mod with KLU. It saves run-time, but also the system at the next iteration may be different.
                     * How do we guarantee that the system is the same at the next iteration? So, the original SPARSE version below sounds like a bug.
                     */
                    if (ft_ngdebug)
                        fprintf (stderr, "Warning: KLU ReFactor failed. Factoring again...\n") ;
                    ckt->CKTniState |= NISHOULDREORDER;
                    ckt->CKTmatrix->SMPkluMatrix->KLUloadDiagGmin = 0 ;
                    error = SMPreorder(ckt->CKTmatrix, ckt->CKTpivotAbsTol, ckt->CKTpivotRelTol, ckt->CKTdiagGmin);
                    ckt->CKTstat->STATreorderTime += SPfrontEnd->IFseconds() - startTime;
                    if (error) {
                        SMPgetError(ckt->CKTmatrix, &i, &j);
                        if (ft_ngdebug || msgcount < 6) {
                            SMPgetError(ckt->CKTmatrix, &i, &j);
                            if (eq(NODENAME(ckt, i), NODENAME(ckt, j)))
                                SPfrontEnd->IFerrorf(ERR_WARNING, "singular matrix:  check node %s\n", NODENAME(ckt, i));
                            else
                                SPfrontEnd->IFerrorf(ERR_WARNING, "singular matrix:  check nodes %s and %s\n", NODENAME(ckt, i), NODENAME(ckt, j));
                            msgcount += 1;
                        }

                        /* CKTload(ckt); */
                        /* SMPprint(ckt->CKTmatrix, stdout); */
                        /* seems to be singular - pass the bad news up */
                        ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                        printf("lufac returned error \n");
#endif
                        FREE(OldCKTstate0);
                        return(error);
                    }
                } else if (error) {
                    if (!(ckt->CKTkluMODE) && (error == E_SINGULAR)) {

                        /* Francesco Lannutti - 25 Aug 2020
                         * If the matrix is numerically singular during ReFactorization, factor it from scratch at the next iteration.
                         * This is the original SPICE3F5 code and uses SPARSE.
                         */

                        ckt->CKTniState |= NISHOULDREORDER;
                        DEBUGMSG(" forced reordering....\n");
                        continue;
                    }
                    /* CKTload(ckt); */
                    /* SMPprint(ckt->CKTmatrix, stdout); */
                    /* seems to be singular - pass the bad news up */
                    ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                    printf("lufac returned error \n");
#endif
                    FREE(OldCKTstate0);
                    return(error);
                }
#else
                if (error) {
                    if (error == E_SINGULAR) {

                        /* Francesco Lannutti - 25 Aug 2020
                         * If the matrix is numerically singular during ReFactorization, factor it from scratch at the next iteration.
                         * This is the original SPICE3F5 code and uses SPARSE.
                         */

                        ckt->CKTniState |= NISHOULDREORDER;
                        DEBUGMSG(" forced reordering....\n");
                        continue;
                    }
                    /* CKTload(ckt); */
                    /* SMPprint(ckt->CKTmatrix, stdout); */
                    /* seems to be singular - pass the bad news up */
                    ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
                    printf("lufac returned error \n");
#endif
                    FREE(OldCKTstate0);
                    return(error);
                }
#endif

            }

            /* moved it to here as if xspice is included then CKTload changes
               CKTnumStates the first time it is run */
            /* fork-lite: OldCKTstate0 is only used by the node damping
               block below, which runs in the operating point only */
            if ((ckt->CKTnodeDamping || dampwin > 0 || dampafter > 0) &&
                (ckt->CKTmode & (MODETRANOP | MODEDCOP | MODEDCTRANCURVE))) {
                if (!OldCKTstate0)
                    OldCKTstate0 = TMALLOC(double, ckt->CKTnumStates + 1);
                if (ckt->CKTstate0)
                    memcpy(OldCKTstate0, ckt->CKTstate0,
                           (size_t) ckt->CKTnumStates * sizeof(double));
            }

            startTime = SPfrontEnd->IFseconds();
            SMPsolve(ckt->CKTmatrix, ckt->CKTrhs, ckt->CKTrhsSpare);
            ckt->CKTstat->STATsolveTime +=
                SPfrontEnd->IFseconds() - startTime;
#ifdef STEPDEBUG
            /*XXXX*/
            if (ckt->CKTrhs[0] != 0.0)
                printf("NIiter: CKTrhs[0] = %g\n", ckt->CKTrhs[0]);
            if (ckt->CKTrhsSpare[0] != 0.0)
                printf("NIiter: CKTrhsSpare[0] = %g\n", ckt->CKTrhsSpare[0]);
            if (ckt->CKTrhsOld[0] != 0.0)
                printf("NIiter: CKTrhsOld[0] = %g\n", ckt->CKTrhsOld[0]);
            /*XXXX*/
#endif
            ckt->CKTrhs[0] = 0;
            ckt->CKTrhsSpare[0] = 0;
            ckt->CKTrhsOld[0] = 0;

            if (iterno > maxIter) {
                ckt->CKTstat->STATnumIter += iterno;
                /* we don't use this info during transient analysis */
                if (ckt->CKTcurrentAnalysis != DOING_TRAN) {
                    FREE(errMsg);
                    errMsg = copy("Too many iterations without convergence");
#ifdef STEPDEBUG
                    fprintf(stderr, "too many iterations without convergence: %d iter's (max iter == %d)\n",
                    iterno, maxIter);
#endif
                    FREE(errMsg);
                }
                FREE(OldCKTstate0);
                return(E_ITERLIM);
            }

            /* fork-lite diagnostic: NGSPICE_NEWTONTRACE=n prints, for the first n
               Newton iterations, the size of the update over the whole solution:
               the largest |dx|/tol (and where), how many unknowns violate their
               tolerance, the largest voltage and current updates, and the init
               mode.  A ratio that drops ~10x per iteration instead of squaring
               means linear (not quadratic) convergence, i.e. an inexact Jacobian
               or a non-smooth model. */
            {
                static long ntmax = -1, ntshown = 0;
                if (ntmax < 0) { const char *e = getenv("NGSPICE_NEWTONTRACE"); ntmax = e ? atol(e) : 0; }
                if (ntshown < ntmax) {
                    CKTnode *nd = ckt->CKTnodes;
                    int k, nviol = 0, kr = 0, kv = 0, ki = 0;
                    int sz = SMPmatSize(ckt->CKTmatrix);
                    double rmax = 0, dvmax = 0, dimax = 0;
                    for (k = 1; k <= sz; k++) {
                        double nw, od, d, tl;
                        nd = nd->next;
                        nw = ckt->CKTrhs[k];
                        od = ckt->CKTrhsOld[k];
                        d = fabs(nw - od);
                        if (nd->type == SP_VOLTAGE) {
                            tl = ckt->CKTreltol * MAX(fabs(od), fabs(nw)) + ckt->CKTvoltTol;
                            if (d > dvmax) { dvmax = d; kv = k; }
                        } else {
                            tl = ckt->CKTreltol * MAX(fabs(od), fabs(nw)) + ckt->CKTabstol;
                            if (d > dimax) { dimax = d; ki = k; }
                        }
                        if (d > tl) nviol++;
                        if (d / tl > rmax) { rmax = d / tl; kr = k; }
                    }
                    ntshown++;
                    fprintf(stderr, "NT it=%d mode=%s loadnc=%d nviol=%d maxratio=%.3e at %s"
                            " | maxdv=%.3e at %s (%.4e -> %.4e) | maxdi=%.3e at %s\n", iterno,
                            (ckt->CKTmode & MODEINITJCT) ? "JCT" :
                            (ckt->CKTmode & MODEINITFIX) ? "FIX" :
                            (ckt->CKTmode & MODEINITFLOAT) ? "FLOAT" : "other",
                            ckt->CKTnoncon, nviol, rmax, kr ? CKTnodName(ckt, kr) : "-",
                            dvmax, kv ? CKTnodName(ckt, kv) : "-",
                            kv ? ckt->CKTrhsOld[kv] : 0.0, kv ? ckt->CKTrhs[kv] : 0.0,
                            dimax, ki ? CKTnodName(ckt, ki) : "-");
                }
            }

            /* fork-lite diagnostic: NGSPICE_VTRACE=t0 tabulates why each iteration
               is declared non-converged (load limiting vs voltage/current tests) */
            {
                static double wt0 = -1; static long wshown = 0;
                if (wt0 < 0) { const char *e = getenv("NGSPICE_VTRACE"); wt0 = e ? atof(e) : 1e300; }
                if (ckt->CKTtime >= wt0 && wshown < 600) {
                    int loadnc = ckt->CKTnoncon, conv = -1;
                    if (loadnc == 0 && iterno != 1)
                        conv = NIconvTest(ckt);
                    wshown++;
                    fprintf(stderr, "ST WHY t=%.6e it=%d loadnoncon=%d convtest=%d node=%s elt=%s\n", ckt->CKTtime, iterno,
                            loadnc, conv, ckt->CKTtroubleNode ? CKTnodName(ckt, ckt->CKTtroubleNode) : "-",
                            ckt->CKTtroubleElt ? ckt->CKTtroubleElt->GENname : "-");
                    if (loadnc == 0 && iterno != 1)
                        ckt->CKTnoncon = conv;
                    else
                        ckt->CKTnoncon = 1;
                    goto why_done;
                }
            }
            if ((ckt->CKTnoncon == 0) && (iterno != 1))
                ckt->CKTnoncon = NIconvTest(ckt);
            else
                ckt->CKTnoncon = 1;
            why_done:;

#ifdef STEPDEBUG
            printf("noncon is %d\n", ckt->CKTnoncon);
#endif
        }

        /* Node damping (.option nodedamping): in the DC operating point and
         * the DC sweep, limit the Newton update so that no node voltage moves
         * more than max(absdv, reldv * |v|) in one iteration.  The whole
         * update (node voltages, branch currents and device states) is scaled
         * by the same factor, so the Newton direction is kept; near the
         * solution the updates are far below the limit and the iteration is
         * unchanged.  Without it, Newton can take steps of many kV from a poor
         * starting point and settle on a non-physical root of a model outside
         * its range of validity.  absdv (0.5 V) and reldv (2.0) are the
         * defaults of the existing options.
         */
        if ((ckt->CKTnodeDamping != 0 || dampwin > 0 || dampafter > 0) &&
            (ckt->CKTnoncon != 0) &&
            (ckt->CKTmode & (MODETRANOP | MODEDCOP | MODEDCTRANCURVE)) &&
            (iterno > 1))
        {
            CKTnode *node;
            double diff, lim, damp_factor = 1.0;
            if (dampwin > 0) {
                /* Window damping: only the nodes whose new value would be
                 * outside |v| <= W (CKTsourceWindows(): per part of the circuit), limit the
                 * update.  Such a node may reach the window edge, and from
                 * outside it may move max(absdv, reldv * |v|), so that a
                 * genuine large voltage is still reached geometrically. */
                double *Wn = CKTsourceWindows(ckt, dampwin);
                for (node = ckt->CKTnodes->next; node; node = node->next)
                    if (node->type == SP_VOLTAGE) {
                        double W = Wn[node->number];
                        double vo = ckt->CKTrhsOld[node->number];
                        double vn = ckt->CKTrhs[node->number];
                        if (fabs(vn) <= W)
                            continue;
                        diff = fabs(vn - vo);
                        if (fabs(vo) < W)   /* inside: stop at the edge */
                            lim = MAX(ckt->CKTabsDv, fabs((vn > 0 ? W : -W) - vo));
                        else                /* outside: grow geometrically */
                            lim = MAX(ckt->CKTabsDv, ckt->CKTrelDv * fabs(vo));
                        if (diff > lim && lim / diff < damp_factor)
                            damp_factor = lim / diff;
                    }
            }
            if ((dampwin <= 0 && dampafter <= 0 && ckt->CKTnodeDamping) ||
                (dampafter > 0 && iterno > dampafter))
            for (node = ckt->CKTnodes->next; node; node = node->next)
                if (node->type == SP_VOLTAGE) {
                    diff = fabs(ckt->CKTrhs[node->number] - ckt->CKTrhsOld[node->number]);
                    lim = MAX(ckt->CKTabsDv,
                              ckt->CKTrelDv * fabs(ckt->CKTrhsOld[node->number]));
                    if (diff > lim && lim / diff < damp_factor)
                        damp_factor = lim / diff;
                }

            if (damp_factor < 1.0) {
                for (node = ckt->CKTnodes->next; node; node = node->next) {
                    diff = ckt->CKTrhs[node->number] - ckt->CKTrhsOld[node->number];
                    ckt->CKTrhs[node->number] =
                        ckt->CKTrhsOld[node->number] + (damp_factor * diff);
                }
                for (i = 0; i < ckt->CKTnumStates; i++) {
                    diff = ckt->CKTstate0[i] - OldCKTstate0[i];
                    ckt->CKTstate0[i] = OldCKTstate0[i] + (damp_factor * diff);
                }
            }
        }

        if (ckt->CKTmode & MODEINITFLOAT) {
            if ((ckt->CKTmode & MODEDC) && ckt->CKThadNodeset) {
                if (ipass)
                    ckt->CKTnoncon = ipass;
                ipass = 0;
            }
            if (ckt->CKTnoncon == 0) {
                ckt->CKTstat->STATnumIter += iterno;
                FREE(OldCKTstate0);
                return(OK);
            }
        } else if (ckt->CKTmode & MODEINITJCT) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFIX;
            ckt->CKTniState |= NISHOULDREORDER;
        } else if (ckt->CKTmode & MODEINITFIX) {
            if (ckt->CKTnoncon == 0)
                ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
            ipass = 1;
        } else if (ckt->CKTmode & MODEINITSMSIG) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else if (ckt->CKTmode & MODEINITTRAN) {
            if (iterno <= 1)
                ckt->CKTniState |= NISHOULDREORDER;
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else if (ckt->CKTmode & MODEINITPRED) {
            ckt->CKTmode = (ckt->CKTmode & ~INITF) | MODEINITFLOAT;
        } else {
            ckt->CKTstat->STATnumIter += iterno;
#ifdef STEPDEBUG
            printf("bad initf state \n");
#endif
            FREE(OldCKTstate0);
            return(E_INTERN);
            /* impossible - no such INITF flag! */
        }

        /* build up the lvnim1 array from the lvn array */
        SWAP(double *, ckt->CKTrhs, ckt->CKTrhsOld);
        /* printf("after loading, after solving\n"); */
        /* CKTdump(ckt); */
    }
    /*NOTREACHED*/
}

void NIresetwarnmsg(void) {
    msgcount = 0;
}
