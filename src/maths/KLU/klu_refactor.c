/* ========================================================================== */
/* === KLU_refactor ========================================================= */
/* ========================================================================== */

/* Factor the matrix, after ordering and analyzing it with KLU_analyze, and
 * factoring it once with KLU_factor.  This routine cannot do any numerical
 * pivoting.  The pattern of the input matrix (Ap, Ai) must be identical to
 * the pattern given to KLU_factor.
 */

#include "klu_internal.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "ngspice/ngpool.h"

#ifndef COMPLEX
/* fork-lite: parallel refactorisation of the largest BTF block.

   klu_refactor on a circuit matrix is memory-latency bound (SPI with
   parasitics: 2.6 Mflop in 13 ms). Column k of the block depends only on the
   columns j with U(j,k) != 0, so the columns are grouped into levels of the
   elimination DAG and each level is processed on the thread pool, each thread
   with its own workspace X. The operations done for a column, and their order,
   are exactly those of the serial loop, so the factors are bit identical.
   Used when the block has at least NGSPICE_KLU_PARMIN columns (4000) and the
   critical-path bound of the DAG is below NGSPICE_KLU_PARCP percent (60);
   levels with fewer than NGSPICE_KLU_PARLEVEL columns (16) run on the master,
   and at most NGSPICE_KLU_PARTHREADS (4) pool threads take part, because the
   work is latency bound and SMT siblings slow it down.
   NGSPICE_KLU_NOPAR=1 disables it; NGSPICE_KLU_PARCHECK=1 runs the serial loop
   as well and compares the factors. */

extern int klu_rf_factor_gen;
static double klu_rf_now(void);

static struct {
    int valid, gen, T, enabled;
    void *numeric;
    Int block, k1, nk, nlev, offtotal, parmin;
    Int *order, *lstart, *offrel;
    Entry *X;                 /* T x nk */
    double cpbound;
    long par_levels, ser_levels, calls;
    double t_par, t_ser;
} ps;

struct rf_ctx {
    Int k1, nzoff, poff0;
    Int *Ap, *Ai, *Q, *Pinv, *Lip, *Llen, *Uip, *Ulen;
    Entry *Az, *Offx, *Udiag;
    double *Rs;
    Unit *LU;
};

static void ps_free(void)
{
    free(ps.order); free(ps.lstart); free(ps.offrel); free(ps.X);
    ps.order = NULL; ps.lstart = NULL; ps.offrel = NULL; ps.X = NULL;
    ps.valid = 0;
}

static int env_int_default(const char *name, int def)
{
    const char *v = getenv(name);
    return v ? atoi(v) : def;
}

static void ps_report(void)
{
    if (ps.calls)
        fprintf(stderr, "ST KLUPAR calls=%ld block=%ld nk=%ld levels=%ld cpbound=%.1f%% threads=%d parallel_levels=%.1f/call (%.2fms) serial_levels=%.1f/call (%.2fms)\n",
                ps.calls, (long) ps.block, (long) ps.nk, (long) ps.nlev, 100 * ps.cpbound, ps.T,
                (double) ps.par_levels / ps.calls, 1e3 * ps.t_par / ps.calls,
                (double) ps.ser_levels / ps.calls, 1e3 * ps.t_ser / ps.calls);
}

/* levels, column order, off-diagonal offsets and workspaces for one block */
static int
ps_build(KLU_numeric *Numeric, Int block, Int k1, Int nk, Int *Ap, Int *Ai, Int *Q, Int *Pinv,
         Int nzoff, Int poff0, int T, int trace)
{
    Unit *LU = (Unit *) Numeric->LUbx [block] ;
    Int *Uip = Numeric->Uip + k1, *Ulen = Numeric->Ulen + k1 ;
    Int *Lip = Numeric->Lip + k1, *Llen = Numeric->Llen + k1 ;
    Int *lev, *cnt, *Ui, *Li, ulen, llen, up, k, p, nlev = 0, cum, l ;
    Entry *Ux, *Lx ;
    double *wl, wtot = 0, cp = 0 ;
    static int reported = 0;

    ps_free();
    lev = (Int *) malloc (sizeof (Int) * (size_t) nk) ;
    for (k = 0 ; k < nk ; k++)
    {
        Int m = 0 ;
        GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, ulen) ;
        for (up = 0 ; up < ulen ; up++)
            if (lev [Ui [up]] > m) m = lev [Ui [up]] ;
        lev [k] = m + 1 ;
        if (lev [k] > nlev) nlev = lev [k] ;
    }
    cnt = (Int *) calloc ((size_t) nlev + 2, sizeof (Int)) ;
    wl = (double *) calloc ((size_t) nlev + 2, sizeof (double)) ;
    for (k = 0 ; k < nk ; k++)
    {
        double w ;
        GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, ulen) ;
        w = (double) ulen ;
        for (up = 0 ; up < ulen ; up++)
        {
            GET_POINTER (LU, Lip, Llen, Li, Lx, Ui [up], llen) ;
            w += (double) llen ;
        }
        GET_POINTER (LU, Lip, Llen, Li, Lx, k, llen) ;
        w += (double) llen + (double)(Ap [Q [k + k1] + 1] - Ap [Q [k + k1]]) ;
        cnt [lev [k]]++ ;
        wl [lev [k]] += w ;
        wtot += w ;
    }
    for (l = 1 ; l <= nlev ; l++)
        if (cnt [l]) cp += wl [l] / cnt [l] ;
    ps.cpbound = wtot > 0 ? cp / wtot : 1.0 ;

    /* counting sort by level, stable in k: deterministic order */
    ps.lstart = (Int *) malloc (sizeof (Int) * ((size_t) nlev + 2)) ;
    ps.lstart [1] = 0 ;
    for (l = 1 ; l <= nlev ; l++)
        ps.lstart [l + 1] = ps.lstart [l] + cnt [l] ;
    ps.order = (Int *) malloc (sizeof (Int) * (size_t) nk) ;
    for (l = 1 ; l <= nlev ; l++) cnt [l] = ps.lstart [l] ;
    for (k = 0 ; k < nk ; k++)
        ps.order [cnt [lev [k]]++] = k ;

    /* off-diagonal entries of each column, counted as the serial loop does */
    ps.offrel = (Int *) malloc (sizeof (Int) * ((size_t) nk + 1)) ;
    cum = 0 ;
    for (k = 0 ; k < nk ; k++)
    {
        Int oldcol = Q [k + k1], pend = Ap [oldcol + 1] ;
        ps.offrel [k] = cum ;
        for (p = Ap [oldcol] ; p < pend ; p++)
            if (Pinv [Ai [p]] - k1 < 0 && poff0 + cum < nzoff)
                cum++ ;
    }
    ps.offrel [nk] = cum ;
    ps.offtotal = cum ;

    ps.X = (Entry *) calloc ((size_t) T * (size_t) nk, sizeof (Entry)) ;
    ps.valid = 1 ; ps.gen = klu_rf_factor_gen ; ps.numeric = Numeric ; ps.T = T ;
    ps.block = block ; ps.k1 = k1 ; ps.nk = nk ; ps.nlev = nlev ;
    ps.parmin = env_int_default ("NGSPICE_KLU_PARLEVEL", 16) ;
    ps.enabled = 100.0 * ps.cpbound < env_int_default ("NGSPICE_KLU_PARCP", 60) ;
    if (trace && !reported) { reported = 1; atexit (ps_report) ; }
    if (trace)
        fprintf (stderr, "ST KLUPAR schedule: block=%ld nk=%ld levels=%ld cpbound=%.1f%% off=%ld threads=%d enabled=%d\n",
                 (long) block, (long) nk, (long) nlev, 100 * ps.cpbound, (long) cum, T, ps.enabled) ;
    free (lev) ; free (cnt) ; free (wl) ;
    return ps.enabled ;
}

/* one column of the block, exactly as the serial loop; returns k when the
   pivot is zero, -1 otherwise */
static inline Int
rf_col(Int k, const struct rf_ctx *c, Entry *X)
{
    Int k1 = c->k1, oldcol = c->Q [k + k1], pend = c->Ap [oldcol + 1] ;
    Int poff = c->poff0 + ps.offrel [k], p, newrow, up, ulen, llen, j, i, *Ui, *Li ;
    Entry *Ux, *Lx, ujk, ukk ;
    Int sing = -1 ;

    if (c->Rs == NULL)
    {
        for (p = c->Ap [oldcol] ; p < pend ; p++)
        {
            newrow = c->Pinv [c->Ai [p]] - k1 ;
            if (newrow < 0 && poff < c->nzoff)
            {
                c->Offx [poff] = c->Az [p] ;
                poff++ ;
            }
            else
            {
                X [newrow] = c->Az [p] ;
            }
        }
    }
    else
    {
        for (p = c->Ap [oldcol] ; p < pend ; p++)
        {
            Int oldrow = c->Ai [p] ;
            newrow = c->Pinv [oldrow] - k1 ;
            if (newrow < 0 && poff < c->nzoff)
            {
                SCALE_DIV_ASSIGN (c->Offx [poff], c->Az [p], c->Rs [oldrow]) ;
                poff++ ;
            }
            else
            {
                SCALE_DIV_ASSIGN (X [newrow], c->Az [p], c->Rs [oldrow]) ;
            }
        }
    }
    GET_POINTER (c->LU, c->Uip, c->Ulen, Ui, Ux, k, ulen) ;
    for (up = 0 ; up < ulen ; up++)
    {
        j = Ui [up] ;
        ujk = X [j] ;
        CLEAR (X [j]) ;
        Ux [up] = ujk ;
        GET_POINTER (c->LU, c->Lip, c->Llen, Li, Lx, j, llen) ;
        for (p = 0 ; p < llen ; p++)
        {
            MULT_SUB (X [Li [p]], Lx [p], ujk) ;
        }
    }
    ukk = X [k] ;
    CLEAR (X [k]) ;
    if (IS_ZERO (ukk))
        sing = k ;
    c->Udiag [k + k1] = ukk ;
    GET_POINTER (c->LU, c->Lip, c->Llen, Li, Lx, k, llen) ;
    for (p = 0 ; p < llen ; p++)
    {
        i = Li [p] ;
        DIV (Lx [p], X [i], ukk) ;
        CLEAR (X [i]) ;
    }
    return sing ;
}

struct rf_job { const struct rf_ctx *c ; Int base ; Int sing [64] ; } ;

static int
rf_item(int i, void *arg)
{
    struct rf_job *j = (struct rf_job *) arg ;
    int t = ngpool_tid () ;
    Int r = rf_col (ps.order [j->base + i], j->c, ps.X + (size_t) t * (size_t) ps.nk) ;
    if (r >= 0 && (j->sing [t] < 0 || r < j->sing [t])) j->sing [t] = r ;
    return 0 ;
}

/* returns 1 when done, 0 when the caller must fall back to the serial loop,
   -1 when singular and halt_if_singular */
static int
ps_run(const struct rf_ctx *c, KLU_common *Common, int T, int trace)
{
    Int l, first_sing = -1 ;
    struct rf_job job ;
    double t0 = 0 ;
    job.c = c ;
    for (l = 1 ; l <= ps.nlev ; l++)
    {
        Int n_l = ps.lstart [l + 1] - ps.lstart [l], i, t ;
        if (trace) t0 = klu_rf_now () ;
        if (n_l < ps.parmin || T < 2)
        {
            for (i = 0 ; i < n_l ; i++)
            {
                Int r = rf_col (ps.order [ps.lstart [l] + i], c, ps.X) ;
                if (r >= 0 && (first_sing < 0 || r < first_sing)) first_sing = r ;
            }
            if (trace) { ps.ser_levels++ ; ps.t_ser += klu_rf_now () - t0 ; }
        }
        else
        {
            int chunk = (int) (n_l / (8 * T)) ;
            if (chunk < 1) chunk = 1 ;
            job.base = ps.lstart [l] ;
            for (t = 0 ; t < T ; t++) job.sing [t] = -1 ;
            ngpool_run_n ((int) n_l, rf_item, &job, chunk, T) ;
            for (t = 0 ; t < T ; t++)
                if (job.sing [t] >= 0 && (first_sing < 0 || job.sing [t] < first_sing)) first_sing = job.sing [t] ;
            if (trace) { ps.par_levels++ ; ps.t_par += klu_rf_now () - t0 ; }
        }
        if (first_sing >= 0 && Common->halt_if_singular)
            break ;
    }
    if (first_sing >= 0)
    {
        Common->status = KLU_SINGULAR ;
        if (Common->numerical_rank == EMPTY)
        {
            Common->numerical_rank = first_sing + c->k1 ;
            Common->singular_col = c->Q [first_sing + c->k1] ;
        }
        if (Common->halt_if_singular)
            return -1 ;
    }
    if (trace) ps.calls++ ;
    return 1 ;
}
#endif /* !COMPLEX */

/* fork-lite diagnostic (NGSPICE_OPTRACE): where the time of a refactorisation
   goes: row scaling, the column loop (scatter + L/U update) and the tail */
static int klu_rf_trace = -1;
static long klu_rf_calls = 0;
static double klu_rf_t_scale = 0, klu_rf_t_loop = 0, klu_rf_t_tail = 0;
static double klu_rf_now(void) { struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts); return ts.tv_sec + 1e-9 * ts.tv_nsec; }
static void klu_rf_report(void)
{
    if (klu_rf_calls)
        fprintf(stderr, "ST KLUREFACTOR calls=%ld scale=%.3fs (%.2fms) loop=%.3fs (%.2fms) tail=%.3fs (%.2fms)\n",
                klu_rf_calls, klu_rf_t_scale, 1e3 * klu_rf_t_scale / klu_rf_calls,
                klu_rf_t_loop, 1e3 * klu_rf_t_loop / klu_rf_calls, klu_rf_t_tail, 1e3 * klu_rf_t_tail / klu_rf_calls);
}


/* ========================================================================== */
/* === KLU_refactor ========================================================= */
/* ========================================================================== */

#ifndef COMPLEX
/* dispatcher: 1 = block done in parallel (poff advanced), 0 = use the serial loop */
static int
klu_rf_par(KLU_numeric *Numeric, KLU_common *Common, Int block, Int k1, Int nk,
           Int *Ap, Int *Ai, Entry *Az, Int *Q, Int *Pinv, Int nzoff, Int *poff, Entry *Offx,
           double *Rs, Entry *Udiag, Int *Lip, Int *Llen, Int *Uip, Int *Ulen, Unit *LU)
{
    static int nopar = -1, parmin = 0, check = 0, trace = 0 ;
    int T = ngpool_threads () ;
    struct rf_ctx c ;
    int r ;

    if (nopar < 0)
    {
        nopar = env_int_default ("NGSPICE_KLU_NOPAR", 0) ;
        parmin = env_int_default ("NGSPICE_KLU_PARMIN", 4000) ;
        check = env_int_default ("NGSPICE_KLU_PARCHECK", 0) ;
        trace = getenv ("NGSPICE_OPTRACE") != NULL ;
    }
    if (nopar || T < 2 || T > 64 || nk < parmin)
        return 0 ;
    {
        static int cap = -1 ;
        if (cap < 0) cap = env_int_default ("NGSPICE_KLU_PARTHREADS", 4) ;
        if (cap > 0 && T > cap) T = cap ;
    }
    if (!ps.valid || ps.gen != klu_rf_factor_gen || ps.numeric != Numeric || ps.block != block ||
        ps.nk != nk || ps.T != T)
    {
        if (!ps_build (Numeric, block, k1, nk, Ap, Ai, Q, Pinv, nzoff, *poff, T, trace))
            return 0 ;
    }
    if (!ps.enabled)
        return 0 ;
    c.k1 = k1 ; c.nzoff = nzoff ; c.poff0 = *poff ;
    c.Ap = Ap ; c.Ai = Ai ; c.Q = Q ; c.Pinv = Pinv ;
    c.Lip = Lip ; c.Llen = Llen ; c.Uip = Uip ; c.Ulen = Ulen ;
    c.Az = Az ; c.Offx = Offx ; c.Udiag = Udiag ; c.Rs = Rs ; c.LU = LU ;
    if (check)
    {
        /* validation: the parallel factors of the previous call were stashed
           and the original serial loop then redid the block; LU/Udiag/Offx
           hold its result now, so compare before doing anything else. Then
           run the parallel path again, stash it, and return 0 so the serial
           loop runs as usual (the simulation uses the serial factors). */
        static Unit *lu_par = NULL ; static Entry *ud_par = NULL, *off_par = NULL ;
        static size_t lus_s = 0 ; static Int nk_s = 0, off_s = 0, poff_s = 0 ; static Entry *udiag_s = NULL, *offx_s = NULL ; static Unit *lu_s = NULL ;
        size_t lus = Numeric->LUsize [block] ;
        if (lu_par)
        {
            Int nbad = 0 ;
            if (memcmp (lu_par, lu_s, lus_s * sizeof (Unit))) nbad++ ;
            if (memcmp (ud_par, udiag_s, sizeof (Entry) * (size_t) nk_s)) nbad += 2 ;
            if (memcmp (off_par, offx_s + poff_s, sizeof (Entry) * (size_t) off_s)) nbad += 4 ;
            fprintf (stderr, "ST KLUPAR check: %s (code %ld)\n", nbad ? "MISMATCH" : "identical", (long) nbad) ;
            free (lu_par) ; free (ud_par) ; free (off_par) ; lu_par = NULL ;
        }
        r = ps_run (&c, Common, T, trace) ;
        if (r <= 0) return 0 ;
        lu_par = (Unit *) malloc (lus * sizeof (Unit)) ; memcpy (lu_par, LU, lus * sizeof (Unit)) ;
        ud_par = (Entry *) malloc (sizeof (Entry) * (size_t) nk) ; memcpy (ud_par, Udiag + k1, sizeof (Entry) * (size_t) nk) ;
        off_par = (Entry *) malloc (sizeof (Entry) * ((size_t) ps.offtotal + 1)) ; memcpy (off_par, Offx + *poff, sizeof (Entry) * (size_t) ps.offtotal) ;
        lus_s = lus ; nk_s = nk ; off_s = ps.offtotal ; poff_s = *poff ; udiag_s = Udiag + k1 ; offx_s = Offx ; lu_s = LU ;
        return 0 ;
    }
    r = ps_run (&c, Common, T, trace) ;
    if (r == 0) return 0 ;
    *poff += ps.offtotal ;
    return 1 ;
}
#endif

Int KLU_refactor        /* returns TRUE if successful, FALSE otherwise */
(
    /* inputs, not modified */
    Int Ap [ ],         /* size n+1, column pointers */
    Int Ai [ ],         /* size nz, row indices */
    double Ax [ ],
    KLU_symbolic *Symbolic,

    /* input/output */
    KLU_numeric *Numeric,
    KLU_common  *Common
)
{
    Entry ukk, ujk, s ;
    Entry *Offx, *Lx, *Ux, *X, *Az, *Udiag ;
    double *Rs ;
    Int *P, *Q, *R, *Pnum, *Offp, *Offi, *Ui, *Li, *Pinv, *Lip, *Uip, *Llen,
        *Ulen ;
    Unit **LUbx ;
    Unit *LU ;
    Int k1, k2, nk, k, block, oldcol, pend, oldrow, n, p, newrow, scale,
        nblocks, poff, i, j, up, ulen, llen, maxblock, nzoff ;

    /* ---------------------------------------------------------------------- */
    /* check inputs */
    /* ---------------------------------------------------------------------- */

    if (Common == NULL)
    {
        return (FALSE) ;
    }
    if (Common->status == KLU_EMPTY_MATRIX)
    {
        return (FALSE) ;
    }
    Common->status = KLU_OK ;
    double klu_rf_t0 = 0, klu_rf_t1 = 0;
    if (klu_rf_trace < 0) { klu_rf_trace = getenv("NGSPICE_OPTRACE") != NULL; if (klu_rf_trace) atexit(klu_rf_report); }

    if (Numeric == NULL)
    {
        /* invalid Numeric object */
        Common->status = KLU_INVALID ;
        return (FALSE) ;
    }

    Common->numerical_rank = EMPTY ;
    Common->singular_col = EMPTY ;

    Az = (Entry *) Ax ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Symbolic object */
    /* ---------------------------------------------------------------------- */

    n = Symbolic->n ;
    P = Symbolic->P ;
    Q = Symbolic->Q ;
    R = Symbolic->R ;
    nblocks = Symbolic->nblocks ;
    maxblock = Symbolic->maxblock ;

    /* ---------------------------------------------------------------------- */
    /* get the contents of the Numeric object */
    /* ---------------------------------------------------------------------- */

    Pnum = Numeric->Pnum ;
    Offp = Numeric->Offp ;
    Offi = Numeric->Offi ;
    Offx = (Entry *) Numeric->Offx ;

    LUbx = (Unit **) Numeric->LUbx ;

    scale = Common->scale ;
    if (scale > 0)
    {
        /* factorization was not scaled, but refactorization is scaled */
        if (Numeric->Rs == NULL)
        {
            Numeric->Rs = KLU_malloc (n, sizeof (double), Common) ;
            if (Common->status < KLU_OK)
            {
                Common->status = KLU_OUT_OF_MEMORY ;
                return (FALSE) ;
            }
        }
    }
    else
    {
        /* no scaling for refactorization; ensure Numeric->Rs is freed.  This
         * does nothing if Numeric->Rs is already NULL. */
        Numeric->Rs = KLU_free (Numeric->Rs, n, sizeof (double), Common) ;
    }
    Rs = Numeric->Rs ;

    Pinv = Numeric->Pinv ;
    X = (Entry *) Numeric->Xwork ;
    Common->nrealloc = 0 ;
    Udiag = Numeric->Udiag ;
    nzoff = Symbolic->nzoff ;

    /* ---------------------------------------------------------------------- */
    /* check the input matrix compute the row scale factors, Rs */
    /* ---------------------------------------------------------------------- */

    /* do no scale, or check the input matrix, if scale < 0 */
    if (scale >= 0)
    {
        /* check for out-of-range indices, but do not check for duplicates */
        if (klu_rf_trace) klu_rf_t0 = klu_rf_now();
        if (!KLU_scale (scale, n, Ap, Ai, Ax, Rs, NULL, Common))
        {
            return (FALSE) ;
        }
    }

    /* ---------------------------------------------------------------------- */
    /* clear workspace X */
    /* ---------------------------------------------------------------------- */

    for (k = 0 ; k < maxblock ; k++)
    {
        /* X [k] = 0 */
        CLEAR (X [k]) ;
    }

    if (klu_rf_trace) { klu_rf_t1 = klu_rf_now(); if (scale > 0) klu_rf_t_scale += klu_rf_t1 - klu_rf_t0; klu_rf_t0 = klu_rf_t1; }
    poff = 0 ;

    /* ---------------------------------------------------------------------- */
    /* factor each block */
    /* ---------------------------------------------------------------------- */

    if (scale <= 0)
    {

        /* ------------------------------------------------------------------ */
        /* no scaling */
        /* ------------------------------------------------------------------ */

        for (block = 0 ; block < nblocks ; block++)
        {

            /* -------------------------------------------------------------- */
            /* the block is from rows/columns k1 to k2-1 */
            /* -------------------------------------------------------------- */

            k1 = R [block] ;
            k2 = R [block+1] ;
            nk = k2 - k1 ;

            if (nk == 1)
            {

                /* ---------------------------------------------------------- */
                /* singleton case */
                /* ---------------------------------------------------------- */

                oldcol = Q [k1] ;
                pend = Ap [oldcol+1] ;
                CLEAR (s) ;
                for (p = Ap [oldcol] ; p < pend ; p++)
                {
                    newrow = Pinv [Ai [p]] - k1 ;
                    if (newrow < 0 && poff < nzoff)
                    {
                        /* entry in off-diagonal block */
                        Offx [poff] = Az [p] ;
                        poff++ ;
                    }
                    else
                    {
                        /* singleton */
                        s = Az [p] ;
                    }
                }
                Udiag [k1] = s ;

            }
            else
            {

                /* ---------------------------------------------------------- */
                /* construct and factor the kth block */
                /* ---------------------------------------------------------- */

                Lip  = Numeric->Lip  + k1 ;
                Llen = Numeric->Llen + k1 ;
                Uip  = Numeric->Uip  + k1 ;
                Ulen = Numeric->Ulen + k1 ;
                LU = LUbx [block] ;
#ifndef COMPLEX
                if (klu_rf_par (Numeric, Common, block, k1, nk, Ap, Ai, Az, Q, Pinv, nzoff, &poff, Offx,
                                scale > 0 ? Rs : NULL, Udiag, Lip, Llen, Uip, Ulen, LU))
                {
                    if (Common->status == KLU_SINGULAR && Common->halt_if_singular)
                        return (FALSE) ;
                    continue ;
                }
#endif

                for (k = 0 ; k < nk ; k++)
                {

                    /* ------------------------------------------------------ */
                    /* scatter kth column of the block into workspace X */
                    /* ------------------------------------------------------ */

                    oldcol = Q [k+k1] ;
                    pend = Ap [oldcol+1] ;
                    for (p = Ap [oldcol] ; p < pend ; p++)
                    {
                        newrow = Pinv [Ai [p]] - k1 ;
                        if (newrow < 0 && poff < nzoff)
                        {
                            /* entry in off-diagonal block */
                            Offx [poff] = Az [p] ;
                            poff++ ;
                        }
                        else
                        {
                            /* (newrow,k) is an entry in the block */
                            X [newrow] = Az [p] ;
                        }
                    }

                    /* ------------------------------------------------------ */
                    /* compute kth column of U, and update kth column of A */
                    /* ------------------------------------------------------ */

                    GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, ulen) ;
                    for (up = 0 ; up < ulen ; up++)
                    {
                        j = Ui [up] ;
                        ujk = X [j] ;
                        /* X [j] = 0 */
                        CLEAR (X [j]) ;
                        Ux [up] = ujk ;
                        GET_POINTER (LU, Lip, Llen, Li, Lx, j, llen) ;
                        for (p = 0 ; p < llen ; p++)
                        {
                            /* X [Li [p]] -= Lx [p] * ujk */
                            MULT_SUB (X [Li [p]], Lx [p], ujk) ;
                        }
                    }
                    /* get the diagonal entry of U */
                    ukk = X [k] ;
                    /* X [k] = 0 */
                    CLEAR (X [k]) ;
                    /* singular case */
                    if (IS_ZERO (ukk))
                    {
                        /* matrix is numerically singular */
                        Common->status = KLU_SINGULAR ;
                        if (Common->numerical_rank == EMPTY)
                        {
                            Common->numerical_rank = k+k1 ;
                            Common->singular_col = Q [k+k1] ;
                        }
                        if (Common->halt_if_singular)
                        {
                            /* do not continue the factorization */
                            return (FALSE) ;
                        }
                    }
                    Udiag [k+k1] = ukk ;
                    /* gather and divide by pivot to get kth column of L */
                    GET_POINTER (LU, Lip, Llen, Li, Lx, k, llen) ;
                    for (p = 0 ; p < llen ; p++)
                    {
                        i = Li [p] ;
                        DIV (Lx [p], X [i], ukk) ;
                        CLEAR (X [i]) ;
                    }

                }
            }
        }

    }
    else
    {

        /* ------------------------------------------------------------------ */
        /* scaling */
        /* ------------------------------------------------------------------ */

        for (block = 0 ; block < nblocks ; block++)
        {

            /* -------------------------------------------------------------- */
            /* the block is from rows/columns k1 to k2-1 */
            /* -------------------------------------------------------------- */

            k1 = R [block] ;
            k2 = R [block+1] ;
            nk = k2 - k1 ;

            if (nk == 1)
            {

                /* ---------------------------------------------------------- */
                /* singleton case */
                /* ---------------------------------------------------------- */

                oldcol = Q [k1] ;
                pend = Ap [oldcol+1] ;
                CLEAR (s) ;
                for (p = Ap [oldcol] ; p < pend ; p++)
                {
                    oldrow = Ai [p] ;
                    newrow = Pinv [oldrow] - k1 ;
                    if (newrow < 0 && poff < nzoff)
                    {
                        /* entry in off-diagonal block */
                        /* Offx [poff] = Az [p] / Rs [oldrow] */
                        SCALE_DIV_ASSIGN (Offx [poff], Az [p], Rs [oldrow]) ;
                        poff++ ;
                    }
                    else
                    {
                        /* singleton */
                        /* s = Az [p] / Rs [oldrow] */
                        SCALE_DIV_ASSIGN (s, Az [p], Rs [oldrow]) ;
                    }
                }
                Udiag [k1] = s ;

            }
            else
            {

                /* ---------------------------------------------------------- */
                /* construct and factor the kth block */
                /* ---------------------------------------------------------- */

                Lip  = Numeric->Lip  + k1 ;
                Llen = Numeric->Llen + k1 ;
                Uip  = Numeric->Uip  + k1 ;
                Ulen = Numeric->Ulen + k1 ;
                LU = LUbx [block] ;
#ifndef COMPLEX
                if (klu_rf_par (Numeric, Common, block, k1, nk, Ap, Ai, Az, Q, Pinv, nzoff, &poff, Offx,
                                scale > 0 ? Rs : NULL, Udiag, Lip, Llen, Uip, Ulen, LU))
                {
                    if (Common->status == KLU_SINGULAR && Common->halt_if_singular)
                        return (FALSE) ;
                    continue ;
                }
#endif

                for (k = 0 ; k < nk ; k++)
                {

                    /* ------------------------------------------------------ */
                    /* scatter kth column of the block into workspace X */
                    /* ------------------------------------------------------ */

                    oldcol = Q [k+k1] ;
                    pend = Ap [oldcol+1] ;
                    for (p = Ap [oldcol] ; p < pend ; p++)
                    {
                        oldrow = Ai [p] ;
                        newrow = Pinv [oldrow] - k1 ;
                        if (newrow < 0 && poff < nzoff)
                        {
                            /* entry in off-diagonal part */
                            /* Offx [poff] = Az [p] / Rs [oldrow] */
                            SCALE_DIV_ASSIGN (Offx [poff], Az [p], Rs [oldrow]);
                            poff++ ;
                        }
                        else
                        {
                            /* (newrow,k) is an entry in the block */
                            /* X [newrow] = Az [p] / Rs [oldrow] */
                            SCALE_DIV_ASSIGN (X [newrow], Az [p], Rs [oldrow]) ;
                        }
                    }

                    /* ------------------------------------------------------ */
                    /* compute kth column of U, and update kth column of A */
                    /* ------------------------------------------------------ */

                    GET_POINTER (LU, Uip, Ulen, Ui, Ux, k, ulen) ;
                    for (up = 0 ; up < ulen ; up++)
                    {
                        j = Ui [up] ;
                        ujk = X [j] ;
                        /* X [j] = 0 */
                        CLEAR (X [j]) ;
                        Ux [up] = ujk ;
                        GET_POINTER (LU, Lip, Llen, Li, Lx, j, llen) ;
                        for (p = 0 ; p < llen ; p++)
                        {
                            /* X [Li [p]] -= Lx [p] * ujk */
                            MULT_SUB (X [Li [p]], Lx [p], ujk) ;
                        }
                    }
                    /* get the diagonal entry of U */
                    ukk = X [k] ;
                    /* X [k] = 0 */
                    CLEAR (X [k]) ;
                    /* singular case */
                    if (IS_ZERO (ukk))
                    {
                        /* matrix is numerically singular */
                        Common->status = KLU_SINGULAR ;
                        if (Common->numerical_rank == EMPTY)
                        {
                            Common->numerical_rank = k+k1 ;
                            Common->singular_col = Q [k+k1] ;
                        }
                        if (Common->halt_if_singular)
                        {
                            /* do not continue the factorization */
                            return (FALSE) ;
                        }
                    }
                    Udiag [k+k1] = ukk ;
                    /* gather and divide by pivot to get kth column of L */
                    GET_POINTER (LU, Lip, Llen, Li, Lx, k, llen) ;
                    for (p = 0 ; p < llen ; p++)
                    {
                        i = Li [p] ;
                        DIV (Lx [p], X [i], ukk) ;
                        CLEAR (X [i]) ;
                    }
                }
            }
        }
    }

    /* ---------------------------------------------------------------------- */
    if (klu_rf_trace) { klu_rf_t1 = klu_rf_now(); klu_rf_t_loop += klu_rf_t1 - klu_rf_t0; klu_rf_t0 = klu_rf_t1; }
    /* fork-lite diagnostic, once: depth and width of the column dependency
       DAG of the largest block (column k depends on every j with U(j,k) != 0).
       Tells whether a level-scheduled parallel refactorisation has enough
       independent columns per level to be worth building. */
    if (klu_rf_trace && klu_rf_calls == 0)
    {
        for (block = 0 ; block < nblocks ; block++)
        {
            Int k1b = R [block], k2b = R [block+1], nkb = k2b - k1b, kk, nlev = 0, lmax = 0, l ;
            Int *lev, *cnt ;
            Unit *LUb ;
            Int *Uipb, *Ulenb, *Uib ;
            Int *LipO = Numeric->Lip + k1b, *LlenO = Numeric->Llen + k1b ;
            Entry *Uxb ;
            double *wl, wmax = 0, wtot = 0 ;
            if (nkb < 2) continue ;
            lev = (Int *) malloc (sizeof (Int) * nkb) ;
            Uipb = Numeric->Uip + k1b ; Ulenb = Numeric->Ulen + k1b ;
            LUb = (Unit *) Numeric->LUbx [block] ;
            {
                Int *Lipb = Numeric->Lip + k1b, *Llenb = Numeric->Llen + k1b, *Lib ;
                Entry *Lxb ;
                Int llenb ;
                for (kk = 0 ; kk < nkb ; kk++)
                {
                    Int m = 0 ;
                    GET_POINTER (LUb, Uipb, Ulenb, Uib, Uxb, kk, ulen) ;
                    for (up = 0 ; up < ulen ; up++)
                        if (lev [Uib [up]] > m) m = lev [Uib [up]] ;
                    lev [kk] = m + 1 ;
                    if (lev [kk] > nlev) nlev = lev [kk] ;
                    (void) Lipb ; (void) Llenb ; (void) Lib ; (void) Lxb ; (void) llenb ;
                }
            }
            cnt = (Int *) calloc ((size_t) nlev + 1, sizeof (Int)) ;
            wl = (double *) calloc ((size_t) nlev + 1, sizeof (double)) ;
            for (kk = 0 ; kk < nkb ; kk++)
            {
                /* work of a column ~ scattered entries + L updates it performs */
                double w = 0 ;
                Int *Lib2 ; Entry *Lxb2 ; Int llen2 ;
                GET_POINTER (LUb, Uipb, Ulenb, Uib, Uxb, kk, ulen) ;
                w += (double) ulen ;
                for (up = 0 ; up < ulen ; up++)
                {
                    GET_POINTER (LUb, LipO, LlenO, Lib2, Lxb2, Uib [up], llen2) ;
                    w += (double) llen2 ;
                }
                cnt [lev [kk]]++ ;
                wl [lev [kk]] += w ;
                wtot += w ;
            }
            for (l = 1 ; l <= nlev ; l++) { if (cnt [l] > lmax) lmax = cnt [l] ; if (wl [l] > wmax) wmax = wl [l] ; }
            fprintf (stderr, "ST KLULEVELS block=%ld n=%ld levels=%ld widest=%ld level1=%ld levels_with_lt8_cols=%ld work_total=%.3g work_widest_level=%.3g (critical path bound %.1f%% of serial)\n",
                     (long) block, (long) nkb, (long) nlev, (long) lmax, (long) cnt [1],
                     (long) ({ Int c = 0 ; for (l = 1 ; l <= nlev ; l++) if (cnt [l] < 8) c++ ; c ; }),
                     wtot, wmax, 100.0 * ({ double cp = 0 ; for (l = 1 ; l <= nlev ; l++) cp += (cnt [l] ? wl [l] / cnt [l] : 0) ; cp ; }) / wtot) ;
            {
                /* level histogram, 12 bins */
                Int b, nb = 12 ;
                fprintf (stderr, "ST KLULEVELS cols per level (level:cols:work):") ;
                for (b = 0 ; b < nb ; b++)
                {
                    Int l0 = 1 + (Int)((long) nlev * b / nb), l1 = 1 + (Int)((long) nlev * (b + 1) / nb), c = 0 ; double w = 0 ;
                    for (l = l0 ; l < l1 ; l++) { c += cnt [l] ; w += wl [l] ; }
                    fprintf (stderr, " %ld-%ld:%ld:%.2g", (long) l0, (long) l1 - 1, (long) c, w) ;
                }
                fprintf (stderr, "\n") ;
            }
            free (lev) ; free (cnt) ; free (wl) ;
        }
    }
    /* permute scale factors Rs according to pivotal row order */
    /* ---------------------------------------------------------------------- */

    if (scale > 0)
    {
        for (k = 0 ; k < n ; k++)
        {
            REAL (X [k]) = Rs [Pnum [k]] ;
        }
        for (k = 0 ; k < n ; k++)
        {
            Rs [k] = REAL (X [k]) ;
        }
    }

#ifndef NDEBUG
    ASSERT (Offp [n] == poff) ;
    ASSERT (Symbolic->nzoff == poff) ;
    PRINTF (("\n------------------- Off diagonal entries, new:\n")) ;
    ASSERT (KLU_valid (n, Offp, Offi, Offx)) ;
    if (Common->status == KLU_OK)
    {
        PRINTF (("\n ########### KLU_BTF_REFACTOR done, nblocks %d\n",nblocks));
        for (block = 0 ; block < nblocks ; block++)
        {
            k1 = R [block] ;
            k2 = R [block+1] ;
            nk = k2 - k1 ;
            PRINTF ((
                "\n================KLU_refactor output: k1 %d k2 %d nk %d\n",
                k1, k2, nk)) ;
            if (nk == 1)
            {
                PRINTF (("singleton  ")) ;
                PRINT_ENTRY (Udiag [k1]) ;
            }
            else
            {
                Lip = Numeric->Lip + k1 ;
                Llen = Numeric->Llen + k1 ;
                LU = (Unit *) Numeric->LUbx [block] ;
                PRINTF (("\n---- L block %d\n", block)) ;
                ASSERT (KLU_valid_LU (nk, TRUE, Lip, Llen, LU)) ;
                Uip = Numeric->Uip + k1 ;
                Ulen = Numeric->Ulen + k1 ;
                PRINTF (("\n---- U block %d\n", block)) ;
                ASSERT (KLU_valid_LU (nk, FALSE, Uip, Ulen, LU)) ;
            }
        }
    }
#endif

    if (klu_rf_trace) { klu_rf_t_tail += klu_rf_now() - klu_rf_t0; klu_rf_calls++; }
    return (TRUE) ;
}
