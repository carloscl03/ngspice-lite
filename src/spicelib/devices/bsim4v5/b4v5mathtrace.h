/* b4v5mathtrace.h -- fork-lite diagnostic, compiled only with -DNG_MATHTRACE.
 *
 * Every exp/log/pow/sqrt in b4v5ld.c goes through mt_call(), which records, per
 * source line:
 *   calls          total calls
 *   same_as_last   calls whose argument(s) are bit-equal to the previous call of
 *                  the same line for the same instance (unchanged between
 *                  Newton iterations)
 *   inv_inst       instances whose argument never changed over the whole run
 *                  (only depends on parameters/temperature: can be precomputed)
 *   dup            calls whose function and argument(s) are bit-equal to a call
 *                  made earlier IN THE SAME EVALUATION from another line
 *                  (and that line), i.e. a recomputation
 * Report at exit on stderr ("MATHTRACE ..."), single-threaded runs only.
 */
#ifndef B4V5MATHTRACE_H
#define B4V5MATHTRACE_H
#include <math.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#define MT_MAXLINE 8192
#define MT_MAXINST 65536
enum { MT_EXP, MT_LOG, MT_POW, MT_SQRT, MT_NF };
static const char *mt_fname[MT_NF] = { "exp", "log", "pow", "sqrt" };

struct mt_line {
    long calls, same, dup, zero;
    int fn, dupwith;
};
static struct mt_line mt_l[MT_MAXLINE];

/* per (line, instance) last argument and "ever changed" flag, in a hash */
struct mt_cell { int line; int inst; double a, b; unsigned char used, changed; };
#define MT_HBITS 21
static struct mt_cell *mt_h;
static const void *mt_inst_ptr[MT_MAXINST];
static int mt_ninst, mt_cur;
static long mt_evals;
/* calls of the current evaluation, for the duplicate test */
struct mt_ev { int fn, line; double a, b; };
static struct mt_ev mt_evbuf[512];
static int mt_nev;
static int mt_on = -1;

static void mt_report(void)
{
    int ln, k;
    long tot[MT_NF] = {0}, invl;
    if (!mt_evals) return;
    for (ln = 0; ln < MT_MAXLINE; ln++)
        if (mt_l[ln].calls) tot[mt_l[ln].fn] += mt_l[ln].calls;
    fprintf(stderr, "MATHTRACE evaluations=%ld instances=%d  per evaluation: exp %.2f log %.2f pow %.2f sqrt %.2f\n",
            mt_evals, mt_ninst, (double) tot[MT_EXP] / mt_evals, (double) tot[MT_LOG] / mt_evals,
            (double) tot[MT_POW] / mt_evals, (double) tot[MT_SQRT] / mt_evals);
    fprintf(stderr, "MATHTRACE   line  fn    calls/eval  same_as_last  invariant_inst  dup(with line)\n");
    for (ln = 0; ln < MT_MAXLINE; ln++) {
        struct mt_line *l = &mt_l[ln];
        if (!l->calls) continue;
        invl = 0;
        for (k = 0; k < (1 << MT_HBITS); k++)
            if (mt_h[k].used && mt_h[k].line == ln && !mt_h[k].changed) invl++;
        fprintf(stderr, "MATHTRACE   %5d %-5s %9.3f  %10.1f%%  %6ld/%-6d zero %5.1f%%  %7.1f%%%s",
                ln, mt_fname[l->fn], (double) l->calls / mt_evals, 100.0 * l->same / l->calls,
                invl, mt_ninst, 100.0 * l->zero / l->calls, 100.0 * l->dup / l->calls, l->dup ? "" : "\n");
        if (l->dup) fprintf(stderr, " (L%d)\n", l->dupwith);
    }
}

static int mt_inst_index(const void *p)
{
    int i;
    for (i = 0; i < mt_ninst; i++) if (mt_inst_ptr[i] == p) return i;
    if (mt_ninst < MT_MAXINST) mt_inst_ptr[mt_ninst] = p;
    return mt_ninst++;
}

static inline void mt_begin(const void *inst)
{
    if (mt_on < 0) {
        mt_on = getenv("NGSPICE_MATHTRACE") != NULL;
        if (mt_on) { mt_h = calloc((size_t) 1 << MT_HBITS, sizeof *mt_h); atexit(mt_report); }
    }
    if (!mt_on) return;
    mt_cur = mt_inst_index(inst);
    mt_nev = 0;
    mt_evals++;
}

static void mt_record(int fn, int line, double a, double b)
{
    struct mt_line *l;
    uint64_t key;
    unsigned h;
    int i;
    if (!mt_on || line <= 0 || line >= MT_MAXLINE) return;
    l = &mt_l[line];
    l->calls++; l->fn = fn;
    if (fn == MT_EXP && a == 0.0) l->zero++;
    /* per instance history */
    key = ((uint64_t) line << 20) ^ (uint64_t) mt_cur;
    h = (unsigned) ((key * 0x9E3779B97F4A7C15ull) >> (64 - MT_HBITS));
    for (;;) {
        struct mt_cell *c = &mt_h[h];
        if (!c->used) { c->used = 1; c->line = line; c->inst = mt_cur; c->a = a; c->b = b; break; }
        if (c->line == line && c->inst == mt_cur) {
            if (memcmp(&c->a, &a, sizeof a) == 0 && memcmp(&c->b, &b, sizeof b) == 0) l->same++;
            else { c->changed = 1; c->a = a; c->b = b; }
            break;
        }
        h = (h + 1) & ((1u << MT_HBITS) - 1);
    }
    /* duplicate within this evaluation, from another line */
    for (i = 0; i < mt_nev; i++)
        if (mt_evbuf[i].fn == fn && mt_evbuf[i].line != line &&
            memcmp(&mt_evbuf[i].a, &a, sizeof a) == 0 && memcmp(&mt_evbuf[i].b, &b, sizeof b) == 0) {
            l->dup++; l->dupwith = mt_evbuf[i].line; break;
        }
    if (mt_nev < 512) { mt_evbuf[mt_nev].fn = fn; mt_evbuf[mt_nev].line = line; mt_evbuf[mt_nev].a = a; mt_evbuf[mt_nev].b = b; mt_nev++; }
}

static inline double mt_exp(int line, double x)  { mt_record(MT_EXP, line, x, 0); return (exp)(x); }
static inline double mt_log(int line, double x)  { mt_record(MT_LOG, line, x, 0); return (log)(x); }
static inline double mt_pow(int line, double x, double y) { mt_record(MT_POW, line, x, y); return (pow)(x, y); }
static inline double mt_sqrt(int line, double x) { mt_record(MT_SQRT, line, x, 0); return (sqrt)(x); }

#define exp(x)    mt_exp(__LINE__, (x))
#define log(x)    mt_log(__LINE__, (x))
#define pow(x, y) mt_pow(__LINE__, (x), (y))
#define sqrt(x)   mt_sqrt(__LINE__, (x))
#endif
