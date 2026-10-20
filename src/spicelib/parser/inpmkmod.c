/**********
Copyright 1990 Regents of the University of California.  All rights reserved.
Author: 1985 Thomas L. Quarles
**********/

#include "ngspice/ngspice.h"
#include <stdio.h>
#include "ngspice/inpdefs.h"
#include "ngspice/iferrmsg.h"
#include "ngspice/hash.h"
#include "inpxx.h"

/*  global input model table.  */
INPmodel *modtab = NULL;
/* Global input model hash table.
   The modelname is the key, the return value is the pointer to the model. */
NGHASHPTR modtabhash = NULL;
/* fork-lite: base model name (binning suffix ".N" stripped) -> newest
   model with that base; the others hang off INPnextBin. INPgetModBin walks
   this chain instead of the whole model table, which with modern PDKs
   (one wrapper subcircuit per transistor, each carrying its binned
   models) made the lookup quadratic in the number of transistors. */
NGHASHPTR modbinhash = NULL;

static char *inp_model_base_name(const char *name)
{
    const char *dot = strrchr(name, '.');
    const char *p;
    if (!dot || !dot[1])
        return copy(name);
    for (p = dot + 1; *p; p++)
        if (!isdigit_c(*p))
            return copy(name);
    return copy_substring(name, dot);
}

INPmodel *INPfindModBins(const char *name)
{
    return modbinhash ? (INPmodel *) nghash_find(modbinhash, (void *) name) : NULL;
}

/*--------------------------------------------------------------
 * This fcn takes the model name and looks to see if it is already
 * in the model table.  If it is, then just return.  Otherwise,
 * stick the model into the model table.
 * Note that the model table INPmodel *modtab is a linked list,
 * in parallel a hash table modtabhash is filled in for faster
 * access to modtab elements by giving the model name.
 *--------------------------------------------------------------*/

int INPmakeMod(char *token, int type, struct card *line)
{
   register INPmodel *newm;
   /* Initialze the hash table. The default key type is string.
      The default comparison function is strcmp.*/
   if (!modtabhash) {
       modtabhash = nghash_init(NGHASH_MIN_SIZE);
       nghash_unique(modtabhash, TRUE);
   }
   /* If the model is already there, just return. */
   else if (nghash_find(modtabhash, token))
       return (OK);

   /* Model name was not already in model table. Therefore stick
      it in the front of the model table, also into the model hash table.
      Then return.  */

#ifdef TRACE
   /* debug statement */
   printf("In INPmakeMod, about to insert new model name = %s . . .\n", token);
#endif

   newm = TMALLOC(INPmodel, 1);
   if (newm == NULL)
      return (E_NOMEM);

   newm->INPmodName = token;                 /* model name */
   newm->INPmodType = type;                  /* model type */
   newm->INPnextModel = modtab;              /* pointer to second model */
   newm->INPmodLine = line;                  /* model line */
   newm->INPmodfast = NULL;

   nghash_insert(modtabhash, token, newm);
   {
       char *base = inp_model_base_name(token);
       if (!modbinhash) {
           modbinhash = nghash_init(NGHASH_MIN_SIZE);
           nghash_unique(modbinhash, TRUE);
       }
       newm->INPnextBin = (INPmodel *) nghash_find(modbinhash, base);
       if (newm->INPnextBin)
           nghash_delete(modbinhash, base);
       nghash_insert(modbinhash, base, newm);   /* the table keeps the key string */
   }

   modtab = newm;

   return (OK);
}

