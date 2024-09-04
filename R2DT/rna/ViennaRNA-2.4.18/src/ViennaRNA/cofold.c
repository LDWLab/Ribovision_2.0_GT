/*
 *                minimum free energy
 *                RNA secondary structure prediction
 *
 *                c Ivo Hofacker, Chrisoph Flamm
 *                original implementation by
 *                Walter Fontana
 *
 *                Vienna RNA package
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <string.h>
#include <limits.h>

#include "ViennaRNA/utils/basic.h"
#include "ViennaRNA/utils/strings.h"
#include "ViennaRNA/utils/structures.h"
#include "ViennaRNA/params/default.h"
#include "ViennaRNA/fold_vars.h"
#include "ViennaRNA/params/basic.h"
#include "ViennaRNA/subopt.h"
#include "ViennaRNA/fold.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/gquad.h"
#include "ViennaRNA/alphabet.h"
#include "ViennaRNA/cofold.h"

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

#ifdef _OPENMP
#include <omp.h>
#endif

#endif

#define MAXSECTORS        500     /* dimension for a backtrack array */

/*
 #################################
 # GLOBAL VARIABLES              #
 #################################
 */


/*
 #################################
 # PRIVATE VARIABLES             #
 #################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/* some backward compatibility stuff */
PRIVATE int                   backward_compat           = 0;
PRIVATE vrna_fold_compound_t  *backward_compat_compound = NULL;

PRIVATE float                 mfe1, mfe2; /* minimum free energies of the monomers */

#ifdef _OPENMP

#pragma omp threadprivate(mfe1, mfe2, backward_compat_compound, backward_compat)

#endif

#endif

/*
 #################################
 # PRIVATE FUNCTION DECLARATIONS #
 #################################
 */

PRIVATE int
backtrack(sect                  bt_stack[],
          vrna_bp_stack_t       *bp_list,
          vrna_fold_compound_t  *vc);


PRIVATE int
fill_arrays(vrna_fold_compound_t  *vc,
            int                   zuker);


PRIVATE void
free_end(int                  *array,
         int                  i,
         int                  start,
         vrna_fold_compound_t *vc);


PRIVATE void
doubleseq(vrna_fold_compound_t *vc);                /* do magic */


PRIVATE void
halfseq(vrna_fold_compound_t *vc);                  /* undo magic */


#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

/* wrappers for old API compatibility */
PRIVATE void
wrap_array_export(int   **f5_p,
                  int   **c_p,
                  int   **fML_p,
                  int   **fM1_p,
                  int   **fc_p,
                  int   **indx_p,
                  char  **ptype_p);


PRIVATE float
wrap_cofold(const char    *string,
            char          *structure,
            vrna_param_t  *parameters,
            int           is_constrained);


PRIVATE SOLUTION *
wrap_zukersubopt(const char   *string,
                 vrna_param_t *parameters);


#endif

/*
 #################################
 # BEGIN OF FUNCTION DEFINITIONS #
 #################################
 */
PUBLIC float
vrna_mfe_dimer(vrna_fold_compound_t *vc,
               char                 *structure)
{
  int             length, energy;
  char            *s;
  sect            bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  vrna_bp_stack_t *bp;

  length = (int)vc->length;

  vc->sequence_encoding[0] = vc->sequence_encoding2[0]; /* store length at pos. 0 in S1 too */

  if (!vrna_fold_compound_prepare(vc, VRNA_OPTION_MFE | VRNA_OPTION_HYBRID)) {
    vrna_message_warning("vrna_mfe_dimer@cofold.c: Failed to prepare vrna_fold_compound");
    return (float)(INF / 100.);
  }

  /* call user-defined recursion status callback function */
  if (vc->stat_cb)
    vc->stat_cb(VRNA_STATUS_MFE_PRE, vc->auxdata);

  energy = fill_arrays(vc, 0);

  /* call user-defined recursion status callback function */
  if (vc->stat_cb)
    vc->stat_cb(VRNA_STATUS_MFE_POST, vc->auxdata);

  if (structure && vc->params->model_details.backtrack) {
    bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + length / 2))); /* add a guess of how many G's may be involved in a G quadruplex */

    backtrack(bt_stack, bp, vc);

    s = vrna_db_from_bp_stack(bp, length);
    strncpy(structure, s, length + 1);
    free(s);
    free(bp);
  }

  if (vc->params->model_details.backtrack_type == 'C')
    return (float)vc->matrices->c[vc->jindx[length] + 1] / 100.;
  else if (vc->params->model_details.backtrack_type == 'M')
    return (float)vc->matrices->fML[vc->jindx[length] + 1] / 100.;
  else
    return (float)energy / 100.;
}


PRIVATE int
fill_arrays(vrna_fold_compound_t  *vc,
            int                   zuker)
{
  /* fill "c", "fML" and "f5" arrays and return  optimal energy */

  unsigned int  strands, *sn, *ss, *se, *so;
  int           i, j, length, energy;
  int           uniq_ML;
  int           no_close, type, maxj, *indx;
  int           *my_f5, *my_c, *my_fML, *my_fM1, *my_fc;
  int           *cc, *cc1;  /* auxilary arrays for canonical structures     */
  int           *Fmi;       /* holds row i of fML (avoids jumps in memory)  */
  int           *DMLi;      /* DMLi[j] holds  MIN(fML[i,k]+fML[k+1,j])      */
  int           *DMLi1;     /*                MIN(fML[i+1,k]+fML[k+1,j])    */
  int           *DMLi2;     /*                MIN(fML[i+2,k]+fML[k+1,j])    */

  int           dangle_model, noGUclosure, noLP, hc_decompose, turn;
  char          *ptype;
  unsigned char *hard_constraints;
  vrna_param_t  *P;
  vrna_mx_mfe_t *matrices;
  vrna_hc_t     *hc;

  length            = (int)vc->length;
  ptype             = vc->ptype;
  indx              = vc->jindx;
  P                 = vc->params;
  dangle_model      = P->model_details.dangles;
  noGUclosure       = P->model_details.noGUclosure;
  noLP              = P->model_details.noLP;
  uniq_ML           = P->model_details.uniq_ML;
  strands           = vc->strands;
  sn                = vc->strand_number;
  ss                = vc->strand_start;
  se                = vc->strand_end;
  so                = vc->strand_order;
  hc                = vc->hc;
  hard_constraints  = hc->mx;
  matrices          = vc->matrices;
  my_f5             = matrices->f5;
  my_c              = matrices->c;
  my_fML            = matrices->fML;
  my_fM1            = matrices->fM1;
  my_fc             = matrices->fc;
  turn              = P->model_details.min_loop_size;

  /* allocate memory for all helper arrays */
  cc    = (int *)vrna_alloc(sizeof(int) * (length + 2));
  cc1   = (int *)vrna_alloc(sizeof(int) * (length + 2));
  Fmi   = (int *)vrna_alloc(sizeof(int) * (length + 1));
  DMLi  = (int *)vrna_alloc(sizeof(int) * (length + 1));
  DMLi1 = (int *)vrna_alloc(sizeof(int) * (length + 1));
  DMLi2 = (int *)vrna_alloc(sizeof(int) * (length + 1));


  /* hard code min_loop_size to 0, since we can not be sure yet that this is already the case */
  turn = 0;

  for (j = 1; j <= length; j++) {
    Fmi[j]    = DMLi[j] = DMLi1[j] = DMLi2[j] = INF;
    my_fc[j]  = 0;
  }

  for (j = 1; j <= length; j++)
    for (i = 1; i <= j; i++) {
      my_c[indx[j] + i] = my_fML[indx[j] + i] = INF;
      if (uniq_ML)
        my_fM1[indx[j] + i] = INF;
    }

  for (i = length - turn - 1; i >= 1; i--) {
    /* i,j in [1..length] */

    maxj = (zuker) ? (MIN2(i + se[so[0]], length)) : length;
    for (j = i + turn + 1; j <= maxj; j++) {
      int ij;
      ij            = indx[j] + i;
      type          = vrna_get_ptype(ij, ptype);
      hc_decompose  = hard_constraints[length * i + j];
      energy        = INF;

      no_close = (((type == 3) || (type == 4)) && noGUclosure);

      if (hc_decompose) {
        /* we have a pair */
        int new_c = INF;

        if (!no_close) {
          /* check for hairpin loop */
          energy  = vrna_E_hp_loop(vc, i, j);
          new_c   = MIN2(new_c, energy);

          /* check for multibranch loops */
          energy  = vrna_E_mb_loop_fast(vc, i, j, DMLi1, DMLi2);
          new_c   = MIN2(new_c, energy);
        }

        if (dangle_model == 3) {
          /* coaxial stacking */
          energy  = vrna_E_mb_loop_stack(vc, i, j);
          new_c   = MIN2(new_c, energy);
        }

        /* check for interior loops */
        energy  = vrna_E_int_loop(vc, i, j);
        new_c   = MIN2(new_c, energy);

        /* remember stack energy for --noLP option */
        if (noLP) {
          if ((sn[i] == sn[i + 1]) && (sn[j - 1] == sn[j])) {
            int stackEnergy = vrna_E_stack(vc, i, j);
            new_c     = MIN2(new_c, cc1[j - 1] + stackEnergy);
            my_c[ij]  = cc1[j - 1] + stackEnergy;
          } else {
            /* currently we don't allow stacking over the cut point */
            my_c[ij] = FORBIDDEN;
          }

          cc[j] = new_c;
        } else {
          my_c[ij] = new_c;
        }
      } /* end >> if (pair) << */
      else {
        my_c[ij] = INF;
      }

      /*
       * done with c[i,j], now compute fML[i,j]
       * free ends ? -----------------------------------------
       */

      my_fML[ij] = vrna_E_ml_stems_fast(vc, i, j, Fmi, DMLi);

      if (uniq_ML)   /* compute fM1 for unique decomposition */
        my_fM1[ij] = E_ml_rightmost_stem(i, j, vc);
    }

    if (i == se[so[0]] + 1)
      for (j = i; j <= maxj; j++)
        free_end(my_fc, j, ss[so[1]], vc);

    if (i <= se[so[0]])
      free_end(my_fc, i, se[so[0]], vc);

    {
      int *FF; /* rotate the auxilliary arrays */
      FF    = DMLi2;
      DMLi2 = DMLi1;
      DMLi1 = DMLi;
      DMLi  = FF;
      FF    = cc1;
      cc1   = cc;
      cc    = FF;
      for (j = 1; j <= maxj; j++)
        cc[j] = Fmi[j] = DMLi[j] = INF;
    }
  }

  /* calculate energies of 5' and 3' fragments */

  for (i = 1; i <= length; i++)
    free_end(my_f5, i, 1, vc);

  if (strands > 1) {
    mfe1  = my_f5[se[so[0]]];
    mfe2  = my_fc[length];
    /* add DuplexInit, check whether duplex*/
    for (i = ss[so[1]]; i <= length; i++) {
      if (my_f5[i] != INF)
        my_f5[i] += P->DuplexInit;

      if ((my_fc[i] != INF) &&
          (my_fc[1] != INF)) {
        energy  = my_fc[i] +
                  my_fc[1];

        if (energy < my_f5[i])
          my_f5[i] = energy;
      }
    }
  }

  energy = my_f5[length];
  if (strands == 1)
    mfe1 = mfe2 = energy;

  /* clean up memory */
  free(cc);
  free(cc1);
  free(Fmi);
  free(DMLi);
  free(DMLi1);
  free(DMLi2);

  return energy;
}


PRIVATE int
backtrack_co(sect                 bt_stack[],
             vrna_bp_stack_t      *bp_list,
             int                  s,
             int                  b, /* b=0: start new structure, b \ne 0: add to existing structure */
             vrna_fold_compound_t *vc)
{
  /*------------------------------------------------------------------
   *  trace back through the "c", "fc", "f5" and "fML" arrays to get the
   *  base pairing list. No search for equivalent structures is done.
   *  This is fast, since only few structure elements are recalculated.
   *  ------------------------------------------------------------------*/

  unsigned int  *se, *so;
  int           i, j, ij, k, length, no_close, type, ret;
  char          *string = vc->sequence;
  vrna_param_t  *P      = vc->params;
  int           *indx   = vc->jindx;
  char          *ptype  = vc->ptype;

  int           noLP            = P->model_details.noLP;
  int           noGUclosure     = P->model_details.noGUclosure;
  char          backtrack_type  = P->model_details.backtrack_type;

  /* the folding matrices */
  int           *my_c;

  ret     = 1;
  length  = vc->length;
  my_c    = vc->matrices->c;
  se      = vc->strand_end;
  so      = vc->strand_order;

  /* int   b=0;*/

  length = strlen(string);
  if (s == 0) {
    bt_stack[++s].i = 1;
    bt_stack[s].j   = length;
    bt_stack[s].ml  = (backtrack_type == 'M') ? 1 : ((backtrack_type == 'C') ? 2 : 0);
  }

  while (s > 0) {
    int ml, cij;
    int canonical = 1;     /* (i,j) closes a canonical structure */

    /* pop one element from stack */
    i   = bt_stack[s].i;
    j   = bt_stack[s].j;
    ml  = bt_stack[s--].ml;

    switch (ml) {
      /* backtrack in f5 */
      case 0:
      {
        int p, q;
        if (vrna_BT_ext_loop_f5(vc, &j, &p, &q, bp_list, &b)) {
          if (j > 0) {
            bt_stack[++s].i = 1;
            bt_stack[s].j   = j;
            bt_stack[s].ml  = 0;
          }

          if (p > 0) {
            i = p;
            j = q;
            goto repeat1;
          }

          continue;
        } else {
          vrna_message_warning("backtracking failed in f5, segment [%d,%d]\n", i, j);
          ret = 0;
          goto backtrack_exit;
        }
      }
      break;

      /* true multi-loop backtrack in fML */
      case 1:
      {
        int p, q, comp1, comp2;
        if (vrna_BT_mb_loop_split(vc, &i, &j, &p, &q, &comp1, &comp2, bp_list, &b)) {
          if (i > 0) {
            bt_stack[++s].i = i;
            bt_stack[s].j   = j;
            bt_stack[s].ml  = comp1;
          }

          if (p > 0) {
            bt_stack[++s].i = p;
            bt_stack[s].j   = q;
            bt_stack[s].ml  = comp2;
          }

          continue;
        } else {
          vrna_message_warning("backtrack failed in fML\n%s", string);
          ret = 0;
          goto backtrack_exit;
        }
      }
      break;

      case 2:
        bp_list[++b].i  = i;
        bp_list[b].j    = j;
        goto repeat1;

      /* backtrack fake-multi loop parts */
      case 3:
      case 4:
      {
        int lower, k, p, q;
        p     = i;
        q     = j;
        lower = (i <= se[so[0]]) ? 1 : 0;

        if (vrna_BT_mb_loop_fake(vc, &k, &i, &j, bp_list, &b)) {
          if (k > 0) {
            bt_stack[++s].i = (lower) ? k : p;
            bt_stack[s].j   = (lower) ? q : k;
            bt_stack[s].ml  = ml;
          }

          if (i > 0)
            goto repeat1;

          continue;
        } else {
          vrna_message_warning("backtrack failed in fc\n%s", string);
          ret = 0;
          goto backtrack_exit;
        }
      }
      break;
    } /* end of switch(ml) */

repeat1:

    /*----- begin of "repeat:" -----*/
    ij = indx[j] + i;

    if (canonical)
      cij = my_c[ij];

    type = vrna_get_ptype(ij, ptype);

    if (noLP) {
      if (vrna_BT_stack(vc, &i, &j, &cij, bp_list, &b)) {
        canonical = 0;
        goto repeat1;
      }
    }

    canonical = 1;

    no_close = (((type == 3) || (type == 4)) && noGUclosure);
    if (no_close) {
      if (cij == FORBIDDEN)
        continue;
    } else {
      if (vrna_BT_hp_loop(vc, i, j, cij, bp_list, &b))
        continue;
    }

    if (vrna_BT_int_loop(vc, &i, &j, cij, bp_list, &b)) {
      if (i < 0)
        continue;
      else
        goto repeat1;
    }

    /* (i.j) must close a fake or true multi-loop */
    int comp1, comp2;

    if (vrna_BT_mb_loop(vc, &i, &j, &k, cij, &comp1, &comp2)) {
      bt_stack[++s].i = i;
      bt_stack[s].j   = k;
      bt_stack[s].ml  = comp1;
      bt_stack[++s].i = k + 1;
      bt_stack[s].j   = j;
      bt_stack[s].ml  = comp2;
    } else {
      vrna_message_warning("backtracking failed in repeat, segment [%d,%d]\n", i, j);
      ret = 0;
      goto backtrack_exit;
    }

    /* end of repeat: --------------------------------------------------*/
  } /* end >> while (s>0) << */

backtrack_exit:

  bp_list[0].i = b;    /* save the total number of base pairs */
  return ret;
}


PRIVATE void
free_end(int                  *array,
         int                  i,
         int                  start,
         vrna_fold_compound_t *vc)
{
  unsigned int  *sn;
  int           inc, type, energy, en, length, j, left, right, dangle_model, with_gquad, *indx,
                *c, *ggg, turn;
  vrna_param_t  *P;
  short         *S1;
  char          *ptype;
  unsigned char *hard_constraints;
  vrna_mx_mfe_t *matrices;
  vrna_hc_t     *hc;
  vrna_sc_t     *sc;

  P                 = vc->params;
  dangle_model      = P->model_details.dangles;
  with_gquad        = P->model_details.gquad;
  turn              = P->model_details.min_loop_size;
  inc               = (i > start) ? 1 : -1;
  length            = (int)vc->length;
  S1                = vc->sequence_encoding;
  ptype             = vc->ptype;
  indx              = vc->jindx;
  sn                = vc->strand_number;
  matrices          = vc->matrices;
  c                 = matrices->c;
  ggg               = matrices->ggg;
  hc                = vc->hc;
  sc                = vc->sc;
  hard_constraints  = hc->mx;

  if (hc->up_ext[i]) {
    if (i == start)
      array[i] = 0;
    else
      array[i] = array[i - inc];

    if (sc) {
      if (sc->energy_up)
        array[i] += sc->energy_up[i][1];

      if (sc->f)
        array[i] += sc->f(start, i, start, i - 1, VRNA_DECOMP_EXT_EXT, sc->data);
    }
  } else {
    array[i] = INF;
  }

  if (inc > 0) {
    left  = start;
    right = i;
  } else {
    left  = i;
    right = start;
  }

  /* hard code min_loop_size to 0, since we can not be sure yet that this is already the case */
  turn = 0;

  for (j = start; inc * (i - j) > turn; j += inc) {
    int   ii, jj;
    short si, sj;
    if (i > j) {
      ii  = j;
      jj  = i;
    }                           /* inc>0 */
    else {
      ii  = i;
      jj  = j;
    }                           /* inc<0 */

    if (hard_constraints[length * ii + jj] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP) {
      type    = vrna_get_ptype(indx[jj] + ii, ptype);
      si      = ((ii > 1) && (sn[ii - 1] == sn[ii])) ? S1[ii - 1] : -1;
      sj      = ((jj < length) && (sn[jj] == sn[jj + 1])) ? S1[jj + 1] : -1;
      energy  = c[indx[jj] + ii];

      if ((sc) && (sc->f))
        energy += sc->f(start, jj, ii - 1, ii, VRNA_DECOMP_EXT_EXT_STEM, sc->data);

      if (energy != INF) {
        switch (dangle_model) {
          case 0:
            if (array[j - inc] != INF) {
              en        = array[j - inc] + energy + vrna_E_ext_stem(type, -1, -1, P);
              array[i]  = MIN2(array[i], en);
            }

            break;
          case 2:
            if (array[j - inc] != INF) {
              en        = array[j - inc] + energy + vrna_E_ext_stem(type, si, sj, P);
              array[i]  = MIN2(array[i], en);
            }

            break;
          default:
            if (array[j - inc] != INF) {
              en        = array[j - inc] + energy + vrna_E_ext_stem(type, -1, -1, P);
              array[i]  = MIN2(array[i], en);
            }

            if (inc > 0) {
              if (j > left) {
                if (hc->up_ext[ii - 1]) {
                  if (array[j - 2] != INF) {
                    en = array[j - 2] + energy + vrna_E_ext_stem(type, si, -1, P);
                    if (sc)
                      if (sc->energy_up)
                        en += sc->energy_up[ii - 1][1];

                    array[i] = MIN2(array[i], en);
                  }
                }
              }
            } else if (j < right) {
              if (hc->up_ext[jj + 1]) {
                if (array[j + 2] != INF) {
                  en = array[j + 2] + energy + vrna_E_ext_stem(type, -1, sj, P);
                  if (sc)
                    if (sc->energy_up)
                      en += sc->energy_up[jj + 1][1];

                  array[i] = MIN2(array[i], en);
                }
              }
            }

            break;
        }
      }
    }

    if (with_gquad) {
      if (sn[ii] == sn[jj])
        if (array[j - inc] != INF)
          array[i] = MIN2(array[i], array[j - inc] + ggg[indx[jj] + ii]);
    }

    if (dangle_model % 2 == 1) {
      /* interval ends in a dangle (i.e. i-inc is paired) */
      if (i > j) {
        ii  = j;
        jj  = i - 1;
      }                             /* inc>0 */
      else {
        ii  = i + 1;
        jj  = j;
      }                             /* inc<0 */

      if (!(hard_constraints[length * ii + jj] & VRNA_CONSTRAINT_CONTEXT_EXT_LOOP))
        continue;

      type    = vrna_get_ptype(indx[jj] + ii, ptype);
      si      = (ii > left) && (sn[ii - 1] == sn[ii]) ? S1[ii - 1] : -1;
      sj      = (jj < right) && (sn[jj] == sn[jj + 1]) ? S1[jj + 1] : -1;
      energy  = c[indx[jj] + ii];
      if (energy != INF) {
        if (inc > 0) {
          if (hc->up_ext[jj - 1]) {
            if (array[j - inc] != INF) {
              en = array[j - inc] + energy + vrna_E_ext_stem(type, -1, sj, P);
              if (sc)
                if (sc->energy_up)
                  en += sc->energy_up[jj + 1][1];

              array[i] = MIN2(array[i], en);
            }
          }
        } else {
          if (hc->up_ext[ii - 1]) {
            if (array[j - inc] != INF) {
              en = array[j - inc] + energy + vrna_E_ext_stem(type, si, -1, P);
              if (sc)
                if (sc->energy_up)
                  en += sc->energy_up[ii - 1][1];

              array[i] = MIN2(array[i], en);
            }
          }
        }

        if (j != start) {
          /* dangle_model on both sides */
          if (hc->up_ext[jj - 1] && hc->up_ext[ii - 1]) {
            if (array[j - 2 * inc] != INF) {
              en = array[j - 2 * inc] + energy + vrna_E_ext_stem(type, si, sj, P);
              if (sc)
                if (sc->energy_up)
                  en += sc->energy_up[ii - 1][1] + sc->energy_up[jj + 1][1];

              array[i] = MIN2(array[i], en);
            }
          }
        }
      }
    }
  }
}


PRIVATE int
backtrack(sect                  bt_stack[],
          vrna_bp_stack_t       *bp_list,
          vrna_fold_compound_t  *vc)
{
  /*routine to call backtrack_co from 1 to n, backtrack type??*/
  return backtrack_co(bt_stack, bp_list, 0, 0, vc);
}


PRIVATE void
doubleseq(vrna_fold_compound_t *vc)
{
  unsigned int length, i, s;

  length = vc->length;

  /* do some magic to re-use cofold code */
  vc->sequence = vrna_realloc(vc->sequence, sizeof(char) * (2 * length + 2));
  memcpy(vc->sequence + length, vc->sequence, sizeof(char) * length);
  vc->sequence[2 * length]  = '\0';
  vc->length                = (unsigned int)strlen(vc->sequence);
  vc->cutpoint              = length + 1;

  vc->strands = 2;

  free(vc->strand_number);
  vc->strand_number = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (vc->length + 1));
  for (s = i = 0; i <= vc->length; i++) {
    if (i == length + 1)
      s++;

    vc->strand_number[i] = s;
  }

  free(vc->strand_order);
  free(vc->strand_start);
  free(vc->strand_end);
  vc->strand_order    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (vc->strands + 1));
  vc->strand_start    = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (vc->strands + 1));
  vc->strand_end      = (unsigned int *)vrna_alloc(sizeof(unsigned int) * (vc->strands + 1));
  vc->strand_order[0] = 0;
  vc->strand_order[1] = 1;
  vc->strand_start[0] = 1;
  vc->strand_end[0]   = vc->strand_start[0] + length - 1;
  vc->strand_start[1] = vc->strand_end[0] + 1;
  vc->strand_end[1]   = vc->strand_start[1] + length - 1;


  vc->sequence_encoding = vrna_realloc(vc->sequence_encoding, sizeof(short) * (vc->length + 2));
  memcpy(vc->sequence_encoding + length + 1, vc->sequence_encoding + 1, sizeof(short) * length);
  vc->sequence_encoding[0]              = vc->sequence_encoding[vc->length];
  vc->sequence_encoding[vc->length + 1] = vc->sequence_encoding[1];

  vc->sequence_encoding2 = vrna_realloc(vc->sequence_encoding2, sizeof(short) * (vc->length + 2));
  memcpy(vc->sequence_encoding2 + length + 1, vc->sequence_encoding2 + 1, sizeof(short) * length);
  vc->sequence_encoding2[0]               = vc->length;
  vc->sequence_encoding2[vc->length + 1]  = 0;

  free(vc->ptype);
  vc->ptype = vrna_ptypes(vc->sequence_encoding2, &(vc->params->model_details));
  free(vc->iindx);
  vc->iindx = vrna_idx_row_wise(vc->length);
  free(vc->jindx);
  vc->jindx = vrna_idx_col_wise(vc->length);

  vrna_hc_init(vc);

  /* add DP matrices */
  vrna_mx_mfe_add(vc, VRNA_MX_DEFAULT, 0);
}


PRIVATE void
halfseq(vrna_fold_compound_t *vc)
{
  unsigned int halflength;

  halflength = vc->length / 2;

  vc->sequence              = vrna_realloc(vc->sequence, sizeof(char) * (halflength + 1));
  vc->sequence[halflength]  = '\0';
  vc->length                = (unsigned int)strlen(vc->sequence);
  vc->cutpoint              = -1;
  vc->strands               = 1;
  vc->strand_number         = (unsigned int *)vrna_realloc(vc->strand_number,
                                                           sizeof(unsigned int) * (vc->length + 1));
  vc->strand_order = (unsigned int *)vrna_realloc(vc->strand_order,
                                                  sizeof(unsigned int) * (vc->strands + 1));
  vc->strand_start = (unsigned int *)vrna_realloc(vc->strand_start,
                                                  sizeof(unsigned int) * (vc->strands + 1));
  vc->strand_end = (unsigned int *)vrna_realloc(vc->strand_end,
                                                sizeof(unsigned int) * (vc->strands + 1));

  vc->sequence_encoding =
    vrna_realloc(vc->sequence_encoding, sizeof(short) * (vc->length + 2));
  vc->sequence_encoding[0]              = vc->sequence_encoding[vc->length];
  vc->sequence_encoding[vc->length + 1] = vc->sequence_encoding[1];

  vc->sequence_encoding2 =
    vrna_realloc(vc->sequence_encoding2, sizeof(short) * (vc->length + 2));
  vc->sequence_encoding2[0]               = vc->length;
  vc->sequence_encoding2[vc->length + 1]  = 0;

  free(vc->ptype);
  vc->ptype = vrna_ptypes(vc->sequence_encoding2, &(vc->params->model_details));
  free(vc->iindx);
  vc->iindx = vrna_idx_row_wise(vc->length);
  free(vc->jindx);
  vc->jindx = vrna_idx_col_wise(vc->length);

  vrna_hc_init(vc);

  /* add DP matrices */
  vrna_mx_mfe_add(vc, VRNA_MX_DEFAULT, 0);
}


typedef struct {
  int i;
  int j;
  int e;
  int idxj;
} zuker_pair;

PRIVATE int
comp_pair(const void  *A,
          const void  *B)
{
  zuker_pair  *x, *y;
  int         ex, ey;

  x   = (zuker_pair *)A;
  y   = (zuker_pair *)B;
  ex  = x->e;
  ey  = y->e;
  if (ex > ey)
    return 1;

  if (ex < ey)
    return -1;

  return x->idxj + x->i - y->idxj + y->i;
}


PUBLIC SOLUTION *
vrna_subopt_zuker(vrna_fold_compound_t *vc)
{
  /* Compute zuker suboptimal. Here, we're abusing the cofold() code
   * "double" sequence, compute dimerarray entries, track back every base pair.
   * This is slightly wasteful compared to the normal solution */

  char            *structure, *mfestructure, **todo, *ptype;
  int             i, j, counter, num_pairs, psize, p, *indx, *c, turn;
  unsigned int    length, doublelength;
  float           energy;
  SOLUTION        *zukresults;
  vrna_bp_stack_t *bp_list;
  zuker_pair      *pairlist;
  sect            bt_stack[MAXSECTORS]; /* stack of partial structures for backtracking */
  vrna_mx_mfe_t   *matrices;
  vrna_md_t       *md;

  md    = &(vc->params->model_details);
  turn  = md->min_loop_size;

  /* do some magic to re-use cofold code although vc is single sequence */
  md->min_loop_size = 0;
  doubleseq(vc);

  if (!vrna_fold_compound_prepare(vc, VRNA_OPTION_MFE | VRNA_OPTION_HYBRID)) {
    vrna_message_warning("vrna_subopt_zuker@cofold.c: Failed to prepare vrna_fold_compound");
    return NULL;
  }

  doublelength    = vc->length;
  length          = doublelength / 2;
  indx            = vc->jindx;
  ptype           = vc->ptype;
  matrices        = vc->matrices;
  c               = matrices->c;
  num_pairs       = counter = 0;
  mfestructure    = (char *)vrna_alloc((unsigned)doublelength + 1);
  structure       = (char *)vrna_alloc((unsigned)doublelength + 1);
  zukresults      = (SOLUTION *)vrna_alloc(((length * (length - 1)) / 2) * sizeof(SOLUTION));
  mfestructure[0] = '\0';

  /* store length at pos. 0 */
  vc->sequence_encoding[0] = vc->sequence_encoding2[0];

  /* get mfe and do forward recursion */
  (void)fill_arrays(vc, 1);

  psize     = length;
  pairlist  = (zuker_pair *)vrna_alloc(sizeof(zuker_pair) * (psize + 1));
  bp_list   = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (1 + length / 2));
  todo      = (char **)vrna_alloc(sizeof(char *) * (length + 1));
  for (i = 1; i < length; i++)
    todo[i] = (char *)vrna_alloc(sizeof(char) * (length + 1));

  /* Make a list of all base pairs */
  for (i = 1; i < length; i++) {
    for (j = i + turn + 1 /*??*/; j <= length; j++) {
      if (ptype[indx[j] + i] == 0)
        continue;

      if (num_pairs >= psize) {
        psize     = 1.2 * psize + 32;
        pairlist  = vrna_realloc(pairlist, sizeof(zuker_pair) * (psize + 1));
      }

      pairlist[num_pairs].i       = i;
      pairlist[num_pairs].j       = j;
      pairlist[num_pairs].e       = c[indx[j] + i] + c[indx[i + length] + j];
      pairlist[num_pairs++].idxj  = indx[j];

      todo[i][j] = 1;
    }
  }

  qsort(pairlist, num_pairs, sizeof(zuker_pair), comp_pair);

  for (p = 0; p < num_pairs; p++) {
    i = pairlist[p].i;
    j = pairlist[p].j;
    if (todo[i][j]) {
      int   k;
      char  *sz;
      bt_stack[1].i   = i;
      bt_stack[1].j   = j;
      bt_stack[1].ml  = 2;
      backtrack_co(bt_stack, bp_list, 1, 0, vc);
      bt_stack[1].i   = j;
      bt_stack[1].j   = i + length;
      bt_stack[1].ml  = 2;
      backtrack_co(bt_stack, bp_list, 1, bp_list[0].i, vc);
      energy                          = pairlist[p].e;
      sz                              = vrna_db_from_bp_stack(bp_list, length);
      zukresults[counter].energy      = energy / 100.;
      zukresults[counter++].structure = sz;
      for (k = 1; k <= bp_list[0].i; k++) {
        /* mark all pairs in structure as done */
        int x, y;
        x = bp_list[k].i;
        y = bp_list[k].j;
        if (x > length)
          x -= length;

        if (y > length)
          y -= length;

        if (x > y) {
          int temp;
          temp  = x;
          x     = y;
          y     = temp;
        }

        todo[x][y] = 0;
      }
    }
  }

  /* clean up */
  free(pairlist);
  for (i = 1; i < length; i++)
    free(todo[i]);
  free(todo);
  free(structure);
  free(mfestructure);
  free(bp_list);

  /* undo magic */
  halfseq(vc);
  md->min_loop_size = turn;

  return zukresults;
}


/*
 *###########################################
 *# deprecated functions below              #
 *###########################################
 */

#ifndef VRNA_DISABLE_BACKWARD_COMPATIBILITY

PRIVATE void
wrap_array_export(int   **f5_p,
                  int   **c_p,
                  int   **fML_p,
                  int   **fM1_p,
                  int   **fc_p,
                  int   **indx_p,
                  char  **ptype_p)
{
  /* make the DP arrays available to routines such as subopt() */
  if (backward_compat_compound) {
    *f5_p     = backward_compat_compound->matrices->f5;
    *c_p      = backward_compat_compound->matrices->c;
    *fML_p    = backward_compat_compound->matrices->fML;
    *fM1_p    = backward_compat_compound->matrices->fM1;
    *fc_p     = backward_compat_compound->matrices->fc;
    *indx_p   = backward_compat_compound->jindx;
    *ptype_p  = backward_compat_compound->ptype;
  }
}


/*--------------------------------------------------------------------------*/

PRIVATE float
wrap_cofold(const char    *string,
            char          *structure,
            vrna_param_t  *parameters,
            int           is_constrained)
{
  unsigned int          length;
  char                  *seq;
  vrna_fold_compound_t  *vc;
  vrna_param_t          *P;
  float                 mfe;

  vc      = NULL;
  length  = strlen(string);

#ifdef _OPENMP
  /* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

  /* we need the parameter structure for hard constraints */
  if (parameters) {
    P = vrna_params_copy(parameters);
  } else {
    vrna_md_t md;
    set_model_details(&md);
    md.temperature  = temperature;
    P               = vrna_params(&md);
  }

  P->model_details.min_loop_size = 0;  /* set min loop length to 0 */

  /* dirty hack to reinsert the '&' according to the global variable 'cut_point' */
  seq = vrna_cut_point_insert(string, cut_point);

  /* get compound structure */
  vc = vrna_fold_compound(seq, &(P->model_details), 0);

  if (parameters) {
    /* replace params if necessary */
    free(vc->params);
    vc->params = P;
  } else {
    free(P);
  }

  /* handle hard constraints in pseudo dot-bracket format if passed via simple interface */
  if (is_constrained && structure) {
    unsigned int constraint_options = 0;
    constraint_options |= VRNA_CONSTRAINT_DB
                          | VRNA_CONSTRAINT_DB_PIPE
                          | VRNA_CONSTRAINT_DB_DOT
                          | VRNA_CONSTRAINT_DB_X
                          | VRNA_CONSTRAINT_DB_ANG_BRACK
                          | VRNA_CONSTRAINT_DB_RND_BRACK
                          | VRNA_CONSTRAINT_DB_INTRAMOL
                          | VRNA_CONSTRAINT_DB_INTERMOL;

    vrna_constraints_add(vc, (const char *)structure, constraint_options);
  }

  if (backward_compat_compound)
    vrna_fold_compound_free(backward_compat_compound);

  backward_compat_compound  = vc;
  backward_compat           = 1;

  /* cleanup */
  free(seq);

  /* call mfe_dimer without backtracing */
  mfe = vrna_mfe_dimer(vc, NULL);

  /* now we backtrace in a backward compatible way */
  if (structure && vc->params->model_details.backtrack) {
    char            *s;
    sect            bt_stack[MAXSECTORS];
    vrna_bp_stack_t *bp;

    bp = (vrna_bp_stack_t *)vrna_alloc(sizeof(vrna_bp_stack_t) * (4 * (1 + length / 2))); /* add a guess of how many G's may be involved in a G quadruplex */

    backtrack(bt_stack, bp, vc);

    s = vrna_db_from_bp_stack(bp, length);
    strncpy(structure, s, length + 1);
    free(s);

    if (base_pair)
      free(base_pair);

    base_pair = bp;
  }

  return mfe;
}


PRIVATE SOLUTION *
wrap_zukersubopt(const char   *string,
                 vrna_param_t *parameters)
{
  vrna_fold_compound_t  *vc;
  vrna_param_t          *P;

  vc = NULL;

#ifdef _OPENMP
  /* Explicitly turn off dynamic threads */
  omp_set_dynamic(0);
#endif

  /* we need the parameter structure for hard constraints */
  if (parameters) {
    P = vrna_params_copy(parameters);
  } else {
    vrna_md_t md;
    set_model_details(&md);
    md.temperature  = temperature;
    P               = vrna_params(&md);
  }

  /* get compound structure */
  vc = vrna_fold_compound(string, &(P->model_details), VRNA_OPTION_DEFAULT);

  if (parameters) {
    /* replace params if necessary */
    free(vc->params);
    vc->params = P;
  } else {
    free(P);
  }

  if (backward_compat_compound)
    vrna_fold_compound_free(backward_compat_compound);

  backward_compat_compound  = vc;
  backward_compat           = 1;

  return vrna_subopt_zuker(vc);
}


PUBLIC void
initialize_cofold(int length)
{
  /* DO NOTHING */ }


PUBLIC void
free_co_arrays(void)
{
  if (backward_compat_compound && backward_compat) {
    vrna_fold_compound_free(backward_compat_compound);
    backward_compat_compound  = NULL;
    backward_compat           = 0;
  }
}


/*--------------------------------------------------------------------------*/

PUBLIC void
export_cofold_arrays_gq(int   **f5_p,
                        int   **c_p,
                        int   **fML_p,
                        int   **fM1_p,
                        int   **fc_p,
                        int   **ggg_p,
                        int   **indx_p,
                        char  **ptype_p)
{
  /* make the DP arrays available to routines such as subopt() */
  wrap_array_export(f5_p, c_p, fML_p, fM1_p, fc_p, indx_p, ptype_p);
  if (backward_compat_compound)
    *ggg_p = backward_compat_compound->matrices->ggg;
}


PUBLIC void
export_cofold_arrays(int  **f5_p,
                     int  **c_p,
                     int  **fML_p,
                     int  **fM1_p,
                     int  **fc_p,
                     int  **indx_p,
                     char **ptype_p)
{
  wrap_array_export(f5_p, c_p, fML_p, fM1_p, fc_p, indx_p, ptype_p);
}


PUBLIC float
cofold(const char *string,
       char       *structure)
{
  return wrap_cofold(string, structure, NULL, fold_constrained);
}


PUBLIC float
cofold_par(const char   *string,
           char         *structure,
           vrna_param_t *parameters,
           int          is_constrained)
{
  return wrap_cofold(string, structure, parameters, is_constrained);
}


PUBLIC SOLUTION *
zukersubopt(const char *string)
{
  return wrap_zukersubopt(string, NULL);
}


PUBLIC SOLUTION *
zukersubopt_par(const char    *string,
                vrna_param_t  *parameters)
{
  return wrap_zukersubopt(string, parameters);
}


PUBLIC void
update_cofold_params(void)
{
  vrna_fold_compound_t *v;

  if (backward_compat_compound && backward_compat) {
    vrna_md_t md;
    v = backward_compat_compound;

    if (v->params)
      free(v->params);

    set_model_details(&md);
    v->params = vrna_params(&md);
  }
}


PUBLIC void
update_cofold_params_par(vrna_param_t *parameters)
{
  vrna_fold_compound_t *v;

  if (backward_compat_compound && backward_compat) {
    v = backward_compat_compound;

    if (v->params)
      free(v->params);

    if (parameters) {
      v->params = vrna_params_copy(parameters);
    } else {
      vrna_md_t md;
      set_model_details(&md);
      md.temperature  = temperature;
      v->params       = vrna_params(&md);
    }
  }
}


PUBLIC void
get_monomere_mfes(float *e1,
                  float *e2)
{
  /*exports monomere free energies*/
  *e1 = mfe1;
  *e2 = mfe2;
}


#endif
