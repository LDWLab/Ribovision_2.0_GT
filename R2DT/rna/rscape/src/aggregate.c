/* aggregate.c */

#include "rscape_config.h"

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <float.h>
#include <limits.h>

#include "easel.h"
#include "esl_alphabet.h"
#include "esl_dmatrix.h"
#include "esl_exponential.h"
#include "esl_gamma.h"
#include "esl_histogram.h"
#include "esl_msa.h"
#include "esl_msafile.h"
#include "esl_sqio.h"
#include "esl_stack.h"
#include "esl_stats.h"
#include "esl_tree.h"
#include "esl_vectorops.h"
#include "esl_wuss.h"

#include "aggregate.h"
#include "contactmap.h"
#include "covariation.h"
#include "correlators.h"
#include "structure.h"

static double *extract_rm_pvals  (RM *rm, SPAIR *spair, char *errbuf, int verbose);
static double *extract_rm_weights(RM *rm, SPAIR *spair, char *errbuf, int verbose);

int
agg_CalculatePvalues(SPAIR *spair, CTLIST *ctlist, RMLIST **ret_rmlist, int helix_unpaired, enum agg_e agg_method, char *errbuf, int verbose)
{
  RMLIST *rmlist = NULL;
  RM     *rm;
  double *pval = NULL;
  double *weights = NULL;
  double  pval_agg;
  int     n;
  int     status;

  rmlist = struct_rmlist_FromCTLIST(helix_unpaired, ctlist, errbuf, verbose);
  if (!rmlist)  ESL_XFAIL(eslFAIL, errbuf, "error int agg_CalculatePvalues()");

  for (n = 0; n < rmlist->nrm; n ++) {
    
    rm   = rmlist->rm[n];
    pval = extract_rm_pvals(rm, spair, errbuf, verbose);
    if (pval[0] < 0) continue; // no covariation analysis
    
    switch(agg_method) {
    case AGG_FISHER:
      status = agg_FISHER(rm->nbp, pval, &pval_agg, errbuf, verbose);
      break;
    case AGG_LANCASTER:
      weights = extract_rm_weights(rm, spair, errbuf, verbose);
     
      status = agg_LANCASTER(rm->nbp, pval, weights, &pval_agg, errbuf, verbose);
      free(weights); weights = NULL;
      break;
    case AGG_SIDAK:
      status = agg_SIDAK(rm->nbp, pval, &pval_agg, errbuf, verbose);
      break;
    case AGG_NONE:
      pval_agg = -1;
      status = eslOK;
      break;
    }
    if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "error in agg_CalculatePvalues()");

    rm->Pval = pval_agg;
    rm->Eval = rmlist->nrm * pval_agg;  // multiple test correction for testing nrm helices

    free(pval); pval = NULL;
  }
    
  rmlist->agg_method = agg_method;
    
  if (verbose) struct_rmlist_Dump(rmlist, NULL, 0);
  
  *ret_rmlist = rmlist;
  if (pval) free(pval);

  return eslOK;

 ERROR:
  if (pval) free(pval);
  return status;
}

// Aggregate p-values with equal weights. Equivalent to the Lancaster method with all p-values weighted at 2.
//
// chisq distribuitio CDF:   CDF(nf,x) = 1/2^{nf/2} * 1/Gamma(nf) * x^{nf/2) * e^{-x/2}
//
// chival = \sim_{i=1}^n -2*log(p_i)
//
// pval_agg  = 1-CDF(nf=n, x=chival)
//
int
agg_FISHER(int n, double *pval, double *ret_pval_agg, char *errbuf, int verbose)
{
  double pval_agg = -1;
  double chival = 0.;
  int    df = 2*n;      // number of degrees of freedom
  int    i;
  int    status;
 
  if (n == 1) {
    *ret_pval_agg = pval[0];
    return eslOK;
  }

  for (i = 0; i < n; i ++)
    chival += -2.0 * log(pval[i]);
  if (isnan(chival)) ESL_XFAIL(eslFAIL, errbuf, "agg_FISHER(): chival is nan");
  
  if (verbose) {
    printf("Fisher chival %f\n", chival);
    for (i = 0; i < n; i ++)
      printf("%f ",  -2.0 * log(pval[i]));
    printf("\n");
  }

  status = esl_stats_ChiSquaredTest(df, chival, &pval_agg);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "agg_FISHER(): esl_stats_ChiSquaredTest failed");
 
  *ret_pval_agg = pval_agg;
  
  return eslOK;

 ERROR:
  return status;
}

// Weighted p-value aggregation.
// doi:10.1111/j.1467-842X.1961.tb00058.x
//
// chisq distribuitio CDF:   CDF(nf,x) = 1/2^{nf/2} * 1/Gamma(nf) * x^{nf/2) * e^{-x/2}
//
// chival = \sim_{i=1}^n CDF^{-1} (nf=w_i, 1-p_i)
//
// pval_agg = 1-CDF(nf=sum_i w_i, x=chival)
//
//
//  if w_i = 2, then pval_agg(lancaster, w_i=2) = pval_agg(fisher)
//
int
agg_LANCASTER(int n, double *pval, double *weights, double *ret_pval_agg, char *errbuf, int verbose)
{
  double *invchi = NULL;
  double  pval_agg = -1;
  double  chival;
  double  df;
  int     i;
  int     status;
  
  if (n == 1) {
    *ret_pval_agg = pval[0];
    return eslOK;
  }

  df = esl_vec_DSum(weights, n);

  ESL_ALLOC(invchi, sizeof(double) * n);
  for (i = 0; i < n; i ++) {
    invchi[i] = esl_gam_invcdf(1.-pval[i], 0.0, 1./2., weights[i]/2.);
  }
  chival = esl_vec_DSum(invchi, n);
  if (isnan(chival)) ESL_XFAIL(eslFAIL, errbuf, "agg_LANCASTER(): chival is nan");

  if (verbose) {
    printf("Lancaster chival %f\n", chival);
    esl_vec_DDump(stdout, invchi, n, NULL);
  }
  
  status = esl_stats_ChiSquaredTest(df, chival, &pval_agg);
  if (status != eslOK) ESL_XFAIL(eslFAIL, errbuf, "agg_FISHER(): esl_stats_ChiSquaredTest failed");

  if (verbose) {
    for (i = 0; i < n; i ++) 
      printf("pval %g weight %g\n", pval[i], weights[i]);
    printf("agg_lancaster %g\n", pval_agg);
  }
  
  *ret_pval_agg = pval_agg;
  
  free(invchi);
  return eslOK;

 ERROR:
  if (invchi) free(invchi);
  return status;
}

// The Sidak method uses the minimum p-value but corrects it for the number of p-values that are aggregated.
// https://www.tandfonline.com/doi/abs/10.1080/01621459.1967.10482935
//
// pval_agg = 1 - (1-min_p) ^ m
//
int
agg_SIDAK(int n, double *pval, double *ret_pval_agg, char *errbuf, int verbose)
{
  double pval_agg = -1;
  double pmin;
  int    status;

  if (n == 1) {
    *ret_pval_agg = pval[0];
    return eslOK;
  }

  pmin = esl_vec_DMin(pval, n);
  if (isnan(pmin)) ESL_XFAIL(eslFAIL, errbuf, "agg_SIDAK(): pmin is nan");
  
  pval_agg = 1. - exp((double)n * log(1. - pmin));

  if (verbose) 
    printf("Sidak pval_min %f pval_agg %f\n", pmin, pval_agg);
  
  *ret_pval_agg = pval_agg;
  
  return eslOK;

 ERROR:
  return status;
}


static double *
extract_rm_pvals(RM *rm, SPAIR *spair, char *errbuf, int verbose)
{
  double *pval = NULL;
  CTLIST *ctlist = rm->ctlist;
  int    *ct;
  int     dim = ctlist->L * (ctlist->L-1) / 2;
  int     b = 0;
  int     c;
  int     idx;
  int     i, j;
  int     status;

  ESL_ALLOC(pval, sizeof(double) * rm->nbp);
  
  for (idx = 0; idx < dim; idx ++) {
    i = spair[idx].i;
    j = spair[idx].j;

    for (c = 0; c < ctlist->nct; c++) {
      ct = ctlist->ct[c];
      if (ct[i+1] == j+1 && ct[j+1] == i+1) {
	pval[b++] = spair[idx].Pval;
	break;
      }
    }
  }
  if (b != rm->nbp) ESL_XFAIL(eslFAIL, errbuf, "error in extract_rm_pvals()");

  if (verbose) esl_vec_DDump(stdout, pval, rm->nbp, NULL);

  return pval;

 ERROR:
  return NULL;
}

static double *
extract_rm_weights(RM *rm, SPAIR *spair, char *errbuf, int verbose)
{
  double *weights = NULL;
  CTLIST *ctlist = rm->ctlist;
  int    *ct;
  int     dim = ctlist->L * (ctlist->L-1) / 2;
  int     b = 0;
  int     c;
  int     idx;
  int     i, j;
  int     status;

  ESL_ALLOC(weights, sizeof(double) * rm->nbp);
  
  for (idx = 0; idx < dim; idx ++) {
    i = spair[idx].i;
    j = spair[idx].j;

    for (c = 0; c < ctlist->nct; c++) {
      ct = ctlist->ct[c];
      if (ct[i+1] == j+1 && ct[j+1] == i+1) {
	//weights[b++] = 2.0; // these weights reproduce the fisher method
	weights[b++] = 2.0 + (double)spair[idx].power;
	break;
      }
    }
  }
  if (b != rm->nbp) ESL_XFAIL(eslFAIL, errbuf, "error in extract_rm_weights()");

  if (verbose) esl_vec_DDump(stdout, weights, rm->nbp, NULL);

  return weights;

 ERROR:
  return NULL;
}


