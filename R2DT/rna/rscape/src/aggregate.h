/* aggregated p-values
 *
 */
#ifndef AGGREGATE_INCLUDED
#define AGGREGATE_INCLUDED


#include <stdio.h>   

#include "easel.h"
#include "correlators.h"
#include "structure.h"

extern int agg_CalculatePvalues(SPAIR *spair, CTLIST *ctlist, RMLIST **ret_rmlist, int helix_unpaired, enum agg_e agg_method, char *errbuf, int verbose);
extern int agg_FISHER   (int n, double *pval,                  double *ret_pval_agg, char *errbuf, int verbose);
extern int agg_LANCASTER(int n, double *pval, double *weights, double *ret_pval_agg, char *errbuf, int verbose);
extern int agg_SIDAK    (int n, double *pval,                  double *ret_pval_agg, char *errbuf, int verbose);

#endif /* AGGREGATE_INCLUDED */

/************************************************************
 * @LICENSE@
 ************************************************************/
