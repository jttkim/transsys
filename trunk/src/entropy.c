/* $Id$ */

/*
 * $Log$
 * Revision 1.1  2005/05/13 15:39:21  jtk
 * added entropy.c
 *
 */

#include <stdio.h>
#include <math.h>

#include "trconfig.h"
#include "transsys.h"


double *expression_entropy(TRANSSYS_INSTANCE **ti, size_t n)
{
  size_t i, f;
  double d, sum, *entropy;
  double log2 = 1.0 / log(2.0);
  const TRANSSYS *transsys = ti[0]->transsys;

  entropy = (double *) malloc(transsys->num_factors * sizeof(double));
  if (entropy == NULL)
  {
    fprintf(stderr, "expression_entropy: malloc failed\n");
    return NULL;
  }
  for (i = 0; i < transsys->num_factors; i++)
  {
    entropy[i] = 0.0;
  }
  for (f = 0; f < transsys->num_factors; f++)
  {
    sum = 0.0;
    for (i = 0; i < n; i++)
    {
      sum += ti[i]->factor_concentration[f];
    }
    for (i = 0; i < n; i++)
    {
      p = ti[i]->factor_concentration[f] / sum;
      entropy[f] -= p * log(p) * log2;
    }
  }
  return (entropy);
}


