/* $Id$ */

/*
 * $Log$
 * Revision 1.2  2005/05/20 10:40:15  jtk
 * differentiated entropy recording for lsys, expression and diffusion phase
 *
 * Revision 1.1  2005/05/13 15:39:21  jtk
 * added entropy.c
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "trconfig.h"
#include "transsys.h"


/*
 * check whether instances in a collection are all instance of the same
 * transsys program. Return 1 if this is the case, 0 otherwise.
 * collection is a NULL-terminated array of pointers to transsys instances
 */

static int collection_is_homogeneous(TRANSSYS_INSTANCE **ti)
{
  const TRANSSYS *transsys;
  size_t i;

  if ((ti == NULL) || (ti[0] == NULL))
  {
    return (1);
  }
  transsys = ti[0]->transsys;
  for (i = 1; ti[i]; i++)
  {
    if (ti[i]->transsys != transsys)
    {
      return (0);
    }
  }
  return (1);
}


/*
 * compute entropy of the distribution of each factor in a collection
 * of transsys instances. Notice all instances must be instances of the
 * same transsys.
 */

double *transsys_collection_factor_entropy(TRANSSYS_INSTANCE **ti)
{
  size_t i, f;
  double p, sum, *entropy;
  double log2 = 1.0 / log(2.0);
  const TRANSSYS *transsys = ti[0]->transsys;

  if (!collection_is_homogeneous(ti))
  {
    fprintf(stderr, "transsys_collection_factor_entropy: inhomogeneous collection\n");
    return NULL;
  }
  entropy = (double *) malloc(transsys->num_factors * sizeof(double));
  if (entropy == NULL)
  {
    fprintf(stderr, "transsys_collection_factor_entropy: malloc failed\n");
    return NULL;
  }
  for (i = 0; i < transsys->num_factors; i++)
  {
    entropy[i] = 0.0;
  }
  for (f = 0; f < transsys->num_factors; f++)
  {
    sum = 0.0;
    for (i = 0; ti[i]; i++)
    {
      sum += ti[i]->factor_concentration[f];
    }
    for (i = 0; ti[i]; i++)
    {
      p = ti[i]->factor_concentration[f] / sum;
      if (p > 0.0)
      {
	entropy[f] -= p * log(p) * log2;
      }
    }
  }
  return (entropy);
}


double *transsys_collection_factor_information(TRANSSYS_INSTANCE **ti)
{
  double *info = transsys_collection_factor_entropy(ti);
  double max_entropy;
  size_t f, n, num_factors;

  if (info == NULL)
  {
    fprintf(stderr, "transsys_collection_factor_information: transsys_collectio_factor_entropy failed\n");
    return (NULL);
  }
  num_factors = ti[0]->transsys->num_factors;
  for (n = 0; ti[n]; n++)
    ;
  max_entropy = log((double) n) / log(2.0);
  for (f = 0; f < num_factors; f++)
  {
    info[f] = max_entropy - info[f];
  }
  return (info);
}
