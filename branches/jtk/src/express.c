/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.1  2005/03/08 17:12:02  jtk
 * Initial revision
 *
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#include <stdio.h>
#include <stdlib.h>

#include "transsys.h"


static double minimal_concentration(int num_binding_factors, int *fi, double *fc)
{
  int i;
  double c;

  c = fc[fi[0]];
  for (i = 1; i < num_binding_factors; i++)
  {
    if (c > fc[fi[i]])
      c = fc[fi[i]];
  }
  return (c);
}


int process_expression(TRANSSYS_INSTANCE *ti)
{
  double *nc;
  int i;
  double d, a, r, cmin, km, max, f_inc;
  PROMOTER_ELEMENT *ae;
  const TRANSSYS_INSTANCE *const_ti = ti;

  if (ti->transsys == NULL)
    return (0);
#ifdef FULL_CHECKING
  if (!ti->transsys->arrayed)
  {
    fprintf(stderr, "process_expression: transsys not arrayed\n");
    return (-1);
  }
#endif
  for (i = 0; i < ti->transsys->num_factors; i++)
  {
    d = evaluate_expression(ti->transsys->factor_list[i].decay_expression, &const_ti);
    ti->new_concentration[i] = ti->factor_concentration[i] * (1.0 - d);
    /* fprintf(stderr, "before decay: [%s] = %g, after decay: [%s] = %g\n", ti->transsys->factor_list[i].name, ti->factor_concentration[i], ti->transsys->factor_list[i].name, ti->new_concentration[i]); */
  }
  for (i = 0; i < ti->transsys->num_genes; i++)
  {
    a = 0.0;
    r = 0.0;
    for (ae = ti->transsys->gene_list[i].promoter_list; ae; ae = ae->next)
    {
      switch (ae->type)
      {
      case ACT_CONSTITUTIVE:
	a += evaluate_expression(ae->expr1, &const_ti);
	break;
      case ACT_ACTIVATE:
	cmin = minimal_concentration(ae->num_binding_factors, ae->factor_index, ti->factor_concentration);
	km = evaluate_expression(ae->expr1, &const_ti);
	max = evaluate_expression(ae->expr2, &const_ti);
	a += max * cmin / (km + cmin);
	break;
      case ACT_REPRESS:
	cmin = minimal_concentration(ae->num_binding_factors, ae->factor_index, ti->factor_concentration);
	km = evaluate_expression(ae->expr1, &const_ti);
	max = evaluate_expression(ae->expr2, &const_ti);
	r += max * cmin / (km + cmin);
	break;
      default:
	fprintf(stderr, "process_expression: unknown activation type %d\n", (int) ae->type);
	break;
      }
    }
    f_inc = a - r;
    if (f_inc < 0.0)
      f_inc = 0.0;
    if ((ti->transsys->gene_list[i].product_index >= 0) && (ti->transsys->gene_list[i].product_index < ti->transsys->num_factors))
      ti->new_concentration[ti->transsys->gene_list[i].product_index] += f_inc;
    else
      fprintf(stderr, "process_expression: transsys \"%s\", gene \"%s\": product index %d out of range\n",
	      ti->transsys->name, ti->transsys->gene_list[i].name, ti->transsys->gene_list[i].product_index);
  }
  nc = ti->new_concentration;
  ti->new_concentration = ti->factor_concentration;
  ti->factor_concentration = nc;
  return (0);
}

