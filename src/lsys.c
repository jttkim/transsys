/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.4  2005/03/31 10:14:05  jtk
 * partial implementation of diffusion
 *
 * Revision 1.3  2005/03/30 18:30:27  jtk
 * progressed transition to arrayred lsys strings
 * introduced lsys string distance matrices
 *
 * Revision 1.2  2005/03/29 17:33:02  jtk
 * introduced arrayed lsys string, with symbol distance matrix.
 *
 * Revision 1.1.1.1  2005/03/08 17:12:02  jtk
 * new cvs after loss at INB
 *
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "trconfig.h"
#include "transsys.h"


static void free_transsys_instance_list(const TRANSSYS_INSTANCE **ti_list)
{
  free(ti_list);
}


static const TRANSSYS_INSTANCE **transsys_instance_list(const SYMBOL_INSTANCE *symbol_list, const RULE_ELEMENT *rule)
{
  int i;
  const TRANSSYS_INSTANCE **ti_list;
  const SYMBOL_INSTANCE *si;

  ti_list = (const TRANSSYS_INSTANCE **) malloc(rule->lhs->num_symbols * sizeof(TRANSSYS_INSTANCE *));
  if (ti_list == NULL)
  {
    fprintf(stderr, "transsys_instance_list: malloc() failed\n");
    return (NULL);
  }
  si = symbol_list;
  for (i = 0; i < rule->lhs->num_symbols; i++)
  {
    if ((si == NULL) || (si->symbol_index != rule->lhs->symbol_list[i].symbol_index))
    {
      fprintf(stderr, "transsys_instance_list: rule mismatch in rule \"%s\"-- expected symbol #%d, found #%d\n", rule->name, rule->lhs->symbol_list[i].symbol_index, si->symbol_index);
      free_transsys_instance_list(ti_list);
      return (NULL);
    }      
    ti_list[i] = &(si->transsys_instance);
    si = si->next;
  }
  return (ti_list);
}


/* now, rewrite evaluate_assignment and evaluate_production routines to process ti_lists rather than single ti contexts */

static int evaluate_assignment(TRANSSYS_INSTANCE *target, const ASSIGNMENT *a, const TRANSSYS_INSTANCE **ti_list)
{
  double x;

  x = evaluate_expression(a->value, ti_list);
  target->factor_concentration[a->factor_index] = x;
  return (0);
}


size_t symbol_strlen(const SYMBOL_INSTANCE *symbol_string)
{
  size_t l = 0;
  const SYMBOL_INSTANCE *symbol;

  for (symbol = symbol_string; symbol; symbol = symbol->next)
    l++;
  return (l);
}


size_t lsys_string_length(const LSYS_STRING *lstr)
{
  if (lstr->arrayed)
  {
    return (lstr->num_symbols);
  }
  else
  {
    return (symbol_strlen(lstr->symbol));
  }
}


static SYMBOL_INSTANCE *evaluate_production(const LSYS_STRING *lstr, const PRODUCTION_ELEMENT *pe, const TRANSSYS_INSTANCE **ti_list, int num_predecessors, int predecessor_index)
{
  SYMBOL_INSTANCE *si;
  ASSIGNMENT *a;
  int return_value;
  int i;

  si = new_symbol_instance(lstr, pe->symbol_index, num_predecessors, predecessor_index, 1);
  if (si == NULL)
    return (NULL);
  if (pe->template_lhs_symbol_index != NO_INDEX)
  {
    if (ti_list[pe->template_lhs_symbol_index]->transsys != si->transsys_instance.transsys)
      fprintf(stderr, "evaluate_production: template transsys \"%s\" and target transsys \"%s\" mismatch\n", ti_list[pe->template_lhs_symbol_index]->transsys->name, si->transsys_instance.transsys->name);
    else
    {
      for (i = 0; i < si->transsys_instance.transsys->num_factors; i++)
	si->transsys_instance.factor_concentration[i] = ti_list[pe->template_lhs_symbol_index]->factor_concentration[i];
    }
  }
  for (a = pe->assignment_list; a; a = a->next)
  {
    return_value = evaluate_assignment(&(si->transsys_instance), a, ti_list);
    if (return_value != 0)
      fprintf(stderr, "evaluate_production: evaluate_assignment() returned %d\n", return_value);
  }
  return (si);
}


static SYMBOL_INSTANCE *evaluate_production_list(const LSYS_STRING *lstr, const PRODUCTION_ELEMENT *plist, const TRANSSYS_INSTANCE **ti_list, int num_predecessors, int predecessor_index)
{
  SYMBOL_INSTANCE *slist, *stail, *si;
  PRODUCTION_ELEMENT *p;

  if (plist == NULL)
    return (NULL);
  slist = evaluate_production(lstr, plist, ti_list, num_predecessors, predecessor_index);
  if (slist == NULL)
    return (NULL);
  stail = slist;
  for (p = plist->next; p; p = p->next)
  {
    si = evaluate_production(lstr, p, ti_list, num_predecessors, predecessor_index);
    if (si == NULL)
    {
      free_symbol_instance_list(slist);
      return (NULL);
    }
    stail->next = si;
    stail = si;
  }
  return (slist);
}


LSYS_STRING *axiom_string(const LSYS *lsys)
{
  LSYS_STRING *axiom;

  axiom = new_lsys_string(lsys);
  axiom->symbol = evaluate_production_list(axiom, lsys->axiom->production_list, NULL, 0, 0);
  arrange_lsys_string_arrays(axiom);
  return(axiom);
}


/*
 * preliminary: use single (i.e. first) symbol for obtaining transsys...
 */

static int rule_match(const SYMBOL_INSTANCE *symbol, const RULE_ELEMENT *rule)
{
  const SYMBOL_INSTANCE *s;
  const TRANSSYS_INSTANCE **ti_list;
  double return_value;
  int i;

  s = symbol;
  for (i = 0; i < rule->lhs->num_symbols; i++)
  {
    if ((s == NULL) || (s->symbol_index != rule->lhs->symbol_list[i].symbol_index))
      return (0);
    s = s->next;
  }
  if (rule->condition == NULL)
    return (1);
  /* fprintf(stderr, "rule_match: rule \"%s\" has %d lhs symbols\n", rule->name, rule->lhs->num_symbols); */
  ti_list = transsys_instance_list(symbol, rule);
  if (ti_list == NULL)
  {
    fprintf(stderr, "rule_match: transsys_instance_list() failed\n");
    return (-1);
  }
  return_value = evaluate_expression(rule->condition, ti_list);
  free_transsys_instance_list(ti_list);
  /* fprintf(stderr, "rule_match: rule \"%s\": condition returns %f, symbol sequence\n", rule->name, return_value); */
/*
  s = symbol;
  for (i = 0; i < rule->lhs->num_symbols; i++)
  {
    fprintf(stderr, "    ");
    fprint_symbol_instance(stderr, s);
    fprintf(stderr, "\n");
    s = s->next;
  }
*/
  if (return_value)
    return (1);
  else
    return (0);
}


int lsys_string_expression(LSYS_STRING *lstr)
{
  SYMBOL_INSTANCE *symbol;
  int return_value = 0;

  for (symbol = lstr->symbol; symbol; symbol = symbol->next)
  {
    return_value = process_expression(&(symbol->transsys_instance));
    if (return_value)
      break;
  }
  return (return_value);
}


/*
 * for all symbol instances, initialise all elements in new_concentration
 * to 0.0.
 * The new_concentration arrays are (mis)used to accumulate the delta
 * concentrations during simulation of diffusion
 */

static void diffusion_init_new_concentration(LSYS_STRING *lstr)
{
  int i, j;
  SYMBOL_INSTANCE *si;

  for (i = 0; i < lstr->num_symbols; i++)
  {
    si = lstr->symbol + i;
    if (si->transsys_instance.transsys)
    {
      for (j = 0; j < si->transsys_instance.transsys->num_factors; j++)
      {
	si->transsys_instance.new_concentration[j] = 0.0;
      }
    }
  }
}


int lsys_string_diffusion(LSYS_STRING *lstr)
{
  int n, i, f;
  double *dx, *w;
  const TRANSSYS *transsys;

  if (!lstr->arrayed)
  {
    fprintf(stderr, "lsys_string_diffusion: cannot process non-arrayed lsys string\n");
    return (-1);
  }
  dx = (double *) malloc(lstr->num_symbols * sizeof(double));
  if (dx == NULL)
  {
    fprintf(stderr, "lsys_string_diffusion: malloc failed\n");
    return (-1);
  }
  w = (double *) malloc(lstr->num_symbols * sizeof(double));
  if (w == NULL)
  {
    free(dx);
    fprintf(stderr, "lsys_string_diffusion: malloc failed\n");
    return (-1);
  }
  diffusion_init_new_concentration(lstr);
  for (n = 0; n < lstr->num_symbols; n++)
  {
    if (lstr->symbol[n].transsys_instance.transsys)
    {
      for (i = 0; i < lstr->num_symbols; i++)
      {
	dx[i] = 0.0;
	w[i] = 0.0;
	if (lstr->symbol[n].transsys_instance.transsys == lstr->symbol[i].transsys_instance.transsys)
	{
	  transsys = lstr->symbol[n].transsys_instance.transsys;
	  for (f = 0; f < transsys->num_factors; f++)
	  {
	  }
	}
      }
    }
  }
  free(dx);
  free(w);
  return (0);
}


int compute_distance_matrix(LSYS_STRING *lstr, const LSYS_STRING *predecessor)
{
  int return_value;
  int i, j, i1, j1, dmin;
  const SYMBOL_INSTANCE *si, *sj;

  if (!lstr->arrayed)
  {
    fprintf(stderr, "compute_distance_matrix: cannot process non-arrayed lsys string\n");
    return (-1);
  }
  if (!predecessor->arrayed)
  {
    fprintf(stderr, "compute_distance_matrix: predecessor not arrayed\n");
    return (-1);
  }
  return_value = alloc_lsys_string_distance(lstr);
  if (return_value != 0)
  {
    return (return_value);
  }
  for (i = 0; i < lstr->num_symbols; i++)
  {
    si = lstr->symbol + i;
    for (j = 0; j < lstr->num_symbols; j++)
    {
      sj = lstr->symbol + j;
      if (i == j)
      {
	lstr->distance[i][j] = 0;
      }
      else if (si->predecessor_index == sj->predecessor_index)
      {
	lstr->distance[i][j] = 1;
      }
      else
      {
	dmin = predecessor->distance[si->predecessor_index][sj->predecessor_index];
	for (i1 = si->predecessor_index; i1 < si->predecessor_index + si->num_predecessors; i1++)
	{
	  for (j1 = sj->predecessor_index; j1 < sj->predecessor_index + sj->num_predecessors; j1++)
	  {
	    if (predecessor->distance[i1][j1] < dmin)
	    {
	      dmin = predecessor->distance[i1][j1];
	    }
	  }
	}
	lstr->distance[i][j] = si->predecessor_distance + sj->predecessor_distance + dmin;
      }
      /* fprintf(stderr, "distance[%d][%d]: %d\n", i, j, lstr->distance[i][j]); */
    }
  }
  return (0);
}


LSYS_STRING *derived_string(const LSYS_STRING *lstr)
{
  LSYS_STRING *dstr;
  SYMBOL_INSTANCE *target_tail = NULL, *target_symbol;
  const RULE_ELEMENT *rule;
  const TRANSSYS_INSTANCE **ti_list;
  int symbol_index;

  if (!lstr->arrayed)
  {
    fprintf(stderr, "derived_string: cannot derive from non-arrayed string\n");
    return (NULL);
  }
  dstr = new_lsys_string(lstr->lsys);
  if (dstr == NULL)
  {
    fprintf(stderr, "derived_string: could not allocate new symbol string\n");
    return (NULL);
  }
  for (symbol_index = 0; symbol_index < lstr->num_symbols; symbol_index++)
  {
    for (rule = lstr->lsys->rule_list; rule; rule = rule->next)
    {
      if (rule_match(lstr->symbol + symbol_index, rule))
      {
	/* fprintf(stderr, "derived_string: applying rule \"%s\"\n", rule->name); */
	ti_list = transsys_instance_list(lstr->symbol + symbol_index, rule);
	target_symbol = evaluate_production_list(dstr, rule->rhs->production_list, ti_list, rule->lhs->num_symbols, symbol_index);
	free_transsys_instance_list(ti_list);
	break;
      }
    }
    if (rule == NULL)
    {
      target_symbol = clone_symbol_instance(lstr->symbol + symbol_index, dstr, symbol_index);
    }
    else
    {
      symbol_index += rule->lhs->num_symbols - 1;
    }
    if (target_symbol)
    {
      if (target_tail == NULL)
	dstr->symbol = target_symbol;
      else
	target_tail->next = target_symbol;
      for (target_tail = target_symbol; target_tail->next; target_tail = target_tail->next)
	;
    }
    else
    {
      if (rule->rhs->num_production_elements > 0)
	fprintf(stderr, "derived_string: failed to produce target symbol\n");
    }
  }
  arrange_lsys_string_arrays(dstr);
  compute_distance_matrix(dstr, lstr);
  /* fprint_lsys_string(stderr, dstr, "***\n"); */
  return (dstr);
}

