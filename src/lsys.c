/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.8  2005/04/13 13:29:30  jtk
 * compute contact_table for diffusion explicitly and only once (optimisation)
 *
 * Revision 1.7  2005/04/05 10:12:39  jtk
 * made diffusion consistent (no oscillation due to overshooting), small fixes
 *
 * Revision 1.6  2005/04/04 21:30:45  jtk
 * fixed diffusion bug (index error)
 *
 * Revision 1.5  2005/03/31 16:07:36  jtk
 * finished (initial) implementation of lsys diffusion
 *
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


int tolerate_negative_concentrations = 0;


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


/*
 * Note on the distance matrix of an axiom string:
 * All symbols in the axiom have distance 1. This reflects
 * the fact that they represent a "germ" structure in which
 * no complex spatial structure has yet differentiated.
 */

LSYS_STRING *axiom_string(const LSYS *lsys)
{
  LSYS_STRING *axiom;
  int i, j;

  axiom = new_lsys_string(lsys);
  axiom->symbol = evaluate_production_list(axiom, lsys->axiom->production_list, NULL, 0, 0);
  arrange_lsys_string_arrays(axiom);
  if(alloc_lsys_string_distance(axiom) != 0)
  {
    free_lsys_string(axiom);
    return (NULL);
  }
  for (i = 0; i < axiom->num_symbols; i++)
  {
    for (j = 0; j < axiom->num_symbols; j++)
    {
      axiom->distance[i][j] = i != j;
    }
  }
  return(axiom);
}


/*
 * FIXME??? (is this already fixed?)
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


static int within_diffusion_range(const LSYS_STRING *lstr, int i, int j)
{
  if (lstr->symbol[i].transsys_instance.transsys != lstr->symbol[j].transsys_instance.transsys)
  {
    return (0);
  }
  return (lstr->distance[i][j] <= lstr->lsys->diffusion_range);
}


static void free_contact_table(const LSYS_STRING *lstr, SYMBOL_INSTANCE ***contact_table)
{
  int i;

  for (i = 0; i < lstr->num_symbols; i++)
  {
    free(contact_table[i]);
  }
  free(contact_table);
}


static SYMBOL_INSTANCE ***symbol_contact_table(const LSYS_STRING *lstr)
{
  SYMBOL_INSTANCE ***contact_table;
  int row_length, i, j;

  contact_table = (SYMBOL_INSTANCE ***) malloc(lstr->num_symbols * sizeof(SYMBOL_INSTANCE **));
  if (contact_table == NULL)
  {
    fprintf(stderr, "symbol_contact_table: malloc failed (contact_table)\n");
    return (NULL);
  }
  for (i = 0; i < lstr->num_symbols; i++)
  {
    contact_table[i] = (SYMBOL_INSTANCE **) malloc((lstr->num_symbols + 1) * sizeof(SYMBOL_INSTANCE *));
    if (contact_table[i] == NULL)
    {
      fprintf(stderr, "symbol_contact_table: malloc failed (row)\n");
      for (j = 0; j < i; j++)
      {
	free(contact_table[j]);
      }
      free(contact_table);
      return (NULL);
    }
    row_length = 0;
    for (j = 0; j < lstr->num_symbols; j++)
    {
      if (within_diffusion_range(lstr, i, j))
      {
	contact_table[i][row_length++] = lstr->symbol + j;
      }
    }
    contact_table[i][row_length++] = NULL;
    realloc(contact_table[i], row_length * sizeof(SYMBOL_INSTANCE *));
  }
  return (contact_table);
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


static void diffusion_weights(const LSYS_STRING *lstr, int i, double *w)
{
  double wsum = 0.0;
  int j;

  for (j = 0; j < lstr->num_symbols; j++)
  {
    if (within_diffusion_range(lstr, i, j))
    {
      w[j] = 1.0;
      wsum += 1.0;
    }
    else
    {
      w[j] = 0.0;
    }
  }
  for (j = 0; j < lstr->num_symbols; j++)
  {
    w[j] /= wsum;
  }
}


/*
 * FIXME: the diffusion algorithm has complexity square of
 * num_symbols in string, could be optimised by initially computing
 * contact graph (num_symbols^2 would remain to be worst case, but
 * one can reasonably assume that contact graph is normally rather
 * sparse)
 */


/*
 * diffusion concept: diffusibility = 1 means that total equilibration
 * *within local neighbourhood* takes place in one time step.
 */

int lsys_string_diffusion(LSYS_STRING *lstr)
{
  int i, j, f;
  double w, dc, dcsum, diffusibility;
  const TRANSSYS *transsys;
  const TRANSSYS_INSTANCE *ti;
  SYMBOL_INSTANCE ***contact_table;

  if (!lstr->arrayed)
  {
    fprintf(stderr, "lsys_string_diffusion: cannot process non-arrayed lsys string\n");
    return (-1);
  }
  contact_table = symbol_contact_table(lstr);
  if (contact_table == NULL)
  {
    fprintf(stderr, "lsys_string_diffusion: symbol_contact_table failed\n");
    return (-1);
  }
  diffusion_init_new_concentration(lstr);
  for (i = 0; i < lstr->num_symbols; i++)
  {
    ti = &(lstr->symbol[i].transsys_instance);
    transsys = ti->transsys;
    if (ti->transsys)
    {
      for (j = 0; contact_table[i][j]; j++)
	;
      w = 1.0 / ((double) j);
      for (f = 0; f < transsys->num_factors; f++)
      {
	diffusibility = evaluate_expression(ti->transsys->factor_list[f].diffusibility_expression, &ti);
	dcsum = 0.0;
	for (j = 0; contact_table[i][j]; j++)
	{
	  dc = (ti->factor_concentration[f] - contact_table[i][j]->transsys_instance.factor_concentration[f]) * diffusibility * w;
	  dcsum += dc;
	  contact_table[i][j]->transsys_instance.new_concentration[f] += dc;
	}
	ti->new_concentration[f] -= dcsum;
      }
    }
  }
  free_contact_table(lstr, contact_table);
  for (i = 0; i < lstr->num_symbols; i++)
  {
    if (lstr->symbol[i].transsys_instance.transsys)
    {
      for (f = 0; f < lstr->symbol[i].transsys_instance.transsys->num_factors; f++)
      {
	lstr->symbol[i].transsys_instance.factor_concentration[f] += lstr->symbol[i].transsys_instance.new_concentration[f];
      }
    }
  }
  return (0);
}


/*
 * IMPLEMENTME: a function to compute the total amount (sum of concentrations)
 * in an lsys string. To be used to check correct implementation of diffusion
 * (conservation)
 */


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

