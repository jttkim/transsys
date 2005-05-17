/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.13  2005/05/17 12:11:30  jtk
 * contact graph works
 *
 * Revision 1.12  2005/05/16 21:03:27  jtk
 * contact graph implementation still buggy
 *
 * Revision 1.11  2005/05/16 12:02:10  jtk
 * in transition from distance matrices to contact graphs
 *
 * Revision 1.10  2005/05/13 20:08:21  jtk
 * implementation information-decreasing diffusion underway
 *
 * Revision 1.9  2005/04/14 18:48:59  jtk
 * diffusion fix: amount transferred must be independent of gradient
 *
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
#include <math.h>

#include "trconfig.h"
#include "transsys.h"


int tolerate_negative_concentrations = 0;


static void free_transsys_instance_list(const TRANSSYS_INSTANCE **ti_list)
{
  free(ti_list);
}


/*
 * Assemble an array of pointers to the transsys instances associated with
 * a substring of symbol instances, starting at symbol_list, and with a length
 * determined by the lhs of the rule.
 */

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


static SYMBOL_INSTANCE *evaluate_production(const LSYS_STRING *lstr, const PRODUCTION_ELEMENT *pe, const TRANSSYS_INSTANCE **ti_list)
{
  SYMBOL_INSTANCE *si;
  ASSIGNMENT *a;
  int return_value;
  int i;

  si = new_symbol_instance(lstr, pe->symbol_index);
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


static SYMBOL_INSTANCE *evaluate_production_list(const LSYS_STRING *lstr, const PRODUCTION_ELEMENT *plist, const TRANSSYS_INSTANCE **ti_list)
{
  SYMBOL_INSTANCE *slist, *stail, *si;
  PRODUCTION_ELEMENT *p;

  if (plist == NULL)
    return (NULL);
  slist = evaluate_production(lstr, plist, ti_list);
  if (slist == NULL)
    return (NULL);
  stail = slist;
  for (p = plist->next; p; p = p->next)
  {
    si = evaluate_production(lstr, p, ti_list);
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


static int connect_symbol_group(LSYS_STRING *lstr, int group_start, int group_length)
{
  int group_end = group_start + group_length, i, j;

  for (i = group_start; i < group_end; i++)
  {
    for (j = group_start; j < i; j++)
    {
      if (connect_lsys_string_symbols(lstr, i, j, 1) != 0)
      {
	fprintf(stderr, "connect_symbol_group: connect_lsys_string_symbols failed\n");
	return (-1);
      }
    }
  }
  return (0);
}


/*
 * Note on the contact graph of an axiom string:
 * All symbols in the axiom have distance 1. This reflects
 * the fact that they represent a "germ" structure in which
 * no complex spatial structure has yet differentiated.
 */

LSYS_STRING *axiom_string(const LSYS *lsys)
{
  LSYS_STRING *axiom;

  axiom = new_lsys_string(lsys);
  axiom->symbol = evaluate_production_list(axiom, lsys->axiom->production_list, NULL);
  arrange_lsys_string_arrays(axiom);
  if (alloc_lsys_string_contact_graph_components(&(axiom->contact_graph), axiom->num_symbols * (axiom->num_symbols - 1) / 2))
  {
    fprintf(stderr, "axiom_string: alloc_lsys_string_contact_graph_components failed\n");
    free_lsys_string(axiom);
    return (NULL);
  }
  if (connect_symbol_group(axiom, 0, axiom->num_symbols) != 0)
  {
    fprintf(stderr, "axiom_string: connect_symbol_group failed\n");
    free_lsys_string(axiom);
    return (NULL);
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
    {
      break;
    }
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


static void diffusion_init_contact_edges(LSYS_STRING *lstr)
{
  size_t i, j;

  for (i = 0; i < lstr->num_symbols; i++)
  {
    for (j = 0; j < lstr->symbol[i].num_contact_edges; j++)
    {
      lstr->symbol[i].contact_edge[j]->amount_diffused = 0.0;
      lstr->symbol[i].contact_edge[j]->amount_valid = 0;
    }
  }
}


static int diffusion_set_transferred_amount(LSYS_STRING_CONTACT_EDGE *edge, int symbol_index, double amount)
{
  if (edge->i2 == symbol_index)
  {
    amount = -amount;
  }
  else if (edge->i1 != symbol_index)
  {
    fprintf(stderr, "diffusion_set_transferred_amount: illegal symbol index\n");
    return (-1);
  }
  if (edge->amount_valid)
  {
    if (fabs(edge->amount_diffused) > fabs(amount))
    {
      edge->amount_diffused = amount;
    }
  }
  else
  {
    edge->amount_diffused = amount;
  }
  edge->amount_valid = 1;
  return (0);
}


static void free_lsys_string_transsys_list(const TRANSSYS **tlist)
{
  free(tlist);
}


static const TRANSSYS **lsys_string_transsys_list(const LSYS_STRING *lstr)
{
  size_t i, j, num_transsys = 0;
  const TRANSSYS *t;
  const TRANSSYS **tlist = (const TRANSSYS **) malloc(sizeof(const TRANSSYS *)), **tl1;

  if (tlist == NULL)
  {
    fprintf(stderr, "lsys_string_transsys_list: malloc failed\n");
    return (NULL);
  }
  for (i = 0; i < lstr->num_symbols; i++)
  {
    t = lstr->symbol[i].transsys_instance.transsys;
    if (t)
    {
      for (j = 0; j < num_transsys; j++)
      {
	if (t == tlist[j])
	{
	  break;
	}
      }
      if (j == num_transsys)
      {
	tl1 = (const TRANSSYS **) realloc(tlist, (num_transsys + 1) * sizeof(const TRANSSYS *));
	if (tl1 == NULL)
	{
	  fprintf(stderr, "lsys_string_transsys_list: realloc failed\n");
	  free(tlist);
	  return (NULL);
	}
	tlist = tl1;
	tlist[num_transsys++] = t;
      }
    }
  }
  tlist[num_transsys] = NULL;
  return (tlist);
}


int other_symbol_instance_index(const LSYS_STRING_CONTACT_EDGE *edge, int si_index)
{
  if (si_index == edge->i1)
  {
    return (edge->i2);
  }
  if (si_index == edge->i2)
  {
    return (edge->i1);
  }
  fprintf(stderr, "other_symbol_instance_index: illegal symbol index %d: not in {%d, %d}\n", si_index, edge->i1, edge->i2);
  return (-1);
}


typedef struct
{
  int num_neighbours;
  int *neighbour_index;
  int *edge_index;
  double *gradient;
  double mean_concentration;
  double gradient_sum;
} NEIGHBOURHOOD;


static void free_neighbourhood(NEIGHBOURHOOD *neighbourhood)
{
  if (neighbourhood->num_neighbours > 0)
  {
    free(neighbourhood->neighbour_index);
    free(neighbourhood->edge_index);
    free(neighbourhood->gradient);
  }
  free(neighbourhood);
}


static NEIGHBOURHOOD *get_neighbourhood(const LSYS_STRING *lstr, int symbol_index, int factor_index)
{
  NEIGHBOURHOOD *neighbourhood = (NEIGHBOURHOOD *) malloc(sizeof(NEIGHBOURHOOD));
  int i, n, other_index;
  const SYMBOL_INSTANCE *si = lstr->symbol + symbol_index;
  const TRANSSYS *transsys = si->transsys_instance.transsys;

  if (neighbourhood == NULL)
  {
    fprintf(stderr, "get_neighbourhood: malloc failed\n");
    return (NULL);
  }
  neighbourhood->num_neighbours = 0;
  for (i = 0; i < si->num_contact_edges; i++)
  {
    other_index = other_symbol_instance_index(si->contact_edge[i], symbol_index);
    if (lstr->symbol[other_index].transsys_instance.transsys == transsys)
    {
      neighbourhood->num_neighbours++;
    }
  }
  if (neighbourhood->num_neighbours == 0)
  {
    neighbourhood->neighbour_index = NULL;
    neighbourhood->edge_index = NULL;
    neighbourhood->gradient = NULL;
    neighbourhood->mean_concentration = si->transsys_instance.factor_concentration[factor_index];
    return (neighbourhood);
  }
  neighbourhood->neighbour_index = (int *) malloc(neighbourhood->num_neighbours * sizeof(int));
  if (neighbourhood->neighbour_index == NULL)
  {
    fprintf(stderr, "get_neighbourhood: malloc for neighbour_index failed\n");
    free(neighbourhood);
    return (NULL);
  }
  neighbourhood->edge_index = (int *) malloc(neighbourhood->num_neighbours * sizeof(int));
  if (neighbourhood->neighbour_index == NULL)
  {
    fprintf(stderr, "get_neighbourhood: malloc for neighbour_index failed\n");
    free(neighbourhood->neighbour_index);
    free(neighbourhood);
    return (NULL);
  }
  neighbourhood->gradient = (double *) malloc(neighbourhood->num_neighbours * sizeof(double));
  if (neighbourhood->gradient == NULL)
  {
    fprintf(stderr, "get_neighbourhood: malloc for gradient failed\n");
    free(neighbourhood->neighbour_index);
    free(neighbourhood->edge_index);
    free(neighbourhood);
    return (NULL);
  }
  neighbourhood->mean_concentration = si->transsys_instance.factor_concentration[factor_index];
  neighbourhood->gradient_sum = 0.0;
  n = 0;
  for (i = 0; i < si->num_contact_edges; i++)
  {
    other_index = other_symbol_instance_index(si->contact_edge[i], symbol_index);
    if (lstr->symbol[other_index].transsys_instance.transsys == transsys)
    {
      neighbourhood->neighbour_index[n] = other_index;
      neighbourhood->edge_index[n] = i;
      neighbourhood->gradient[n] = lstr->symbol[other_index].transsys_instance.factor_concentration[factor_index] - si->transsys_instance.factor_concentration[factor_index];
      neighbourhood->gradient_sum += neighbourhood->gradient[n];
      neighbourhood->mean_concentration += lstr->symbol[other_index].transsys_instance.factor_concentration[factor_index];
      n++;
    }
  }
  neighbourhood->mean_concentration /= (neighbourhood->num_neighbours + 1);
  return (neighbourhood);
}


static int diffuse_along_contact_edges(LSYS_STRING *lstr, const TRANSSYS *transsys, int factor_index)
{
  int e;
  SYMBOL_INSTANCE *s1, *s2;

  for (e = 0; e < lstr->contact_graph.num_edges; e++)
  {
    s1 = lstr->symbol + lstr->contact_graph.edge[e].i1;
    s2 = lstr->symbol + lstr->contact_graph.edge[e].i2;
    if ((s1->transsys_instance.transsys == transsys) && (s2->transsys_instance.transsys == transsys))
    {
      s1->transsys_instance.factor_concentration[factor_index] += lstr->contact_graph.edge[e].amount_diffused;
      s2->transsys_instance.factor_concentration[factor_index] -= lstr->contact_graph.edge[e].amount_diffused;
    }
  }
  return (0);
}


/*
 * diffusion concept: diffusibility = 1 means that total equilibration
 * *within local neighbourhood* takes place in one time step (to the extent
 * this is feasible).
 */

int lsys_string_diffusion(LSYS_STRING *lstr)
{
  int i, j, f, t;
  double diffusibility, d;
  const TRANSSYS *transsys, **tlist;
  const TRANSSYS_INSTANCE *ti;
  SYMBOL_INSTANCE *si;
  NEIGHBOURHOOD *neighbourhood;

  if (!lstr->arrayed)
  {
    fprintf(stderr, "lsys_string_diffusion: cannot process non-arrayed lsys string\n");
    return (-1);
  }
  tlist = lsys_string_transsys_list(lstr);
  if (tlist == NULL)
  {
    return (-1);
  }
  diffusion_init_new_concentration(lstr);
  for (t = 0; tlist[t]; t++)
  {
    for (f = 0; f < tlist[t]->num_factors; f++)
    {
      diffusion_init_contact_edges(lstr);
      for (i = 0; i < lstr->num_symbols; i++)
      {
	si = lstr->symbol + i;
	ti = &(si->transsys_instance);
	transsys = ti->transsys;
	if (transsys == tlist[t])
	{
	  neighbourhood = get_neighbourhood(lstr, i, f);
	  if (neighbourhood == NULL)
	  {
	    return (-1);
	  }
	  /* FIXME (?) should diffusibility be clipped to [0, 1] ?? Warning? */
	  diffusibility = evaluate_expression(transsys->factor_list[f].diffusibility_expression, &ti);
	  for (j = 0; j < neighbourhood->num_neighbours; j++)
	  {
	    /* FIXME: this is numerically very unstable. Some more maths may help fixing this... */
	    if (neighbourhood->gradient_sum == 0.0)
	    {
	      if (si->transsys_instance.factor_concentration[f] != neighbourhood->mean_concentration)
	      {
		double relative_error = (neighbourhood->mean_concentration - si->transsys_instance.factor_concentration[f]) / (neighbourhood->mean_concentration + si->transsys_instance.factor_concentration[f]) * 0.5;
		if (relative_error > 1e-15)
		{
		  fprintf(stderr, "lsys_string_diffusion: gradient sum 0 but local mean - local concentration = %g, ratio = %g\n", neighbourhood->mean_concentration - si->transsys_instance.factor_concentration[f], relative_error);
		}
	      }
	      d = diffusibility;
	    }
	    else
	    {
	      d = (neighbourhood->mean_concentration - si->transsys_instance.factor_concentration[f]) / neighbourhood->gradient_sum * diffusibility;
	    }
	    diffusion_set_transferred_amount(si->contact_edge[neighbourhood->edge_index[j]], i, d * neighbourhood->gradient[j]);
	  }
	  free_neighbourhood(neighbourhood);
	}
      }
      /* fprintf(stderr, "lsys_string_diffusion: diffusing #%d along edges\n", f); */
      diffuse_along_contact_edges(lstr, tlist[t], f);
    }
  }
  free_lsys_string_transsys_list(tlist);
  return (0);
}


/*
 * IMPLEMENTME: a function to compute the total amount (sum of concentrations)
 * in an lsys string. To be used to check correct implementation of diffusion
 * (conservation)
 */


LSYS_STRING_CONTACT_GRAPH *group_contact_graph(const LSYS_STRING *lstr)
{
  LSYS_STRING_CONTACT_GRAPH *gcgraph = (LSYS_STRING_CONTACT_GRAPH *) malloc(sizeof(LSYS_STRING_CONTACT_GRAPH));
  int e, ge, i1, i2, return_value, ls1, ls2, ltmp;

  if (gcgraph == NULL)
  {
    fprintf(stderr, "group_contact_graph: malloc failed\n");
    return (NULL);
  }
  init_lsys_string_contact_graph_components(gcgraph);
  return_value = alloc_lsys_string_contact_graph_components(gcgraph, lstr->contact_graph.num_edges);
  if (return_value != 0)
  {
    fprintf(stderr, "group_contact_graph: alloc_lsys_string_contact_graph_components failed\n");
    free(gcgraph);
    return (NULL);
  }
  for (e = 0; e < lstr->contact_graph.num_edges; e++)
  {
    i1 = lstr->contact_graph.edge[e].i1;
    i2 = lstr->contact_graph.edge[e].i2;
    ls1 = lstr->symbol[i1].lhs_group_start;
    ls2 = lstr->symbol[i2].lhs_group_start;
    /* always have l1 < l2 in group contact graph helps to control time complexity */
    if (ls1 > ls2)
    {
      l_tmp = ls1;
      ls1 = ls2;
      ls2 = l_tmp;
    }
    /* FIXME: this linear search is a major time burner */
    for (ge = 0; ge < gcgraph->num_edges; ge++)
    {
      if ((gcgraph->edge[ge].i1 == ls1) && (gcgraph->edge[ge].i2 == ls2))
      {
	break;
      }
    }
    if (ge == gcgraph->num_edges)
    {
      if (add_lsys_string_contact_edge(gcgraph, ls1, ls2, lstr->contact_graph.edge[e].distance) != 0)
      {
	fprintf(stderr, "group_contact_graph: add_lsys_string_contact_edge failed\n");
      }
    }
    else
    {
      if (gcgraph->edge[ge].distance > lstr->contact_graph.edge[e].distance)
      {
	gcgraph->edge[ge].distance = lstr->contact_graph.edge[e].distance;
      }
    }
  }
  return (gcgraph);
}


static int compute_contact_graph(LSYS_STRING *lstr, const LSYS_STRING *predecessor)
{
  int e, n, p_i1, p_i2, i, j, distance;
  const SYMBOL_INSTANCE *p1, *p2;
  LSYS_STRING_CONTACT_GRAPH *gcgraph;

  if (!lstr->arrayed)
  {
    fprintf(stderr, "compute_contact_graph: cannot process non-arrayed lsys string\n");
    return (-1);
  }
  if (!predecessor->arrayed)
  {
    fprintf(stderr, "compute_contact_graph: predecessor not arrayed\n");
    return (-1);
  }
  gcgraph = group_contact_graph(predecessor);
  if (gcgraph == NULL)
  {
    fprintf(stderr, "compute_contact_graph: group_contact_graph failed\n");
    return (-1);
  }
  /* FIXME: loop is repeatedly done for finding number of edges and for adding edges */
  n = 0;
  /* edges between symbols with same predecessors */
  i = 0;
  while (i < predecessor->num_symbols)
  {
    if (predecessor->symbol[i].lhs_group_length == 0)
    {
      fprintf(stderr, "compute_contact_graph: lhs group info corrupted\n");
    }
    else
    {
      if (predecessor->symbol[i].num_successors > 0)
      {
	n += predecessor->symbol[i].num_successors * (predecessor->symbol[i].num_successors - 1) / 2;
      }
      i += predecessor->symbol[i].lhs_group_length;
    }
  }
  /* edges between symbols with different predecessors */
  for (e = 0; e < gcgraph->num_edges; e++)
  {
    p_i1 = gcgraph->edge[e].i1;
    p_i2 = gcgraph->edge[e].i2;
    p1 = predecessor->symbol + p_i1;
    p2 = predecessor->symbol + p_i2;
    for (i = p1->successor_index; i < p1->successor_index + p1->num_successors; i++)
    {
      /* FIXME: this could probably be improved in speed */
      for (j = p2->successor_index; j < p2->successor_index + p2->num_successors; j++)
      {
	distance = gcgraph->edge[e].distance + p1->successor_distance + p2->successor_distance;
	if (distance <= lstr->lsys->diffusion_range)
	{
	  n++;
	}
      }
    }
  }
  if (alloc_lsys_string_contact_graph_components(&(lstr->contact_graph), n) != 0)
  {
    fprintf(stderr, "compute_contact_graph: alloc_lsys_string_contact_graph_components failed\n");
    free_lsys_string_contact_graph(gcgraph);
    return (-1);
  }
  /* edges between symbols with same predecessors */
  i = 0;
  while (i < predecessor->num_symbols)
  {
    if (predecessor->symbol[i].lhs_group_length == 0)
    {
      fprintf(stderr, "compute_contact_graph: lhs group info corrupted\n");
    }
    else
    {
      if (predecessor->symbol[i].num_successors > 0)
      {
	if (connect_symbol_group(lstr, predecessor->symbol[i].successor_index, predecessor->symbol[i].num_successors) != 0)
	{
	  /* FIXME: should bail out here */
	  fprintf(stderr, "compute_contact_graph: connect_symbol_group failed\n");
	}
      }
      i += predecessor->symbol[i].lhs_group_length;
    }
  }
  /* edges between symbols with different predecessors */
  for (e = 0; e < gcgraph->num_edges; e++)
  {
    p_i1 = gcgraph->edge[e].i1;
    p_i2 = gcgraph->edge[e].i2;
    p1 = predecessor->symbol + p_i1;
    p2 = predecessor->symbol + p_i2;
    for (i = p1->successor_index; i < p1->successor_index + p1->num_successors; i++)
    {
      /* FIXME: this could probably be improved in speed */
      for (j = p2->successor_index; j < p2->successor_index + p2->num_successors; j++)
      {
	distance = gcgraph->edge[e].distance + p1->successor_distance + p2->successor_distance;
	if (distance <= lstr->lsys->diffusion_range)
	{
	  if (connect_lsys_string_symbols(lstr, i, j, distance) != 0)
	  {
	    fprintf(stderr, "compute_contact_graph: connect_lsys_string_symbols failed\n");
	  }
	}
      }
    }
  }
  free_lsys_string_contact_graph(gcgraph);
  return (0);
}


/* incomplete: successor info not stored in symbol instances (?) */

LSYS_STRING *derived_string(LSYS_STRING *lstr)
{
  LSYS_STRING *dstr;
  SYMBOL_INSTANCE *target_tail = NULL, *target_symbol;
  const RULE_ELEMENT *rule;
  const TRANSSYS_INSTANCE **ti_list;
  int si_index, i, lhs_length;
  int num_successors, successor_index;

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
  successor_index = 0;
  si_index = 0;
  while (si_index < lstr->num_symbols)
  {
    for (rule = lstr->lsys->rule_list; rule; rule = rule->next)
    {
      if (rule_match(lstr->symbol + si_index, rule))
      {
	/* fprintf(stderr, "derived_string: applying rule \"%s\"\n", rule->name); */
	break;
      }
    }
    if (rule)
    {
      /* FIXME: should bail out if there's a problem here... */
      ti_list = transsys_instance_list(lstr->symbol + si_index, rule);
      if (ti_list == NULL)
      {
	fprintf(stderr, "derived_string: transsys_instance_list failed\n");
	free_lsys_string(dstr);
	return (NULL);
      }
      lhs_length = rule->lhs->num_symbols;
      target_symbol = evaluate_production_list(dstr, rule->rhs->production_list, ti_list);
      free_transsys_instance_list(ti_list);
      if (target_symbol == NULL)
      {
	fprintf(stderr, "derived_string: failed to produce target symbol\n");
	free_lsys_string(dstr);
	return (NULL);
      }
    }
    else
    {
      target_symbol = clone_symbol_instance(lstr->symbol + si_index, dstr);
      lhs_length = 1;
    }
    lstr->symbol[si_index].lhs_group_length = lhs_length;
    lstr->symbol[si_index].lhs_group_start = si_index;
    for (i = 1; i < lhs_length; i++)
    {
      /* an lhs group length of 0 indicates that this symbol is not the first one in an lhs group */
      lstr->symbol[si_index + i].lhs_group_length = 0;
      lstr->symbol[si_index + i].lhs_group_start = si_index;
    }
    if (target_tail == NULL)
    {
      dstr->symbol = target_symbol;
    }
    else
    {
      target_tail->next = target_symbol;
    }
    num_successors = 0;
    /* advance target_tail to new tail and count number of successors in the process */
    if (target_symbol != NULL)
    {
      num_successors = 1;
      for (target_tail = target_symbol; target_tail->next; target_tail = target_tail->next)
      {
	num_successors++;
      }
    }
    for (i = 0; i < lhs_length; i++)
    {
      lstr->symbol[si_index].num_successors = num_successors;
      lstr->symbol[si_index].successor_index = successor_index;
      if  (rule != NULL)
      {
	lstr->symbol[si_index].successor_distance = 1;
      }
      else
      {
	lstr->symbol[si_index].successor_distance = 0;
      }
      si_index++;
    }
    successor_index += num_successors;
  }
  if (arrange_lsys_string_arrays(dstr) != 0)
  {
    fprintf(stderr, "derived_string: arrange_lsys_string_arrays failed\n");
    free_lsys_string(dstr);
    return (NULL);
  }
  if (compute_contact_graph(dstr, lstr) != 0)
  {
    fprintf(stderr, "derived_string: compute_contact_graph failed\n");
    free_lsys_string(dstr);
    return (NULL);
  }
  /* fprint_lsys_string(stderr, dstr, "***\n"); */
  return (dstr);
}

