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


static SYMBOL_INSTANCE *evaluate_production(const LSYS *lsys, const PRODUCTION_ELEMENT *pe, const TRANSSYS_INSTANCE **ti_list)
{
  SYMBOL_INSTANCE *si;
  ASSIGNMENT *a;
  int return_value;
  int i;

/*   fprintf(stderr, "evaluate_production: producing symbol \"%s\"\n", lsys->symbol_list[pe->symbol_index].name); */
/*   if (context) */
/*     fprintf(stderr, "evaluate_production: context transsys is \"%s\"\n", context->transsys->name); */
/*   else */
/*     fprintf(stderr, "evaluate_production: no context transsys\n"); */
  si = new_symbol_instance(lsys, pe->symbol_index);
  if (si == NULL)
    return (NULL);
  /* copy transsys instance if target transsys is source transsys -- how to extend this for multi-symbol lhs? */
/*   if (si->transsys_instance.transsys && ti_list && (si->transsys_instance.transsys == context->transsys)) */
/*   { */
/*     for (i = 0; i < context->transsys->num_factors; i++) */
/*       si->transsys_instance.factor_concentration[i] = context->factor_concentration[i]; */
/*   } */
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


static SYMBOL_INSTANCE *evaluate_production_list(const LSYS *lsys, const PRODUCTION_ELEMENT *plist, const TRANSSYS_INSTANCE **ti_list)
{
  SYMBOL_INSTANCE *slist, *stail, *si;
  PRODUCTION_ELEMENT *p;

  if (plist == NULL)
    return (NULL);
  slist = evaluate_production(lsys, plist, ti_list);
  if (slist == NULL)
    return (NULL);
  stail = slist;
  for (p = plist->next; p; p = p->next)
  {
    si = evaluate_production(lsys, p, ti_list);
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


SYMBOL_INSTANCE *axiom_string(const LSYS *lsys)
{
  return (evaluate_production_list(lsys, lsys->axiom->production_list, NULL));
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


int string_transsys_expression(SYMBOL_INSTANCE *symbol_string)
{
  SYMBOL_INSTANCE *symbol;
  int return_value;

  for (symbol = symbol_string; symbol; symbol = symbol->next)
  {
    return_value = process_expression(&(symbol->transsys_instance));
    if (return_value)
      break;
  }
  return (0);
}


SYMBOL_INSTANCE *derived_string(const LSYS *lsys, const SYMBOL_INSTANCE *source_string)
{
  const SYMBOL_INSTANCE *source_symbol;
  SYMBOL_INSTANCE *target_string = NULL, *target_tail = NULL, *target_symbol;
  const RULE_ELEMENT *rule;
  const TRANSSYS_INSTANCE **ti_list;
  int i;

  for (source_symbol = source_string; source_symbol; source_symbol = source_symbol->next)
  {
    for (rule = lsys->rule_list; rule; rule = rule->next)
    {
      if (rule_match(source_symbol, rule))
      {
	/* fprintf(stderr, "derived_string: applying rule \"%s\"\n", rule->name); */
	ti_list = transsys_instance_list(source_symbol, rule);
	target_symbol = evaluate_production_list(lsys, rule->rhs->production_list, ti_list);
	free_transsys_instance_list(ti_list);
	break;
      }
    }
    if (rule == NULL)
      target_symbol = clone_symbol_instance(source_symbol);
    else
    {
      for (i = 1; i < rule->lhs->num_symbols; i++)
	source_symbol = source_symbol->next;
    }
    if (target_symbol)
    {
      if (target_tail == NULL)
	target_string = target_symbol;
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
  return (target_string);
}

