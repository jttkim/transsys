/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 */

/*
 * collection of functions to support parsing and conversion
 * from other languages such as Python.
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "trconfig.h"
#include "transsys.h"


void add_factor_definition(TRANSSYS *ts, FACTOR_ELEMENT *fe)
{
  fe->next = ts->factor_list;
  ts->factor_list = fe;
  fe->index = ts->num_factors++;
  /* fprintf(stderr, "factor \"%s\" has index %d\n", fe->name, fe->index); */
}


void add_gene_definition(TRANSSYS *ts, GENE_ELEMENT *ge)
{
  ge->next = ts->gene_list;
  ts->gene_list = ge;
  ge->index = ts->num_genes++;
}


void add_factordef_decay(EXPRESSION_NODE *decay_expression, FACTOR_ELEMENT *fe)
{
  free_expression_tree(fe->decay_expression);
  fe->decay_expression = decay_expression;
}


void add_factordef_diffusibility(EXPRESSION_NODE *diffusibility_expression, FACTOR_ELEMENT *fe)
{
  free_expression_tree(fe->diffusibility_expression);
  fe->diffusibility_expression = diffusibility_expression;
}


PROMOTER_ELEMENT *extend_promoter_list(PROMOTER_ELEMENT *alist, PROMOTER_ELEMENT *a)
{
  PROMOTER_ELEMENT *a1;

  if (alist == NULL)
  {
    return (a);
  }
  else
  {
    for (a1 = alist; a1->next; a1 = a1->next)
      ;
    a1->next = a;
    return (alist);
  }
}


int find_factor_index(const TRANSSYS *ts, const char *name)
{
  FACTOR_ELEMENT *fe;

  if (ts == NULL)
    return (-1);
  for (fe = ts->factor_list; fe; fe = fe->next)
  {
    if (!strcmp(name, fe->name))
      return (fe->index);
  }
  fprintf(stderr, "unknown factor \"%s\"\n", name);
  return (-1);
}


int find_lhs_symbol_index(const RULE_ELEMENT *re, const char *transsys_label)
{
  const LHS_SYMBOL *lhs_symbol;

  if (re == NULL)
  {
    return (-1);
  }
  for (lhs_symbol = re->lhs->symbol_list; lhs_symbol; lhs_symbol = lhs_symbol->next)
  {
    if (!strcmp(transsys_label, lhs_symbol->transsys_label))
      return (lhs_symbol->index);
  }
  fprintf(stderr, "unknown transsys label \"%s\" in rule \"%s\"\n", transsys_label, re->name);
  return (-1);
}


int resolve_simple_identifier(EXPRESSION_NODE *node, const TRANSSYS *transsys)
{
  int factor_index;

  factor_index = find_factor_index(transsys, node->content.raw_identifier.factor_name);
  if (factor_index < 0)
  {
    if (transsys)
    {
      fprintf(stderr, "resolve_simple_identifier: no factor \"%s\" in transsys \"%s\"\n", node->content.raw_identifier.factor_name, transsys->name);
    }
    else
    {
      fprintf(stderr, "resolve_simple_identifier: no transsys, no factor \"%s\"\n", node->content.raw_identifier.factor_name);
    }
    free(node->content.raw_identifier.factor_name);
    node->content.identifier.factor_index = factor_index;
    return (-1);
  }
  free(node->content.raw_identifier.factor_name);
  node->content.identifier.factor_index = factor_index;
  node->content.identifier.lhs_symbol_index = NO_INDEX;
  return (0);
}


int resolve_complex_identifier(const LSYS *lsys, EXPRESSION_NODE *node, const RULE_ELEMENT *rule_element)
{
  int lhs_symbol_index, return_value;
  int si;
  const SYMBOL_ELEMENT *sym;

  lhs_symbol_index = find_lhs_symbol_index(rule_element, node->content.raw_identifier.transsys_label);
  /* fprintf(stderr, "resolve_complex_identifier: resolved transsys label \"%s\" to %d\n", node->content.raw_identifier.transsys_label, lhs_symbol_index); */
  if (lhs_symbol_index < 0)
  {
    if (rule_element)
      fprintf(stderr, "resolve_complex_identifier: no transsys \"%s\" in rule \"%s\"\n", node->content.raw_identifier.transsys_label, rule_element->name);
    else
      fprintf(stderr, "resolve_complex_identifier: no rule, no transsys \"%s\"\n", node->content.raw_identifier.transsys_label);
    free(node->content.raw_identifier.transsys_label);
    free(node->content.raw_identifier.factor_name);
    node->content.identifier.lhs_symbol_index = NO_INDEX;
    node->content.identifier.factor_index = NO_INDEX;
    return (-1);
  }
  free(node->content.raw_identifier.transsys_label);
  si = rule_element->lhs->symbol_list[lhs_symbol_index].symbol_index;
  /* fprintf(stderr, "resolve_complex_identifier: symbol index is %d\n", si); */
  for (sym = lsys->symbol_list; sym; sym = sym->next)
  {
    if (sym->index == si)
      break;
  }
  if (sym == NULL)
  {
    fprintf(stderr, "resolve_complex_identifier: internal error: no symbol with index %d\n", si);
    node->content.identifier.factor_index = NO_INDEX;
    return (-1);
  }
  /* fprintf(stderr, "resolve_complex_identifier: symbol \"%s\" (index %d) has transsys %p\n", sym->name, sym->index, sym->transsys); */
  return_value = resolve_simple_identifier(node, sym->transsys);
  node->content.identifier.lhs_symbol_index = lhs_symbol_index;
  return (return_value);
}


int resolve_identifiers(EXPRESSION_NODE *node, const TRANSSYS *transsys, const LSYS *lsys, const RULE_ELEMENT *rule_element)
{
  int return_value;

  switch(node->type)
  {
  case NT_VALUE:
  case NT_IDENTIFIER:
    return (0);
  case NT_RAW_IDENTIFIER:
    node->type = NT_IDENTIFIER;
    if (node->content.raw_identifier.transsys_label)
    {
      /* fprintf(stderr, "resolving \"%s.%s\"\n", node->content.raw_identifier.transsys_label, node->content.raw_identifier.factor_name); */
      return (resolve_complex_identifier(lsys, node, rule_element));
    }
    else
    {
      /* fprintf(stderr, "resolving \"%s\"\n", node->content.raw_identifier.factor_name); */
      return (resolve_simple_identifier(node, transsys));
    }
  case NT_NOT:
  case NT_ATAN:
    return (resolve_identifiers(node->content.argument[0], transsys, lsys, rule_element));
  case NT_LOGICAL_OR:
  case NT_LOGICAL_AND:
  case NT_LOWER:
  case NT_LOWER_EQUAL:
  case NT_GREATER:
  case NT_GREATER_EQUAL:
  case NT_EQUAL:
  case NT_UNEQUAL:
  case NT_ADD:
  case NT_SUBTRACT:
  case NT_MULT:
  case NT_DIV:
  case NT_RANDOM:
  case NT_GAUSS:
  case NT_POW:
  case NT_LOG:
    /* CHECKME: is there a reason for this strange return strategy? */
    return_value = resolve_identifiers(node->content.argument[0], transsys, lsys, rule_element);
    if (resolve_identifiers(node->content.argument[1], transsys, lsys, rule_element) < 0)
    {
      return (-1);
    }
    else
    {
      return (return_value);
    }
  default:
    fprintf(stderr, "resolve_identifiers: unknown expression type %d\n", (int) node->type);
    return (-1);
  }
}


INTEGER_ARRAY *extend_factor_combination(INTEGER_ARRAY *ia, const char *name, TRANSSYS *transsys)
{
  int fi;

  fi = find_factor_index(transsys, name);
  if (fi < 0)
  {
    if (ia)
    {
      if (ia->length)
      {
	free(ia->array);
      }
      free(ia);
    }
    return (NULL);
  }
  ia = extend_integer_array(ia, fi);
  return (ia);
}


PROMOTER_ELEMENT *create_promoter(PROMOTERELEMENT_TYPE atype, INTEGER_ARRAY *ia, EXPRESSION_NODE *expr1, EXPRESSION_NODE *expr2)
{
  PROMOTER_ELEMENT *a;

  if (ia)
  {
    a = new_promoter_element(atype, ia->length, ia->array, expr1, expr2);
    free(ia);
  }
  else
  {
    a = new_promoter_element(atype, 0, NULL, expr1, expr2);
  }
  if (a == NULL)
  {
    fprintf(stderr, "create_promoter failed");
    return (NULL);
  }
  return (a);
}


GENE_ELEMENT *create_gene(const char *name, PROMOTER_ELEMENT *alist, int gene_product)
{
  GENE_ELEMENT *ge;

  ge = new_gene_element(name, alist, gene_product);
  if (ge == NULL)
  {
    fprintf(stderr, "create_gene failed");
    return (NULL);
  }
  return (ge);
}


/*
 * resolves identifier nodes. should only be called once after
 * assembly by parsing (or conversion from Python etc.)
 */
int resolve_transsys(TRANSSYS *tr)
{
  PROMOTER_ELEMENT *pe;
  FACTOR_ELEMENT *fe;
  GENE_ELEMENT *ge;
  int return_value;

  if (!tr->arrayed)
  {
    return_value = arrange_transsys_arrays(tr);
    if (return_value != 0)
    {
      fprintf(stderr, "resolve_transsys: arrange_transsys_arrays() returned %d\n", return_value);
      return (-1);
    }
  }
  for (fe = tr->factor_list; fe; fe = fe->next)
  {
    resolve_identifiers(fe->diffusibility_expression, tr, NULL, NULL);
    resolve_identifiers(fe->decay_expression, tr, NULL, NULL);
  }
  for (ge = tr->gene_list; ge; ge = ge->next)
  {
    for (pe = ge->promoter_list; pe; pe = pe->next)
    {
      resolve_identifiers(pe->expr1, tr, NULL, NULL);
      if (pe->type != PROMOTERELEMENT_CONSTITUTIVE)
      {
	resolve_identifiers(pe->expr2, tr, NULL, NULL);
      }
    }
  }
  return (0);
}


TRANSSYS *add_transsys(TRANSSYS *trlist, TRANSSYS *tr)
{
  TRANSSYS *tr1;

  if (resolve_transsys(tr) != 0)
  {
    free_transsys_list(trlist);
    free_transsys_list(tr);
    return (NULL);
  }
  if (trlist == NULL)
  {
    return (tr);
  }
  else
  {
    for (tr1 = trlist; tr1->next; tr1 = tr1->next)
      ;
    tr1->next = tr;
    return (trlist);
  }
}


GRAPHICS_PRIMITIVE *add_graphics_primitive(GRAPHICS_PRIMITIVE *gplist, GRAPHICS_PRIMITIVE *gp)
{
  gp->next = gplist;
  return (gp);
}


SYMBOL_ELEMENT *find_symbol(const LSYS *ls, const char *name)
{
  SYMBOL_ELEMENT *se;

  for (se = ls->symbol_list; se; se = se->next)
  {
    if (!strcmp(name, se->name))
    {
      return (se);
    }
  }
  fprintf(stderr, "unknown symbol \"%s\"\n", name);
  return (NULL);
}


int find_symbol_index(const LSYS *ls, const char *name)
{
  SYMBOL_ELEMENT *se;

  se = find_symbol(ls, name);
  if (se == NULL)
  {
    return (-1);
  }
  return (se->index);
}


const TRANSSYS *find_transsys(const TRANSSYS *trlist, const char *name)
{
  while (trlist)
  {
    if (!strcmp(name, trlist->name))
    {
      return (trlist);
    }
    trlist = trlist->next;
  }
  /* fprintf(stderr, "unknown transsys \"%s\"\n", name); */
  return (NULL);
}


SYMBOL_ELEMENT *create_symbol_element(const char *name, const char *transsys_name, const TRANSSYS *transsys_list)
{
  const TRANSSYS *tr = NULL;

  if (transsys_name)
  {
    tr = find_transsys(transsys_list, transsys_name);
    if (tr == NULL)
    {
      fprintf(stderr, "create_symbol_element: unknown transsys \"%s\"\n", transsys_name);
      return (NULL);
    }
  }
  return (new_symbol_element(name, tr));
}


int find_lhs_symbol_by_transsys_label(const LHS_SYMBOL *symbol_list, const char *transsys_label)
{
  const LHS_SYMBOL *sym;

/*
  {
    int i = 0;
    for (sym = symbol_list; sym; sym = sym->next)
    {
      fprintf(stderr, "  %d\n", i);
      i++;
    }
    fprintf(stderr, "find_lhs_symbol_by_transsys_label: symbol_list = %p, %d lhs_symbols\n", (void *) symbol_list, i);
  }
*/
  for (sym = symbol_list; sym; sym = sym->next)
  {
    if (!strcmp(transsys_label, sym->transsys_label))
    {
      return (sym->index);
    }
  }
  return (NO_INDEX);
}


ASSIGNMENT_RESOLUTION_RESULT resolve_raw_assignments(const SYMBOL_ELEMENT *se, LHS_SYMBOL *lhs_symbol_list, RAW_ASSIGNMENT *ra_list)
{
  ASSIGNMENT_RESOLUTION_RESULT a_result = { NULL, NO_INDEX, 0 };
  RAW_ASSIGNMENT *ra;
  ASSIGNMENT *a;
  int factor_index, lhs_symbol_index;

/*
  fprintf(stderr, "resolve_raw_assignments: symbol is \"%s\"", se->name);
  if (se->transsys)
    fprintf(stderr, ", transsys is \"%s\"", se->transsys->name);
  fprintf(stderr, "\n");
*/
  while (ra_list)
  {
    ra = ra_list;
    if (ra->identifier[0])
    {
      factor_index = find_factor_index(se->transsys, ra->identifier);
      if (factor_index >= 0)
      {
	/* fprintf(stderr, "resolve_raw_assignments: found factor index %d for factor \"%s\"\n", factor_index, ra->identifier); */
	a = new_assignment(se->transsys, factor_index, ra->value);
	if (a)
	{
	  a->next = a_result.assignment_list;
	  a_result.assignment_list = a;
	}
	else
	{
	  fprintf(stderr, "resolve_raw_assignments: new_assignment() failed");
	  free_expression_tree(ra->value);
	}
      }
      else
      {
	a_result.error = 1;
	if (se->transsys)
	{
	  /* FIXME: this should be an error resulting in termination of parsing */
	  fprintf(stderr, "resolve_raw_assignments: symbol \"%s\": no factor \"%s\" in transsys \"%s\"\n", se->name, ra->identifier, se->transsys->name);
	}
	else
	{
	  fprintf(stderr, "resolve_raw_assignments: symbol \"%s\": no transsys, no factor \"%s\"\n", se->name, ra->identifier);
	}
      }
    }
    else if (ra->transsys_label[0])
    {
      lhs_symbol_index = find_lhs_symbol_by_transsys_label(lhs_symbol_list, ra->transsys_label);
      /* fprintf(stderr, "resolve_raw_assignments: transsys label \"%s\": resolves to index %d", ra->transsys_label, lhs_symbol_index); */
      if (a_result.source_lhs_symbol_index != NO_INDEX)
      {
        fprintf(stderr, "resolve_raw_assignments: multiple transsys labels specified (inconsistent parser execution?)");
      }
      a_result.source_lhs_symbol_index = lhs_symbol_index;
    }
    ra_list = ra_list->next;
    free(ra);
  }
  /* fprintf(stderr, "resolve_raw_assignments: done\n\n"); */
  return (a_result);
}


LHS_DESCRIPTOR *create_lhs_descriptor(const LSYS *lsys, LHS_SYMBOL *symlist)
{
  LHS_DESCRIPTOR *lhs_descriptor;

  lhs_descriptor = new_lhs_descriptor(symlist);
  /* fprintf(stderr, "create_lhs_descriptor: single lhs_symbol \"%s\" has symbol_index %d = %d\n", symbol_name, lhs_descriptor->symbol_list[0].symbol_index, symbol->index); */
  arrange_lhs_descriptor_arrays(lhs_descriptor);
  /* this nonsense is hopefully now finally outdated for good. */
  /* provide global pointer to current list so RHS routines can resolve
     transsys labels */
  /* current_lhs_symbol_list = lhs_descriptor->symbol_list; */
  return (lhs_descriptor);
}


PRODUCTION_ELEMENT *create_production_element(const char *symbol_name, RAW_ASSIGNMENT *ra_list, const LSYS *lsys, LHS_SYMBOL *lhs_symbol_list)
{
  SYMBOL_ELEMENT *se;
  PRODUCTION_ELEMENT *sp;
  ASSIGNMENT_RESOLUTION_RESULT a_result;

  se = find_symbol(lsys, symbol_name);
  if (se == NULL)
  {
    fprintf(stderr, "create_production_element: unknown symbol \"%s\"\n", symbol_name);
    return (NULL);
  }
  if ((se->transsys == NULL) && (ra_list != NULL))
  {
    fprintf(stderr, "create_production_element: symbol \"%s\" has no transsys, but assignment list given\n", se->name);
    return (NULL);
  }
  a_result = resolve_raw_assignments(se, lhs_symbol_list, ra_list);
  if (a_result.error)
  {
    free_assignment_list(a_result.assignment_list);
    return (NULL);
  }
  sp = new_production_element(se->index, a_result.source_lhs_symbol_index, a_result.assignment_list);
  if (sp == NULL)
  {
    fprintf(stderr, "create_production_element: new_production_element() failed");
  }
  return (sp);
}


PRODUCTION_ELEMENT *add_production_element(PRODUCTION_ELEMENT *splist, PRODUCTION_ELEMENT *sp)
{
  if (sp == NULL)
    return (splist);
  else
  {
    sp->next = splist;
    return (sp);
  }
}


LHS_SYMBOL *create_lhs_symbol(const char *symbol_name, const char *transsys_label, const LSYS *lsys)
{
  SYMBOL_ELEMENT *se;
  LHS_SYMBOL *lhs_symbol;

  se = find_symbol(lsys, symbol_name);
  if (se == NULL)
  {
    return (NULL);
  }
  if (transsys_label)
  {
    lhs_symbol = new_lhs_symbol(transsys_label, se->transsys, se->index);
  }
  else
  {
    lhs_symbol = new_lhs_symbol(NULL, NULL, se->index);
  }
  if (lhs_symbol == NULL)
  {
    fprintf(stderr, "create_lhs_symbol: new_lhs_symbol() failed");
    return (NULL);
  }
  return (lhs_symbol);
}


/*
 * add LHS symbol to symlist. Differently from other add...() functions, also
 * manages indices
 */

LHS_SYMBOL *add_lhs_symbol(LHS_SYMBOL *symlist, LHS_SYMBOL *symbol)
{
  if (symlist == NULL)
  {
    symbol->index = 0;
    symbol->next = NULL;
  }
  else
  {
    symbol->index = symlist->index + 1;
    symbol->next = symlist;
  }
  /* fprintf(stderr, "add_lhs_symbol: added lhs_symbol with index %d\n", symbol->index); */
  return (symbol);
}


RAW_ASSIGNMENT *create_raw_assignment(const char *factor_name, const char *transsys_label, EXPRESSION_NODE *value)
{
  RAW_ASSIGNMENT *a;

  /* fprintf(stderr, "create_raw_assignment: identifier \"%s\"\n", factor_name); */
  a = (RAW_ASSIGNMENT *) malloc(sizeof(RAW_ASSIGNMENT));
  if (a == NULL)
  {
    fprintf(stderr, "create_raw_assignment: malloc failed");
    return (NULL);
  }
  a->next = NULL;
  if (factor_name)
  {
    strncpy(a->identifier, factor_name, IDENTIFIER_MAX);
    a->identifier[IDENTIFIER_MAX] = '\0';
    a->transsys_label[0] = '\0';
    a->value = value;
  }
  else if (transsys_label)
  {
    strncpy(a->transsys_label, transsys_label, IDENTIFIER_MAX);
    a->transsys_label[IDENTIFIER_MAX] = '\0';
    a->identifier[0] = '\0';
  }
  else
  {
    fprintf(stderr, "create_raw_assignment: neither factor name nor transsys label specified\n");
    a->transsys_label[0] = '\0';
  }
  return (a);
}


RAW_ASSIGNMENT *add_raw_assignment(RAW_ASSIGNMENT *ra_list, RAW_ASSIGNMENT *ra)
{
  ra->next = ra_list;
  return (ra);
}


int resolve_graphics_identifiers(GRAPHICS_PRIMITIVE *gp, const TRANSSYS *transsys)
{
  int return_value, i;

  /* fprintf(stderr, "resolve_graphics_identifiers: start\n"); */
  for (i = 0; i < gp->num_arguments; i++)
  {
    return_value = resolve_identifiers(gp->argument[i], transsys, NULL, NULL);
    if (return_value != 0)
    {
      fprintf(stderr, "resolve_graphics_identifiers: resolve_identifiers() returned %d", return_value);
      return (-1);
    }
  }
  /* fprintf(stderr, "resolve_graphics_identifiers: resolved %d expressions\n", i); */
  return (0);
}


int resolve_lsys(LSYS *ls)
{
  int return_value;
  SYMBOL_ELEMENT *se;
  GRAPHICS_PRIMITIVE *gp;
  int error = 0;

  for (se = ls->symbol_list; se; se = se->next)
  {
    /* fprintf(stderr, "resolve_lsys: resolving graphics expressions for symbol \"%s\"\n", se->name); */
    for (gp = se->graphics_primitive_list; gp; gp = gp->next)
    {
      return_value = resolve_graphics_identifiers(gp, se->transsys);
      if (return_value != 0)
      {
	error = 1;
      }
    }
    return_value = arrange_symbol_element_arrays(se);
    if (return_value != 0)
    {
      fprintf(stderr, "add_lsys: arrange_symbol_element_arrays() returned %d\n", return_value);
    }
  }
  if (error)
  {
    return (-1);
  }
/*
  fprintf(stderr, "resolve_lsys: resolved graphics identifiers\n");
  if (lsyslist)
    fprint_lsys(stderr, 4, lsyslist);
*/
  return_value = arrange_lsys_arrays(ls);
/*
  fprintf(stderr, "resolve_lsys: arranged arrays\n");
  if (lsyslist)
    fprint_lsys(stderr, 4, lsyslist);
*/
  if (return_value != 0)
  {
    fprintf(stderr, "resolve_lsys: arrange_lsys_arrays() returned %d\n", return_value);
  }
  return (0);
}


LSYS *add_lsys(LSYS *lsyslist, LSYS *ls)
{
  LSYS *ls1;

  if (resolve_lsys(ls) != 0)
  {
    free_lsys_list(ls);
    free_lsys_list(lsyslist);
    return (NULL);
  }
  /* fprintf(stderr, "add_lsys: added lsys \"%s\"\n", ls->name); */
  ls->next = NULL;
  if (lsyslist == NULL)
  {
    lsyslist = ls;
  }
  else
  {
    for (ls1 = lsyslist; ls1->next; ls1 = ls1->next)
      ;
    ls1->next = ls;
  }
/*
  {
    LSYS *l;

    for (l = lsyslist; l; l = l->next)
    {
      fprintf(stderr, "add_lsys: lsys \"%s\"\n", l->name);
      fprint_lsys(stderr, 4, l);
    }
  }
*/
  return (lsyslist);
}


void add_symbol_definition(LSYS *ls, SYMBOL_ELEMENT *se)
{
  se->next = ls->symbol_list;
  ls->symbol_list = se;
  se->index = ls->num_symbols++;
  /* fprintf(stderr, "add_symbol_definition: adding symbol \"%s\" to lsys \"%s\"\n", se->name, ls->name); */
}


void add_axiom_definition(LSYS *ls, SYMBOL_PRODUCTION *sp)
{
  if (ls->axiom)
  {
    free_symbol_production_components(ls->axiom);
    free(ls->axiom);
  }
  if (arrange_symbol_production_arrays(sp) != 0)
    fprintf(stderr, "add_axiom_definition: arrange_symbol_production_arrays() failed");
  ls->axiom = sp;
  /* ls->axiom->transsys = NULL; */
  /* fprintf(stderr, "add_axiom_definition: added axiom\n"); */
}


void add_diffusionrange_definition(LSYS *ls, double diffusionrange)
{
  ls->diffusion_range = (int) floor(diffusionrange);
}


int resolve_rule_identifiers(RULE_ELEMENT *re, LSYS *lsys)
{
  PRODUCTION_ELEMENT *pe;
  ASSIGNMENT *a;
  int return_value;

  for (pe = re->rhs->production_list; pe; pe = pe->next)
  {
    for (a = pe->assignment_list; a; a = a-> next)
    {
      return_value = resolve_identifiers(a->value, NULL, lsys, re);
      if (return_value != 0)
      {
	return (return_value);
      }
    }
  }
  return (0);
}


int add_rule_definition(LSYS *ls, RULE_ELEMENT *re)
{
  SYMBOL_ELEMENT *se;
  int return_value;

  /* fprintf(stderr, "add_rule_definition: starting, ls = %p, re = %p\n", (void *) ls, (void *) re); */
  for (se = ls->symbol_list; se; se = se->next)
  {
    /* fprintf(stderr, "add_rule_definition: se = %p, se->index = %d, re->lhs->symbol_list = %p\n", (void *) se, se->index, (void *) re->lhs->symbol_list); */
    if (se->index == re->lhs->symbol_list->symbol_index)
    {
      /* re->rhs->transsys = se->transsys; */
      /* fprintf(stderr, "add_rule_definition: resolving rule identifiers\n"); */
      return_value = resolve_rule_identifiers(re, ls);
      if (return_value != 0)
      {
	fprintf(stderr, "add_rule_definition: resolve_rule_identifiers returned %d", return_value);
	return (return_value);
      }
      if (re->condition)
      {
	/* fprintf(stderr, "add_rule_definition: resolving condition identifiers\n"); */
	return_value = resolve_identifiers(re->condition, NULL, ls, re);
	if (return_value != 0)
	{
	  fprintf(stderr, "add_rule_definition: resolve_identifiers() returned %d", return_value);
	  return (return_value);
	}
      }
      break;
    }
  }
  re->next = ls->rule_list;
  ls->rule_list = re;
  ls->num_rules++;
  return (0);
}


void add_graphics_to_symbol(SYMBOL_ELEMENT *se, GRAPHICS_PRIMITIVE *gp_list)
{
  if (se->graphics_primitive_list)
  {
    fprintf(stderr, "add_graphics_to_named_symbol: graphics of symbol \"%s\" superseded", se->name);
    free_graphics_primitive_list(se->graphics_primitive_list);
  }
  se->graphics_primitive_list = gp_list;
}


void add_graphics_to_named_symbol(LSYS *ls, const char *symbol_name, GRAPHICS_PRIMITIVE *gp_list)
{
  SYMBOL_ELEMENT *se;

  se = find_symbol(ls, symbol_name);
  if (se == NULL)
  {
    fprintf(stderr, "add_graphics_to_named_symbol: discarding graphics for unknown symbol %s", symbol_name);
    free_graphics_primitive_list(gp_list);
    return;
  }
  add_graphics_to_symbol(se, gp_list);
}
