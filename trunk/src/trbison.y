/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

%{

/*
 * $Id$
 *
 * $Log$
 * Revision 1.6  2005/06/22 09:58:36  jtk
 * prevented unknown variables from resulting in core-dump eliciting parsing results
 *
 * Revision 1.5  2005/06/15 22:17:13  jtk
 * counting number of transsys programs in lsys (deprecating multiples)
 *
 * Revision 1.4  2005/04/14 18:44:46  jtk
 * fixed parser slipthrough of symbols declared with nonexistent transsys programs
 *
 * Revision 1.3  2005/04/01 13:37:05  jtk
 * Made parser somewhat more strict
 *
 * Revision 1.2  2005/03/31 16:07:36  jtk
 * finished (initial) implementation of lsys diffusion
 *
 * Revision 1.1.1.1  2005/03/08 17:12:02  jtk
 * new cvs after loss at INB
 *
 * Revision 1.2  2003/01/22 11:39:24  kim
 * added dot.c
 *
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#include <ctype.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "trconfig.h"
#include "transsys.h"


typedef struct
{
  ASSIGNMENT *assignment_list;
  int source_lhs_symbol_index;
  int error;
} ASSIGNMENT_RESOLUTION_RESULT;


long lineno = 1;
char *yyin_name = NULL;
extern int yychar;


TRANSSYS *parsed_transsys = NULL;
LSYS *parsed_lsys = NULL;

static LSYS *current_lsys = NULL;
static TRANSSYS *current_transsys = NULL;
static FACTOR_ELEMENT *current_factor = NULL;
static LHS_SYMBOL *current_lhs_symbol_list;


static void list_parsed_stuff(const char *msg)
{
  TRANSSYS *t;
  LSYS *l;
  int indent = 0;

  if (msg)
  {
    fprintf(stderr, "list_parsed_stuff: %s\n", msg);
    indent = 4;
  }
  return;
/*
  for (t = parsed_transsys; t; t = t->next)
    fprint_transsys(stderr, indent, t);
*/
  for (l = parsed_lsys; l; l = l->next)
    fprint_lsys(stderr, indent, l);
}


void yyerror(const char *s, ...)
{
  va_list arglist;

  va_start(arglist, s);
  fprintf(stderr, "%s:%ld: ", yyin_name, lineno);
  vfprintf(stderr, s, arglist);
  fprintf(stderr, "\n");
  switch (yychar)
  {
  case REALVALUE:
    fprintf(stderr, "    token: REALVALUE:\n");
    break;
  case IDENTIFIER:
    fprintf(stderr, "    token: IDENTIFIER:\n");
    break;
  case FACTOR_DEF:
    fprintf(stderr, "    token: FACTOR_DEF:\n");
    break;
  case GENE_DEF:
    fprintf(stderr, "    token: GENE_DEF:\n");
    break;
  case PROMOTER_DEF:
    fprintf(stderr, "    token: PROMOTER_DEF:\n");
    break;
  case PRODUCT_DEF:
    fprintf(stderr, "    token: PRODUCT_DEF:\n");
    break;
  case CONSTITUTIVE:
    fprintf(stderr, "    token: CONSTITUTIVE:\n");
    break;
  case ACTIVATE:
    fprintf(stderr, "    token: ACTIVATE:\n");
    break;
  case REPRESS:
    fprintf(stderr, "    token: REPRESS:\n");
    break;
  case DEFAULT:
    fprintf(stderr, "    token: DEFAULT:\n");
    break;
  case RANDOM:
    fprintf(stderr, "    token: RANDOM:\n");
    break;
  case GAUSS:
    fprintf(stderr, "    token: GAUSS:\n");
    break;
  case TRANSSYS_DEF:
    fprintf(stderr, "    token: TRANSSYS_DEF:\n");
    break;
  case DECAY_DEF:
    fprintf(stderr, "    token: DECAY_DEF:\n");
    break;
  case DIFFUSIBILITY_DEF:
    fprintf(stderr, "    token: DIFFUSIBILITY_DEF:\n");
    break;
  case DIFFUSIONRANGE_DEF:
    fprintf(stderr, "    token: DIFFUSIONRANGE_DEF:\n");
    break;
  case LSYS_DEF:
    fprintf(stderr, "    token: LSYS_DEF:\n");
    break;
  case SYMBOL_DEF:
    fprintf(stderr, "    token: SYMBOL_DEF:\n");
    break;
  case RULE_DEF:
    fprintf(stderr, "    token: RULE_DEF:\n");
    break;
  case AXIOM_DEF:
    fprintf(stderr, "    token: AXIOM_DEF:\n");
    break;
  case GRAPHICS_DEF:
    fprintf(stderr, "    token: GRAPHICS_DEF:\n");
    break;
  case ARROW:
    fprintf(stderr, "    token: ARROW:\n");
    break;
  case MOVE:
    fprintf(stderr, "    token: MOVE:\n");
    break;
  case COLOR:
    fprintf(stderr, "    token: COLOR:\n");
    break;
  case SPHERE:
    fprintf(stderr, "    token: SPHERE:\n");
    break;
  case CYLINDER:
    fprintf(stderr, "    token: CYLINDER:\n");
    break;
  case BOX:
    fprintf(stderr, "    token: BOX:\n");
    break;
  case TURN:
    fprintf(stderr, "    token: TURN:\n");
    break;
  case ROLL:
    fprintf(stderr, "    token: ROLL:\n");
    break;
  case BANK:
    fprintf(stderr, "    token: BANK:\n");
    break;
  case LOWER_EQUAL:
    fprintf(stderr, "    token: LOWER_EQUAL:\n");
    break;
  case GREATER_EQUAL:
    fprintf(stderr, "    token: GREATER_EQUAL:\n");
    break;
  case EQUAL:
    fprintf(stderr, "    token: EQUAL:\n");
    break;
  case UNEQUAL:
    fprintf(stderr, "    token: UNEQUAL:\n");
    break;
  case LOGICAL_AND:
    fprintf(stderr, "    token: LOGICAL_AND:\n");
    break;
  case LOGICAL_OR:
    fprintf(stderr, "    token: LOGICAL_OR:\n");
    break;
  default:
    if (yychar < 256)
    {
      if (isprint(yychar))
	fprintf(stderr, "    token: \'%c\'\n", yychar);
      else
	fprintf(stderr, "    token: #%d\n", yychar);
    }
    else
      fprintf(stderr, "    token: #%d\n", yychar);
    break;
  }
}


static void ignore_message(char *format, ...)
{
  va_list ap;

  va_start(ap, format);
  fprintf(stderr, "%s:%ld: ignoring ", yyin_name, lineno);
  vfprintf(stderr, format, ap);
  fprintf(stderr, "\n");
}


/*** original transsys functions ***/

static void add_factor_definition(TRANSSYS *ts, FACTOR_ELEMENT *fe)
{
  fe->next = ts->factor_list;
  ts->factor_list = fe;
  fe->index = ts->num_factors++;
  /* fprintf(stderr, "factor \"%s\" has index %d\n", fe->name, fe->index); */
}


static void add_gene_definition(TRANSSYS *ts, GENE_ELEMENT *ge)
{
  ge->next = ts->gene_list;
  ts->gene_list = ge;
  ge->index = ts->num_genes++;
}


static int setup_current_factor(const char *name)
{
  EXPRESSION_NODE *de, *di;

  de = new_expression_node(NT_VALUE, 1.0);
  if (de == NULL)
  {
    yyerror("new_expression_node (de) failed in setup_current_factor");
    return (-1);
  }
  di = new_expression_node(NT_VALUE, 1.0);
  if (di == NULL)
  {
    yyerror("new_expression_node (di) failed in setup_current_factor");
    return (-1);
  }
  current_factor = new_factor_element(name, de, di);
  if (current_factor == NULL)
  {
    yyerror("setup_current_factor failed");
    return (-1);
  }
  /* fprintf(stderr, "created factor \"%s\"\n", current_factor->name); */
  return (0);
}


static void add_factordef_decay(EXPRESSION_NODE *decay_expression)
{
  free_expression_tree(current_factor->decay_expression);
  current_factor->decay_expression = decay_expression;
}


static void add_factordef_diffusibility(EXPRESSION_NODE *diffusibility_expression)
{
  free_expression_tree(current_factor->diffusibility_expression);
  current_factor->diffusibility_expression = diffusibility_expression;
}


static PROMOTER_ELEMENT *extend_promoter_list(PROMOTER_ELEMENT *alist, PROMOTER_ELEMENT *a)
{
  PROMOTER_ELEMENT *a1;

  if (alist == NULL)
    return (a);
  else
  {
    for (a1 = alist; a1->next; a1 = a1->next)
      ;
    a1->next = a;
    return (alist);
  }
}


static int find_factor_index(const TRANSSYS *ts, const char *name)
{
  FACTOR_ELEMENT *fe;

  if (ts == NULL)
    return (-1);
  for (fe = ts->factor_list; fe; fe = fe->next)
  {
    if (!strcmp(name, fe->name))
      return (fe->index);
  }
  yyerror("unknown factor \"%s\"", name);
  return (-1);
}


static int find_lhs_symbol_index(const RULE_ELEMENT *re, const char *transsys_label)
{
  const LHS_SYMBOL *lhs_symbol;

  if (re == NULL)
    return (-1);
  for (lhs_symbol = re->lhs->symbol_list; lhs_symbol; lhs_symbol = lhs_symbol->next)
  {
    if (!strcmp(transsys_label, lhs_symbol->transsys_label))
      return (lhs_symbol->index);
  }
  yyerror("unknown transsys label \"%s\" in rule \"%s\"\n", transsys_label, re->name);
  return (-1);
}


static int resolve_simple_identifier(EXPRESSION_NODE *node, const TRANSSYS *transsys)
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


static int resolve_complex_identifier(EXPRESSION_NODE *node, const RULE_ELEMENT *rule_element)
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
  for (sym = current_lsys->symbol_list; sym; sym = sym->next)
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


static int resolve_identifiers(EXPRESSION_NODE *node, const TRANSSYS *transsys, const RULE_ELEMENT *rule_element)
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
      return (resolve_complex_identifier(node, rule_element));
    }
    else
    {
      /* fprintf(stderr, "resolving \"%s\"\n", node->content.raw_identifier.factor_name); */
      return (resolve_simple_identifier(node, transsys));
    }
  case NT_NOT:
    return (resolve_identifiers(node->content.argument[0], transsys, rule_element));
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
    return_value = resolve_identifiers(node->content.argument[0], transsys, rule_element);
    if (resolve_identifiers(node->content.argument[1], transsys, rule_element) < 0)
      return (-1);
    else
      return (return_value);
  default:
    fprintf(stderr, "resolve_identifiers: unknown expression type %d\n", (int) node->type);
    return (-1);
  }
}


static INTEGER_ARRAY *extend_factor_combination(INTEGER_ARRAY *ia, const char *name)
{
  int fi;

  fi = find_factor_index(current_transsys, name);
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


static PROMOTER_ELEMENT *create_promoter(ACTIVATION_TYPE atype, INTEGER_ARRAY *ia, EXPRESSION_NODE *expr1, EXPRESSION_NODE *expr2)
{
  PROMOTER_ELEMENT *a;

  if (ia)
  {
    a = new_promoter_element(atype, ia->length, ia->array, expr1, expr2);
    free(ia);
  }
  else
    a = new_promoter_element(atype, 0, NULL, expr1, expr2);
  if (a == NULL)
  {
    yyerror("create_promoter failed");
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
    yyerror("create_gene failed");
    return (NULL);
  }
  return (ge);
}


static int setup_current_transsys(const char *name)
{
  current_transsys = new_transsys(name);
  if (current_transsys == NULL)
  {
    yyerror("failed to set up new transsys");
    return (-1);
  }
  return (0);
}


static TRANSSYS *add_transsys(TRANSSYS *trlist, TRANSSYS *tr)
{
  int return_value;
  PROMOTER_ELEMENT *pe;
  FACTOR_ELEMENT *fe;
  GENE_ELEMENT *ge;
  TRANSSYS *tr1;

  return_value = arrange_transsys_arrays(tr);
  if (return_value != 0)
    fprintf(stderr, "add_transsys: arrange_transsys_arrays() returned %d\n", return_value);
  for (fe = tr->factor_list; fe; fe = fe->next)
  {
    resolve_identifiers(fe->diffusibility_expression, tr, NULL);
    resolve_identifiers(fe->decay_expression, tr, NULL);
  }
  for (ge = tr->gene_list; ge; ge = ge->next)
  {
    for (pe = ge->promoter_list; pe; pe = pe->next)
    {
      resolve_identifiers(pe->expr1, tr, NULL);
      if (pe->type != ACT_CONSTITUTIVE)
	resolve_identifiers(pe->expr2, tr, NULL);
    }
  }
  tr->next = NULL;
  if (trlist == NULL)
    return (tr);
  else
  {
    for (tr1 = trlist; tr1->next; tr1 = tr1->next)
      ;
    tr1->next = tr;
    return (trlist);
  }
}


static GRAPHICS_PRIMITIVE *add_graphics_primitive(GRAPHICS_PRIMITIVE *gplist, GRAPHICS_PRIMITIVE *gp)
{
  gp->next = gplist;
  return (gp);
}


static SYMBOL_ELEMENT *find_symbol(const LSYS *ls, const char *name)
{
  SYMBOL_ELEMENT *se;

  for (se = ls->symbol_list; se; se = se->next)
  {
    if (!strcmp(name, se->name))
      return (se);
  }
  yyerror("unknown symbol \"%s\"", name);
  return (NULL);
}


static int find_symbol_index(const LSYS *ls, const char *name)
{
  SYMBOL_ELEMENT *se;

  se = find_symbol(ls, name);
  if (se == NULL)
    return (-1);
  return (se->index);
}


static TRANSSYS *find_transsys(TRANSSYS *trlist, const char *name)
{
  while (trlist)
  {
    if (!strcmp(name, trlist->name))
      return (trlist);
    trlist = trlist->next;
  }
  yyerror("unknown transsys \"%s\"", name);
  return (NULL);
}


static SYMBOL_ELEMENT *create_symbol_element(const char *name, const char *transsys_name)
{
  TRANSSYS *tr = NULL;

  if (transsys_name)
  {
    tr = find_transsys(parsed_transsys, transsys_name);
    if (tr == NULL)
      return (NULL);
  }
  return (new_symbol_element(name, tr));
}


static int find_lhs_symbol_by_transsys_label(const LHS_SYMBOL *symbol_list, const char *transsys_label)
{
  const LHS_SYMBOL *sym;

  for (sym = symbol_list; sym; sym = sym->next)
  {
    if (!strcmp(transsys_label, sym->transsys_label))
      return (sym->index);
  }
  return (-1);
}


static ASSIGNMENT_RESOLUTION_RESULT resolve_raw_assignments(const SYMBOL_ELEMENT *se, RAW_ASSIGNMENT *ra_list)
{
  ASSIGNMENT_RESOLUTION_RESULT a_result = { NULL, NO_INDEX, 0 };
  RAW_ASSIGNMENT *ra;
  ASSIGNMENT *a;
  int factor_index, lhs_symbol_index;

/*
  fprintf(stderr, "resolve_raw_assignment: symbol is \"%s\"", se->name);
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
	  /* FIXME: is this really an yyerror? */
	  yyerror("resolve_raw_assignments: new_assignment() failed");
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
      lhs_symbol_index = find_lhs_symbol_by_transsys_label(current_lhs_symbol_list, ra->transsys_label);
      /* yyerror("resolve_raw_assignments: transsys label \"%s\": resolves to index %d", ra->transsys_label, lhs_symbol_index); */
      if (a_result.source_lhs_symbol_index != NO_INDEX)
        yyerror("resolve_raw_assignments: multiple transsys labels specified (inconsistent parser execution?)");
      a_result.source_lhs_symbol_index = lhs_symbol_index;
    }
    ra_list = ra_list->next;
    free(ra);
  }
  /* fprintf(stderr, "resolve_raw_assignments: done\n\n"); */
  return (a_result);
}


/*
 * preliminary: make lhs descriptor with just one lhs symbol
 */

static LHS_DESCRIPTOR *create_lhs_descriptor(const LSYS *lsys, LHS_SYMBOL *symlist)
{
  LHS_DESCRIPTOR *lhs_descriptor;

  lhs_descriptor = new_lhs_descriptor(symlist);
  /* fprintf(stderr, "create_lhs_descriptor: single lhs_symbol \"%s\" has symbol_index %d = %d\n", symbol_name, lhs_descriptor->symbol_list[0].symbol_index, symbol->index); */
  arrange_lhs_descriptor_arrays(lhs_descriptor);
  
  /* provide global pointer to current list so RHS routines can resolve
     transsys labels */
  current_lhs_symbol_list = lhs_descriptor->symbol_list;
  return (lhs_descriptor);
}


static PRODUCTION_ELEMENT *create_production_element(const char *symbol_name, RAW_ASSIGNMENT *ra_list)
{
  SYMBOL_ELEMENT *se;
  PRODUCTION_ELEMENT *sp;
  ASSIGNMENT_RESOLUTION_RESULT a_result;

  se = find_symbol(current_lsys, symbol_name);
  if (se == NULL)
    return (NULL);
  a_result = resolve_raw_assignments(se, ra_list);
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


static PRODUCTION_ELEMENT *add_production_element(PRODUCTION_ELEMENT *splist, PRODUCTION_ELEMENT *sp)
{
  if (sp == NULL)
    return (splist);
  else
  {
    sp->next = splist;
    return (sp);
  }
}


static LHS_SYMBOL *create_lhs_symbol(const char *symbol_name, const char *transsys_label)
{
  SYMBOL_ELEMENT *se;
  LHS_SYMBOL *lhs_symbol;

  se = find_symbol(current_lsys, symbol_name);
  if (se == NULL)
    return (NULL);
  if (transsys_label)
    lhs_symbol = new_lhs_symbol(transsys_label, se->transsys, se->index);
  else
    lhs_symbol = new_lhs_symbol(NULL, NULL, se->index);
  if (lhs_symbol == NULL)
  {
    yyerror("create_lhs_descriptor: new_lhs_symbol() failed");
    return (NULL);
  }
  return (lhs_symbol);
}


/*
 * add LHS symbol to symlist. Differently from other add...() functions, also
 * manages indices
 */

static LHS_SYMBOL *add_lhs_symbol(LHS_SYMBOL *symlist, LHS_SYMBOL *symbol)
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


static RAW_ASSIGNMENT *create_raw_assignment(const char *factor_name, const char *transsys_label, EXPRESSION_NODE *value)
{
  RAW_ASSIGNMENT *a;

  /* fprintf(stderr, "create_raw_assignment: identifier \"%s\"\n", factor_name); */
  a = (RAW_ASSIGNMENT *) malloc(sizeof(RAW_ASSIGNMENT));
  if (a == NULL)
    yyerror("create_raw_assignment: malloc failed");
  else
  {
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
      a->transsys_label[0] = '\0';
    }
  }
  return (a);
}


static RAW_ASSIGNMENT *add_raw_assignment(RAW_ASSIGNMENT *ra_list, RAW_ASSIGNMENT *ra)
{
  ra->next = ra_list;
  return (ra);
}


static int setup_current_lsys(const char *name)
{
  current_lsys = new_lsys(name);
  if (current_lsys == NULL)
  {
    yyerror("failed to set up new lsys");
    return (-1);
  }
  return (0);
}


static int resolve_graphics_identifiers(GRAPHICS_PRIMITIVE *gp, const TRANSSYS *transsys)
{
  int return_value, i;

  /* fprintf(stderr, "resolve_graphics_identifiers: start\n"); */
  for (i = 0; i < gp->num_arguments; i++)
  {
    return_value = resolve_identifiers(gp->argument[i], transsys, NULL);
    if (return_value != 0)
    {
      yyerror("resolve_graphics_identifiers: resolve_identifiers() returned %d", return_value);
      return (-1);
    }
  }
  /* fprintf(stderr, "resolve_graphics_identifiers: resolved %d expressions\n", i); */
  return (0);
}


static LSYS *add_lsys(LSYS *lsyslist, LSYS *ls)
{
  int return_value;
  SYMBOL_ELEMENT *se;
  GRAPHICS_PRIMITIVE *gp;
  LSYS *ls1;
  int error = 0;

/*
  fprintf(stderr, "add_lsys: starting\n");
  if (lsyslist)
    fprint_lsys(stderr, 4, lsyslist);
*/
  for (se = ls->symbol_list; se; se = se->next)
  {
    /* fprintf(stderr, "add_lsys: resolving graphics expressions for symbol \"%s\"\n", se->name); */
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
      yyerror("add_lsys: arrange_symbol_element_arrays() returned %d\n", return_value);
    }
  }
  if (error)
  {
    free_lsys_list(ls);
    free_lsys_list(lsyslist);
    return (NULL);
  }
/*
  fprintf(stderr, "add_lsys: resolved graphics identifiers\n");
  if (lsyslist)
    fprint_lsys(stderr, 4, lsyslist);
*/
  return_value = arrange_lsys_arrays(ls);
/*
  fprintf(stderr, "add_lsys: arranged arrays\n");
  if (lsyslist)
    fprint_lsys(stderr, 4, lsyslist);
*/
  if (return_value != 0)
  {
    yyerror("add_lsys: arrange_lsys_arrays() returned %d\n", return_value);
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


static void add_symbol_definition(LSYS *ls, SYMBOL_ELEMENT *se)
{
  se->next = ls->symbol_list;
  ls->symbol_list = se;
  se->index = ls->num_symbols++;
  /* fprintf(stderr, "add_symbol_definition: adding symbol \"%s\" to lsys \"%s\"\n", se->name, ls->name); */
}


static void add_axiom_definition(LSYS *ls, SYMBOL_PRODUCTION *sp)
{
  if (ls->axiom)
  {
    free_symbol_production_components(ls->axiom);
    free(ls->axiom);
  }
  if (arrange_symbol_production_arrays(sp) != 0)
    yyerror("add_axiom_definition: arrange_symbol_production_arrays() failed");
  ls->axiom = sp;
  ls->axiom->transsys = NULL;
  /* fprintf(stderr, "add_axiom_definition: added axiom\n"); */
}


static void add_diffusionrange_definition(LSYS *ls, double diffusionrange)
{
  ls->diffusion_range = (int) floor(diffusionrange);
}


static int resolve_rule_identifiers(RULE_ELEMENT *re)
{
  PRODUCTION_ELEMENT *pe;
  ASSIGNMENT *a;
  int return_value;

  for (pe = re->rhs->production_list; pe; pe = pe->next)
  {
    for (a = pe->assignment_list; a; a = a-> next)
    {
      /* obsolescent call: resolve_identifiers(a->value, re->rhs->transsys); */
      return_value = resolve_identifiers(a->value, NULL, re);
      if (return_value != 0)
      {
	return (return_value);
      }
    }
  }
  return (0);
}


static RULE_ELEMENT *complete_rule(const char *name, RULE_ELEMENT *re)
{
  strncpy(re->name, name, IDENTIFIER_MAX);
  re->name[IDENTIFIER_MAX - 1] = '\0';
  if (arrange_symbol_production_arrays(re->rhs) != 0)
    yyerror("complete_rule: arrange_symbol_production_arrays() failed");
  return (re);
}


/*
 * preliminary: use single (i.e. first) lhs symbol to obtain "governing"
 * transsys for rhs -- this concept needs to be replaced anyway...
 */ 

static int add_rule_definition(LSYS *ls, RULE_ELEMENT *re)
{
  SYMBOL_ELEMENT *se;
  int return_value;

  for (se = ls->symbol_list; se; se = se->next)
  {
    if (se->index == re->lhs->symbol_list->symbol_index)
    {
      re->rhs->transsys = se->transsys;
      return_value = resolve_rule_identifiers(re);
      if (return_value != 0)
      {
	yyerror("add_rule_definition: resolve_rule_identifiers returned %d", return_value);
	return (return_value);
      }
      if (re->condition)
      {
	return_value = resolve_identifiers(re->condition, NULL, re);
	if (return_value != 0)
	{
	  yyerror("add_rule_definition: resolve_identifiers() returned %d", return_value);
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


static void add_graphics_to_symbol(LSYS *ls, const char *symbol_name, GRAPHICS_PRIMITIVE *gp_list)
{
  SYMBOL_ELEMENT *se;

  se = find_symbol(ls, symbol_name);
  if (se == NULL)
  {
    yyerror("add_graphics_to_symbol: discarding graphics for unknown symbol %s", symbol_name);
    free_graphics_primitive_list(gp_list);
    return;
  }
  if (se->graphics_primitive_list)
  {
    yyerror("add_graphics_to_symbol: graphics of symbol \"%s\" superseded", symbol_name);
    free_graphics_primitive_list(se->graphics_primitive_list);
  }
  se->graphics_primitive_list = gp_list;
}


%}


%union {
  int integer;
  INTEGER_ARRAY *intarr;
  double real;
  char identifier[IDENTIFIER_MAX];
  EXPRESSION_NODE *expression;
  FACTOR_ELEMENT *factor;
  PROMOTER_ELEMENT *promoter;
  GENE_ELEMENT *gene;
  TRANSSYS *transsys;
  SYMBOL_ELEMENT *symbol_element;
  RULE_ELEMENT *rule_element;
  LHS_DESCRIPTOR *lhs_descriptor;
  LHS_SYMBOL *lhs_symbol;
  PRODUCTION_ELEMENT *production_element;
  SYMBOL_PRODUCTION *symbol_production;
  RAW_ASSIGNMENT *raw_assignment;
  GRAPHICS_PRIMITIVE *graphics_primitive;
}


%left '+' '-'
%left '*' '/'
%left LOGICAL_AND LOGICAL_OR
%token <real> REALVALUE
%token <identifier> IDENTIFIER
%token FACTOR_DEF GENE_DEF PROMOTER_DEF PRODUCT_DEF CONSTITUTIVE ACTIVATE REPRESS DEFAULT RANDOM GAUSS TRANSSYS_DEF DECAY_DEF DIFFUSIBILITY_DEF
%token LSYS_DEF SYMBOL_DEF RULE_DEF AXIOM_DEF DIFFUSIONRANGE_DEF GRAPHICS_DEF ARROW MOVE COLOR SPHERE CYLINDER BOX TURN ROLL BANK PUSH POP
%token LOWER_EQUAL GREATER_EQUAL EQUAL UNEQUAL LOGICAL_AND LOGICAL_OR

%type <factor> factor_definition
%type <gene> gene_definition
%type <integer> product_component product_statements
%type <lhs_descriptor> rule_lhs
%type <intarr> factor_combination
%type <promoter> promoter_component promoter_statements promoter_statement
%type <expression> expr and_expr not_expr cmp_expr arithmetic_expr term value
%type <symbol_element> symbol_definition
%type <real> diffusionrange_definition
%type <production_element> production_element_string production_element
%type <symbol_production> axiom_definition rule_rhs
%type <rule_element> rule_definition rule_components
%type <graphics_primitive> graphcmd graphcmd_list
%type <raw_assignment> assignment assignment_list source_transsys_specifier transsys_initializer
%type <lhs_symbol> lhs_element lhs_element_string

%%

ltrfile
	: /* empty */
	| ltrfile  lsys { parsed_lsys = add_lsys(parsed_lsys, current_lsys); if (parsed_lsys == NULL) YYABORT; }
	| ltrfile transsys  { parsed_transsys  = add_transsys(parsed_transsys, current_transsys); }
	| error ';' { yyerrok; }
	;

lsys
	: LSYS_DEF IDENTIFIER { setup_current_lsys($2); } '{' lsys_element_list '}' {}
	;

lsys_element_list
	: lsys_element {}
	| lsys_element_list lsys_element {}
	;

lsys_element
	: symbol_definition { add_symbol_definition(current_lsys, $1); }
	| axiom_definition { add_axiom_definition(current_lsys, $1); }
	| diffusionrange_definition { add_diffusionrange_definition(current_lsys, $1); }
	| rule_definition { if (add_rule_definition(current_lsys, $1) != 0) YYABORT; }
	| graphics_definition {}
	;

symbol_definition
	: SYMBOL_DEF IDENTIFIER ';' { $$ = create_symbol_element($2, NULL); if ($$ == NULL) YYABORT; }
	| SYMBOL_DEF IDENTIFIER '(' IDENTIFIER ')' ';' { $$ = create_symbol_element($2, $4); if ($$ == NULL) YYABORT; }
	| SYMBOL_DEF '[' ';' { $$ = create_symbol_element("[", NULL); if ($$ == NULL) YYABORT; }
	| SYMBOL_DEF '[' '(' IDENTIFIER ')' ';' { $$ = create_symbol_element("[", $4); if ($$ == NULL) YYABORT; }
	| SYMBOL_DEF ']' ';' { $$ = create_symbol_element("]", NULL); if ($$ == NULL) YYABORT; }
	| SYMBOL_DEF ']' '(' IDENTIFIER ')' ';' { $$ = create_symbol_element("]", $4); if ($$ == NULL) YYABORT; }
	;

axiom_definition
	: AXIOM_DEF production_element_string ';' { $$ = new_symbol_production(NULL, $2); }
	;

diffusionrange_definition
	: DIFFUSIONRANGE_DEF ':' REALVALUE ';' { $$ = $3; }
	;

rule_definition
	: RULE_DEF IDENTIFIER '{' rule_components '}' { $$ = complete_rule($2, $4); }
	;

rule_components
	: rule_lhs ':' expr ARROW rule_rhs { $$ = new_rule_element(NULL, $1, $3, $5); }
	| rule_lhs ARROW rule_rhs { $$ = new_rule_element(NULL, $1, NULL, $3); }
	;

rule_lhs
	: lhs_element_string { $$ = create_lhs_descriptor(current_lsys, $1); }
	;

lhs_element_string
	: /* empty */ { $$ = (LHS_SYMBOL *) NULL; }
	| lhs_element_string lhs_element { $$ = add_lhs_symbol($1, $2); }
	;

lhs_element
	: IDENTIFIER  { $$ = create_lhs_symbol($1, NULL); }
	| IDENTIFIER '(' IDENTIFIER ')' { $$ = create_lhs_symbol($1, $3); }
	;

rule_rhs
	: production_element_string { $$ = new_symbol_production(NULL, $1); }
	;

production_element_string
	: /* empty */ { $$ = (PRODUCTION_ELEMENT *) NULL; }
	| production_element_string production_element { $$ = add_production_element($1, $2); }
	;

/*
 * brackets are allowed as symbols too, this implementation is hacky and
 * probably induces anomalies -- should be fixed
 */

production_element
	: IDENTIFIER { $$ = create_production_element($1, NULL); if ($$ == NULL) YYABORT; }
	| IDENTIFIER '(' transsys_initializer ')' { $$ = create_production_element($1, $3); if ($$ == NULL) YYABORT; }
	| '[' { $$ = create_production_element("[", NULL); if ($$ == NULL) YYABORT; }
	| ']' { $$ = create_production_element("]", NULL); if ($$ == NULL) YYABORT; }
	;

transsys_initializer
	: source_transsys_specifier assignment_list { $$ = add_raw_assignment($2, $1); }
	| assignment_list { $$ = $1; }
	;

source_transsys_specifier
	: TRANSSYS_DEF IDENTIFIER ':' { $$ = create_raw_assignment(NULL, $2, NULL); } 
	;

assignment_list
	: /* empty */ { $$ = (RAW_ASSIGNMENT *) NULL; }
	| assignment { $$ = $1; }
	| assignment_list ',' assignment { $$ = add_raw_assignment($1, $3); }
	;

assignment
	: IDENTIFIER '=' expr { $$ = create_raw_assignment($1, NULL, $3); }
	;

graphics_definition
	: GRAPHICS_DEF '{' symgraph_list '}' {}
	;

symgraph_list
	: /* empty */ {}
	| symgraph_list symgraph {}
	;

symgraph
	: IDENTIFIER '{' graphcmd_list '}' { add_graphics_to_symbol(current_lsys, $1, $3); }
	| '[' '{' graphcmd_list '}' { add_graphics_to_symbol(current_lsys, "[", $3); }
	| ']' '{' graphcmd_list '}' { add_graphics_to_symbol(current_lsys, "]", $3); }
	;

graphcmd_list
	: /* empty */ { $$ = (GRAPHICS_PRIMITIVE *) NULL; }
	| graphcmd_list graphcmd { $$ = add_graphics_primitive($1, $2); }
	;

graphcmd
	: MOVE '(' expr ')' ';' { $$ = new_graphics_primitive(GRAPHICS_MOVE, $3); }
	| PUSH '(' ')' ';' { $$ = new_graphics_primitive(GRAPHICS_PUSH, NULL); }
	| POP '(' ')' ';' { $$ = new_graphics_primitive(GRAPHICS_POP, NULL); }
	| TURN '(' expr ')' ';' { $$ = new_graphics_primitive(GRAPHICS_TURN, $3); }
	| ROLL '(' expr ')' ';' { $$ = new_graphics_primitive(GRAPHICS_ROLL, $3); }
	| BANK '(' expr ')' ';' { $$ = new_graphics_primitive(GRAPHICS_BANK, $3); }
	| SPHERE '(' expr ')' ';' { $$ = new_graphics_primitive(GRAPHICS_SPHERE, $3); }
	| CYLINDER '(' expr ',' expr ')' ';' { $$ = new_graphics_primitive(GRAPHICS_CYLINDER, $3, $5); }
	| BOX '(' expr ',' expr ',' expr ')' ';' { $$ = new_graphics_primitive(GRAPHICS_BOX, $3, $5, $7); }
	| COLOR '(' expr ',' expr ',' expr ')' ';' { $$ = new_graphics_primitive(GRAPHICS_COLOR, $3, $5, $7); }
	;

/*** original transsys spec ***/

transsys
	: TRANSSYS_DEF IDENTIFIER { setup_current_transsys($2); } '{' transsys_element_list '}' {}
	;

transsys_element_list
	: /* empty */
	| transsys_element_list factor_definition { add_factor_definition(current_transsys, $2); }
	| transsys_element_list gene_definition { add_gene_definition(current_transsys, $2); }
	| error ';' { yyerrok; }
	;

factor_definition
	: FACTOR_DEF IDENTIFIER { setup_current_factor($2); } '{' factordef_components '}' { $$ = current_factor; }
	;


factordef_components
	: /* empty */
	| factordef_components factordef_component
	;

factordef_component
	: DECAY_DEF ':' expr ';' { add_factordef_decay($3); }
	| DIFFUSIBILITY_DEF ':' expr ';' { add_factordef_diffusibility($3); }
	;

gene_definition
	: GENE_DEF IDENTIFIER '{' promoter_component product_component '}' { $$ = create_gene($2, $4, $5); }
	;

promoter_component
	: PROMOTER_DEF '{' promoter_statements '}' { $$ = $3; }
	;

product_component
	: PRODUCT_DEF '{' product_statements '}' { $$ = $3; }
	;

promoter_statements
	: promoter_statement { $$ = $1; }
	| promoter_statements promoter_statement { $$ = extend_promoter_list($1, $2); }
	;

promoter_statement
	: CONSTITUTIVE ':' expr ';' { $$ = create_promoter(ACT_CONSTITUTIVE, NULL, $3, NULL); }
	| factor_combination ':' ACTIVATE '(' expr ',' expr ')' ';' { $$ = create_promoter(ACT_ACTIVATE, $1, $5, $7); }
	| factor_combination ':' REPRESS '(' expr ',' expr ')' ';' { $$ = create_promoter(ACT_REPRESS, $1, $5, $7); }
	;

factor_combination
	: IDENTIFIER { $$ = extend_factor_combination(NULL, $1); if ($$ == NULL) YYABORT; }
	| factor_combination '+' IDENTIFIER { $$ = extend_factor_combination($1, $3); if ($$ == NULL) YYABORT; }
	;

product_statements
	: DEFAULT ':' IDENTIFIER ';' { $$ = find_factor_index(current_transsys, $3); }
	;

expr
	: expr LOGICAL_OR and_expr { $$ = new_expression_node(NT_LOGICAL_OR, $1, $3); }
	| and_expr { $$ = $1; }
	;

and_expr
	: and_expr LOGICAL_AND not_expr { $$ = new_expression_node(NT_LOGICAL_AND, $1, $3); }
	| not_expr { $$ = $1; }
	;

not_expr
	: '!' not_expr { $$ = new_expression_node(NT_NOT, $2); }
	| cmp_expr { $$ = $1; }
	;

cmp_expr
	: cmp_expr '<' arithmetic_expr { $$ = new_expression_node(NT_LOWER, $1, $3); }
	| cmp_expr '>' arithmetic_expr { $$ = new_expression_node(NT_GREATER, $1, $3); }
	| cmp_expr LOWER_EQUAL arithmetic_expr { $$ = new_expression_node(NT_LOWER_EQUAL, $1, $3); }
	| cmp_expr GREATER_EQUAL arithmetic_expr { $$ = new_expression_node(NT_GREATER_EQUAL, $1, $3); }
	| cmp_expr EQUAL arithmetic_expr { $$ = new_expression_node(NT_EQUAL, $1, $3); }
	| cmp_expr UNEQUAL arithmetic_expr { $$ = new_expression_node(NT_UNEQUAL, $1, $3); }
	| arithmetic_expr { $$ = $1; }
	;

arithmetic_expr
	: arithmetic_expr '+' term { $$ = new_expression_node(NT_ADD, $1, $3); }
	| arithmetic_expr '-' term { $$ = new_expression_node(NT_SUBTRACT, $1, $3); }
	| term { $$ = $1; } 
	;

term
	: term '*' value { $$ = new_expression_node(NT_MULT, $1, $3); }
	| term '/' value { $$ = new_expression_node(NT_DIV, $1, $3); }
	| value { $$ = $1; }
	;

value
	: REALVALUE { $$ = new_expression_node(NT_VALUE, $1); }
	| '(' expr ')' { $$ = $2; }
	| RANDOM '(' expr ',' expr ')' { $$ = new_expression_node(NT_RANDOM, $3, $5); }
	| GAUSS '(' expr ',' expr ')' { $$ = new_expression_node(NT_GAUSS, $3, $5); }
	| IDENTIFIER { $$ = new_expression_node(NT_RAW_IDENTIFIER, NULL, $1); }
	| IDENTIFIER '.' IDENTIFIER { $$ = new_expression_node(NT_RAW_IDENTIFIER,$1, $3); }
	;

%%


#ifdef TEST_PARSER

char *yyin_name = "stdin";
TRANSSYS *parsed_transsys = NULL;


int main(int argc, char *argv[])
{
  int retval;

  retval = yyparse();
  printf("\n*** yyparse returned %d ***\n", retval);
  return (0);
}

#endif

