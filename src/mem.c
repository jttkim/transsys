/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.11  2005/06/22 09:58:36  jtk
 * prevented unknown variables from resulting in core-dump eliciting parsing results
 *
 * Revision 1.10  2005/06/16 09:36:26  jtk
 * implemented rule statistics gathering
 *
 * Revision 1.9  2005/06/15 22:17:13  jtk
 * counting number of transsys programs in lsys (deprecating multiples)
 *
 * Revision 1.8  2005/05/17 12:11:30  jtk
 * contact graph works
 *
 * Revision 1.7  2005/05/16 12:02:10  jtk
 * in transition from distance matrices to contact graphs
 *
 * Revision 1.6  2005/04/05 10:12:39  jtk
 * made diffusion consistent (no oscillation due to overshooting), small fixes
 *
 * Revision 1.5  2005/03/31 10:14:05  jtk
 * partial implementation of diffusion
 *
 * Revision 1.4  2005/03/30 18:30:27  jtk
 * progressed transition to arrayred lsys strings
 * introduced lsys string distance matrices
 *
 * Revision 1.3  2005/03/30 09:51:02  jtk
 * fixed free_lsys_string to actually free the string itself
 *
 * Revision 1.2  2005/03/29 17:33:02  jtk
 * introduced arrayed lsys string, with symbol distance matrix.
 *
 * Revision 1.1.1.1  2005/03/08 17:12:02  jtk
 * new cvs after loss at INB
 *
 * Revision 1.3  2002/01/25 03:35:03  kim
 * Added gnuplot link functionality to transexpr, transscatter, improved
 *     PostScript stuff
 *
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "trconfig.h"
#include "transsys.h"

/* #define FREE_DEADBEEF */

void free_deadbeef(void *p)
{
  int *ip = p, i;

  i = *ip;
  *ip = 0xdeadbeef;
#ifdef FREE_DEADBEEF
#  undef free
#endif
  fprintf(stderr, "free(%p): %08x\n", p, i);
  free(p);
#ifdef FREE_DEADBEEF
#  define free(x) free_deadbeef(x)
#endif

}


/*
 * this is just for flushing out problems in the malloc
 * arena...
 */

void malloc_testloop(size_t n, const char *fname, int lineno)
{
  size_t i;
  void *p;

  fprintf(stderr, "malloc_testloop: %s, %d\n", fname, lineno);
  for (i = 1; i <= n; i++)
  {
    p = malloc(i);
    free(p);
  }
}


void dump_memory(void *p, size_t n, const char *fname, int lineno)
{
  size_t i;
  unsigned char *cp = p;

  fprintf(stderr, "dump_memory: %s, %d:", fname, lineno);
  for (i = 0; i < n; i++)
  {
    fprintf(stderr, " %02x", *cp++);
  }
  fprintf(stderr, "\n");
}


void free_integer_array(INTEGER_ARRAY *ia)
{
  if (ia != NULL)
  {
    if (ia->array != NULL)
    {
      free(ia->array);
    }
    free(ia);
  }
}


INTEGER_ARRAY *extend_integer_array(INTEGER_ARRAY *ia, int v)
{
  int *new_array;

  if (ia == NULL)
  {
    ia = (INTEGER_ARRAY *) malloc(sizeof(INTEGER_ARRAY));
    if (ia == NULL)
    {
      fprintf(stderr, "extend_integer_array: failed to create integer array\n");
      return (NULL);
    }
    ia->length = 0;
    ia->array = NULL;
  }
  if (ia->length > 0)
  {
    new_array = (int *) realloc(ia->array, (ia->length + 1) * sizeof(int));
  }
  else
  {
    new_array = (int *) malloc((ia->length + 1) * sizeof(int));
  }
  if (new_array == NULL)
  {
    fprintf(stderr, "extend_integer_array: failed to extend integer array\n");
    return (ia);
  }
  ia->array = new_array;
  ia->array[ia->length] = v;
  ia->length++;
  return (ia);
}


void free_expression_tree(EXPRESSION_NODE *node)
{
  if (node == NULL)
  {
    /* allow attempts to free nonexisting nodes to support optional expressions in factors etc. */
    return;
  }
  switch(node->type)
  {
  case NT_VALUE:
  case NT_IDENTIFIER:
    break;
  case NT_RAW_IDENTIFIER:
    if (node->content.raw_identifier.factor_name)
      free(node->content.raw_identifier.factor_name);
    if (node->content.raw_identifier.transsys_label)
      free(node->content.raw_identifier.transsys_label);
    break;
  case NT_NOT:
  case NT_ATAN:
    free_expression_tree(node->content.argument[0]);
    free(node->content.argument);
    break;
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
    free_expression_tree(node->content.argument[0]);
    free_expression_tree(node->content.argument[1]);
    free(node->content.argument);
    break;
  default:
    fprintf(stderr, "free_expression_tree: unknown expression type %d\n", (int) node->type);
    break;
  }
  free(node);
}


EXPRESSION_NODE *new_expression_node(EXPR_NODE_TYPE type, ...)
{
  char *identifier_name;
  EXPRESSION_NODE *node;
  va_list arglist;

  node = (EXPRESSION_NODE *) malloc(sizeof(EXPRESSION_NODE));
  if (node == NULL)
    return (NULL);
  node->type = type;
  node->content.argument = NULL;
  va_start(arglist, type);
  switch(type)
  {
  case NT_VALUE:
    node->content.value = va_arg(arglist, double);
    break;
  case NT_IDENTIFIER:
    node->content.identifier.lhs_symbol_index = NO_INDEX; /*preliminary */
    node->content.identifier.factor_index = va_arg(arglist, int);
    break;
  case NT_RAW_IDENTIFIER:
    identifier_name = va_arg(arglist, char *);
    if (identifier_name)
    {
      node->content.raw_identifier.transsys_label = (char *) malloc((strlen(identifier_name) + 1) * sizeof(char));
      if (node->content.raw_identifier.transsys_label == NULL)
      {
	free(node);
	return (NULL);
      }
      strcpy(node->content.raw_identifier.transsys_label, identifier_name);
    }
    else
    {
      node->content.raw_identifier.transsys_label = NULL;
    }
    identifier_name = va_arg(arglist, char *);
    node->content.raw_identifier.factor_name = (char *) malloc((strlen(identifier_name) + 1) * sizeof(char));
    if (node->content.raw_identifier.factor_name == NULL)
    {
      free(node->content.raw_identifier.transsys_label);
      free(node);
      return (NULL);
    }
    strcpy(node->content.raw_identifier.factor_name, identifier_name);
    break;
  case NT_NOT:
  case NT_ATAN:
    node->content.argument = (EXPRESSION_NODE **) malloc(sizeof(EXPRESSION_NODE *));
    if (node->content.argument == NULL)
    {
      free(node);
      return (NULL);
    }
    node->content.argument[0] = va_arg(arglist, EXPRESSION_NODE *);
    break;
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
    node->content.argument = (EXPRESSION_NODE **) malloc(2 * sizeof(EXPRESSION_NODE *));
    if (node->content.argument == NULL)
    {
      free(node);
      return (NULL);
    }
    node->content.argument[0] = va_arg(arglist, EXPRESSION_NODE *);
    node->content.argument[1] = va_arg(arglist, EXPRESSION_NODE *);
    break;
  default:
    fprintf(stderr, "new_expression_node: unknown expression type %d\n", (int) type);
    break;
  }
  va_end(arglist);
  return (node);
}


void free_expression_context_entry_components(EXPRESSION_CONTEXT_ENTRY *e)
{
}


void free_expression_context_entry_list(EXPRESSION_CONTEXT_ENTRY *e)
{
  EXPRESSION_CONTEXT_ENTRY *e1;

  while (e)
  {
    e1 = e;
    e = e->next;
    free_expression_context_entry_components(e);
    free(e1);
  }
}


EXPRESSION_CONTEXT_ENTRY *new_expression_context_entry(const TRANSSYS *transsys, const char *label)
{
  EXPRESSION_CONTEXT_ENTRY *e;

  e = (EXPRESSION_CONTEXT_ENTRY *) malloc(sizeof(EXPRESSION_CONTEXT_ENTRY));
  if (e == NULL)
    return (NULL);
  e->next = NULL;
  e->transsys = transsys;
  strncpy(e->label, label, IDENTIFIER_MAX);
  e->label[IDENTIFIER_MAX - 1] = '\0';
  return (e);
}


static void free_promoter_components(PROMOTER_ELEMENT *a)
{
  if (a->num_binding_factors)
  {
    free(a->factor_index);
  }
  if (a->type == PROMOTERELEMENT_CONSTITUTIVE)
  {
    free_expression_tree(a->expr1);
  }
  else
  {
    free_expression_tree(a->expr1);
    free_expression_tree(a->expr2);
  }
}


void free_promoter_list(PROMOTER_ELEMENT *alist)
{
  PROMOTER_ELEMENT *a;

  while (alist)
  {
    a = alist;
    alist = alist->next;
    free_promoter_components(a);
    free(a);
  }
}


PROMOTER_ELEMENT *new_promoter_element(PROMOTERELEMENT_TYPE type, int num_binding_factors, int *factors, EXPRESSION_NODE *expr1, EXPRESSION_NODE *expr2)
{
  PROMOTER_ELEMENT *a;

  a = (PROMOTER_ELEMENT *) malloc(sizeof(PROMOTER_ELEMENT));
  if (a == NULL)
    return (NULL);
  a->next = NULL;
  a->type = type;
  a->num_binding_factors = num_binding_factors;
  if (num_binding_factors)
  {
    a->factor_index = factors;
  }
  else
  {
    a->factor_index = NULL;
  }
  if (type == PROMOTERELEMENT_CONSTITUTIVE)
  {
    a->expr1 = expr1;
  }
  else
  {
    a->expr1 = expr1;
    a->expr2 = expr2;
  }
  return (a);
}


static void free_factor_components(FACTOR_ELEMENT *fe)
{
  free_expression_tree(fe->decay_expression);
  free_expression_tree(fe->diffusibility_expression);
  free_expression_tree(fe->synthesis_expression);
  if (fe->num_producing_genes)
  {
    free(fe->gene_index);
  }
}


static void free_factor_list(FACTOR_ELEMENT *flist)
{
  FACTOR_ELEMENT *fe;

  while (flist)
  {
    fe = flist;
    flist = flist->next;
    free_factor_components(fe);
    free(fe);
  }
}


FACTOR_ELEMENT *new_factor_element(const char *name, EXPRESSION_NODE *decay_expression, EXPRESSION_NODE *diffusibility_expression, EXPRESSION_NODE *synthesis_expression)
{
  FACTOR_ELEMENT *fe;

  fe = (FACTOR_ELEMENT *) malloc(sizeof(FACTOR_ELEMENT));
  if (fe == NULL)
    return (NULL);
  fe->next = NULL;
  fe->index = NO_INDEX;
  fe->num_producing_genes = 0;
  fe->gene_index = NULL;
  strncpy(fe->name, name, IDENTIFIER_MAX);
  fe->name[IDENTIFIER_MAX - 1] = '\0';
  fe->decay_expression = decay_expression;
  fe->diffusibility_expression = diffusibility_expression;
  fe->synthesis_expression = synthesis_expression;
  return (fe);
}


static void free_gene_components(GENE_ELEMENT *ge)
{
  free_promoter_list(ge->promoter_list);
}


static void free_gene_list(GENE_ELEMENT *glist)
{
  GENE_ELEMENT *ge;

  while (glist)
  {
    ge = glist;
    glist = glist->next;
    free_gene_components(ge);
    free(ge);
  }
}


GENE_ELEMENT *new_gene_element(const char *name, PROMOTER_ELEMENT *promoter_list, int product_index)
{
  GENE_ELEMENT *ge;

  ge = (GENE_ELEMENT *) malloc(sizeof(GENE_ELEMENT));
  if (ge == NULL)
    return (NULL);
  ge->next = NULL;
  ge->index = NO_INDEX;
  strncpy(ge->name, name, IDENTIFIER_MAX);
  ge->name[IDENTIFIER_MAX - 1] = '\0';
  ge->promoter_list = promoter_list;
  ge->product_index = product_index;
  return (ge);
}


void free_transsys_components(TRANSSYS *tr)
{
  int i;

  if (tr->arrayed)
  {
    for (i = 0; i < tr->num_factors; i++)
      free_factor_components(tr->factor_list + i);
    if (tr->factor_list)
      free(tr->factor_list);
    for (i = 0; i < tr->num_genes; i++)
      free_gene_components(tr->gene_list + i);
    if (tr->gene_list)
      free(tr->gene_list);
  }
  else
  {
    free_factor_list(tr->factor_list);
    free_gene_list(tr->gene_list);
  }
}


void free_transsys_list(TRANSSYS *tr)
{
  TRANSSYS *tr1;

  while (tr)
  {
    free_transsys_components(tr);
    tr1 = tr;
    tr = tr->next;
    free(tr1);
  }
}


TRANSSYS *new_transsys(const char *name)
{
  TRANSSYS *tr;

  tr = (TRANSSYS *) malloc(sizeof(TRANSSYS));
  if (tr == NULL)
  {
    return (NULL);
  }
  tr->next = NULL;
  tr->arrayed = 0;
  tr->num_factors = 0;
  tr->num_genes = 0;
  tr->factor_list = NULL;
  tr->gene_list = NULL;
  strncpy(tr->name, name, IDENTIFIER_MAX);
  tr->name[IDENTIFIER_MAX - 1] = '\0';
  return (tr);
}


static int cmp_index_factor(const void *p1, const void *p2)
{
  FACTOR_ELEMENT *f1 = (FACTOR_ELEMENT *) p1,  *f2 = (FACTOR_ELEMENT *) p2;

  if (f1->index < f2->index)
    return (-1);
  else if (f1->index > f2->index)
    return (1);
  else
    return (0);
}


static int cmp_index_gene(const void *p1, const void *p2)
{
  GENE_ELEMENT *g1 = (GENE_ELEMENT *) p1,  *g2 = (GENE_ELEMENT *) p2;

  if (g1->index < g2->index)
    return (-1);
  else if (g1->index > g2->index)
    return (1);
  else
    return (0);
}


int arrange_transsys_arrays(TRANSSYS *transsys)
{
  int num_factors, num_genes, i, j;
  GENE_ELEMENT *ge, *ge1, *ge_arr = NULL;
  FACTOR_ELEMENT *fe, *fe1, *fe_arr;
  INTEGER_ARRAY *ia;

  if (transsys->arrayed)
  {
    fprintf(stderr, "arrange_transsys_arrays: attempt to multiply arrange transsys \"%s\"\n", transsys->name);
    return (-1);
  }
  num_genes = 0;
  for (ge = transsys->gene_list; ge; ge = ge->next)
  {
    num_genes++;
  }
  num_factors = 0;
  for (fe = transsys->factor_list; fe; fe = fe->next)
  {
    num_factors++;
  }
  /* fprintf(stderr, "arrange_transsys_arrays: %d factors, %d genes\n", num_factors, num_genes); */
  if (num_genes > 0)
  {
    ge_arr = (GENE_ELEMENT *) malloc(num_genes * sizeof(GENE_ELEMENT));
    if (ge_arr == NULL)
    {
      return (-1);
    }
    ge = transsys->gene_list;
    for (i = 0; i < num_genes; i++)
    {
      ge_arr[i] = *ge;
      ge_arr[i].next = NULL;
      ge1 = ge;
      ge = ge->next;
      free(ge1);
    }
    qsort(ge_arr, num_genes, sizeof(GENE_ELEMENT), cmp_index_gene);
    for (i = 0; i < num_genes - 1; i++)
    {
      ge_arr[i].next = ge_arr + i + 1;
    }
    ge_arr[i].next = NULL;
    transsys->gene_list = ge_arr;
  }
  if (num_factors > 0)
  {
    fe_arr = (FACTOR_ELEMENT *) malloc(num_factors * sizeof(FACTOR_ELEMENT));
    if (fe_arr == NULL)
    {
      if (num_genes > 0)
        free(ge_arr);
      return (-1);
    }
    fe = transsys->factor_list;
    for (i = 0; i < num_factors; i++)
    {
      fe_arr[i] = *fe;
      fe_arr[i].next = NULL;
      fe1 = fe;
      fe = fe->next;
      free(fe1);
    }
    qsort(fe_arr, num_factors, sizeof(FACTOR_ELEMENT), cmp_index_factor);
    for (i = 0; i < num_factors - 1; i++)
    {
      fe_arr[i].next = fe_arr + i + 1;
    }
    fe_arr[i].next = NULL;
    transsys->factor_list = fe_arr;
    for (i = 0; i < transsys->num_factors; i++)
    {
      ia = NULL;
      for (j = 0; j < transsys->num_genes; j++)
      {
	if (transsys->gene_list[j].product_index == i)
	{
	  ia = extend_integer_array(ia, j);
	}
      }
      if (ia)
      {
	transsys->factor_list[i].num_producing_genes = ia->length;
	transsys->factor_list[i].gene_index = ia->array;
	free(ia);
      }
/*
      else
      {
	fprintf(stderr, "warning: factor \"%s\" is not produced by any gene\n", transsys->factor_list[i].name);
      }
*/
    }
  }
  transsys->arrayed = 1;
  return (0);
}


void free_graphics_primitive_components(GRAPHICS_PRIMITIVE *gp)
{
  int i;

  for (i = 0; i < gp->num_arguments; i++)
  {
    free_expression_tree(gp->argument[i]);
  }
  if (gp->argument)
  {
    free(gp->argument);
  }
}


void free_graphics_primitive_list(GRAPHICS_PRIMITIVE *glist)
{
  GRAPHICS_PRIMITIVE *gp;

  /* fprintf(stderr, "free_graphics_primitive_list: start\n"); */
  while (glist)
  {
    gp = glist;
    glist = glist->next;
    free_graphics_primitive_components(gp);
    free(gp);
  }
}


GRAPHICS_PRIMITIVE *new_graphics_primitive(GRAPHICS_PRIMITIVE_TYPE type, ...)
{
  GRAPHICS_PRIMITIVE *gp;
  va_list arglist;

  gp = (GRAPHICS_PRIMITIVE *) malloc(sizeof(GRAPHICS_PRIMITIVE));
  if (gp == NULL)
    return (NULL);
  gp->next = NULL;
  gp->type = type;
  gp->argument = NULL;
  va_start(arglist, type);
  switch(type)
  {
  case GRAPHICS_PUSH:
  case GRAPHICS_POP:
    gp->num_arguments = 0;
    break;
  case GRAPHICS_MOVE:
  case GRAPHICS_TURN:
  case GRAPHICS_ROLL:
  case GRAPHICS_BANK:
  case GRAPHICS_SPHERE:
    gp->argument = (EXPRESSION_NODE **) malloc(sizeof(EXPRESSION_NODE *));
    if (gp->argument == NULL)
    {
      free(gp);
      return (NULL);
    }
    gp->num_arguments = 1;
    gp->argument[0] = va_arg(arglist, EXPRESSION_NODE *);
    break;
  case GRAPHICS_CYLINDER:
    gp->argument = (EXPRESSION_NODE **) malloc(2 * sizeof(EXPRESSION_NODE *));
    if (gp->argument == NULL)
    {
      free(gp);
      return (NULL);
    }
    gp->num_arguments = 2;
    gp->argument[0] = va_arg(arglist, EXPRESSION_NODE *);
    gp->argument[1] = va_arg(arglist, EXPRESSION_NODE *);
    break;
  case GRAPHICS_BOX:
  case GRAPHICS_COLOR:
    gp->argument = (EXPRESSION_NODE **) malloc(3 * sizeof(EXPRESSION_NODE *));
    if (gp->argument == NULL)
    {
      free(gp);
      return (NULL);
    }
    gp->num_arguments = 3;
    gp->argument[0] = va_arg(arglist, EXPRESSION_NODE *);
    gp->argument[1] = va_arg(arglist, EXPRESSION_NODE *);
    gp->argument[2] = va_arg(arglist, EXPRESSION_NODE *);
    break;
  default:
    fprintf(stderr, "new_graphics_primitive: unknown type %d\n", (int) type);
    break;
  }
  return (gp);
}


void free_symbol_element_components(SYMBOL_ELEMENT *sym)
{
  int i;

  if (sym->arrayed)
  {
    if (sym->graphics_primitive_list)
    {
      for (i = 0; i < sym->num_graphics_primitives; i++)
      {
	free_graphics_primitive_components(sym->graphics_primitive_list + i);
      }
      free(sym->graphics_primitive_list);
    }
  }
  else
  {
    if (sym->graphics_primitive_list)
    {
      free_graphics_primitive_list(sym->graphics_primitive_list);
    }
  }
}


void free_symbol_element_list(SYMBOL_ELEMENT *sym)
{
  SYMBOL_ELEMENT *s;

  while (sym)
  {
    s = sym;
    sym = sym->next;
    free_symbol_element_components(s);
    free(s);
  }
}


SYMBOL_ELEMENT *new_symbol_element(const char *name, const TRANSSYS *transsys)
{
  SYMBOL_ELEMENT *s;

  s = (SYMBOL_ELEMENT *) malloc(sizeof(SYMBOL_ELEMENT));
  if (s == NULL)
    return (NULL);
  s->next = NULL;
  s->index = 0;
  s->arrayed = 0;
  s->transsys = transsys;
  s->num_graphics_primitives = 0;
  s->graphics_primitive_list = NULL;
  strncpy(s->name, name, IDENTIFIER_MAX);
  s->name[IDENTIFIER_MAX - 1] = '\0';
  return (s);
}


int arrange_symbol_element_arrays(SYMBOL_ELEMENT *se)
{
  GRAPHICS_PRIMITIVE *gp_arr = NULL, *gp, *gp1;
  int num_graphics_primitives, i;

  if (se->arrayed)
    return (0);
  for (num_graphics_primitives = 0, gp = se->graphics_primitive_list; gp; gp = gp->next)
    num_graphics_primitives++;
  if (num_graphics_primitives > 0)
  {
    gp_arr = (GRAPHICS_PRIMITIVE *) malloc(num_graphics_primitives * sizeof(GRAPHICS_PRIMITIVE));
    if (gp_arr == NULL)
      return (-1);
    i = num_graphics_primitives;
    gp = se->graphics_primitive_list;
    while (gp)
    {
      gp_arr[--i] = *gp;
      gp_arr[i].next = NULL;
      gp1 = gp;
      gp = gp->next;
      free(gp1);
    }
    se->num_graphics_primitives = num_graphics_primitives;
    se->graphics_primitive_list = gp_arr;
    for (i = 1; i < num_graphics_primitives; i++)
      se->graphics_primitive_list[i - 1].next = se->graphics_primitive_list + i;
    se->graphics_primitive_list[i - 1].next = NULL;
  }
  se->arrayed = 1;
  return (0);
}


void free_lhs_symbol_components(LHS_SYMBOL *ls)
{
}


void free_lhs_symbol_list(LHS_SYMBOL *ls_list)
{
  LHS_SYMBOL *l;

  while (ls_list)
  {
    l = ls_list;
    ls_list = ls_list->next;
    free_lhs_symbol_components(l);
    free(l);
  }
}


LHS_SYMBOL *new_lhs_symbol(const char *transsys_label, const TRANSSYS *transsys, int symbol_index)
{
  LHS_SYMBOL *l;

  l = (LHS_SYMBOL *) malloc(sizeof(LHS_SYMBOL));
  if (l == NULL)
    return (NULL);
  l->next = NULL;
  l->index = NO_INDEX;
  l->symbol_index = symbol_index;
  if (transsys_label)
  {
    strncpy(l->transsys_label, transsys_label, IDENTIFIER_MAX);
    l->transsys_label[IDENTIFIER_MAX - 1] = '\0';
  }
  else
    l->transsys_label[0] = '\0';
  l->transsys = transsys;
  return (l);
}


void free_assignment_components(ASSIGNMENT *a)
{
  free_expression_tree(a->value);
}


void free_assignment_list(ASSIGNMENT *alist)
{
  ASSIGNMENT *a;

  while (alist)
  {
    a = alist;
    alist = alist->next;
    free_assignment_components(a);
    free(a);
  }
}


void free_lhs_descriptor_components(LHS_DESCRIPTOR *l)
{
  int i;

  if (l->arrayed)
  {
    for (i = 0; i < l->num_symbols; i++)
      free_lhs_symbol_components(l->symbol_list + i);
    free(l->symbol_list);
  }
  else
    free_lhs_symbol_list(l->symbol_list);
}


LHS_DESCRIPTOR *new_lhs_descriptor(LHS_SYMBOL *symbol_list)
{
  LHS_DESCRIPTOR *l;

  l = (LHS_DESCRIPTOR *) malloc(sizeof(LHS_DESCRIPTOR));
  if (l == NULL)
  {
    return (NULL);
  }
  l->arrayed = 0;
  l->num_symbols = 0;
  l->symbol_list = symbol_list;
  return (l);
}


int arrange_lhs_descriptor_arrays(LHS_DESCRIPTOR *l)
{
  int num_symbols, i;
  LHS_SYMBOL *s_arr, *s, *s1;

  if (l->arrayed)
  {
    return (0);
  }
  for (num_symbols = 0, s = l->symbol_list; s; s = s->next)
  {
    num_symbols++;
  }
  if (num_symbols > 0)
  {
    s_arr = (LHS_SYMBOL *) malloc(num_symbols * sizeof(LHS_SYMBOL));
    if (s_arr == NULL)
      return (-1);
    i = num_symbols;
    s = l->symbol_list;
    while (s)
    {
      s_arr[--i] = *s;
      s_arr[i].next = NULL;
      s1 = s;
      s = s->next;
      free(s1);
    }
    l->symbol_list = s_arr;
  }
  l->num_symbols = num_symbols;
  if (num_symbols > 0)
  {
    for (i = 1; i < num_symbols; i++)
    {
      l->symbol_list[i - 1].next = l->symbol_list + i;
    }
    l->symbol_list[i - 1].next = NULL;
  }
  l->arrayed = 1;
  return (0);
}


ASSIGNMENT *new_assignment(const TRANSSYS *target_transsys, int factor_index, EXPRESSION_NODE *value)
{
  ASSIGNMENT *a;

  a = (ASSIGNMENT *) malloc(sizeof(ASSIGNMENT));
  if (a == NULL)
    return (NULL);
  a->next = NULL;
  a->target_transsys = target_transsys;
  a->factor_index = factor_index;
  a->value = value;
  return (a);
}


void free_production_element_components(PRODUCTION_ELEMENT *sp)
{
  int i;

  if (sp->arrayed)
  {
    for (i = 0; i < sp->num_assignments; i++)
      free_assignment_components(sp->assignment_list + i);
    free(sp->assignment_list);
  }
  else
  {
    free_assignment_list(sp->assignment_list);
  }
}


void free_production_element_list(PRODUCTION_ELEMENT *slist)
{
  PRODUCTION_ELEMENT *sp;

  while (slist)
  {
    sp = slist;
    slist = slist->next;
    free_production_element_components(sp);
    free(sp);
  }
}


PRODUCTION_ELEMENT *new_production_element(int symbol_index, int template_lhs_symbol_index, ASSIGNMENT *assignment_list)
{
  PRODUCTION_ELEMENT *sp;

  sp = (PRODUCTION_ELEMENT *) malloc(sizeof(PRODUCTION_ELEMENT));
  if (sp == NULL)
  {
    return (NULL);
  }
  sp->next = NULL;
  sp->arrayed = 0;
  sp->symbol_index = symbol_index;
  sp->template_lhs_symbol_index = template_lhs_symbol_index;
  sp->assignment_list = assignment_list;
  return (sp);
}


int arrange_production_element_arrays(PRODUCTION_ELEMENT *sp)
{
  ASSIGNMENT *a_arr, *a, *a1;
  int num_assignments, i;

  if (sp->arrayed)
    return (0);
  for (num_assignments = 0, a = sp->assignment_list; a; a = a->next)
    num_assignments++;
  if (num_assignments > 0)
  {
    a_arr = (ASSIGNMENT *) malloc(num_assignments * sizeof(ASSIGNMENT));
    if (a_arr == NULL)
      return (-1);
    i = num_assignments;
    a = sp->assignment_list;
    while (a)
    {
      a_arr[--i] = *a;
      a_arr[i].next = NULL;
      a1 = a;
      a = a->next;
      free(a1);
    }
    sp->assignment_list = a_arr;
  }
  sp->num_assignments = num_assignments;
  if (num_assignments > 0)
  {
    for (i = 1; i < num_assignments; i++)
      sp->assignment_list[i - 1].next = sp->assignment_list + i;
    sp->assignment_list[i - 1].next = NULL;
  }
  sp->arrayed = 1;
  return (0);
}


void free_symbol_production_components(SYMBOL_PRODUCTION *sp)
{
  int i;

  if (sp->arrayed)
  {
    for (i = 0; i < sp->num_production_elements; i++)
      free_production_element_components(sp->production_list + i);
    free(sp->production_list);
  }
  else
    free_production_element_list(sp->production_list);
}


SYMBOL_PRODUCTION *new_symbol_production(PRODUCTION_ELEMENT *production_list)
{
  SYMBOL_PRODUCTION *sp;

  sp = (SYMBOL_PRODUCTION *) malloc(sizeof(SYMBOL_PRODUCTION));
  if (sp == NULL)
    return (NULL);
  sp->arrayed = 0;
  sp->num_production_elements = 0;
  sp->production_list = production_list;
  return (sp);
}


int arrange_symbol_production_arrays(SYMBOL_PRODUCTION *sp)
{
  PRODUCTION_ELEMENT *pe_arr, *pe, *pe1;
  int num_production_elements, i;

  for (num_production_elements = 0, pe = sp->production_list; pe; pe = pe->next)
    num_production_elements++;
  sp->num_production_elements = num_production_elements;
  if (num_production_elements > 0)
  {
    pe_arr = (PRODUCTION_ELEMENT *) malloc(num_production_elements * sizeof(PRODUCTION_ELEMENT));
    if (pe_arr == NULL)
      return (-1);
    pe = sp->production_list;
    i = num_production_elements;
    while (pe)
    {
      pe_arr[--i] = *pe;
      pe_arr[i].next = NULL;
      pe1 = pe;
      pe = pe->next;
      free(pe1);
    }
    sp->production_list = pe_arr;
  }
  sp->num_production_elements = num_production_elements;
  if (num_production_elements > 0)
  {
    for (i = 1; i < num_production_elements; i++)
      sp->production_list[i - 1].next = sp->production_list + i;
    sp->production_list[i - 1].next = NULL;
  }
  sp->arrayed = 1;
  return (0);
}


void free_rule_element_components(RULE_ELEMENT *r)
{
  if (r->condition)
  {
    free_expression_tree(r->condition);
  }
  if (r->lhs)
  {
    free_lhs_descriptor_components(r->lhs);
    free(r->lhs);
  }
  if (r->rhs)
  {
    free_symbol_production_components(r->rhs);
    free(r->rhs);
  }
}


void free_rule_element_list(RULE_ELEMENT *rlist)
{
  RULE_ELEMENT *r;

  while (rlist)
  {
    r = rlist;
    rlist = rlist->next;
    free_rule_element_components(r);
    free(r);
  }
}


RULE_ELEMENT *new_rule_element(const char *name, LHS_DESCRIPTOR *lhs, EXPRESSION_NODE *condition, SYMBOL_PRODUCTION *rhs)
{
  RULE_ELEMENT *r;

  r = (RULE_ELEMENT *) malloc(sizeof(RULE_ELEMENT));
  if (r == NULL)
  {
    return (NULL);
  }
  r->next = NULL;
  if (name)
  {
    strncpy(r->name, name, IDENTIFIER_MAX);
    r->name[IDENTIFIER_MAX - 1] = '\0';
  }
  else
  {
    r->name[0] = '\0';
  }
  r->lhs = lhs;
  r->condition = condition;
  r->rhs = rhs;
  return (r);
}


static void free_lsys_components(LSYS *lsys)
{
  int i;

  if (lsys->arrayed)
  {
    if (lsys->num_rules > 0)
    {
      for (i = 0; i < lsys->num_rules; i++)
      {
	free_rule_element_components(lsys->rule_list + i);
      }
      free(lsys->rule_list);
    }
    if (lsys->num_symbols > 0)
    {
      for (i = 0; i < lsys->num_symbols; i++)
      {
	free_symbol_element_components(lsys->symbol_list + i);
      }
      free(lsys->symbol_list);
    }
    if (lsys->num_transsys > 0)
    {
      free(lsys->transsys_list);
    }
  }
  else
  {
    free_rule_element_list(lsys->rule_list);
    free_symbol_element_list(lsys->symbol_list);
    if (lsys->transsys_list)
    {
      free(lsys->transsys_list);
    }
  }
  if (lsys->axiom)
  {
    free_symbol_production_components(lsys->axiom);
    free(lsys->axiom);
  }
}


void free_lsys_list(LSYS *ls)
{
  LSYS *ls1;

  while (ls)
  {
    free_lsys_components(ls);
    ls1 = ls;
    ls = ls->next;
    free(ls1);
  }
}


/*
 * Returns a NULL terminated array of pointers to all transsys instances in lsys.
 *
 * Returns NULL if malloc for the array fails or if lsys is not arrayed (error).
 *
 * If successful, an array is returned. The array may consist of a NULL pointer
 * only, and must be freed by caller after use.
 */
static const TRANSSYS **lsys_transsys_list(const LSYS *lsys)
{
  const TRANSSYS **tl = NULL, **tl1;
  int n, i, j;

  if (!lsys->arrayed)
  {
    fprintf(stderr, "lsys_transsys_list: lsys not arrayed\n");
    return (NULL);
  }
  tl = (const TRANSSYS **) malloc(sizeof(TRANSSYS *));
  if (tl == NULL)
  {
    fprintf(stderr, "lsys_transsys_list: malloc failed\n");
    return (NULL);
  }
  n = 0;
  for (i = 0; i < lsys->num_symbols; i++)
  {
    if (lsys->symbol_list[i].transsys != NULL)
    {
      for (j = 0; (j < n) && (tl[j] != lsys->symbol_list[i].transsys); j++)
	;
      if (j == n)
      {
	tl1 = (const TRANSSYS **) realloc(tl, (n + 1) * sizeof(TRANSSYS *));
	if (tl1 == NULL)
	{
	  fprintf(stderr, "lsys_transsys_list: realloc failed\n");
	  free(tl);
	  return (NULL);
	}
	tl = tl1;
	tl[n] = lsys->symbol_list[i].transsys;
	n++;
      }
    }
  }
  tl[n] = NULL;
  return (tl);
}


/*
 * free one lsys and all transsys programs it uses in its symbols.
 */
void free_lsys_with_transsys(LSYS *ls)
{
  const TRANSSYS **transsys_list;
  int i;

  transsys_list = lsys_transsys_list(ls);
  if (transsys_list != NULL)
  {
    for (i = 0; transsys_list[i] != NULL; i++)
    {
      /* cast from const TRANSSYS * to TRANSSYS * so we can alter and free transsys */
      TRANSSYS *transsys = (TRANSSYS *) transsys_list[i];

      if (transsys->next != NULL)
      {
	fprintf(stderr, "free_lsys_with_transsys: transsys \"%s\" has a successor (next)\n", transsys->name);
      }
      free_transsys_list(transsys);
    }
  }
  free(transsys_list);
  free_lsys_components(ls);
  free(ls);
}


LSYS *new_lsys(const char *name)
{
  LSYS *lsys;

  lsys = (LSYS *) malloc(sizeof(LSYS));
  if (lsys == NULL)
  {
    return (NULL);
  }
  lsys->next = NULL;
  lsys->arrayed = 0;
  lsys->num_symbols = 0;
  lsys->symbol_list = NULL;
  lsys->axiom = NULL;
  lsys->num_rules = 0;
  lsys->rule_list = NULL;
  lsys->num_transsys = 0;
  lsys->transsys_list = NULL;
  strncpy(lsys->name, name, IDENTIFIER_MAX);
  lsys->name[IDENTIFIER_MAX - 1] = '\0';
  return (lsys);
}


static int set_lsys_transsys_list(LSYS *lsys)
{
  int i, j;
  const TRANSSYS **tl;

  lsys->num_transsys = 0;
  for (i = 0; i < lsys->num_symbols; i++)
  {
    if (lsys->symbol_list[i].transsys != NULL)
    {
      for (j = 0; j < lsys->num_transsys; j++)
      {
	if (lsys->transsys_list[j] == lsys->symbol_list[i].transsys)
	{
	  break;
	}
      }
      if (j == lsys->num_transsys)
      {
	if (lsys->num_transsys == 0)
	{
	  lsys->transsys_list = (const TRANSSYS **) malloc(sizeof(TRANSSYS *));
	  if (lsys->transsys_list == NULL)
	  {
	    fprintf(stderr, "set_lsys_transsys_list: malloc failed\n");
	    return (-1);
	  }
	}
	else
	{
	  tl = (const TRANSSYS **) realloc(lsys->transsys_list, (lsys->num_transsys + 1) * sizeof(const TRANSSYS *));
	  if (tl == NULL)
	  {
	    fprintf(stderr, "set_lsys_transsys_list: realloc failed\n");
	    return(-1);
	  }
	  lsys->transsys_list = tl;
	}
	lsys->transsys_list[lsys->num_transsys++] = lsys->symbol_list[i].transsys;
      }
    }
  }
  if (lsys->num_transsys > 1)
  {
    fprintf(stderr, "warning: lsys %s comprises multiple (%d) transsys programs, this is deprecated\n", lsys->name, lsys->num_transsys);
  }
  return (0);
}


int arrange_lsys_arrays(LSYS *lsys)
{
  SYMBOL_ELEMENT *se_arr = NULL, *se, *se1;
  RULE_ELEMENT *re_arr = NULL, *re, *re1;
  int num_symbols, num_rules, i;

  if (lsys->arrayed)
  {
    return (0);
  }
  for (num_symbols = 0, se = lsys->symbol_list; se; se = se->next)
  {
    num_symbols++;
  }
  for (num_rules = 0, re = lsys->rule_list; re; re = re->next)
  {
    num_rules++;
  }
  if (num_symbols > 0)
  {
    se_arr = (SYMBOL_ELEMENT *) malloc(num_symbols * sizeof(SYMBOL_ELEMENT));
    if (se_arr == NULL)
    {
      return (-1);
    }
  }
  if (num_rules > 0)
  {
    re_arr = (RULE_ELEMENT *) malloc(num_rules * sizeof(RULE_ELEMENT));
    if (re_arr == NULL)
    {
      if (se_arr)
      {
	free(se_arr);
      }
      return (-1);
    }
  }
  lsys->num_symbols = num_symbols;
  if (num_symbols > 0)
  {
    i = num_symbols;
    se = lsys->symbol_list;
    while (se)
    {
      se_arr[--i] = *se;
      se_arr[i].next = NULL;
      se1 = se;
      se = se->next;
      free(se1);
    }
    lsys->symbol_list = se_arr;
    for (i = 1; i < num_symbols; i++)
    {
      lsys->symbol_list[i - 1].next = lsys->symbol_list + i;
    }
    lsys->symbol_list[i - 1].next = NULL;
  }
  lsys->num_rules = num_rules;
  if (num_rules > 0)
  {
    i = num_rules;
    re = lsys->rule_list;
    while (re)
    {
      re_arr[--i] = *re;
      re_arr[i].next = NULL;
      re1 = re;
      re = re->next;
      free(re1);
    }
    lsys->rule_list = re_arr;
    for (i = 1; i < num_rules; i++)
      lsys->rule_list[i - 1].next = lsys->rule_list + i;
    lsys->rule_list[i - 1].next = NULL;
  }
  lsys->arrayed = 1;
  if (set_lsys_transsys_list(lsys) != 0)
  {
    return (-1);
  }
  return (0);
}


void init_transsys_instance_components(TRANSSYS_INSTANCE *ti)
{
  ti->factor_concentration = NULL;
  ti->new_concentration = NULL;
  ti->transsys = NULL;
}


void free_transsys_instance_components(TRANSSYS_INSTANCE *ti)
{
  if (ti->factor_concentration)
  {
    free(ti->factor_concentration);
  }
  if (ti->new_concentration)
  {
    free(ti->new_concentration);
  }
  init_transsys_instance_components(ti);
}


/*
 * attempt to allocate components for transsys = NULL is legal
 * and results in returning 0 (success) without any action taking
 * place. This behaviour is used in symbol instance construction.
 */

int alloc_transsys_instance_components(TRANSSYS_INSTANCE *ti, const TRANSSYS *transsys)
{
  int i;

  free_transsys_instance_components(ti);
  if (transsys == NULL)
  {
    return (0);
  }
  if (transsys->num_factors > 0)
  {
    ti->factor_concentration = (double *) malloc(transsys->num_factors * sizeof(double));
    if (ti->factor_concentration == NULL)
    {
      return (-1);
    }
    ti->new_concentration = (double *) malloc(transsys->num_factors * sizeof(double));
    if (ti->new_concentration == NULL)
    {
      free(ti->factor_concentration);
      ti->factor_concentration = NULL;
      return (-1);
    }
    for (i = 0; i < transsys->num_factors; i++)
    {
      ti->factor_concentration[i] = 0.0;
      ti->new_concentration[i] = 0.0;
    }
  }
  else
  {
    ti->factor_concentration = NULL;
    ti->new_concentration = NULL;
  }
  ti->transsys = transsys;
  /* fprintf(stderr, "alloc_transsys_instance_components: transsys %s, %d factors, factor_concentration = %p, new_concentration = %p\n", transsys->name, transsys->num_factors, ti->factor_concentration, ti->new_concentration); */
  return (0);
}


int clone_transsys_instance(TRANSSYS_INSTANCE *ti, const TRANSSYS_INSTANCE *source)
{
  int i, return_value;

  return_value = alloc_transsys_instance_components(ti, source->transsys);
  if (return_value != 0)
  {
    return (return_value);
  }
  for (i = 0; i < ti->transsys->num_factors; i++)
  {
    ti->factor_concentration[i] = source->factor_concentration[i];
  }
  return (0);
}


void free_transsys_instance(TRANSSYS_INSTANCE *ti)
{
  free_transsys_instance_components(ti);
  free(ti);
}


TRANSSYS_INSTANCE *new_transsys_instance(const TRANSSYS *transsys)
{
  TRANSSYS_INSTANCE *ti = (TRANSSYS_INSTANCE *) malloc(sizeof(TRANSSYS_INSTANCE));

  if (ti == NULL)
  {
    return (NULL);
  }
  init_transsys_instance_components(ti);
  if (alloc_transsys_instance_components(ti, transsys) != 0)
  {
    free(ti);
    return (NULL);
  }
  return (ti);
}


void free_cell_components(CELL *cell)
{
  free_transsys_instance_components(&(cell->transsys_instance));
  if (cell->neighbor)
    free(cell->neighbor);
  if (cell->contact_weight)
    free(cell->contact_weight);
}


void free_cells(int num_cells, CELL *cell)
{
  int i;

  for (i = 0; i < num_cells; i++)
    free_cell_components(cell + i);
  free(cell);
}


int alloc_cell_components(CELL *c, const TRANSSYS *transsys)
{
  int return_value;

  return_value = alloc_transsys_instance_components(&(c->transsys_instance), transsys);
  if (return_value != 0)
    return (return_value);
  c->alive = 0;
  c->num_neighbors = 0;
  c->neighbor = NULL;
  c->contact_weight = NULL;
  return (0);
}


CELL *new_cells(int num_cells, const TRANSSYS *transsys)
{
  CELL *c;
  int i;

  c = (CELL *) malloc(num_cells * sizeof(CELL));
  if (c == NULL)
    return (NULL);
  for (i = 0; i < num_cells; i++)
  {
    init_transsys_instance_components(&(c[i].transsys_instance));
    c[i].existing = 0;
    c[i].alive = 0;
    c[i].num_neighbors = 0;
    c[i].neighbor = NULL;
    c[i].contact_weight = NULL;
  }
  for (i = 0; i < num_cells; i++)
  {
    if (alloc_cell_components(c + i, transsys))
      break;
  }
  if (i < num_cells)
  {
    free_cells(num_cells, c);
    return (NULL);
  }
  return (c);
}


void free_symbol_instance_components(SYMBOL_INSTANCE *si)
{
  if (si->num_contact_edges > 0)
  {
    free(si->contact_edge);
  }
  free_transsys_instance_components(&(si->transsys_instance));
}


void free_symbol_instance_list(SYMBOL_INSTANCE *slist)
{
  SYMBOL_INSTANCE *si;

  while (slist)
  {
    free_symbol_instance_components(slist);
    si = slist;
    slist = slist->next;
    free(si);
  }
}


int extend_symbol_contact_edge_array(SYMBOL_INSTANCE *si, LSYS_STRING_CONTACT_EDGE *edge)
{
  LSYS_STRING_CONTACT_EDGE **e;

  if (si->num_contact_edges > 0)
  {
    e = (LSYS_STRING_CONTACT_EDGE **) realloc(si->contact_edge, (si->num_contact_edges + 1) * sizeof(LSYS_STRING_CONTACT_EDGE *));
  }
  else
  {
    e = (LSYS_STRING_CONTACT_EDGE **) malloc(sizeof(LSYS_STRING_CONTACT_EDGE *));
  }
  if (e == NULL)
  {
    fprintf(stderr, "add_symbol_contact_edge: malloc / realloc failed\n");
    return (-1);
  }
  si->contact_edge = e;
  si->contact_edge[si->num_contact_edges] = edge;
  si->num_contact_edges++;
  return (0);
}


SYMBOL_INSTANCE *new_symbol_instance(const LSYS_STRING *lsys_string, int symbol_index)
{
  SYMBOL_INSTANCE *si;

  si = (SYMBOL_INSTANCE *) malloc(sizeof(SYMBOL_INSTANCE));
  if (si == NULL)
    return (NULL);
  si->next = NULL;
  si->lsys_string = lsys_string;
  si->symbol_index = symbol_index;
  si->rule_index = NO_INDEX;
  si->num_successors = 0;
  si->successor_index = NO_INDEX;
  si->successor_distance = -1;
  si->num_contact_edges = 0;
  si->contact_edge = NULL;
  init_transsys_instance_components(&(si->transsys_instance));
  if (alloc_transsys_instance_components(&(si->transsys_instance), lsys_string->lsys->symbol_list[symbol_index].transsys) != 0)
  {
    free(si);
    return (NULL);
  }
  return (si);
}


SYMBOL_INSTANCE *clone_symbol_instance(const SYMBOL_INSTANCE *source, const LSYS_STRING *lsys_string)
{
  SYMBOL_INSTANCE *si;
  int i;

  si = (SYMBOL_INSTANCE *) malloc(sizeof(SYMBOL_INSTANCE));
  if (si == NULL)
    return (NULL);
  si->next = NULL;
  si->lsys_string = lsys_string;
  si->symbol_index = source->symbol_index;
  si->num_successors = 0;
  si->successor_index = -NO_INDEX;
  si->successor_distance = -1;
  si->num_contact_edges = 0;    /* contact graph must be established by caller */
  si->contact_edge = NULL;
  init_transsys_instance_components(&(si->transsys_instance));
  if (alloc_transsys_instance_components(&(si->transsys_instance), source->transsys_instance.transsys) != 0)
  {
    free(si);
    return (NULL);
  }
  si->transsys_instance.transsys = source->transsys_instance.transsys;
  if (source->transsys_instance.transsys)
  {
    for (i = 0; i < source->transsys_instance.transsys->num_factors; i++)
      si->transsys_instance.factor_concentration[i] = source->transsys_instance.factor_concentration[i];
  }
  /* fprint_symbol_instance(stderr, source); */
  /* fprintf(stderr, "\n"); */
  /* fprint_symbol_instance(stderr, si); */
  /* fprintf(stderr, "\n"); */
  return (si);
}


void init_lsys_string_contact_graph_components(LSYS_STRING_CONTACT_GRAPH *g)
{
  g->array_size = 0;
  g->num_edges = 0;
  g->edge = NULL;
}


void free_lsys_string_contact_graph_components(LSYS_STRING_CONTACT_GRAPH *g)
{
  if (g->array_size > 0)
  {
    free(g->edge);
  }
  init_lsys_string_contact_graph_components(g);
}


void free_lsys_string_contact_graph(LSYS_STRING_CONTACT_GRAPH *g)
{
  free_lsys_string_contact_graph_components(g);
  free(g);
}


int alloc_lsys_string_contact_graph_components(LSYS_STRING_CONTACT_GRAPH *g, size_t array_size)
{
  if (g->array_size > 0)
  {
    fprintf(stderr, "alloc_lsys_string_contact_graph_components: array already present\n");
    return (-1);
  }
  if (array_size == 0)
  {
    free_lsys_string_contact_graph_components(g);
    return (0);
  }
  g->edge = (LSYS_STRING_CONTACT_EDGE *) malloc(array_size * sizeof(LSYS_STRING_CONTACT_EDGE));
  if (g->edge == NULL)
  {
    fprintf(stderr, "alloc_lsys_string_contact_graph_components: malloc failed\n");
    return (-1);
  }
  g->array_size = array_size;
  g->num_edges = 0;
  return (0);
}


LSYS_STRING_CONTACT_GRAPH *new_lsys_string_contact_graph(size_t num_edges)
{
  LSYS_STRING_CONTACT_GRAPH *g = (LSYS_STRING_CONTACT_GRAPH *) malloc(sizeof(LSYS_STRING_CONTACT_GRAPH));

  if (g == NULL)
  {
    fprintf(stderr, "new_lsys_string_contact_graph: malloc failed\n");
    return (NULL);
  }
  init_lsys_string_contact_graph_components(g);
  if (num_edges > 0)
  {
    if (alloc_lsys_string_contact_graph_components(g, num_edges) != 0)
    {
      fprintf(stderr, "new_lsys_string_contact_graph: alloc_lsys_string_graph_components failed\n");
      free(g);
      return (NULL);
    }
  }
  return (g);
}


int add_lsys_string_contact_edge(LSYS_STRING_CONTACT_GRAPH *g, int i1, int i2, int distance)
{
  if (g->num_edges >= g->array_size)
  {
    fprintf(stderr, "add_lsys_string_contact_edge: contact graph full\n");
    return (-1);
  }
  g->edge[g->num_edges].i1 = i1;
  g->edge[g->num_edges].i2 = i2;
  g->edge[g->num_edges].distance = distance;
  g->num_edges++;
  return (0);
}


int connect_lsys_string_symbols(LSYS_STRING *lstr, int i1, int i2, int distance)
{
  SYMBOL_INSTANCE *s1 = lstr->symbol + i1, *s2 = lstr->symbol + i2;
  LSYS_STRING_CONTACT_EDGE *edge;

  if (add_lsys_string_contact_edge(&(lstr->contact_graph), i1, i2, distance) != 0)
  {
    return (-1);
  }
  edge = lstr->contact_graph.edge + (lstr->contact_graph.num_edges - 1);
  if (extend_symbol_contact_edge_array(s1, edge) != 0)
  {
    fprintf(stderr, "connect_contacting_symbols: failed to add edge to symbol\n");
    return (-1);
  }
  if (extend_symbol_contact_edge_array(s2, edge) != 0)
  {
    fprintf(stderr, "connect_contacting_symbols: failed to add edge to symbol\n");
    return (-1);
  }
  s1->contact_edge[s1->num_contact_edges - 1] = edge;
  s2->contact_edge[s2->num_contact_edges - 1] = edge;
  return (0);
}


void free_lsys_string(LSYS_STRING *lstr)
{
  size_t i;

  free_lsys_string_contact_graph_components(&(lstr->contact_graph));
  if (lstr->arrayed)
  {
    for (i = 0; i < lstr->num_symbols; i++)
    {
      free_symbol_instance_components(lstr->symbol + i);
    }
    free(lstr->symbol);
  }
  else
  {
    free_symbol_instance_list(lstr->symbol);
  }
  free(lstr);
}


int arrange_lsys_string_arrays(LSYS_STRING *lstr)
{
  SYMBOL_INSTANCE *si_arr, *si, *si1;
  size_t i;

  /*
  fprintf(stderr, "arrange_lsys_string_arrays: start\n");
  verify_Watchdogs();
  fprintf(stderr, "arrange_lsys_string_arrays: start\n");
  */
  if (lstr->arrayed)
  {
    fprintf(stderr, "multiple array arrangement for LSYS_STRING of \"%s\"\n", lstr->lsys->name);
    return (-1);
  }
  lstr->num_symbols = 0;
  for (si = lstr->symbol; si; si = si->next)
  {
    lstr->num_symbols++;
  }
  if (lstr->num_symbols > 0)
  {
    si_arr = (SYMBOL_INSTANCE *) malloc(lstr->num_symbols * sizeof(SYMBOL_INSTANCE));
    if (si_arr == NULL)
    {
      return (-1);
    }
    si = lstr->symbol;
    for (i = 0; i < lstr->num_symbols; i++)
    {
      si_arr[i] = *si;
      si_arr[i].next = NULL;
      si1 = si;
      si = si->next;
      free(si1);
    }
    for (i = 1; i < lstr->num_symbols; i++)
    {
      si_arr[i - 1].next = si_arr + i;
    }
    lstr->symbol = si_arr;
  }
  lstr->arrayed = 1;
  /*
  fprintf(stderr, "arrange_lsys_string_arrays: end\n");
  verify_Watchdogs();
  fprintf(stderr, "arrange_lsys_string_arrays: end\n");
  */
  return (0);
}


LSYS_STRING *new_lsys_string(const LSYS *lsys)
{
  LSYS_STRING *lstr;

  lstr = (LSYS_STRING *) malloc(sizeof(LSYS_STRING));
  if (lstr == NULL)
  {
    return (NULL);
  }
  lstr->lsys = lsys;
  lstr->num_symbols = 0;
  lstr->arrayed = 0;
  lstr->symbol = NULL;
  init_lsys_string_contact_graph_components(&(lstr->contact_graph));
  return (lstr);
}

