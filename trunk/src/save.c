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
 * Revision 1.8  2005/05/16 21:03:27  jtk
 * contact graph implementation still buggy
 *
 * Revision 1.7  2005/05/16 12:02:10  jtk
 * in transition from distance matrices to contact graphs
 *
 * Revision 1.6  2005/04/04 21:30:07  jtk
 * differentiated fprint_lsys_string and fprint_lsys_string_distances
 *
 * Revision 1.5  2005/04/04 09:39:54  jtk
 * added lsys capabilities to transexpr, various small changes
 *
 * Revision 1.4  2005/03/31 16:07:36  jtk
 * finished (initial) implementation of lsys diffusion
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
 * Revision 1.3  2003/02/04 23:48:12  kim
 * hack-fixed bug regarding ! by parentheses, no real fix!!
 *
 * Revision 1.2  2002/05/28 08:52:50  kim
 * added some discrete network stuff to sources
 *
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#include <stdio.h>
#include <stdarg.h>

#include "trconfig.h"
#include "transsys.h"


static void fprint_expression_tree(FILE *f, const EXPRESSION_NODE *node, const FACTOR_ELEMENT *factor, const RULE_ELEMENT *rule);


static void fprint_binary_expression(FILE *f, const char *prefix, const char *infix, const char *postfix, const EXPRESSION_NODE *arg1, const EXPRESSION_NODE *arg2, const FACTOR_ELEMENT *factor, const RULE_ELEMENT *rule)
{
  fprintf(f, "%s", prefix);
  fprint_expression_tree(f, arg1, factor, rule);
  fprintf(f, "%s", infix);
  fprint_expression_tree(f, arg2, factor, rule);
  fprintf(f, "%s", postfix);
}


static void fprint_simple_identifier(FILE *f, const EXPRESSION_NODE *node, const FACTOR_ELEMENT *factor)
{
  if ((factor == NULL) || (node->content.identifier.factor_index < 0))
    fprintf(f, "<unknown factor #%d>", node->content.identifier.factor_index);
  else
    fprintf(f, "%s", factor[node->content.identifier.factor_index].name);
}


static void fprint_complex_identifier(FILE *f, const EXPRESSION_NODE *node, const RULE_ELEMENT *rule)
{
  const TRANSSYS *transsys;

  if ((rule == NULL) || (node->content.identifier.lhs_symbol_index < 0))
    fprintf(f, "<unknown lhs #%d>.", node->content.identifier.lhs_symbol_index);
  else
    fprintf(f, "%s.", rule->lhs->symbol_list[node->content.identifier.lhs_symbol_index].transsys_label);
  transsys = rule->lhs->symbol_list[node->content.identifier.lhs_symbol_index].transsys;
  if (transsys)
    fprint_simple_identifier(f, node, transsys->factor_list);
  else
    fprint_simple_identifier(f, node, NULL);
}


static void fprint_expression_tree(FILE *f, const EXPRESSION_NODE *node, const FACTOR_ELEMENT *factor, const RULE_ELEMENT *rule)
{
  switch(node->type)
  {
  case NT_VALUE:
    fprintf(f, "%g", node->content.value);
    break;
  case NT_IDENTIFIER:
    if (node->content.identifier.lhs_symbol_index >= 0)
      fprint_complex_identifier(f, node, rule);
    else
      fprint_simple_identifier(f, node, factor);
    break;
  case NT_RAW_IDENTIFIER:
    fprintf(f, "raw(");
    if (node->content.raw_identifier.transsys_label)
      fprintf(f, "%s.", node->content.raw_identifier.transsys_label);
    fprintf(f, "%s)", node->content.raw_identifier.factor_name);
    break;
  case NT_NOT:
    fprintf(f, "(!(");
    fprint_expression_tree(f, node->content.argument[0], factor, rule);
    fprintf(f, "))");
    break;
  case NT_LOGICAL_OR:
    fprint_binary_expression(f, "(", " || ", ")", node->content.argument[0], node->content.argument[1], factor, rule);
    break;
  case NT_LOGICAL_AND:
    fprint_binary_expression(f, "(", " && ", ")", node->content.argument[0], node->content.argument[1], factor, rule);
    break;
  case NT_LOWER:
    fprint_binary_expression(f, "(", " < ", ")", node->content.argument[0], node->content.argument[1], factor, rule);
    break;
  case NT_LOWER_EQUAL:
    fprint_binary_expression(f, "(", " <= ", ")", node->content.argument[0], node->content.argument[1], factor, rule);
    break;
  case NT_GREATER:
    fprint_binary_expression(f, "(", " > ", ")", node->content.argument[0], node->content.argument[1], factor, rule);
    break;
  case NT_GREATER_EQUAL:
    fprint_binary_expression(f, "(", " >= ", ")", node->content.argument[0], node->content.argument[1], factor, rule);
    break;
  case NT_EQUAL:
    fprint_binary_expression(f, "(", " == ", ")", node->content.argument[0], node->content.argument[1], factor, rule);
    break;
  case NT_UNEQUAL:
    fprint_binary_expression(f, "(", " != ", ")", node->content.argument[0], node->content.argument[1], factor, rule);
    break;
  case NT_ADD:
    fprint_binary_expression(f, "(", " + ", ")", node->content.argument[0], node->content.argument[1], factor, rule);
    break;
  case NT_SUBTRACT:
    fprint_binary_expression(f, "(", " - ", ")", node->content.argument[0], node->content.argument[1], factor, rule);
    break;
  case NT_MULT:
    fprint_binary_expression(f, "(", " * ", ")", node->content.argument[0], node->content.argument[1], factor, rule);
    break;
  case NT_DIV:
    fprint_binary_expression(f, "(", " / ", ")", node->content.argument[0], node->content.argument[1], factor, rule);
    break;
  case NT_RANDOM:
    fprint_binary_expression(f, "random(", ", ", ")", node->content.argument[0], node->content.argument[1], factor, rule);
    break;
  case NT_GAUSS:
    fprint_binary_expression(f, "gauss(", ", ", ")", node->content.argument[0], node->content.argument[1], factor, rule);
    break;
  default:
    fprintf(stderr, "fprint_expression_tree: unknown expression type %d\n", (int) node->type);
    break;
  }
}


static void fprint_indented(FILE *f, int indent_depth, char *format, ...)
{
  int i;
  va_list arglist;

  for (i = 0; i < indent_depth; i++)
    fputc(' ', f);
  va_start(arglist, format);
  vfprintf(f, format, arglist);
}


static void fprint_promoter(FILE *f, int indent_depth, const PROMOTER_ELEMENT *a, const FACTOR_ELEMENT *factor)
{
  int i;

  fprint_indented(f, indent_depth, "");
  switch(a->type)
  {
  case PROMOTERELEMENT_CONSTITUTIVE:
    fprintf(f, "constitutive: ");
    fprint_expression_tree(f, a->expr1, factor, NULL);
    break;
  case PROMOTERELEMENT_ACTIVATE:
    for (i = 0; i < a->num_binding_factors; i++)
    {
      if (i == 0)
	fprintf(f, "%s", factor[a->factor_index[i]].name);
      else
	fprintf(f, " + %s", factor[a->factor_index[i]].name);
    }
    fprintf(f, ": activate(");
    fprint_expression_tree(f, a->expr1, factor, NULL);
    fprintf(f, ", ");
    fprint_expression_tree(f, a->expr2, factor, NULL);
    fprintf(f, ")");
    break;
  case PROMOTERELEMENT_REPRESS:
    for (i = 0; i < a->num_binding_factors; i++)
    {
      if (i == 0)
	fprintf(f, "%s", factor[a->factor_index[i]].name);
      else
	fprintf(f, " + %s", factor[a->factor_index[i]].name);
    }
    fprintf(f, ": repress(");
    fprint_expression_tree(f, a->expr1, factor, NULL);
    fprintf(f, ", ");
    fprint_expression_tree(f, a->expr2, factor, NULL);
    fprintf(f, ")");
    break;
  default:
    fprintf(stderr, "fprint_promoter: unknown activation type %d\n", (int) a->type);
    break;
  }
  fprintf(f, ";\n");
}


static void fprint_factor(FILE *f, int indent_depth, const FACTOR_ELEMENT *factor, const FACTOR_ELEMENT *factor_array)
{
  fprint_indented(f, indent_depth, "factor %s\n", factor->name);
  fprint_indented(f, indent_depth, "{\n");
  fprint_indented(f, indent_depth + 2, "decay: ");
  fprint_expression_tree(f, factor->decay_expression, factor_array, NULL);
  fprintf(f, ";\n");
  fprint_indented(f, indent_depth + 2, "diffusibility: ");
  fprint_expression_tree(f, factor->diffusibility_expression, factor_array, NULL);
  fprintf(f, ";\n");
  fprint_indented(f, indent_depth, "}\n");
}


static void fprint_gene(FILE *f, int indent_depth, const GENE_ELEMENT *gene, const FACTOR_ELEMENT *factor)
{
  const PROMOTER_ELEMENT *a;

  fprint_indented(f, indent_depth, "gene %s\n", gene->name);
  fprint_indented(f, indent_depth, "{\n");
  fprint_indented(f, indent_depth + 2, "promoter\n");
  fprint_indented(f, indent_depth + 2, "{\n");
  for (a = gene->promoter_list; a; a = a->next)
    fprint_promoter(f, indent_depth + 4, a, factor);
  fprint_indented(f, indent_depth + 2, "}\n");
  fprint_indented(f, indent_depth + 2, "product\n");
  fprint_indented(f, indent_depth + 2, "{\n");
  if (gene->product_index < 0)
    fprint_indented(f, indent_depth + 4, "default: <unknown product %d>;\n", gene->product_index);
  else
    fprint_indented(f, indent_depth + 4, "default: %s;\n", factor[gene->product_index].name);
  fprint_indented(f, indent_depth + 2, "}\n");
  fprint_indented(f, indent_depth, "}\n");
}


void fprint_transsys(FILE *f, int indent_depth, const TRANSSYS *transsys)
{
  const FACTOR_ELEMENT *fe;
  const GENE_ELEMENT *ge;

  if (!transsys->arrayed)
  {
    fprintf(stderr, "fprint_transsys: printing non-arrayed transsys \"%s\" (may not work)\n", transsys->name);
    fprint_indented(f, indent_depth, "#non-arrayed transsys\n");
    /*
    fprintf(stderr, "fprint_transsys: printing non-arrayed transcription system, arraying first\n");
    if (arrange_transsys_arrays((TRANSSYS *) transsys) != 0)
    {
      return;
    }
    */
  }
  fprint_indented(f, indent_depth, "transsys %s\n", transsys->name);
  fprint_indented(f, indent_depth, "{\n");
  for (fe = transsys->factor_list; fe; fe = fe->next)
  {
    fprint_factor(f, indent_depth + 2, fe, transsys->factor_list);
  }
  for (ge = transsys->gene_list; ge; ge = ge->next)
  {
    fprintf(f, "\n");
    fprint_gene(f, indent_depth + 2, ge, transsys->factor_list);
  }
  fprint_indented(f, indent_depth, "}\n\n");
}


static void fprint_gene_as_discretenet_node(FILE *f, const TRANSSYS *transsys, int gene_index)
{
  const PROMOTER_ELEMENT *a;
  const GENE_ELEMENT *ge;
  int i;
  const char *glue = "";

  fprintf(f, "  node %s(", transsys->gene_list[gene_index].name);
  for (a = transsys->gene_list[gene_index].promoter_list; a; a = a->next)
  {
    switch (a->type)
    {
    case PROMOTERELEMENT_ACTIVATE:
    case PROMOTERELEMENT_REPRESS:
      for (i = 0; i < a->num_binding_factors; i++)
      {
	for (ge = transsys->gene_list; ge; ge = ge->next)
	{
	  if (ge->product_index == a->factor_index[i])
	  {
	    fprintf(f, "%s%s", glue, ge->name);
	    glue = ", ";
	  }
	}
      }
      break;
    default:
      fprintf(stderr, "fprint_gene_as_discretenet_node: unhandled promoter element %d\n", (int) a->type);
      break;
    }
  }
  fprintf(f, ") {}\n");
}


void fprint_transsys_as_discretenet(FILE *f, const TRANSSYS *transsys)
{
  int gene_index;

  if (!transsys->arrayed)
  {
    fprintf(stderr, "fprint_transsys_as_discretenet: non-arrayed transsys %s\n", transsys->name);
    return;
  }
  fprintf(f, "network %s\n", transsys->name);
  fprintf(f, "{\n");
  fprintf(f, "  {\n");
  fprintf(f, "    comment: \"discrete network with topology of transsys %s\"\n", transsys->name);
  fprintf(f, "  }\n");
  for (gene_index = 0; gene_index < transsys->num_genes; gene_index++)
  {
    fprint_gene_as_discretenet_node(f, transsys, gene_index);
  }
  fprintf(f, "}\n");
}


static void fprint_graphics_primitive(FILE *f, int indent_depth, const GRAPHICS_PRIMITIVE *gp, const FACTOR_ELEMENT *factor)
{
  if (gp == NULL)
    return;
  switch (gp->type)
  {
  case GRAPHICS_PUSH:
    fprint_indented(f, indent_depth, "push();\n");
    break;
  case GRAPHICS_POP:
    fprint_indented(f, indent_depth, "pop();\n");
    break;
  case GRAPHICS_MOVE:
    fprint_indented(f, indent_depth, "move(");
    fprint_expression_tree(f, gp->argument[0], factor, NULL);
    fprintf(f, ");\n");
    break;
  case GRAPHICS_TURN:
    fprint_indented(f, indent_depth, "turn(");
    fprint_expression_tree(f, gp->argument[0], factor, NULL);
    fprintf(f, ");\n");
    break;
  case GRAPHICS_ROLL:
    fprint_indented(f, indent_depth, "roll(");
    fprint_expression_tree(f, gp->argument[0], factor, NULL);
    fprintf(f, ");\n");
    break;
  case GRAPHICS_BANK:
    fprint_indented(f, indent_depth, "bank(");
    fprint_expression_tree(f, gp->argument[0], factor, NULL);
    fprintf(f, ");\n");
    break;
  case GRAPHICS_SPHERE:
    fprint_indented(f, indent_depth, "sphere(");
    fprint_expression_tree(f, gp->argument[0], factor, NULL);
    fprintf(f, ");\n");
    break;
  case GRAPHICS_CYLINDER:
    fprint_indented(f, indent_depth, "cylinder(");
    fprint_expression_tree(f, gp->argument[0], factor, NULL);
    fprintf(f, ", ");
    fprint_expression_tree(f, gp->argument[1], factor, NULL);
    fprintf(f, ");\n");
    break;
  case GRAPHICS_BOX:
    fprint_indented(f, indent_depth, "box(");
    fprint_expression_tree(f, gp->argument[0], factor, NULL);
    fprintf(f, ", ");
    fprint_expression_tree(f, gp->argument[1], factor, NULL);
    fprintf(f, ", ");
    fprint_expression_tree(f, gp->argument[2], factor, NULL);
    fprintf(f, ");\n");
    break;
  case GRAPHICS_COLOR:
    fprint_indented(f, indent_depth, "color(");
    fprint_expression_tree(f, gp->argument[0], factor, NULL);
    fprintf(f, ", ");
    fprint_expression_tree(f, gp->argument[1], factor, NULL);
    fprintf(f, ", ");
    fprint_expression_tree(f, gp->argument[2], factor, NULL);
    fprintf(f, ");\n");
    break;
  default:
    fprintf(stderr, "fprint_graphics_primitive: unknown primitive type %d\n", (int) gp->type);
    break;
  }
}


static void fprint_graphics_primitive_list(FILE *f, int indent_depth, const GRAPHICS_PRIMITIVE *gp_list, const FACTOR_ELEMENT *factor)
{
  const GRAPHICS_PRIMITIVE *gp;

  for (gp = gp_list; gp; gp = gp->next)
    fprint_graphics_primitive(f, indent_depth, gp, factor);
}


static void fprint_symbol(FILE *f, int indent_depth, const SYMBOL_ELEMENT *se)
{
  fprint_indented(f, indent_depth, "symbol %s", se->name);
  if (se->transsys)
    fprintf(f, "(%s)", se->transsys->name);
  fprintf(f, ";\n");
}


static void fprint_assignment(FILE *f, const ASSIGNMENT *a, const RULE_ELEMENT *rule)
{
  fprintf(f, "%s = ", a->target_transsys->factor_list[a->factor_index].name);
  fprint_expression_tree(f, a->value, NULL, rule);
}


static void fprint_production_element(FILE *f, const PRODUCTION_ELEMENT *sp, const SYMBOL_ELEMENT *symbol_list, const RULE_ELEMENT *rule)
{
  ASSIGNMENT *a;

  fprintf(f, "%s", symbol_list[sp->symbol_index].name);
  if (sp->assignment_list)
  {
    fprintf(f, "(");
    if (sp->template_lhs_symbol_index != NO_INDEX)
      fprintf(f, "transsys %s: ", rule->lhs->symbol_list[sp->template_lhs_symbol_index].transsys_label);
    for (a = sp->assignment_list; a; a = a->next)
    {
      fprint_assignment(f, a, rule);
      if (a->next)
	fprintf(f, ", ");
    }
    fprintf(f, ")");
  }
  else if (sp->template_lhs_symbol_index != NO_INDEX)
    fprintf(f, "(transsys %s: )", rule->lhs->symbol_list[sp->template_lhs_symbol_index].transsys_label);
}


static void fprint_production_element_list(FILE *f, const PRODUCTION_ELEMENT *sp_list, const SYMBOL_ELEMENT *symbol_list, const RULE_ELEMENT *rule)
{
  const PRODUCTION_ELEMENT *sp;

  for (sp = sp_list; sp; sp = sp->next)
  {
    fprint_production_element(f, sp, symbol_list, rule);
    if (sp->next)
      fprintf(f, " ");
  }

}


static void fprint_symbol_production(FILE *f, const SYMBOL_PRODUCTION *sp, const SYMBOL_ELEMENT *symbol_list, const RULE_ELEMENT *rule)
{
  fprint_production_element_list(f, sp->production_list, symbol_list, rule);
}


static void fprint_lhs_symbol(FILE *f, const LHS_SYMBOL *lhs_symbol, const SYMBOL_ELEMENT *symbol_list)
{
  fprintf(f, "%s(%s)", symbol_list[lhs_symbol->symbol_index].name, lhs_symbol->transsys_label);
  /* fprintf(f, "\n# lhs_symbol->index: %d\n", lhs_symbol->index); */
}


static void fprint_lhs(FILE *f, int indent_depth, const LHS_DESCRIPTOR *lhs, const SYMBOL_ELEMENT *symbol_list)
{
  const LHS_SYMBOL *lhs_symbol;

  fprint_indented(f, indent_depth, "");
  for (lhs_symbol = lhs->symbol_list; lhs_symbol; lhs_symbol = lhs_symbol->next)
  {
    fprint_lhs_symbol(f, lhs_symbol, symbol_list);
    if (lhs_symbol->next)
      fprintf(f, " ");
  }
}


/*
 * preliminary: use single (i.e. first) lhs symbol to obtain transsys...
 */

static void fprint_rule(FILE *f, int indent_depth, const RULE_ELEMENT *re, const SYMBOL_ELEMENT *symbol_list)
{
  fprint_indented(f, indent_depth, "rule %s\n", re->name);
  fprint_indented(f, indent_depth, "{\n");
  fprint_lhs(f, indent_depth + 2, re->lhs, symbol_list);
  if (re->condition)
  {
    fprintf(f, ": ");
    fprint_expression_tree(f, re->condition, NULL, re);
  }
  fprintf(f, "\n");
  fprint_indented(f, indent_depth + 2, "--> ");
  fprint_symbol_production(f, re->rhs, symbol_list, re);
  fprintf(f, "\n");
  fprint_indented(f, indent_depth, "}\n");
}


void fprint_lsys(FILE *f, int indent_depth, const LSYS *lsys)
{
  const SYMBOL_ELEMENT *se;
  const RULE_ELEMENT *re;

  if (!lsys->arrayed)
  {
    fprintf(stderr, "fprint_lsys: printing non-arrayed lsys \"%s\" (may  not work)\n", lsys->name);
    /* arrange_lsys_arrays((LSYS *) lsys); */
    fprint_indented(f, indent_depth, "# non arrayed lsys\n");
  }
  fprint_indented(f, indent_depth, "lsys %s\n", lsys->name);
  fprint_indented(f, indent_depth, "{\n");
  fprint_indented(f, indent_depth + 2, "diffusionrange: %d;\n\n", lsys->diffusion_range);
  for (se = lsys->symbol_list; se; se = se->next)
  {
    fprint_symbol(f, indent_depth + 2, se);
  }
/*
  fprint_indented(f, indent_depth + 2, "// list of transsys programs used:");
  glue = "";
  for (i = 0; i < lsys->num_transsys; i++)
  {
    fprintf(f, "%s %s", glue, lsys->transsys_list[i]->name);
    glue = ",";
  }
*/
  if (lsys->axiom)
  {
    fprintf(f, "\n");
    fprint_indented(f, indent_depth + 2, "axiom ");
    fprint_symbol_production(f, lsys->axiom, lsys->symbol_list, NULL);
    fprintf(f, ";\n");
  }
  else
    fprint_indented(f, indent_depth, "# no axiom\n");
  for(re = lsys->rule_list; re; re = re->next)
  {
    fprintf(f, "\n");
    fprint_rule(f, indent_depth + 2, re, lsys->symbol_list);
  }
  fprintf(f, "\n");
  fprint_indented(f, indent_depth + 2, "graphics\n");
  fprint_indented(f, indent_depth + 2, "{\n");
  for (se = lsys->symbol_list; se; se = se->next)
  {
    fprint_indented(f, indent_depth + 4, "%s\n", se->name);
    fprint_indented(f, indent_depth + 4, "{\n");
    /* fprintf(f, "# %d graphics primitives\n", se->num_graphics_primitives); */
    if (se->transsys)
      fprint_graphics_primitive_list(f, indent_depth + 6, se->graphics_primitive_list, se->transsys->factor_list);
    else
      fprint_graphics_primitive_list(f, indent_depth + 6, se->graphics_primitive_list, NULL);
    fprint_indented(f, indent_depth + 4, "}\n");
  }
  fprint_indented(f, indent_depth + 2, "}\n");
  fprint_indented(f, indent_depth, "}\n\n");
}


void fprint_symbol_instance(FILE *f, const SYMBOL_INSTANCE *si)
{
  int i;

  fprintf(f, "%s", si->lsys_string->lsys->symbol_list[si->symbol_index].name);
  if (si->transsys_instance.transsys)
  {
    fprintf(f, "(");
    for (i = 0; i < si->transsys_instance.transsys->num_factors; i++)
    {
      if (i > 0)
	fprintf(f, ", ");
      fprintf(f, "%s = %f", si->transsys_instance.transsys->factor_list[i].name, si->transsys_instance.factor_concentration[i]);
    }
    fprintf(f, ")");
  }
  if (si->num_successors == 0)
  {
    fprintf(f, " [no successors]");
  }
  else
  {
    fprintf(f, " [successors:%d:%d, pdist:%d]", si->successor_index, si->successor_index + si->num_successors - 1, si->successor_distance);
  }
}


void fprint_symbol_instance_list(FILE *f, const SYMBOL_INSTANCE *si, const char *sep)
{
  for ( ; si; si = si->next)
  {
    fprint_symbol_instance(f, si);
    fprintf(f, "%s", sep);
  }
}


/* FIXME: no functions for printing contact graphs yet */


void fprint_lsys_string(FILE *f, const LSYS_STRING *lstr, const char *sep)
{
  fprintf(f, "symbol string of lsys %s\n", lstr->lsys->name);
  fprint_symbol_instance_list(f, lstr->symbol, sep);
}


void fprint_lsys_string_contact_graph(FILE *f, const LSYS_STRING *lstr)
{
  int i, e, other_index;

  if (!lstr->arrayed)
  {
    fprintf(stderr, "fprint_lsys_string_contact_graph: cannot operate with non-arrayed lsys string\n");
    return;
  }
  for (i = 0; i < lstr->num_symbols; i++)
  {
    fprintf(f, "%d (%s):", i, lstr->lsys->symbol_list[lstr->symbol[i].symbol_index].name);
    for (e = 0; e < lstr->symbol[i].num_contact_edges; e++)
    {
      other_index = other_symbol_instance_index(lstr->symbol[i].contact_edge[e], i);
      fprintf(f, " %d(d=%d)", other_index, lstr->symbol[i].contact_edge[e]->distance);
    }
    fprintf(f, "\n");
  }
}


void fprint_transsys_instance(FILE *f, const TRANSSYS_INSTANCE *ti)
{
  int i;

  fprintf(f, "instance of transsys \"%s\"\n", ti->transsys->name);
  for (i = 0; i < ti->transsys->num_factors; i++)
  {
    fprintf(f, "[%s] = %g\n", ti->transsys->factor_list[i].name, ti->factor_concentration[i]);
  }
}


void fprint_transsys_instance_values(FILE *f, const TRANSSYS_INSTANCE *ti)
{
  int i;

  for (i = 0; i < ti->transsys->num_factors; i++)
  {
    fprintf(f, " %1.12e", ti->factor_concentration[i]);
  }
}


void fprint_cell(FILE *f, const CELL *cell)
{
  int i;

  fprintf(f, "existing: %d, alive: %d\n", cell->existing, cell->alive);
  fprintf(f, "neighbors: %d\n", cell->num_neighbors);
  fprint_transsys_instance(f, &(cell->transsys_instance));
  if (cell->num_neighbors)
  {
    fprintf(f, "contact weights:");
    for (i = 0; i < cell->num_neighbors; i++)
      fprintf(f, " %g", cell->contact_weight[i]);
  }
  fprintf(f, "\n");
}

