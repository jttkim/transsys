/* Copyright (C) 2002 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.1  2005/03/08 17:12:02  jtk
 * Initial revision
 *
 * Revision 1.1  2003/01/22 11:39:24  kim
 * added dot.c
 *
 */


#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "transsys.h"


#define DOT_ID_MAX 1024


static char *dot_expression_identifier(const EXPRESSION_NODE *expression, char buf[DOT_ID_MAX])
{
  sprintf(buf, "expr%p", expression);
  return (buf);
}


static char *dot_constitutive_identifier(const PROMOTER_ELEMENT *a, char buf[DOT_ID_MAX])
{
  sprintf(buf, "const%p", a);
  return (buf);
}


static int dot_expression_node(FILE *f, const TRANSSYS *transsys, const EXPRESSION_NODE *expression)
{
  char expr_id1[DOT_ID_MAX];

  fprintf(f, "        %s", dot_expression_identifier(expression, expr_id1));
  switch (expression->type)
  {
  case NT_NONE:
    fprintf(f, " [label=\"???\"]");
    break;
  case NT_VALUE:
    fprintf(f, " [label=\"%g\"]", expression->content.value);
    break;
  case NT_IDENTIFIER:
    fprintf(f, " [label=\"%s\"]", transsys->factor_list[expression->content.identifier.factor_index].name);
    break;
  case NT_RAW_IDENTIFIER:
    fprintf(f, " [label=\"(raw) %s.%s\"]", expression->content.raw_identifier.transsys_label, expression->content.raw_identifier.factor_name);
    break;
  case NT_ADD:
    fprintf(f, " [label=\"+\"]");
    break;
  case NT_SUBTRACT:
    fprintf(f, " [label=\"-\"]");
    break;
  case NT_MULT:
    fprintf(f, " [label=\"*\"]");
    break;
  case NT_DIV:
    fprintf(f, " [label=\"/\"]");
    break;
  case NT_LOWER:
    fprintf(f, " [label=\"<\"]");
    break;
  case NT_LOWER_EQUAL:
    fprintf(f, " [label=\"<=\"]");
    break;
  case NT_GREATER:
    fprintf(f, " [label=\">\"]");
    break;
  case NT_GREATER_EQUAL:
    fprintf(f, " [label=\">=\"]");
    break;
  case NT_EQUAL:
    fprintf(f, " [label=\"==\"]");
    break;
  case NT_UNEQUAL:
    fprintf(f, " [label=\"!=\"]");
    break;
  case NT_NOT:
    fprintf(f, " [label=\"!\"]");
    break;
  case NT_LOGICAL_AND:
    fprintf(f, " [label=\"&&\"]");
    break;
  case NT_LOGICAL_OR:
    fprintf(f, " [label=\"||\"]");
    break;
  case NT_RANDOM:
    fprintf(f, " [label=\"rnd\"]");
    break;
  case NT_GAUSS:
    fprintf(f, " [label=\"gauss\"]");
    break;
  default:
    fprintf(stderr, "dot_expression_node: unknown expression type %d\n", (int) expression->type);
    break;
  }
  fprintf(f, ";\n");
  switch (expression->type)
  {
  case NT_NONE:
  case NT_VALUE:
  case NT_IDENTIFIER:
  case NT_RAW_IDENTIFIER:
    break;
  case NT_NOT:
    dot_expression_node(f, transsys, expression->content.argument[0]);
    break;
  case NT_ADD:
  case NT_SUBTRACT:
  case NT_MULT:
  case NT_DIV:
  case NT_LOWER:
  case NT_LOWER_EQUAL:
  case NT_GREATER:
  case NT_GREATER_EQUAL:
  case NT_EQUAL:
  case NT_UNEQUAL:
  case NT_LOGICAL_AND:
  case NT_LOGICAL_OR:
  case NT_RANDOM:
  case NT_GAUSS:
  default:
    dot_expression_node(f, transsys, expression->content.argument[0]);
    dot_expression_node(f, transsys, expression->content.argument[1]);
    break;
  }
}


static int dot_expression_links(FILE *f, const TRANSSYS *transsys, const EXPRESSION_NODE *expression)
{
  char expr_id1[DOT_ID_MAX], expr_id2[DOT_ID_MAX];

  switch (expression->type)
  {
  case NT_IDENTIFIER:
    fprintf(f, "        %s -> %s [arrowhead=diamond];\n", transsys->factor_list[expression->content.identifier.factor_index].name, dot_expression_identifier(expression, expr_id1));
    break;
  case NT_NONE:
  case NT_VALUE:
  case NT_RAW_IDENTIFIER:
    break;
  case NT_NOT:
    fprintf(f, "        %s -> %s [arrowhead=none];\n", dot_expression_identifier(expression, expr_id1), dot_expression_identifier(expression->content.argument[0], expr_id2));
    dot_expression_links(f, transsys, expression->content.argument[0]);
    break;
  case NT_ADD:
  case NT_SUBTRACT:
  case NT_MULT:
  case NT_DIV:
  case NT_LOWER:
  case NT_LOWER_EQUAL:
  case NT_GREATER:
  case NT_GREATER_EQUAL:
  case NT_EQUAL:
  case NT_UNEQUAL:
  case NT_LOGICAL_AND:
  case NT_LOGICAL_OR:
  case NT_RANDOM:
  case NT_GAUSS:
  default:
    fprintf(f, "        %s -> %s [arrowhead=none];\n", dot_expression_identifier(expression, expr_id1), dot_expression_identifier(expression->content.argument[0], expr_id2));
    dot_expression_links(f, transsys, expression->content.argument[0]);
    fprintf(f, "        %s -> %s [arrowhead=none];\n", dot_expression_identifier(expression, expr_id1), dot_expression_identifier(expression->content.argument[1], expr_id2));
    dot_expression_links(f, transsys, expression->content.argument[1]);
    break;
  }
}




static int dot_gene_node(FILE *f, const TRANSSYS *transsys, const GENE_ELEMENT *ge)
{
  const PROMOTER_ELEMENT *a;
  char s[DOT_ID_MAX];

  fprintf(f, "    subgraph cluster_%s\n", ge->name);
  fprintf(f, "    {\n");
  fprintf(f, "      node [shape=box];\n");
  fprintf(f, "      %s [shape=box, peripheries=2];\n", ge->name);
  for (a = ge->promoter_list; a; a = a->next)
  {
    switch (a->type)
    {
    case ACT_CONSTITUTIVE:
      fprintf(f, "      subgraph cluster_%s;\n", dot_constitutive_identifier(a, s));
      fprintf(f, "      {\n");
      dot_expression_node(f, transsys, a->expr1);
      dot_expression_links(f, transsys, a->expr1);
      fprintf(f, "      }\n");
      fprintf(f, "      %s -> %s [arrowhead=dot,arrowtail=dot];\n", ge->name, dot_expression_identifier(a->expr1, s));
      break;
    }
  }
  fprintf(f, "    }\n");
  return (0);
}


static int dot_gene_product_link(FILE *f, const TRANSSYS *transsys, const GENE_ELEMENT *ge)
{
  fprintf(f, "  %s -> %s;\n", ge->name, transsys->factor_list[ge->product_index].name);
  return (0);
}


static int dot_factor(FILE *f, const TRANSSYS *transsys, const FACTOR_ELEMENT *fe)
{
  int i, j;
  const GENE_ELEMENT *ge;
  const PROMOTER_ELEMENT *a;

  for (i = 0; i < transsys->num_genes; i++)
  {
    ge = transsys->gene_list + i;
    for (a = ge->promoter_list; a; a = a->next)
    {
      switch (a->type)
      {
      case ACT_CONSTITUTIVE:
	break;
      case ACT_ACTIVATE:
	for (j = 0; j < a->num_binding_factors; j++)
	{
	  if (a->factor_index[j] == fe->index)
	    fprintf(f, "  %s -> %s [arrowhead=normal];\n", fe->name, ge->name);
	}
	break;
      case ACT_REPRESS:
	for (j = 0; j < a->num_binding_factors; j++)
	{
	  if (a->factor_index[j] == fe->index)
	    fprintf(f, "  %s -> %s [arrowhead=tee];\n", fe->name, ge->name);
	}
	break;
      default:
	fprintf(stderr, "dot_factor: cannot draw activation type %d\n", (int) a->type);
	break;
      }
    }
  }
  return (0);
}


int dot_transsys(FILE *f, const TRANSSYS *transsys)
{
  int i;

  fprintf(f, "digraph %s\n", transsys->name);
  fprintf(f, "{\n");
  fprintf(f, "  subgraph cluster_factors\n");
  fprintf(f, "  {\n");
  fprintf(f, "    style=invis;\n");
  fprintf(f, "    node [shape=ellipse, peripheries=2];\n");
  for (i = 0; i < transsys->num_factors; i++)
    fprintf(f, "    %s;\n", transsys->factor_list[i].name);
  fprintf(f, "  }\n");
  fprintf(f, "  subgraph cluster_genes\n");
  fprintf(f, "  {\n");
  fprintf(f, "    style=invis;\n");
  fprintf(f, "    node [shape=box];\n");
  for (i = 0; i < transsys->num_genes; i++)
  {
    dot_gene_node(f, transsys, transsys->gene_list + i);
  }
  fprintf(f, "  }\n");
  for (i = 0; i < transsys->num_genes; i++)
  {
    dot_gene_product_link(f, transsys, transsys->gene_list +i);
  }
  for (i = 0; i < transsys->num_factors; i++)
  {
    dot_factor(f, transsys, transsys->factor_list + i);
  }
  fprintf(f, "}\n");
}
