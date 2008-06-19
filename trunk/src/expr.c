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

#include "trconfig.h"
#include "trtypes.h"
#include "transsys.h"


double evaluate_expression(const EXPRESSION_NODE *expr, const TRANSSYS_INSTANCE **ti_list)
{
  double arg1, arg2, return_value;
  int transsys_index;

  switch(expr->type)
  {
  case NT_VALUE:
    /* fprintf(stderr, "returning constant %g\n", expr->content.value); */
    return(expr->content.value);
  case NT_IDENTIFIER:
    if (expr->content.identifier.lhs_symbol_index < 0)
      transsys_index = 0;
    else
      /* ti is expected to contain pointers to transsys instances in lhs symbol order */
      transsys_index = expr->content.identifier.lhs_symbol_index;
    if ((expr->content.identifier.factor_index < 0) || (expr->content.identifier.factor_index >= ti_list[transsys_index]->transsys->num_factors))
    {
      fprintf(stderr, "Identifier index %d out of range [0, %d] for transsys \"%s\"\n", expr->content.identifier.factor_index, ti_list[transsys_index]->transsys->num_factors - 1, ti_list[transsys_index]->transsys->name);
      return (0.0);
    }
    /* fprintf(stderr, "returning value(%s) = %g\n", ti_list[transsys_index]->transsys->factor_list[expr->content.identifier.factor_index].name, concentration[expr->content.identifier.factor_index]); */
    return (ti_list[transsys_index]->factor_concentration[expr->content.identifier.factor_index]);
  case NT_LOGICAL_OR:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    return ((arg1 || arg2) ? 1.0 : 0.0);
  case NT_LOGICAL_AND:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    return ((arg1 && arg2) ? 1.0 : 0.0);
  case NT_LOWER:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    return ((arg1 < arg2) ? 1.0 : 0.0);
  case NT_LOWER_EQUAL:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    return ((arg1 <= arg2) ? 1.0 : 0.0);
  case NT_GREATER:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    return ((arg1 > arg2) ? 1.0 : 0.0);
  case NT_GREATER_EQUAL:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    return ((arg1 >= arg2) ? 1.0 : 0.0);
  case NT_EQUAL:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    return ((arg1 == arg2) ? 1.0 : 0.0);
  case NT_UNEQUAL:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    return ((arg1 != arg2) ? 1.0 : 0.0);
  case NT_NOT:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    return ((arg1 != 0.0) ? 1.0 : 0.0);
  case NT_ADD:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    /* fprintf(stderr, "returning %g + %g\n", arg1, arg2); */
    return (arg1 + arg2);
  case NT_SUBTRACT:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    /* fprintf(stderr, "returning %g - %g\n", arg1, arg2); */
    return (arg1 - arg2);
  case NT_MULT:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    /* fprintf(stderr, "returning %g * %g\n", arg1, arg2); */
    return (arg1 * arg2);
  case NT_DIV:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    /* fprintf(stderr, "returning %g / %g\n", arg1, arg2); */
    return (arg1 / arg2);
  case NT_RANDOM:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    if (arg1 < arg2)
      return_value = arg1 + (arg2 - arg1) * urandom_double();
    else
      return_value = arg2 + (arg1 - arg2) * urandom_double();
/*     fprintf(stderr, "evaluate_expression: random(%g, %g) = %g\n", arg1, arg2, return_value); */
    return (return_value);
  case NT_GAUSS:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    return (arg1 + arg2 * urandom_gauss());
  case NT_POW:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    /* FIXME: should check for valid arguments (positive base etc.) */
    return (pow(arg1, arg2));
  case NT_LOG:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    arg2 = evaluate_expression(expr->content.argument[1], ti_list);
    /* FIXME: should check for valid arguments */
    return (log(arg1) / log(arg2));
  case NT_ATAN:
    arg1 = evaluate_expression(expr->content.argument[0], ti_list);
    return (atan(arg1));
  default:
    fprintf(stderr, "evaluate_expression: unknown expression type %d\n", (int) expr->type);
    return (0.0);
  }
}
