/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

%{

/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  2005/04/05 20:26:35  jtk
 * fixed realvalue regexp
 *
 * Revision 1.2  2005/03/31 16:07:36  jtk
 * finished (initial) implementation of lsys diffusion
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
#include <stdarg.h>

#include "trconfig.h"
#include "transsys.h"

#include "trbison.h"


%}

%%
\/\/.*			;
factor			{ return (FACTOR_DEF); }
gene			{ return (GENE_DEF); }
promoter		{ return (PROMOTER_DEF); }
product			{ return (PRODUCT_DEF); }
constitutive		{ return (CONSTITUTIVE); }
activate		{ return (ACTIVATE); }
repress			{ return (REPRESS); }
default			{ return (DEFAULT); }
transsys		{ return (TRANSSYS_DEF); }
decay			{ return (DECAY_DEF); }
diffusibility		{ return (DIFFUSIBILITY_DEF); }
synthesis		{ return (SYNTHESIS_DEF); }
diffusionrange		{ return (DIFFUSIONRANGE_DEF); }
lsys			{ return (LSYS_DEF); }
symbol			{ return (SYMBOL_DEF); }
axiom			{ return (AXIOM_DEF); }
rule			{ return (RULE_DEF); }
-->			{ return (ARROW); }
graphics		{ return (GRAPHICS_DEF); }
move			{ return (MOVE); }
sphere			{ return (SPHERE); }
cylinder		{ return (CYLINDER); }
box			{ return (BOX); }
turn			{ return (TURN); }
roll			{ return (ROLL); }
bank			{ return (BANK); }
color			{ return (COLOR); }
push			{ return (PUSH); }
pop			{ return (POP); }
pow			{ return (NUMERIC_POW); }
log			{ return (NUMERIC_LOG); }
atan			{ return (NUMERIC_ATAN); }
gauss			{ return (NUMERIC_GAUSS); }
random			{ return (NUMERIC_RANDOM); }
\<=			{ return (LOWER_EQUAL); }
\>=			{ return (GREATER_EQUAL); }
==			{ return (EQUAL); }
!=			{ return (UNEQUAL); }
&&			{ return (LOGICAL_AND); }
\|\|			{ return (LOGICAL_OR); }
#.*			{ /* fprintf(stderr, "line %ld: comment\n", lineno); */ }
[A-Za-z_][A-Za-z0-9_]* 	{
			  strncpy(yylval.identifier, yytext, IDENTIFIER_MAX - 1);
			  if (yyleng > IDENTIFIER_MAX - 1)
			  {
			    fprintf(stderr, "%ld: excessively long identifier\n", lineno);
			    yylval.identifier[IDENTIFIER_MAX - 1] = '\0';
			  }
			  return (IDENTIFIER);
			}
[+-]?[0-9]+(\.[0-9]*)?([Ee][+-]?[0-9]+)?	{
			  yylval.real = strtod(yytext, NULL);
			  return (REALVALUE);
			}
[+-]?\.[0-9]+([Ee][+-]?[0-9]+)?	{
			  yylval.real = strtod(yytext, NULL);
			  return (REALVALUE);
			}
[\t ]			;
.			{ return (yytext[0]); }
\n			{ lineno++; }

%%

/*** for testing the scanner ***/

#ifdef TEST_SCANNER

YYSTYPE yylval;
long lineno = 0;


void printf_token(char *fmt, ...)
{
  va_list ap;

  va_start(ap, fmt);
  printf("%ld: ", lineno);
  vprintf(fmt, ap);
  printf("\n");
}


int main(int argc, char *argv[])
{
  int retval;

  while (1)
  {
    retval = yylex();
    if (retval == 0)
      break;
    switch (retval)
    {
    case IDENTIFIER:
      printf_token("<IDENTIFIER> %s ", yylval.identifier);
      break;
    case REALVALUE:
      printf_token("<REALVALUE> %f ", yylval.real);
      break;
    case FACTOR_DEF:
      printf_token("<FACTOR_DEF> ");
      break;
    case GENE_DEF:
      printf_token("<GENE_DEF> ");
      break;
    case PROMOTER_DEF:
      printf_token("<PROMOTER_DEF> ");
      break;
    case PRODUCT_DEF:
      printf_token("<PRODUCT_DEF> ");
      break;
    case CONSTITUTIVE:
      printf_token("<CONSTITUTIVE> ");
      break;
    case ACTIVATE:
      printf_token("<ACTIVATE> ");
      break;
    case REPRESS:
      printf_token("<REPRESS> ");
      break;
    case DEFAULT:
      printf_token("<DEFAULT> ");
      break;
    case RANDOM:
      printf_token("<RANDOM> ");
      break;
    case GAUSS:
      printf_token("<GAUSS> ");
      break;
    case DECAY_DEF:
      printf_token("<DECAY_DEF> ");
      break;
    case DIFFUSIBILITY_DEF:
      printf_token("<DIFFUSIBILITY_DEF> ");
      break;
    case SYNTHESIS_DEF:
      printf_token("<SYNTHESIS_DEF> ");
      break;
    case TRANSSYS_DEF:
      printf_token("<TRANSSYS_DEF> ");
      break;
    case LSYS_DEF:
      printf_token("<LSYS_DEF> ");
      break;
    case SYMBOL_DEF:
      printf_token("<SYMBOL_DEF> ");
      break;
    case AXIOM_DEF:
      printf_token("<AXIOM_DEF> ");
      break;
    case RULE_DEF:
      printf_token("<RULE_DEF> ");
      break;
    case ARROW:
      printf_token("<ARROW> ");
      break;
    case GRAPHICS_DEF:
      printf_token("<GRAPHICS_DEF> ");
      break;
    case MOVE:
      printf_token("<MOVE> ");
      break;
    case SPHERE:
      printf_token("<SPHERE> ");
      break;
    case CYLINDER:
      printf_token("<CYLINDER> ");
      break;
    case BOX:
      printf_token("<BOX> ");
      break;
    case TURN:
      printf_token("<TURN> ");
      break;
    case ROLL:
      printf_token("<ROLL> ");
      break;
    case BANK:
      printf_token("<BANK> ");
      break;
    case LOWER_EQUAL:
      printf_token("<LOWER_EQUAL> ");
      break;
    case GREATER_EQUAL:
      printf_token("<GREATER_EQUAL> ");
      break;
    case EQUAL:
      printf_token("<EQUAL> ");
      break;
    case UNEQUAL:
      printf_token("<UNEQUAL> ");
      break;
    case LOGICAL_AND:
      printf_token("<LOGICAL_AND> ");
      break;
    case LOGICAL_OR:
      printf_token("<LOGICAL_OR> ");
      break;
    case '\n':
      printf_token("\a<<newline>>");
      break;
    default:
      if (retval < 256)
	printf_token("\'%c\' ", retval);
      else
        printf_token("<TOKEN %d> ", retval);
      break;
    }
  }
  printf("\n");
  printf("feof(stdin) = %d\n", feof(stdin));
  return (EXIT_SUCCESS);
}

#endif /* TEST_SCANNER */

