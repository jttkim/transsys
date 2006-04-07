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
  /* TRANSSYS *t; */
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


static RULE_ELEMENT *complete_rule(const char *name, RULE_ELEMENT *re)
{
  strncpy(re->name, name, IDENTIFIER_MAX);
  re->name[IDENTIFIER_MAX - 1] = '\0';
  if (arrange_symbol_production_arrays(re->rhs) != 0)
    yyerror("complete_rule: arrange_symbol_production_arrays() failed");
  return (re);
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
	| ltrfile transsys  { parsed_transsys  = add_transsys(parsed_transsys, current_transsys); if (parsed_transsys == NULL) YYABORT; }
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
	: SYMBOL_DEF IDENTIFIER ';' { $$ = create_symbol_element($2, NULL, current_transsys); if ($$ == NULL) YYABORT; }
	| SYMBOL_DEF IDENTIFIER '(' IDENTIFIER ')' ';' { $$ = create_symbol_element($2, $4, current_transsys); if ($$ == NULL) YYABORT; }
	| SYMBOL_DEF '[' ';' { $$ = create_symbol_element("[", NULL, current_transsys); if ($$ == NULL) YYABORT; }
	| SYMBOL_DEF '[' '(' IDENTIFIER ')' ';' { $$ = create_symbol_element("[", $4, current_transsys); if ($$ == NULL) YYABORT; }
	| SYMBOL_DEF ']' ';' { $$ = create_symbol_element("]", NULL, current_transsys); if ($$ == NULL) YYABORT; }
	| SYMBOL_DEF ']' '(' IDENTIFIER ')' ';' { $$ = create_symbol_element("]", $4, current_transsys); if ($$ == NULL) YYABORT; }
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
	: IDENTIFIER  { $$ = create_lhs_symbol($1, NULL, current_lsys); }
	| IDENTIFIER '(' IDENTIFIER ')' { $$ = create_lhs_symbol($1, $3, current_lsys); }
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
	: IDENTIFIER { $$ = create_production_element($1, NULL, current_lsys, current_lhs_symbol_list); if ($$ == NULL) YYABORT; }
	| IDENTIFIER '(' transsys_initializer ')' { $$ = create_production_element($1, $3, current_lsys, current_lhs_symbol_list); if ($$ == NULL) YYABORT; }
	| '[' { $$ = create_production_element("[", NULL, current_lsys, current_lhs_symbol_list); if ($$ == NULL) YYABORT; }
	| ']' { $$ = create_production_element("]", NULL, current_lsys, current_lhs_symbol_list); if ($$ == NULL) YYABORT; }
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
	: DECAY_DEF ':' expr ';' { add_factordef_decay($3, current_factor); }
	| DIFFUSIBILITY_DEF ':' expr ';' { add_factordef_diffusibility($3, current_factor); }
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
	: CONSTITUTIVE ':' expr ';' { $$ = create_promoter(PROMOTERELEMENT_CONSTITUTIVE, NULL, $3, NULL); }
	| factor_combination ':' ACTIVATE '(' expr ',' expr ')' ';' { $$ = create_promoter(PROMOTERELEMENT_ACTIVATE, $1, $5, $7); }
	| factor_combination ':' REPRESS '(' expr ',' expr ')' ';' { $$ = create_promoter(PROMOTERELEMENT_REPRESS, $1, $5, $7); }
	;

factor_combination
	: IDENTIFIER { $$ = extend_factor_combination(NULL, $1, current_transsys); if ($$ == NULL) YYABORT; }
	| factor_combination '+' IDENTIFIER { $$ = extend_factor_combination($1, $3, current_transsys); if ($$ == NULL) YYABORT; }
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

