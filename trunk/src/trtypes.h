/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.5  2005/05/16 12:02:10  jtk
 * in transition from distance matrices to contact graphs
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
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#ifndef TRTYPES_H
#define TRTYPES_H

#include <stddef.h>  /* for size_t */

#include "trconfig.h"


typedef enum
{
  NT_NONE,
  NT_VALUE,
  NT_IDENTIFIER,
  NT_RAW_IDENTIFIER,
  NT_ADD,
  NT_SUBTRACT,
  NT_MULT,
  NT_DIV,
  NT_LOWER,
  NT_LOWER_EQUAL,
  NT_GREATER,
  NT_GREATER_EQUAL,
  NT_EQUAL,
  NT_UNEQUAL,
  NT_NOT,
  NT_LOGICAL_AND,
  NT_LOGICAL_OR,
  NT_RANDOM,
  NT_GAUSS
} EXPR_NODE_TYPE;

typedef enum
{
  ACT_NONE,
  ACT_CONSTITUTIVE,
  ACT_ACTIVATE,
  ACT_REPRESS
} ACTIVATION_TYPE;

typedef enum
{
  GRAPHICS_NONE,
  GRAPHICS_MOVE,
  GRAPHICS_PUSH,
  GRAPHICS_POP,
  GRAPHICS_TURN,
  GRAPHICS_ROLL,
  GRAPHICS_BANK,
  GRAPHICS_SPHERE,
  GRAPHICS_CYLINDER,
  GRAPHICS_BOX,
  GRAPHICS_COLOR
} GRAPHICS_PRIMITIVE_TYPE;

typedef enum
{
  EC_SINGLE_TRANSSYS,
  EC_MULTI_TRANSSYS
} EXPRESSION_CONTEXT_TYPE;


struct tag_transsys;
typedef struct tag_transsys TRANSSYS;


typedef struct
{
  size_t length;
  int *array;
} INTEGER_ARRAY;

typedef struct
{
  double x, y, w, h;
} POSTSCRIPT_BOX;

typedef struct
{
  POSTSCRIPT_BOX pagebox;
  int page_num;
  double lnd_x0, lnd_y0;
  char eps;
  double linewidth;
  double box_linewidth;
  double arrow_width, arrow_length, arrow_linewidth;
  double fontheight;
  char fontname[IDENTIFIER_MAX];
} POSTSCRIPT_STYLE;

typedef struct tag_expression_node
{
  EXPR_NODE_TYPE type;
  union
  {
    struct
    {
      int lhs_symbol_index;
      int factor_index;
    } identifier;
    double value;
    struct tag_expression_node **argument;
    struct
    {
      char *factor_name;
      char *transsys_label;
    } raw_identifier;
  } content;
} EXPRESSION_NODE;

typedef struct tag_promoter
{
  struct tag_promoter *next;
  ACTIVATION_TYPE type;
  int num_binding_factors;
  int *factor_index;
  EXPRESSION_NODE *expr1, *expr2;
} PROMOTER_ELEMENT;

typedef struct tag_factor
{
  struct tag_factor *next;
  int index;
  int num_producing_genes;
  int *gene_index;
  EXPRESSION_NODE *diffusibility_expression;
  EXPRESSION_NODE *decay_expression;
  char name[IDENTIFIER_MAX];
} FACTOR_ELEMENT;

typedef struct tag_gene
{
  struct tag_gene *next;
  int index;
  int product_index;
  PROMOTER_ELEMENT *promoter_list;
  POSTSCRIPT_BOX box;
  char name [IDENTIFIER_MAX];
} GENE_ELEMENT;

struct tag_transsys
{
  struct tag_transsys *next;
  int arrayed;
  int num_factors, num_genes;
  FACTOR_ELEMENT *factor_list;
  GENE_ELEMENT *gene_list;
  char name[IDENTIFIER_MAX];
};

typedef struct tag_transsys_instance
{
  const TRANSSYS *transsys;
  double *factor_concentration, *new_concentration;
} TRANSSYS_INSTANCE;  

typedef struct tag_cell
{
  int existing, alive;
  TRANSSYS_INSTANCE transsys_instance;
  int num_neighbors;
  struct tag_cell **neighbor;
  double *contact_weight;
} CELL;

typedef struct tag_graphics
{
  struct tag_graphics *next;
  GRAPHICS_PRIMITIVE_TYPE type;
  int num_arguments;
  EXPRESSION_NODE **argument;
} GRAPHICS_PRIMITIVE;

typedef struct tag_symbol
{
  struct tag_symbol *next;
  char name[IDENTIFIER_MAX];
  int index;
  int arrayed;
  const TRANSSYS *transsys;
  int num_graphics_primitives;
  GRAPHICS_PRIMITIVE *graphics_primitive_list;
} SYMBOL_ELEMENT;

typedef struct tag_raw_assignment
{
  struct tag_raw_assignment *next;
  char identifier[IDENTIFIER_MAX];
  char transsys_label[IDENTIFIER_MAX];
  EXPRESSION_NODE *value;
} RAW_ASSIGNMENT;

typedef struct tag_assignment
{
  struct tag_assignment *next;
  const TRANSSYS *target_transsys;
  int factor_index;
  EXPRESSION_NODE *value;
} ASSIGNMENT;

typedef struct tag_production_element
{
  struct tag_production_element *next;
  int arrayed;
  int symbol_index;
  int template_lhs_symbol_index;
  int num_assignments;
  ASSIGNMENT *assignment_list;
} PRODUCTION_ELEMENT;

typedef struct tag_symbol_production
{
  int arrayed;
  const TRANSSYS *transsys; /* transsys governing production -- obsolescent */
  int num_production_elements;
  PRODUCTION_ELEMENT *production_list;
} SYMBOL_PRODUCTION;

typedef struct tag_lhs_symbol
{
  struct tag_lhs_symbol *next;
  int index;
  char transsys_label[IDENTIFIER_MAX];
  const TRANSSYS *transsys;
  int symbol_index;
} LHS_SYMBOL;

typedef struct tag_lhs_descriptor
{
  int arrayed;
  int num_symbols;
  LHS_SYMBOL *symbol_list;
} LHS_DESCRIPTOR;

/* not used yet -- will this ever be useful?? */

typedef struct tag_expression_context_entry
{
  struct tag_expression_context_entry *next;
  const TRANSSYS *transsys;
  char label[IDENTIFIER_MAX];
} EXPRESSION_CONTEXT_ENTRY;

typedef struct tag_expression_context
{
  int num_transsys;
  EXPRESSION_CONTEXT_ENTRY *entry;
} EXPRESSION_CONTEXT;

typedef struct tag_rule
{
  struct tag_rule *next;
  char name[IDENTIFIER_MAX];
  LHS_DESCRIPTOR *lhs;
  EXPRESSION_NODE *condition;
  SYMBOL_PRODUCTION *rhs;
} RULE_ELEMENT;

typedef struct tag_lsys
{
  struct tag_lsys *next;
  char name[IDENTIFIER_MAX];
  int arrayed;
  int diffusion_range;
  int num_symbols, num_rules;
  SYMBOL_ELEMENT *symbol_list;
  SYMBOL_PRODUCTION *axiom;
  RULE_ELEMENT *rule_list;
} LSYS;


struct tag_lsys_string;
typedef struct tag_lsys_string LSYS_STRING;
struct tag_lsys_string_contact_edge;
typedef struct tag_lsys_string_contact_edge LSYS_STRING_CONTACT_EDGE;

typedef struct tag_symbol_instance
{
  struct tag_symbol_instance *next;
  const LSYS_STRING *lsys_string;
  int symbol_index;
  int lhs_group_start, lhs_group_length;
  int successor_distance;
  int num_successors, successor_index;
  int num_contact_edges;
  LSYS_STRING_CONTACT_EDGE **contact_edge;
  TRANSSYS_INSTANCE transsys_instance;
} SYMBOL_INSTANCE;

/*
 * Implementation notices: Even though the contact graph is undirected
 * conceptually, the order of i1 and i2 does matter: Flow is directed
 * from symbol i1 to symbol i2, a net flow from symbol i2 to symbol i1
 * is represented by negative values of amount_diffused.
 * The type of i1 and i2 is int, not size_t, so that -1 can be used to
 * represent "not an index".
 */

struct tag_lsys_string_contact_edge
{
  int i1, i2;  /* indices of contacting symbols */
  int distance;
  int amount_valid;
  double amount_diffused;
}; /* LSYS_STRING_CONTACT_EDGE */

/*
 * Implementation notice:
 * the edge component of a contact graph may not be reallocated or
 * otherwise moved around in memory, unless it is certain that there
 * are no references to edges from symbol instances that could be
 * rendered invalid by such relocations.
 */

typedef struct
{
  size_t array_size;
  size_t num_edges;
  LSYS_STRING_CONTACT_EDGE *edge;
} LSYS_STRING_CONTACT_GRAPH;

struct tag_lsys_string
{
  const LSYS *lsys;
  int arrayed;
  size_t num_symbols;
  SYMBOL_INSTANCE *symbol;
  LSYS_STRING_CONTACT_GRAPH contact_graph;
};

#endif /* TRTYPES_H */

