/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.4  2005/03/30 18:30:27  jtk
 * progressed transition to arrayred lsys strings
 * introduced lsys string distance matrices
 *
 * Revision 1.3  2005/03/30 09:51:56  jtk
 * added include for dmalloc (contitional upon DMALLOC macro)
 *
 * Revision 1.2  2005/03/29 17:33:02  jtk
 * introduced arrayed lsys string, with symbol distance matrix.
 *
 * Revision 1.1.1.1  2005/03/08 17:12:02  jtk
 * new cvs after loss at INB
 *
 * Revision 1.4  2002/12/09 19:37:55  kim
 * added dot rendering of transsys network graphs (initial version)
 *
 * Revision 1.3  2002/05/28 08:52:51  kim
 * added some discrete network stuff to sources
 *
 * Revision 1.2  2002/01/25 03:35:04  kim
 * Added gnuplot link functionality to transexpr, transscatter, improved
 *     PostScript stuff
 *
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#ifndef TRANSSYS_H
#define TRANSSYS_H

#ifdef MEMDEBUG
#  include <memdebug.h>
#endif
#ifdef DMALLOC
#  include <dmalloc.h>
#endif

#include "trtypes.h"

#define NO_INDEX -1
#define GRAPHICS_PRIMITIVE_ARGMAX 3


extern TRANSSYS *parsed_transsys;
extern long lineno;
extern FILE *yyin;
extern char *yyin_name;
extern int yylex(void);
extern int yyparse(void);
extern void yyerror(const char *s, ...);

extern void ulong_srandom( unsigned int x);
extern unsigned long *ulong_initstate(unsigned int seed, unsigned long *arg_state);
extern unsigned long *ulong_setstate(unsigned long *arg_state);
extern long int ulong_random(void);
extern unsigned long random_long(unsigned long range);
extern double urandom_double(void);
extern double urandom_gauss(void);
extern int write_rndgenerator_state(FILE *f);
extern int read_rndgenerator_state(FILE *f);
extern void seed_ulong_random(int seed);
extern void random_shuffle(long num, long *a);

void free_deadbeef(void *p);
extern INTEGER_ARRAY *extend_integer_array(INTEGER_ARRAY *ia, int v);
extern EXPRESSION_NODE *new_expression_node(EXPR_NODE_TYPE type, ...);
extern PROMOTER_ELEMENT *new_promoter_element(ACTIVATION_TYPE type, int num_binding_factors, int *factors, EXPRESSION_NODE *expr1, EXPRESSION_NODE *expr2);
extern FACTOR_ELEMENT *new_factor_element(const char *name, EXPRESSION_NODE *decay_expression, EXPRESSION_NODE *diffusibility_expression);
extern GENE_ELEMENT *new_gene_element(const char *name, PROMOTER_ELEMENT *promoter_list, int product_index);
extern void free_expression_tree(EXPRESSION_NODE *node);
extern void free_transsys_list(TRANSSYS *tr);
extern TRANSSYS *new_transsys(const char *name);
extern int arrange_transsys_arrays(TRANSSYS *transsys);
extern void free_graphics_primitive_components(GRAPHICS_PRIMITIVE *gp);
extern void free_graphics_primitive_list(GRAPHICS_PRIMITIVE *glist);
extern GRAPHICS_PRIMITIVE *new_graphics_primitive(GRAPHICS_PRIMITIVE_TYPE type, ...);
extern void free_lhs_symbol_components(LHS_SYMBOL *ls);
extern void free_lhs_symbol_list(LHS_SYMBOL *ls_list);
extern LHS_SYMBOL *new_lhs_symbol(const char *transsys_label, const TRANSSYS *transsys, int symbol_index);
extern void free_lhs_descriptor_components(LHS_DESCRIPTOR *l);
extern LHS_DESCRIPTOR *new_lhs_descriptor(LHS_SYMBOL *symbol_list);
extern int arrange_lhs_descriptor_arrays(LHS_DESCRIPTOR *l);
extern void free_symbol_element_components(SYMBOL_ELEMENT *sym);
extern void free_symbol_element_list(SYMBOL_ELEMENT *sym);
extern SYMBOL_ELEMENT *new_symbol_element(const char *name, const TRANSSYS *transsys);
extern int arrange_symbol_element_arrays(SYMBOL_ELEMENT *se);
extern void free_assignment_components(ASSIGNMENT *a);
extern void free_assignment_list(ASSIGNMENT *alist);
extern ASSIGNMENT *new_assignment(const TRANSSYS *target_transsys, int factor_index, EXPRESSION_NODE *value);
extern void free_production_element_components(PRODUCTION_ELEMENT *sp);
extern void free_production_element_list(PRODUCTION_ELEMENT *slist);
extern PRODUCTION_ELEMENT *new_production_element(int symbol_index, int template_lhs_symbol_index, ASSIGNMENT *assignment_list);
extern void free_symbol_production_components(SYMBOL_PRODUCTION *sp);
extern SYMBOL_PRODUCTION *new_symbol_production(const TRANSSYS *transsys, PRODUCTION_ELEMENT *production_list);
extern int arrange_symbol_production_arrays(SYMBOL_PRODUCTION *sp);
extern void free_rule_element_components(RULE_ELEMENT *r);
extern void free_rule_element_list(RULE_ELEMENT *rlist);
extern RULE_ELEMENT *new_rule_element(const char *name, LHS_DESCRIPTOR *lhs, EXPRESSION_NODE *condition, SYMBOL_PRODUCTION *rhs);
extern void free_lsys_list(LSYS *ls);
extern LSYS *new_lsys(const char *name);
extern int arrange_lsys_arrays(LSYS *lsys);

extern void malloc_testloop(size_t n, const char *fname, int lineno);
extern void dump_memory(void *p, size_t n, const char *fname, int lineno);
extern void init_transsys_instance_components(TRANSSYS_INSTANCE *ti);
extern void free_transsys_instance_components(TRANSSYS_INSTANCE *ti);
extern int alloc_transsys_instance_components(TRANSSYS_INSTANCE *ti, const TRANSSYS *transsys);
extern int clone_transsys_instance(TRANSSYS_INSTANCE *ti, const TRANSSYS_INSTANCE *source);
extern void free_cell_components(CELL *cell);
extern void free_cells(int num_cells, CELL *cell);
extern int alloc_cell_components(CELL *c, const TRANSSYS *transsys);
extern CELL *new_cells(int num_cells, const TRANSSYS *transsys);
extern void free_symbol_instance_components(SYMBOL_INSTANCE *si);
extern void free_symbol_instance_list(SYMBOL_INSTANCE *slist);
extern SYMBOL_INSTANCE *new_symbol_instance(const LSYS_STRING *lsys_string, int symbol_index, int num_predecessors, int predecessor_index, int predecessor_distance);
extern SYMBOL_INSTANCE *clone_symbol_instance(const SYMBOL_INSTANCE *source, const LSYS_STRING *lsys_string, int predecessor_index);
extern void free_lsys_string_distance(LSYS_STRING *lstr);
extern void free_lsys_string(LSYS_STRING *lstr);
extern int alloc_lsys_string_distance(LSYS_STRING *lstr);
extern int arrange_lsys_string_arrays(LSYS_STRING *lstr);
extern LSYS_STRING *new_lsys_string(const LSYS *lsys);

extern int add_neighbor(CELL *cell1, CELL *cell2, double weight);
extern int set_contact_weight(CELL *cell1, CELL *cell2, double weight);
extern int cell_expression(CELL *cell);
extern int diffuse(int num_cells, CELL *cell);
extern double total_concentration(int num_cells, const CELL *cell, int factor_no);
extern double max_concentration(int num_cells, const CELL *cell, int factor_no);

extern size_t symbol_strlen(const SYMBOL_INSTANCE *symbol_string);
extern size_t lsys_string_length(const LSYS_STRING *lstr);
extern int lsys_string_expression(LSYS_STRING *lstr);
extern int lsys_string_diffusion(LSYS_STRING *lstr);
extern LSYS_STRING *axiom_string(const LSYS *lsys);
extern LSYS_STRING *derived_string(const LSYS_STRING *lstr);

extern void fprint_transsys(FILE *f, int indent_depth, const TRANSSYS *transsys);

void fprint_transsys_as_discretenet(FILE *f, const TRANSSYS *transsys);

extern void fprint_lsys(FILE *f, int indent_depth, const LSYS *lsys);
extern void fprint_factorconc_commentline(FILE *f, const TRANSSYS_INSTANCE *ti);
extern void fprint_factorconc_line(FILE *f, const TRANSSYS_INSTANCE *ti, unsigned long time_step);
extern void fprint_transsys_instance(FILE *f, const TRANSSYS_INSTANCE *ti);
extern void fprint_cell(FILE *f, const CELL *cell);
extern void fprint_symbol_instance(FILE *f, const SYMBOL_INSTANCE *si);
extern void fprint_symbol_instance_list(FILE *f, const SYMBOL_INSTANCE *si, const char *sep);
extern void fprint_lsys_string(FILE *f, const LSYS_STRING *lstr, const char *sep); 
extern double evaluate_expression(const EXPRESSION_NODE *expr, const TRANSSYS_INSTANCE **ti_list);

extern int process_expression(TRANSSYS_INSTANCE *ti);

extern int set_default_postscript_style(int eps, const TRANSSYS *transsys, POSTSCRIPT_STYLE *style);
extern int postscript_transsys_init(FILE *f, const char *title, const POSTSCRIPT_STYLE *style);
extern int postscript_transsys_finish(FILE *f, const POSTSCRIPT_STYLE *style);
extern int postscript_transsys(FILE *f, const TRANSSYS *transsys, const POSTSCRIPT_STYLE *style);
extern int postscript_symbol_string(FILE *f, const SYMBOL_INSTANCE *symbol_string, const POSTSCRIPT_STYLE *style, double scale);
extern int postscript_lsys_prolog(FILE *f, const LSYS *lsys, const POSTSCRIPT_STYLE *style);

extern int dot_transsys(FILE *f, const TRANSSYS *transsys);

extern char *prgname;
extern char *yyin_name;
extern TRANSSYS *parsed_transsys;
extern LSYS *parsed_lsys;

#endif /* TRANSSYS_H */

