/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.10  2005/06/16 09:36:26  jtk
 * implemented rule statistics gathering
 *
 * Revision 1.9  2005/05/20 10:40:15  jtk
 * differentiated entropy recording for lsys, expression and diffusion phase
 *
 * Revision 1.8  2005/05/17 12:11:30  jtk
 * contact graph works
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
extern void free_integer_array(INTEGER_ARRAY *ia);
extern INTEGER_ARRAY *extend_integer_array(INTEGER_ARRAY *ia, int v);
extern EXPRESSION_NODE *new_expression_node(EXPR_NODE_TYPE type, ...);
extern void free_promoter_list(PROMOTER_ELEMENT *alist);
extern PROMOTER_ELEMENT *new_promoter_element(PROMOTERELEMENT_TYPE type, int num_binding_factors, int *factors, EXPRESSION_NODE *expr1, EXPRESSION_NODE *expr2);
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
extern int extend_symbol_contact_edge_array(SYMBOL_INSTANCE *si, LSYS_STRING_CONTACT_EDGE *edge);
extern SYMBOL_ELEMENT *new_symbol_element(const char *name, const TRANSSYS *transsys);
extern int arrange_symbol_element_arrays(SYMBOL_ELEMENT *se);
extern void free_assignment_components(ASSIGNMENT *a);
extern void free_assignment_list(ASSIGNMENT *alist);
extern ASSIGNMENT *new_assignment(const TRANSSYS *target_transsys, int factor_index, EXPRESSION_NODE *value);
extern void free_production_element_components(PRODUCTION_ELEMENT *sp);
extern void free_production_element_list(PRODUCTION_ELEMENT *slist);
extern PRODUCTION_ELEMENT *new_production_element(int symbol_index, int template_lhs_symbol_index, ASSIGNMENT *assignment_list);
extern void free_symbol_production_components(SYMBOL_PRODUCTION *sp);
extern SYMBOL_PRODUCTION *new_symbol_production(PRODUCTION_ELEMENT *production_list);
extern int arrange_symbol_production_arrays(SYMBOL_PRODUCTION *sp);
extern void free_rule_element_components(RULE_ELEMENT *r);
extern void free_rule_element_list(RULE_ELEMENT *rlist);
extern RULE_ELEMENT *new_rule_element(const char *name, LHS_DESCRIPTOR *lhs, EXPRESSION_NODE *condition, SYMBOL_PRODUCTION *rhs);
extern void free_lsys_list(LSYS *ls);
extern void free_lsys_with_transsys(LSYS *ls);
extern LSYS *new_lsys(const char *name);
extern int arrange_lsys_arrays(LSYS *lsys);

extern void malloc_testloop(size_t n, const char *fname, int lineno);
extern void dump_memory(void *p, size_t n, const char *fname, int lineno);
extern void init_transsys_instance_components(TRANSSYS_INSTANCE *ti);
extern void free_transsys_instance_components(TRANSSYS_INSTANCE *ti);
extern int alloc_transsys_instance_components(TRANSSYS_INSTANCE *ti, const TRANSSYS *transsys);
extern int clone_transsys_instance(TRANSSYS_INSTANCE *ti, const TRANSSYS_INSTANCE *source);
extern void free_transsys_instance(TRANSSYS_INSTANCE *ti);
extern TRANSSYS_INSTANCE *new_transsys_instance(const TRANSSYS *transsys);
extern void free_cell_components(CELL *cell);
extern void free_cells(int num_cells, CELL *cell);
extern int alloc_cell_components(CELL *c, const TRANSSYS *transsys);
extern CELL *new_cells(int num_cells, const TRANSSYS *transsys);
extern void free_symbol_instance_components(SYMBOL_INSTANCE *si);
extern void free_symbol_instance_list(SYMBOL_INSTANCE *slist);
extern SYMBOL_INSTANCE *new_symbol_instance(const LSYS_STRING *lsys_string, int symbol_index);
extern SYMBOL_INSTANCE *clone_symbol_instance(const SYMBOL_INSTANCE *source, const LSYS_STRING *lsys_string);
extern void init_lsys_string_contact_graph_components(LSYS_STRING_CONTACT_GRAPH *g);
extern void free_lsys_string_contact_graph_components(LSYS_STRING_CONTACT_GRAPH *g);
extern void free_lsys_string_contact_graph(LSYS_STRING_CONTACT_GRAPH *g);
extern int alloc_lsys_string_contact_graph_components(LSYS_STRING_CONTACT_GRAPH *g, size_t array_size);
extern LSYS_STRING_CONTACT_GRAPH *new_lsys_string_contact_graph(size_t num_edges);
extern int add_lsys_string_contact_edge(LSYS_STRING_CONTACT_GRAPH *g, int i1, int i2, int distance);
extern int connect_lsys_string_symbols(LSYS_STRING *lstr, int i1, int i2, int distance);
extern void free_lsys_string(LSYS_STRING *lstr);
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
extern int other_symbol_instance_index(const LSYS_STRING_CONTACT_EDGE *edge, int si_index);
extern int lsys_string_diffusion(LSYS_STRING *lstr);
extern LSYS_STRING *axiom_string(const LSYS *lsys);
extern LSYS_STRING *derived_string(LSYS_STRING *lstr);

extern void fprint_transsys(FILE *f, int indent_depth, const TRANSSYS *transsys);

void fprint_transsys_as_discretenet(FILE *f, const TRANSSYS *transsys);

extern void fprint_lsys(FILE *f, int indent_depth, const LSYS *lsys);
extern void fprint_transsys_instance(FILE *f, const TRANSSYS_INSTANCE *ti);
extern void fprint_transsys_instance_values(FILE *f, const TRANSSYS_INSTANCE *ti);
extern void fprint_cell(FILE *f, const CELL *cell);
extern void fprint_symbol_instance(FILE *f, const SYMBOL_INSTANCE *si);
extern void fprint_symbol_instance_list(FILE *f, const SYMBOL_INSTANCE *si, const char *sep);
extern void fprint_lsys_string(FILE *f, const LSYS_STRING *lstr, const char *sep); 
extern void fprint_lsys_string_contact_graph(FILE *f, const LSYS_STRING *lstr); 
extern double evaluate_expression(const EXPRESSION_NODE *expr, const TRANSSYS_INSTANCE **ti_list);

extern int process_expression(TRANSSYS_INSTANCE *ti);

extern int set_default_postscript_style(int eps, const TRANSSYS *transsys, POSTSCRIPT_STYLE *style);
extern int postscript_transsys_init(FILE *f, const char *title, const POSTSCRIPT_STYLE *style);
extern int postscript_transsys_finish(FILE *f, const POSTSCRIPT_STYLE *style);
extern int postscript_transsys(FILE *f, const TRANSSYS *transsys, const POSTSCRIPT_STYLE *style);
extern int postscript_symbol_string(FILE *f, const SYMBOL_INSTANCE *symbol_string, const POSTSCRIPT_STYLE *style, double scale);
extern int postscript_lsys_prolog(FILE *f, const LSYS *lsys, const POSTSCRIPT_STYLE *style);

extern int dot_transsys(FILE *f, const TRANSSYS *transsys);

extern double *transsys_collection_factor_entropy(TRANSSYS_INSTANCE **ti);
extern double *transsys_collection_factor_information(TRANSSYS_INSTANCE **ti);

/* parse.c */
extern void add_factor_definition(TRANSSYS *ts, FACTOR_ELEMENT *fe);
extern void add_gene_definition(TRANSSYS *ts, GENE_ELEMENT *ge);
extern void add_factordef_decay(EXPRESSION_NODE *decay_expression, FACTOR_ELEMENT *fe);
extern void add_factordef_diffusibility(EXPRESSION_NODE *diffusibility_expression, FACTOR_ELEMENT *fe);
extern PROMOTER_ELEMENT *extend_promoter_list(PROMOTER_ELEMENT *alist, PROMOTER_ELEMENT *a);
extern int find_factor_index(const TRANSSYS *ts, const char *name);
extern int find_lhs_symbol_index(const RULE_ELEMENT *re, const char *transsys_label);
extern int resolve_simple_identifier(EXPRESSION_NODE *node, const TRANSSYS *transsys);
extern int resolve_complex_identifier(const LSYS *lsys, EXPRESSION_NODE *node, const RULE_ELEMENT *rule_element);
extern int resolve_identifiers(EXPRESSION_NODE *node, const TRANSSYS *transsys, const LSYS *lsys, const RULE_ELEMENT *rule_element);
extern INTEGER_ARRAY *extend_factor_combination(INTEGER_ARRAY *ia, const char *name, TRANSSYS *transsys);
extern PROMOTER_ELEMENT *create_promoter(PROMOTERELEMENT_TYPE atype, INTEGER_ARRAY *ia, EXPRESSION_NODE *expr1, EXPRESSION_NODE *expr2);
extern GENE_ELEMENT *create_gene(const char *name, PROMOTER_ELEMENT *alist, int gene_product);
extern int resolve_transsys(TRANSSYS *tr);
extern TRANSSYS *add_transsys(TRANSSYS *trlist, TRANSSYS *tr);
extern GRAPHICS_PRIMITIVE *add_graphics_primitive(GRAPHICS_PRIMITIVE *gplist, GRAPHICS_PRIMITIVE *gp);
extern SYMBOL_ELEMENT *find_symbol(const LSYS *ls, const char *name);
extern int find_symbol_index(const LSYS *ls, const char *name);
extern const TRANSSYS *find_transsys(const TRANSSYS *trlist, const char *name);
extern SYMBOL_ELEMENT *create_symbol_element(const char *name, const char *transsys_name, const TRANSSYS *transsys_list);
extern int find_lhs_symbol_by_transsys_label(const LHS_SYMBOL *symbol_list, const char *transsys_label);

extern ASSIGNMENT_RESOLUTION_RESULT resolve_raw_assignments(const SYMBOL_ELEMENT *se, LHS_SYMBOL *lhs_symbol_list, RAW_ASSIGNMENT *ra_list);
extern LHS_DESCRIPTOR *create_lhs_descriptor(const LSYS *lsys, LHS_SYMBOL *symlist);
extern PRODUCTION_ELEMENT *create_production_element(const char *symbol_name, RAW_ASSIGNMENT *ra_list, const LSYS *lsys, LHS_SYMBOL *lhs_symbol_list);
extern PRODUCTION_ELEMENT *add_production_element(PRODUCTION_ELEMENT *splist, PRODUCTION_ELEMENT *sp);
extern LHS_SYMBOL *create_lhs_symbol(const char *symbol_name, const char *transsys_label, const LSYS *lsys);
extern LHS_SYMBOL *add_lhs_symbol(LHS_SYMBOL *symlist, LHS_SYMBOL *symbol);
extern RAW_ASSIGNMENT *create_raw_assignment(const char *factor_name, const char *transsys_label, EXPRESSION_NODE *value);
extern RAW_ASSIGNMENT *add_raw_assignment(RAW_ASSIGNMENT *ra_list, RAW_ASSIGNMENT *ra);
extern int resolve_graphics_identifiers(GRAPHICS_PRIMITIVE *gp, const TRANSSYS *transsys);
extern int resolve_lsys(LSYS *ls);
extern LSYS *add_lsys(LSYS *lsyslist, LSYS *ls);
extern void add_symbol_definition(LSYS *ls, SYMBOL_ELEMENT *se);
extern void add_axiom_definition(LSYS *ls, SYMBOL_PRODUCTION *sp);
extern void add_diffusionrange_definition(LSYS *ls, double diffusionrange);
extern int resolve_rule_identifiers(RULE_ELEMENT *re, LSYS *lsys);
extern int add_rule_definition(LSYS *ls, RULE_ELEMENT *re);
extern void add_graphics_to_symbol(SYMBOL_ELEMENT *se, GRAPHICS_PRIMITIVE *gp_list);
extern void add_graphics_to_named_symbol(LSYS *ls, const char *symbol_name, GRAPHICS_PRIMITIVE *gp_list);

extern char *prgname;
extern char *yyin_name;
extern TRANSSYS *parsed_transsys;
extern LSYS *parsed_lsys;

extern const char transsys_revision[];

#endif /* TRANSSYS_H */

