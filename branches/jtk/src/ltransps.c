/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.1  2005/03/08 17:12:02  jtk
 * Initial revision
 *
 * Revision 1.2  2003/02/26 17:50:04  kim
 * fixed bug of returning success when parse errors occurred
 *
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include "trconfig.h"
#include "trtypes.h"
#include "transsys.h"


char *prgname = "transps";


int main(int argc, char **argv)
{
  int oc;
  extern char *optarg;
  extern int optind;
  char *outfile_name = NULL;
  FILE *outfile;
  int verbose = 0;
  int yyreturn;
  int num_derivations = 3;
  double scale = 1.0;
  POSTSCRIPT_STYLE style;
  SYMBOL_INSTANCE *symbol_string, *dstring;

  set_default_postscript_style(0, NULL, &style);
  while ((oc = getopt(argc, argv, "n:s:w:vh")) != -1)
  {
    switch(oc)
    {
    case 'n':
      num_derivations = strtol(optarg, NULL, 10);
      break;
    case 's':
      scale = strtod(optarg, NULL);
      break;
    case 'w':
      style.linewidth = strtod(optarg, NULL);
      break;
    case 'v':
      verbose = 1;
      fprintf(stderr, "verbose mode activated\n");
      break;
    case 'h':
      printf("ltransps -- make PostScript ltranssys derivation series\n");
      printf("usage: %s [options] <infile> <outfile>\n", argv[0]);
      printf("options:\n");
      printf("-n <num>: set number of derivations\n");
      printf("-h: print this help and exit\n");
      exit(EXIT_SUCCESS);
    }
  }
  if (optind < argc)
    yyin_name = argv[optind++];
  if (optind < argc)
    outfile_name = argv[optind++];
  if (yyin_name)
  {
    if ((yyin = fopen(yyin_name, "r")) == NULL)
    {
      fprintf(stderr, "Failed to open \"%s\" for input -- exit\n", yyin_name);
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    yyin = stdin;
    yyin_name = "stdin";
  }
  if (outfile_name)
  {
    if ((outfile = fopen(outfile_name, "w")) == NULL)
    {
      fprintf(stderr, "Failed to open \"%s\" for output -- exit\n", outfile_name);
      if (yyin != stdin)
        fclose(yyin);
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    outfile = stdout;
    outfile_name = "stdout";
  }
  yyreturn = yyparse();
  if (yyreturn)
  {
    fprintf(stderr, "Failed to parse from \"%s\" -- exit\n", yyin_name);
    if (yyin != stdin)
      fclose(yyin);
    if (outfile != stdout)
      fclose(outfile);
    exit(EXIT_FAILURE);
  }
  symbol_string = axiom_string(parsed_lsys);
  if (symbol_string == NULL)
  {
    fprintf(stderr, "ltransps: axiom_string() returned NULL\n");
    exit(EXIT_FAILURE);
  }
  postscript_lsys_prolog(outfile, parsed_lsys, &style);
  for (style.page_num = 0; style.page_num < num_derivations; style.page_num++)
  {
    if (verbose)
    {
      fprintf(stderr, "***** step %d *****\n", style.page_num);
      fprint_symbol_instance_list(stderr, symbol_string, "\n");
      fprintf(stderr, "\n");
    }
    postscript_symbol_string(outfile, symbol_string, &style, scale);
    string_transsys_expression(symbol_string);
    dstring = derived_string(parsed_lsys, symbol_string);
    free_symbol_instance_list(symbol_string);
    symbol_string = dstring;
    if (dstring == NULL)
    {
      fprintf(stderr, "ltransps: derived_string() returned NULL\n");
      break;
    }
  }
  if (symbol_string)
  {
    postscript_symbol_string(outfile, symbol_string, &style, scale);
    if (verbose)
    {
      fprintf(stderr, "***** step %d *****\n", style.page_num);
      fprint_symbol_instance_list(stderr, symbol_string, "\n");
      fprintf(stderr, "\n");
    }
  }
  free_lsys_list(parsed_lsys);
  free_transsys_list(parsed_transsys);
  if (yyin != stdin)
    fclose(yyin);
  if (outfile != stdout)
    fclose(outfile);
  return (EXIT_SUCCESS);
}

