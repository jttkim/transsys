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


extern int yydebug;


int main(int argc, char **argv)
{
  int oc;
  extern char *optarg;
  extern int optind;
  char *outfile_name = NULL;
  FILE *outfile;
  int yyreturn;
  int num_derivations = 10, i;
  int process_transsys = 0;
  SYMBOL_INSTANCE *symbol_string, *dstring;
  LSYS *ls;

  while ((oc = getopt(argc, argv, "n:th")) != -1)
  {
    switch(oc)
    {
    case 't':
      process_transsys = 1;
      break;
    case 'n':
      num_derivations = strtol(optarg, NULL, 10);
      break;
    case 'h':
      printf("-n <num>: set number of derivations\n");
      printf("-t: activate transsys between derivations\n");
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
  for (ls = parsed_lsys; ls; ls = ls->next)
  {
    symbol_string = axiom_string(ls);
    if (symbol_string == NULL)
      fprintf(stderr, "axiom_string() returned NULL\n");
    else
    {
      for (i = 0; i < num_derivations; i++)
      {
	fprintf(outfile, "***** step %d *****\n", i);
	fprint_symbol_instance_list(outfile, symbol_string, "\n");
	fprintf(outfile, "\n");
	if (process_transsys)
	  string_transsys_expression(symbol_string);
	dstring = derived_string(ls, symbol_string);
	free_symbol_instance_list(symbol_string);
	symbol_string = dstring;
	if (dstring == NULL)
	{
	  fprintf(stderr, "ltrcheck: derived_string() returned NULL\n");
	  break;
	}
      }
      if (symbol_string)
      {
	fprintf(outfile, "***** step %d *****\n", i);
	fprint_symbol_instance_list(outfile, symbol_string, "\n");
	fprintf(outfile, "\n");
	free_symbol_instance_list(symbol_string);
      }
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

