/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
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
  LSYS_STRING *lstr, *dstring;
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
    lstr = axiom_string(ls);
    if (lstr == NULL)
      fprintf(stderr, "axiom_string() returned NULL\n");
    else
    {
      for (i = 0; i < num_derivations; i++)
      {
	/* print_MemdebugStatistics(); */
	/* verify_Watchdogs(); */
	fprintf(outfile, "***** step %d *****\n", i);
	fprint_lsys_string(outfile, lstr, "\n");
	fprintf(outfile, "\n");
	if (process_transsys)
	  lsys_string_expression(lstr);
	dstring = derived_string(lstr);
	free_lsys_string(lstr);
	lstr = dstring;
	if (dstring == NULL)
	{
	  fprintf(stderr, "ltrcheck: derived_string() returned NULL\n");
	  break;
	}
      }
      if (lstr)
      {
	/* verify_Watchdogs(); */
	fprintf(outfile, "***** step %d *****\n", i);
	fprint_lsys_string(outfile, lstr, "\n");
	fprintf(outfile, "\n");
	free_lsys_string(lstr);
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

