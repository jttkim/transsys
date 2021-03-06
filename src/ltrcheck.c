/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.6  2005/06/16 09:36:26  jtk
 * implemented rule statistics gathering
 *
 * Revision 1.5  2005/04/04 09:39:54  jtk
 * added lsys capabilities to transexpr, various small changes
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


static int fprint_rulestats(FILE *f, int n, const LSYS_STRING *lstr)
{
  int i;

  for (i = 0; i < lstr->num_symbols; i++)
  {
    if (lstr->lsys->symbol_list[lstr->symbol[i].symbol_index].transsys)
    {
      fprintf(f, "%d %s", n, lstr->lsys->symbol_list[lstr->symbol[i].symbol_index].name);
      if (lstr->symbol[i].rule_index == NO_INDEX)
      {
	fprintf(f, " clone-rule");
      }
      else
      {
	fprintf(f, " %s", lstr->lsys->rule_list[lstr->symbol[i].rule_index].name);
      }
      fprintf(f, " %s", lstr->lsys->symbol_list[lstr->symbol[i].symbol_index].transsys->name);
      fprint_transsys_instance_values(f, &(lstr->symbol[i].transsys_instance));
      fprintf(f, "\n");
    }
    else
    {
      fprintf(f, "# bare symbol %s\n", lstr->lsys->symbol_list[lstr->symbol[i].symbol_index].name);
    }
  }
  return (0);
}


int main(int argc, char **argv)
{
  int oc;
  extern char *optarg;
  extern int optind;
  char *outfile_name = NULL, *rulestatsfile_name = NULL;
  FILE *outfile, *rulestatsfile = NULL;
  int yyreturn;
  int num_derivations = 10, i;
  int process_transsys = 1, process_diffusion = 1;
  unsigned int rndseed = 1;
  LSYS_STRING *lstr, *dstring;
  LSYS *ls;

  while ((oc = getopt(argc, argv, "s:d:r:DXh")) != -1)
  {
    switch(oc)
    {
    case 'X':
      process_transsys = 0;
      break;
    case 'D':
      process_diffusion = 0;
      break;
    case 'd':
      num_derivations = strtol(optarg, NULL, 10);
      break;
    case 'r':
      rulestatsfile_name = optarg;
      break;
    case 's':
      rndseed = strtoul(optarg, NULL, 10);
      break;
    case 'h':
      printf("-d <num>: set number of derivations\n");
      printf("-s <num>: set random seed\n");
      printf("-r <filename>: specify rule stats file\n");
      printf("-X: suppress gene expression (and factor decay) between derivations\n");
      printf("-D: suppress diffusion of factors between derivations\n");
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
      {
        fclose(yyin);
      }
      exit(EXIT_FAILURE);
    }
  }
  else
  {
    outfile = stdout;
    outfile_name = "stdout";
  }
  if (rulestatsfile_name)
  {
    rulestatsfile = fopen(rulestatsfile_name, "w");
    if (rulestatsfile == NULL)
    {
      fprintf(stderr, "Failed to open \"%s\" for rule stats -- exit\n", rulestatsfile_name);
      if (outfile != stdout)
      {
	fclose(outfile);
      }
      if (yyin != stdin)
      {
        fclose(yyin);
      }
      exit(EXIT_FAILURE);
    }
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
    ulong_srandom(rndseed);
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
	if (process_diffusion)
	  lsys_string_diffusion(lstr);
	dstring = derived_string(lstr);
	if (rulestatsfile)
	{
	  fprint_rulestats(rulestatsfile, i, lstr);
	}
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
  if (rulestatsfile)
  {
    fclose(rulestatsfile);
  }
  if (yyin != stdin)
  {
    fclose(yyin);
  }
  if (outfile != stdout)
  {
    fclose(outfile);
  }
  return (EXIT_SUCCESS);
}

