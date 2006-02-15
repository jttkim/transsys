/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 */

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "trconfig.h"
#include "trtypes.h"
#include "transsys.h"


static const TRANSSYS *find_transsys_by_name(const TRANSSYS *trlist, const char *trname)
{
  const TRANSSYS *tr;

  for (tr = trlist; tr; tr = tr->next)
  {
    if (!strcmp(trname, tr->name))
    {
      return (tr);
    }
  }
  return (NULL);
}


static int find_symbol_index_by_name(const LSYS *lsys, const char *symname)
{
  int s;

  for (s = 0; s < lsys->num_symbols; s++)
  {
    if (!strcmp(lsys->symbol_list[s].name, symname))
    {
      return (s);
    }
  }
  return (NO_INDEX);
}


/*
 * print a table line containing the symbol's index in the lsys string,
 * the symbol's name, and, where present, the transsys name and factor
 * concentrations.
 */
static int fprint_symbol_line(FILE *outfile, const LSYS_STRING *lstr, size_t symbol_instance_index)
{
  const SYMBOL_INSTANCE *si = lstr->symbol + symbol_instance_index;
  const SYMBOL_ELEMENT *symbol = si->lsys_string->lsys->symbol_list +si->symbol_index;
  fprintf(outfile, "%lu %s", (unsigned long) symbol_instance_index, symbol->name);
  if (symbol->transsys)
  {
    const TRANSSYS *transsys = symbol->transsys;
    const TRANSSYS_INSTANCE *ti = &(si->transsys_instance);
    size_t f;

    fprintf(outfile, " %s", transsys->name);
    for (f = 0; f < transsys->num_factors; f++)
    {
      fprintf(outfile, " %g", ti->factor_concentration[f]);
    }
  }
  fprintf(outfile, "\n");
  return (0);
}


/*
 * determine whether a the i-th symbolof lstr passes the filter for printing.
 */
static int symbol_filter(const LSYS_STRING *lstr, size_t i, const TRANSSYS *transsys, int symbol_index)
{
  if (transsys && (transsys != lstr->lsys->symbol_list[lstr->symbol[i].symbol_index].transsys))
  {
    return (0);
  }
  if ((symbol_index != NO_INDEX) && (lstr->symbol[i].symbol_index))
  {
    return (0);
  }
  return (1);
}


/*
 * print a table of lines, each line describing a symbol of the lsys string
 */
static int fprint_lstr_table(FILE *outfile, const LSYS_STRING *lstr, const TRANSSYS *transsys, int symbol_index)
{
  size_t i;

  for (i = 0; i < lstr->num_symbols; i++)
  {
    if (symbol_filter(lstr, i, transsys, symbol_index))
    {
      fprint_symbol_line(outfile, lstr, i);
    }
  }
  return (0);
}


static int ltransexpr(FILE *outfile, const LSYS *lsys, const TRANSSYS *transsys, const char *symbol_name, unsigned int rndseed, unsigned long num_timesteps, unsigned long output_period)
{
  LSYS_STRING *lstr, *dstr;
  unsigned long t;
  int symbol_index = -2;  /* initialiser pacifies -Wall */
  int return_value;

  if (symbol_name)
  {
    symbol_index = find_symbol_index_by_name(lsys, symbol_name);
    if (symbol_index == NO_INDEX)
    {
      fprintf(stderr, "no symbol \"%s\" in lsys \"%s\"\n", symbol_name, lsys->name);
      return (-1);
    }
  }
  else
  {
    symbol_index = NO_INDEX;
  }
  ulong_srandom(rndseed);
  lstr = axiom_string(lsys);
  for (t = 0; t < num_timesteps; t++)
  {
    return_value = lsys_string_expression(lstr);
    if (return_value != 0)
    {
      free_lsys_string(lstr);
      return (-1);
    }
    return_value = lsys_string_diffusion(lstr);
    if (return_value != 0)
    {
      free_lsys_string(lstr);
      return (-1);
    } 
    dstr = derived_string(lstr);
    if ((t % output_period) == 0)
    {
      fprint_lstr_table(outfile, lstr, transsys, symbol_index);
    }
    free_lsys_string(lstr);
    lstr = dstr;
  }
  return (0);
}


int main(int argc, char **argv)
{
  int oc;
  extern char *optarg;
  extern int optind;
  char *outfile_name = NULL;
  FILE *outfile = NULL;
  int yyreturn;
  unsigned long num_timesteps = 0, output_period = 1;
  const TRANSSYS *tr = NULL;
  LSYS *lsys;
  const char *transsys_name = NULL;
  int return_value = 0;
  unsigned int rndseed = 1;
  const char *symbol_name = NULL;

  while ((oc = getopt(argc, argv, "d:n:t:s:r:y:h")) != -1)
  {
    switch(oc)
    {
    case 'd':
      output_period = strtoul(optarg, NULL, 10);
      break;
    case 't':
      transsys_name = optarg;
      break;
    case 'n':
      num_timesteps = strtoul(optarg, NULL, 10);
      break;
    case 's':
      rndseed = strtoul(optarg, NULL, 10);
      break;
    case 'y':
      symbol_name = optarg;
      break;
    case 'h':
      printf("-c <filename>: specify gnuplot command file\n");
      printf("-d <intnum>: specify period of output (i.e. 0, d, 2*d etc. will be printed)\n");
      printf("-n <num>: specify number of time steps\n");
      printf("-t <name>: specify name of transsys to process\n");
      printf("-s <num>: specify random seed\n");
      printf("-y <name>: specify name of symbol to be monitored (lsys mode)\n");
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
  /* fprintf(stderr, "parsing from %s\n", yyin_name); */
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
  if (transsys_name)
  {
    tr = find_transsys_by_name(parsed_transsys, transsys_name);
    if (tr == NULL)
    {
      fprintf(stderr, "ltransexpr: no transsys \"%s\" found in file \"%s\"\n", transsys_name, yyin_name);
      return_value = -1;
    }
  }
  for (lsys = parsed_lsys; lsys; lsys = lsys->next)
  {
    return_value = ltransexpr(outfile, lsys, tr, symbol_name, rndseed, num_timesteps, output_period);
    if (return_value != 0)
    {
      break;
    }
  }
  free_lsys_list(parsed_lsys);
  free_transsys_list(parsed_transsys);
  if (yyin != stdin)
    fclose(yyin);
  if (outfile != stdout)
    fclose(outfile);
  if (return_value != 0)
  {
    exit(EXIT_FAILURE);
  }
  return (EXIT_SUCCESS);
}

