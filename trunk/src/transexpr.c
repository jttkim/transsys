/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.3  2005/04/01 18:16:34  jtk
 * proper exit codes for transexpr
 *
 * Revision 1.2  2005/03/29 17:33:02  jtk
 * introduced arrayed lsys string, with symbol distance matrix.
 *
 * Revision 1.1.1.1  2005/03/08 17:12:02  jtk
 * new cvs after loss at INB
 *
 * Revision 1.6  2003/03/05 13:58:35  kim
 * added missing error messages upon failure exit
 *
 * Revision 1.5  2003/02/26 17:50:04  kim
 * fixed bug of returning success when parse errors occurred
 *
 * Revision 1.4  2002/05/28 08:52:51  kim
 * added some discrete network stuff to sources
 *
 * Revision 1.3  2002/01/25 03:35:03  kim
 * Added gnuplot link functionality to transexpr, transscatter, improved
 *     PostScript stuff
 *
 * Revision 1.2  2001/04/05 17:24:50  kim
 * Added initialization of transsys_name = NULL
 *
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

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


int fprint_plotcommands(FILE *f, const TRANSSYS *tr, const char *outfile_name, int extensiveness)
{
  int i, j;

  fprintf(f, "# transsys \"%s\"\n", tr->name);
  for (i = 0; i < tr->num_factors; i++)
  {
    fprintf(f, "plot \'%s\' using 1:%d title \'%s:%s\' with lines\n",
	    outfile_name, i + 2, tr->name, tr->factor_list[i].name);
    fprintf(f, "pause -1 \'Hit return\'\n");
  }
  if (extensiveness == 0)
    return (0);
  for (i = 0; i < tr->num_factors; i++)
  {
    for (j = i + 1; j < tr->num_factors; j++)
    {
      fprintf(f, "plot \'%s\' using %d:%d title \'%s:%s/%s\'\n",
	      outfile_name, i + 2, j + 2, tr->name, tr->factor_list[i].name, tr->factor_list[j].name);
      fprintf(f, "pause -1 \'Hit return\'\n");
    }
  }
  return (0);
}


int factor_concentrations_from_string(TRANSSYS_INSTANCE *ti, const char *str)
{
  char *s;
  int i;

  for (i = 0; i < ti->transsys->num_factors; i++)
  {
    while (isspace(*str))
      str++;
    if (!*str)
      break;
    ti->factor_concentration[i] = strtod(str, &s);
    if (str == s)
      return (-1);
    str = s;
  }
  return (0);
}


static int transexpr(FILE *outfile, const TRANSSYS *transsys, unsigned long num_timesteps, unsigned long output_period, double factorconc_init, const char *factorconc_init_string, FILE *plotcmdfile, const char *outfile_name, int plot_extensiveness)
{
  TRANSSYS_INSTANCE ti;
  unsigned long t;
  int i;

  init_transsys_instance_components(&ti);
  alloc_transsys_instance_components(&ti, transsys);
  if (plotcmdfile)
  {
    fprint_plotcommands(plotcmdfile, transsys, outfile_name, plot_extensiveness);
  }
  for (i = 0; i < ti.transsys->num_factors; i++)
  {
    ti.factor_concentration[i] = factorconc_init;
  }
  if (factorconc_init_string)
  {
    if (factor_concentrations_from_string(&ti, factorconc_init_string))
    {
      fprintf(stderr, "factor_concentrations_from_string: error\n");
      fprintf(stderr, "    \"%s\"\n", factorconc_init_string);
    }
  }
  fprint_factorconc_commentline(outfile, &ti);
  for (t = 0; t < num_timesteps; t++)
  {
    if ((t % output_period) == 0)
      fprint_factorconc_line(outfile, &ti, t);
    process_expression(&ti);
  }
  /* fprint_transsys_instance(stderr, &ti); */
  free_transsys_instance_components(&ti);
  return (0);
}


int fprint_lsys_average_expression_line(FILE *outfile, const LSYS_STRING *lstr, const TRANSSYS *transsys, unsigned long t)
{
  /*
   * this transsys instance is (ab-) used to accumulate the averages
   * in factor_concentration and standard deviations in new_concentration
   */
  TRANSSYS_INSTANCE ti;
  int i, f, n;
  double d;
  
  init_transsys_instance_components(&ti);
  alloc_transsys_instance_components(&ti, transsys);
  n = 0;
  for (i = 0; i < lstr->num_symbols; i++)
  {
    if (lstr->symbol[i].transsys_instance.transsys == transsys)
    {
      n++;
      for (f = 0; f < transsys->num_factors; f++)
      {
	ti.factor_concentration[f] += lstr->symbol[i].transsys_instance.factor_concentration[f];
      }
    }
  }
  for (f = 0; f < transsys->num_factors; f++)
  {
    ti.factor_concentration[f] /= n;
  }
  for (i = 0; i < lstr->num_symbols; i++)
  {
    if (lstr->symbol[i].transsys_instance.transsys == transsys)
    {
      for (f = 0; f < transsys->num_factors; f++)
      {
	d = lstr->symbol[i].transsys_instance.factor_concentration[f] - ti.factor_concentration[f];
	ti.new_concentration[f] += d * d;
      }
    }
  }
  for (f = 0; f < transsys->num_factors; f++)
  {
    ti.new_concentration[f] = sqrt(ti.new_concentration[f]) / (n - 1);
  }
  fprintf(outfile, "%lu", t);
  for (f = 0; f < transsys->num_factors; f++)
  {
    fprintf(outfile, "  %g %g", ti.factor_concentration[f], ti.new_concentration[f]);
  }
  fprintf(outfile, "\n");
  free_transsys_instance_components(&ti);
  return (0);
}


static int transexpr_lsys(FILE *outfile, const LSYS *lsys, const TRANSSYS *transsys, unsigned long num_timesteps, unsigned long output_period, FILE *plotcmdfile, const char *outfile_name, int plot_extensiveness)
{
  LSYS_STRING *lstr, *dstr;
  unsigned long t;
  int return_value;

  lstr = axiom_string(lsys);
  for (t = 0; t < num_timesteps; t++)
  {
    if ((t % output_period) == 0)
    {
      fprint_lsys_average_expression_line(outfile, lstr, transsys, t);
    }
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
  char *outfile_name = NULL, *plotcmdfile_name = NULL;
  FILE *outfile = NULL, *plotcmdfile = NULL;
  int yyreturn;
  unsigned long num_timesteps = 100, output_period = 1;
  const TRANSSYS *tr;
  LSYS *lsys;
  const char *transsys_name = NULL, *factorconc_init_string = NULL;
  double factorconc_init = 0.0;
  int plot_extensiveness = 0;
  int lsys_mode = 0;
  int return_value = 0;

  while ((oc = getopt(argc, argv, "c:d:F:f:n:t:lh")) != -1)
  {
    switch(oc)
    {
    case 'l':
      lsys_mode = 1;
      break;
    case 'c':
      plotcmdfile_name = optarg;
      break;
    case 'd':
      output_period = strtoul(optarg, NULL, 10);
      break;
    case 'f' :
      factorconc_init = strtod(optarg, NULL);
      break;
    case 'F' :
      factorconc_init_string = optarg;
      break;
    case 't':
      transsys_name = optarg;
      break;
    case 'n':
      num_timesteps = strtoul(optarg, NULL, 10);
      break;
    case 'h':
      printf("-l: run in lsys mode\n");
      printf("-c <filename>: specify gnuplot command file\n");
      printf("-f <num>: specify uniform initial factor concentration\n");
      printf("-d <intnum>: specify period of output (i.e. 0, d, 2*d etc. will be printed)\n");
      printf("-F <string>: specify initial factor concentrations (ws separated)\n");
      printf("-n <num>: specify number of time steps\n");
      printf("-t <name>: specify name of transsys to process\n");
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
  if (plotcmdfile_name)
  {
    plotcmdfile = fopen(plotcmdfile_name, "w");
    if (plotcmdfile == NULL)
    {
      fprintf(stderr, "failed to open plot command file \"%s\" -- exit\n", plotcmdfile_name);
      if (yyin != stdin)
	fclose(yyin);
      if (outfile != stdout)
	fclose(outfile);
      exit(EXIT_FAILURE);
    }
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
  if (lsys_mode)
  {
    if (transsys_name)
    {
      tr = find_transsys_by_name(parsed_transsys, transsys_name);
      if (tr)
      {
	for (lsys = parsed_lsys; lsys; lsys = lsys->next)
	{
	  return_value = transexpr_lsys(outfile, lsys, tr, num_timesteps, output_period, plotcmdfile, outfile_name, plot_extensiveness);
	  if (return_value != 0)
	  {
	    break;
	  }
	}
      }
      else
      {
	fprintf(stderr, "transexpr: no transsys \"%s\" found in file \"%s\"\n", transsys_name, yyin_name);
	return_value = -1;
      }
    }
    else
    {
      fprintf(stderr, "transexpr: a transsys must be specified by name in lsys mode\n");
      return_value = -1;
    }
  }
  else
  {
    if (transsys_name)
    {
      tr = find_transsys_by_name(parsed_transsys, transsys_name);
      if (tr)
      {
	transexpr(outfile, tr, num_timesteps, output_period, factorconc_init, factorconc_init_string, plotcmdfile, outfile_name, plot_extensiveness);
      }
      else
      {
	fprintf(stderr, "transexpr: no transsys \"%s\" found in file \"%s\"\n", transsys_name, yyin_name);
	return_value = -1;
      }
    }
    else
    {
      for (tr = parsed_transsys; tr; tr = tr->next)
      {
	return_value = transexpr(outfile, tr, num_timesteps, output_period, factorconc_init, factorconc_init_string, plotcmdfile, outfile_name, plot_extensiveness);
	if (return_value != 0)
	{
	  break;
	}
      }
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

