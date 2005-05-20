/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.9  2005/05/20 10:40:15  jtk
 * differentiated entropy recording for lsys, expression and diffusion phase
 *
 * Revision 1.8  2005/04/14 20:43:15  jtk
 * entropy bugfix: log(2) -> 1.0 / log(2), changed entropy plot command to lines
 *
 * Revision 1.7  2005/04/14 18:50:15  jtk
 * added entropy computation, introduced magic line for transexpr
 *
 * Revision 1.6  2005/04/08 20:04:28  jtk
 * added missing options in help output
 *
 * Revision 1.5  2005/04/05 10:12:39  jtk
 * made diffusion consistent (no oscillation due to overshooting), small fixes
 *
 * Revision 1.4  2005/04/04 09:39:54  jtk
 * added lsys capabilities to transexpr, various small changes
 *
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
  return (-1);
}


int fprint_plotcommands(FILE *f, const TRANSSYS *tr, const char *outfile_name, int extensiveness)
{
  int i, j;

  fprintf(f, "# transsys \"%s\"\n", tr->name);
  for (i = 0; i < tr->num_factors; i++)
  {
    fprintf(f, "plot \'%s\' using 1:%d:%d title \'%s:%s\' with errorbars\n",
	    outfile_name, i * 3 + 3, i * 3 + 4, tr->name, tr->factor_list[i].name);
    fprintf(f, "pause -1 \'Hit return\'\n");
    fprintf(f, "plot \'%s\' using 1:%d title \'%s:%s, entropy\' with lines\n", outfile_name, i * 3 + 5, tr->name, tr->factor_list[i].name);
    fprintf(f, "pause -1 \'Hit return\'\n");
  }
  if (extensiveness == 0)
    return (0);
  for (i = 0; i < tr->num_factors; i++)
  {
    for (j = i + 1; j < tr->num_factors; j++)
    {
      fprintf(f, "plot \'%s\' using %d:%d title \'%s:%s/%s\'\n",
	      outfile_name, i * 2 + 3, j * 2 + 3, tr->name, tr->factor_list[i].name, tr->factor_list[j].name);
      fprintf(f, "pause -1 \'Hit return\'\n");
    }
  }
  return (0);
}


static int factor_concentrations_from_string(TRANSSYS_INSTANCE *ti, const char *str)
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


static int fprint_null_expression_line(FILE *outfile, const TRANSSYS *transsys, unsigned long t)
{
  int i;

  fprintf(outfile, "%lu 0", t);
  for (i = 0; i < transsys->num_factors; i++)
  {
    fprintf(outfile, "  0.0 0.0");
  }
  fprintf(outfile, "\n");
  return (0);
}


/*
 * expression_record_version *must* be incremented each time the record
 * format changes *semantically*. It should not be changed on other
 * occasions.
 * The main use of this versioning is to provide tools like transsys.py
 * (TranssysInstance.time_series) with a means to detect whether they
 * are working on a compatible numeric output produced by transexpr.
 *
 * History:
 * 1.0: Introduction of versioning. Record format is
 *     - 0: time step
 *     - 1: number of transsys instances
 *     - 2 + 3 * i: average concentration of factor i
 *     - 3 + 3 * i: standard deviation of concentration of factor i
 *     - 4 + 3 * i: Shannon entropy of factor i.
 */

static const char *expression_record_version = "1.0";

static void fprint_expression_header(FILE *f, const TRANSSYS *transsys, const LSYS *lsys)
{
  int i;
  time_t t;

  fprintf(f, "# transsys expression records %s\n", expression_record_version);
  fprintf(f, "# time n.instances");
  for (i = 0; i < transsys->num_factors; i++)
  {
    fprintf(f, "  %s.avg %s.stddev %s.entropy", transsys->factor_list[i].name, transsys->factor_list[i].name, transsys->factor_list[i].name);
  }
  fprintf(f, "\n");
  t = time(NULL);
  fprintf(f, "# transsys %s\n", transsys->name);
  if (lsys)
  {
    fprintf(f, "# lsys %s\n", lsys->name);
  }
  else
  {
    fprintf(f, "# no lsys\n");
  }
  fprintf(f, "# run on %s", ctime(&t));
}


/*
 * Remember to change the transsys_record_version on *each* change of
 * the record format.
 */

int fprint_expression_record(FILE *outfile, TRANSSYS_INSTANCE **ti, unsigned long t)
{
  size_t i, f, n;
  double d, *average, *stddev, *entropy;
  const TRANSSYS *transsys = ti[0]->transsys;

  average = (double *) malloc(2 * transsys->num_factors * sizeof(double));
  if (average == NULL)
  {
    fprintf(stderr, "fprint_expression_record: malloc failed\n");
    return (-1);
  }
  stddev = average + transsys->num_factors;
  entropy = transsys_collection_factor_entropy(ti);
  if (entropy == NULL)
  {
    fprintf(stderr, "fprint_expression_record: transsys_collection_factor_entropy failed\n");
    free(average);
    return (-1);
  }
  for (i = 0; i < transsys->num_factors; i++)
  {
    average[i] = 0.0;
    stddev[i] = 0.0;
  }
  for (n = 0; ti[n]; n++)
  {
    for (f = 0; f < transsys->num_factors; f++)
    {
      average[f] += ti[n]->factor_concentration[f];
    }
  }
  for (f = 0; f < transsys->num_factors; f++)
  {
    average[f] /= n;
  }
  if (n > 1)
  {
    for (i = 0; i < n; i++)
    {
      for (f = 0; f < transsys->num_factors; f++)
      {
	d = ti[i]->factor_concentration[f] - average[f];
	stddev[f] += d * d;
      }
    }
    for (f = 0; f < transsys->num_factors; f++)
    {
      stddev[f] = sqrt(stddev[f]) / (n - 1);
    }
  }
  fprintf(outfile, "%lu %lu", t, (unsigned long) n);
  for (f = 0; f < transsys->num_factors; f++)
  {
    fprintf(outfile, "  %g %g %g", average[f], stddev[f], entropy[f]);
  }
  fprintf(outfile, "\n");
  free(average);
  free(entropy);
  return (0);
}


void free_ti_array(TRANSSYS_INSTANCE **ti)
{
  size_t i;

  for (i = 0; ti[i]; i++)
  {
    free_transsys_instance_components(ti[i]);
  }
  free(ti[0]);
  free(ti);
}


TRANSSYS_INSTANCE **alloc_ti_array(const TRANSSYS *transsys, size_t n)
{
  TRANSSYS_INSTANCE **ti;
  size_t i, j;

  ti = (TRANSSYS_INSTANCE **) malloc((n + 1) * sizeof(TRANSSYS_INSTANCE *));
  if (ti == NULL)
  {
    return (NULL);
  }
  ti[0] = (TRANSSYS_INSTANCE *) malloc(n * sizeof(TRANSSYS_INSTANCE));
  if (ti[0] == NULL)
  {
    free(ti);
    return (NULL);
  }
  for (i = 0; i < n; i++)
  {
    ti[i] = ti[0] + i;
    init_transsys_instance_components(ti[i]);
    if (alloc_transsys_instance_components(ti[i], transsys) != 0)
    {
      for (j = 0; j < i; j++)
      {
	free_transsys_instance_components(ti[j]);
	free(ti[0]);
	free(ti);
	return (NULL);
      }
    }
  }
  ti[n] = NULL;
  return (ti);
}


int init_ti_array(TRANSSYS_INSTANCE **ti, double factorconc_init, const char *factorconc_init_string)
{
  size_t r, f;

  for (r = 0; ti[r]; r++)
  {
    for (f = 0; f < ti[r]->transsys->num_factors; f++)
    {
      ti[r]->factor_concentration[f] = factorconc_init;
    }
  }
  if (factorconc_init_string)
  {
    for (r = 0; ti[r]; r++)
    {
      if (factor_concentrations_from_string(ti[r], factorconc_init_string))
      {
	fprintf(stderr, "factor_concentrations_from_string: error\n");
	fprintf(stderr, "    \"%s\"\n", factorconc_init_string);
	return (-1);
      }
    }
  }
  return (0);
}


static int transexpr(FILE *outfile, const TRANSSYS *transsys, unsigned int rndseed, int num_repeats, unsigned long num_timesteps, unsigned long output_period, double factorconc_init, const char *factorconc_init_string, FILE *plotcmdfile, const char *outfile_name, int plot_extensiveness)
{
  TRANSSYS_INSTANCE **ti;
  unsigned long t;
  int r;

  if (num_repeats == 0)
  {
    fprintf(stderr, "transexpr: pointless to run with 0 repeats\n");
    return (-1);
  }
  ulong_srandom(rndseed);
  ti = alloc_ti_array(transsys, num_repeats);
  if (ti == NULL)
  {
    fprintf(stderr, "transexpr: alloc_ti_array failed\n");
    return (-1);
  }
  if (init_ti_array(ti, factorconc_init, factorconc_init_string) != 0)
  {
    free_ti_array(ti);
    return (-1);
  }
  if (plotcmdfile)
  {
    fprint_plotcommands(plotcmdfile, transsys, outfile_name, plot_extensiveness);
  }
  fprint_expression_header(outfile, transsys, NULL);
  for (t = 0; t < num_timesteps; t++)
  {
    if ((t % output_period) == 0)
    {
      fprint_expression_record(outfile, ti, t);
    }
    for (r = 0; r < num_repeats; r++)
    {
      process_expression(ti[r]);
    }
  }
  /* fprint_transsys_instance(stderr, &ti); */
  free_ti_array(ti);
  return (0);
}


/*
 * compute a NULL-terminated array of pointers to transsys instances
 * in an lsys string. All transsys instances are instances of transsys.
 * if symbol_index is non-negative, only instances of symbol specified
 * by symbol_index are included.
 */

static TRANSSYS_INSTANCE **lsys_transsys_instance_array(const LSYS_STRING *lstr, const TRANSSYS *transsys, int symbol_index)
{
  TRANSSYS_INSTANCE **ti;
  size_t n, i;

  /* FIXME: allocating this large an array may be inefficient*/
  /* ... but with the perspective to require all symbols to have the same
     transsys instance, it's not that bad after all... */
  ti = (TRANSSYS_INSTANCE **) malloc((lstr->num_symbols + 1) * sizeof(TRANSSYS_INSTANCE *));
  if (ti == NULL)
  {
    fprintf(stderr, "lsys_transsys_instance_list: malloc failed\n");
    return (NULL);
  }
  n = 0;
  for (i = 0; i < lstr->num_symbols; i++)
  {
    if (lstr->symbol[i].transsys_instance.transsys == transsys)
    {
      if ((symbol_index < 0) || (lstr->symbol[i].symbol_index == symbol_index))
      {
	ti[n++] = &(lstr->symbol[i].transsys_instance);
      }
    }
  }
  ti[n] = NULL;
  return (ti);
}


static int transexpr_lsys(FILE *outfile, FILE *entropyfile, const LSYS *lsys, const TRANSSYS *transsys, const char *symbol_name, unsigned int rndseed, unsigned long num_timesteps, unsigned long output_period, FILE *plotcmdfile, const char *outfile_name, int plot_extensiveness)
{
  LSYS_STRING *lstr, *dstr;
  TRANSSYS_INSTANCE **ti = NULL;
  double *entropy_lsys = NULL, *entropy_expression = NULL, *entropy_diffusion = NULL;
  unsigned long t;
  size_t n, i, f;
  int symbol_index = -2;  /* initialiser pacifies -Wall */
  int return_value;

  if (symbol_name)
  {
    symbol_index = find_symbol_index_by_name(lsys, symbol_name);
    if (symbol_index == -1)
    {
      fprintf(stderr, "no symbol \"%s\" in lsys \"%s\"\n", symbol_name, lsys->name);
      return (-1);
    }
  }
  if (plotcmdfile)
  {
    fprint_plotcommands(plotcmdfile, transsys, outfile_name, plot_extensiveness);
  }
  ulong_srandom(rndseed);
  lstr = axiom_string(lsys);
  fprint_expression_header(outfile, transsys, lsys);
  for (t = 0; t < num_timesteps; t++)
  {
    if ((t % output_period) == 0)
    {
      ti = lsys_transsys_instance_array(lstr, transsys, symbol_index);
      if (ti == NULL)
      {
	fprintf(stderr, "transexpr_lsys: lsys_transsys_instance_list failed\n");
	free_lsys_string(lstr);
	return (-1);
      }
      n = 0;
      for (i = 0; i < lstr->num_symbols; i++)
      {
	if (lstr->symbol[i].transsys_instance.transsys == transsys)
	{
	  if ((symbol_name == NULL) || (lstr->symbol[i].symbol_index == symbol_index))
	  {
	    ti[n++] = &(lstr->symbol[i].transsys_instance);
	  }

	}
      }
      if (n > 0)
      {
	fprint_expression_record(outfile, ti, t);
      }
      else
      {
	fprint_null_expression_line(outfile, transsys, t);
      }
      if (entropyfile)
      {
	entropy_lsys = transsys_collection_factor_entropy(ti);
      }
    }
    return_value = lsys_string_expression(lstr);
    if (return_value != 0)
    {
      if (entropy_lsys)
      {
	free(entropy_lsys);
      }
      free_lsys_string(lstr);
      return (-1);
    }
    if (entropyfile && (t % output_period) == 0)
    {
      entropy_expression = transsys_collection_factor_entropy(ti);
    }
    return_value = lsys_string_diffusion(lstr);
    if (return_value != 0)
    {
      if (entropy_lsys)
      {
	free(entropy_lsys);
      }
      if (entropy_expression)
      {
	free(entropy_expression);
      }
      free_lsys_string(lstr);
      return (-1);
    } 
    if (entropyfile && (t % output_period) == 0)
    {
      entropy_diffusion = transsys_collection_factor_entropy(ti);
      fprintf(entropyfile, "%lu", t);
      for (f = 0; f < transsys->num_factors; f++)
      {
	fprintf(entropyfile, " %g %g %g", entropy_lsys[f], entropy_expression[f], entropy_diffusion[f]);
      }
      fprintf(entropyfile, "\n");
      free(entropy_lsys);
      entropy_lsys = NULL;
      free(entropy_expression);
      entropy_expression = NULL;
      free(entropy_diffusion);
      entropy_diffusion = NULL;
    }
    dstr = derived_string(lstr);
    free_lsys_string(lstr);
    if (ti)
    {
      free(ti);
      ti = NULL;
    }
    lstr = dstr;
  }
  return (0);
}


int main(int argc, char **argv)
{
  int oc;
  extern char *optarg;
  extern int optind;
  char *outfile_name = NULL, *plotcmdfile_name = NULL, *entropyfile_name = NULL;
  FILE *outfile = NULL, *plotcmdfile = NULL, *entropyfile = NULL;
  int yyreturn;
  unsigned long num_timesteps = 100, output_period = 1;
  const TRANSSYS *tr;
  LSYS *lsys;
  const char *transsys_name = NULL, *factorconc_init_string = NULL;
  double factorconc_init = 0.0;
  int plot_extensiveness = 0;
  int lsys_mode = 0;
  int return_value = 0;
  unsigned int rndseed = 1;
  int num_repeats = 1;
  const char *symbol_name = NULL;

  while ((oc = getopt(argc, argv, "c:e:d:F:f:n:t:s:r:y:lh")) != -1)
  {
    switch(oc)
    {
    case 'l':
      lsys_mode = 1;
      break;
    case 'c':
      plotcmdfile_name = optarg;
      break;
    case 'e':
      entropyfile_name = optarg;
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
    case 's':
      rndseed = strtoul(optarg, NULL, 10);
      break;
    case 'r':
      num_repeats = strtol(optarg, NULL, 10);
      break;
    case 'y':
      symbol_name = optarg;
      break;
    case 'h':
      printf("-l: run in lsys mode\n");
      printf("-c <filename>: specify gnuplot command file\n");
      printf("-f <num>: specify uniform initial factor concentration\n");
      printf("-d <intnum>: specify period of output (i.e. 0, d, 2*d etc. will be printed)\n");
      printf("-F <string>: specify initial factor concentrations (ws separated)\n");
      printf("-n <num>: specify number of time steps\n");
      printf("-r <num>: specify number of repeats\n");
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
  if (entropyfile_name)
  {
    if (lsys_mode)
    {
      entropyfile = fopen(entropyfile_name, "w");
      if (entropyfile == NULL)
      {
	fprintf(stderr, "failed to open entropy file \"%s\" -- exit\n", entropyfile_name);
	if (yyin != stdin)
	  fclose(yyin);
	if (outfile != stdout)
	  fclose(outfile);
	exit(EXIT_FAILURE);
      }
    }
    else
    {
      fprintf(stderr, "Not running in lsys mode: detailed entropy tracing (-e) ignored\n");
    }
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
      if (entropyfile)
	fclose(entropyfile);
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
	  return_value = transexpr_lsys(outfile, entropyfile, lsys, tr, symbol_name, rndseed, num_timesteps, output_period, plotcmdfile, outfile_name, plot_extensiveness);
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
	transexpr(outfile, tr, rndseed, num_repeats, num_timesteps, output_period, factorconc_init, factorconc_init_string, plotcmdfile, outfile_name, plot_extensiveness);
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
	return_value = transexpr(outfile, tr, rndseed, num_repeats, num_timesteps, output_period, factorconc_init, factorconc_init_string, plotcmdfile, outfile_name, plot_extensiveness);
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
  if (entropyfile)
    fclose(entropyfile);
  if (return_value != 0)
  {
    exit(EXIT_FAILURE);
  }
  return (EXIT_SUCCESS);
}

