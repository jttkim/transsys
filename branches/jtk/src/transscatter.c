/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.1  2005/03/08 17:12:02  jtk
 * Initial revision
 *
 * Revision 1.4  2003/03/05 13:58:35  kim
 * added missing error messages upon failure exit
 *
 * Revision 1.3  2003/02/26 17:50:04  kim
 * fixed bug of returning success when parse errors occurred
 *
 * Revision 1.2  2002/01/25 03:35:04  kim
 * Added gnuplot link functionality to transexpr, transscatter, improved
 *     PostScript stuff
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


int fprint_plotcommands(FILE *f, const TRANSSYS *tr, const char *outfile_name, int extensiveness)
{
  int i, j;

  fprintf(f, "# transsys \"%s\"\n", tr->name);
  for (i = 0; i < tr->num_factors; i++)
  {
    fprintf(f, "plot \'%s\' using 1:%d title \'%s:%s\'\n",
	    outfile_name, i + 3, tr->name, tr->factor_list[i].name);
    fprintf(f, "pause -1 \'Hit return\'\n");
  }
  if (extensiveness == 0)
    return (0);
  for (i = 0; i < tr->num_factors; i++)
  {
    for (j = i + 1; j < tr->num_factors; j++)
    {
      fprintf(f, "plot \'%s\' using %d:%d title \'%s:%s/%s\'\n",
	      outfile_name, i + 3, j + 3, tr->name, tr->factor_list[i].name, tr->factor_list[j].name);
      fprintf(f, "pause -1 \'Hit return\'\n");
    }
  }
  return (0);
}


/*
 * create a value that can be represented by printf -- lest we get lost
 * in some Lorenz sensitivity...
 */

double printable_random_value(double maxval)
{
  char buf[50];

  sprintf(buf, "%1.12g", urandom_double() * maxval);
  return (strtod(buf, NULL));
}


int main(int argc, char **argv)
{
  int oc;
  extern char *optarg;
  extern int optind;
  char *outfile_name = NULL, *plotcmdfile_name = NULL;
  FILE *outfile = NULL, *plotcmdfile = NULL;
  int yyreturn;
  int num_timesteps = 100, num_restarts = 100, r, i;
  TRANSSYS *tr;
  TRANSSYS_INSTANCE ti;
  const char *transsys_name = NULL;
  unsigned long t;
  double max_concentration = 1.0;
  int plot_extensiveness = 0;

  while ((oc = getopt(argc, argv, "c:m:n:r:t:h")) != -1)
  {
    switch(oc)
    {
    case 'c':
      plotcmdfile_name = optarg;
      break;
    case 'm':
      max_concentration = strtod(optarg, NULL);
      break;
    case 'n':
      num_timesteps = strtol(optarg, NULL, 10);
      break;
    case 'r':
      num_restarts = strtol(optarg, NULL, 10);
      break;
    case 't':
      transsys_name = optarg;
      break;
    case 'h':
      printf("-c <filename>: specify gnuplot command file\n");
      printf("-m <num>: specify maximum of randomly generated initial concenctration\n");
      printf("-n <num>: specify number of time steps\n");
      printf("-r <num>: specify number of restarts\n");
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
  init_transsys_instance_components(&ti);
  for (tr = parsed_transsys; tr; tr = tr->next)
  {
    if (transsys_name && strcmp(tr->name, transsys_name))
      continue;
    fprintf(outfile, "# transsys \"%s\"\n", tr->name);
    if (plotcmdfile)
      fprint_plotcommands(plotcmdfile, tr, outfile_name, plot_extensiveness);
    alloc_transsys_instance_components(&ti, tr);
    for (r = 0; r < num_restarts; r++)
    {
      fprintf(outfile ,"# init:");
      for (i = 0; i < ti.transsys->num_factors; i++)
      {
	ti.factor_concentration[i] = printable_random_value(max_concentration);
	fprintf(outfile, " %1.12g", ti.factor_concentration[i]);
      }
      fprintf(outfile, "\n");
      for (t = 0; t < num_timesteps; t++)
      {
	/* fprint_factorconc_line(outfile, cell); */
	process_expression(&ti);
      }
      fprintf(outfile, "%d ", r);
      fprint_factorconc_line(outfile, &ti, t);
      /* fprint_cell(stderr, cell); */
    }
    free_transsys_instance_components(&ti);
    if (transsys_name && !strcmp(tr->name, transsys_name))
      break;
  }
  free_transsys_list(parsed_transsys);
  if (yyin != stdin)
    fclose(yyin);
  if (outfile != stdout)
    fclose(outfile);
  return (EXIT_SUCCESS);
}

