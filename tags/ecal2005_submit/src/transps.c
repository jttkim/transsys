/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.2  2005/03/29 17:33:02  jtk
 * introduced arrayed lsys string, with symbol distance matrix.
 *
 * Revision 1.1.1.1  2005/03/08 17:12:02  jtk
 * new cvs after loss at INB
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
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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
  int eps_mode = 0;
  int yyreturn;
  TRANSSYS *trlist, *tr;
  POSTSCRIPT_STYLE style;
  const char *transsys_name = NULL;
  double fontheight = -1.0;

  while ((oc = getopt(argc, argv, "t:f:eh")) != -1)
  {
    switch(oc)
    {
    case 't':
      transsys_name = optarg;
      break;
    case 'e':
      eps_mode = 1;
      break;
    case 'f' :
      fontheight = strtod(optarg, NULL);
      break;
    case 'h':
      printf("transps -- make PostScript representation of a transsys\n");
      printf("usage: %s [options] <infile> <outfile>\n", argv[0]);
      printf("options:\n");
      printf("-t <name>: specify transsys by name\n");
      printf("-f <num>: specify font height\n");
      printf("-e: output encapsulated PostScript (EPS)\n");
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
  if (transsys_name)
  {
    for (trlist = parsed_transsys; trlist; trlist = trlist->next)
    {
      if (!strcmp(trlist->name, transsys_name))
	break;
    }
    if (trlist == NULL)
    {
      fprintf(stderr, "no transsys \"%s\" in file \"%s\" -- exit\n", transsys_name, yyin_name);
      if (yyin != stdin)
	fclose(yyin);
      if (outfile != stdout)
	fclose(outfile);
      exit(EXIT_FAILURE);
    }
  }
  else
    trlist = parsed_transsys;
  if (eps_mode)
  {
    if (!transsys_name && trlist->next)
      fprintf(stderr, "EPS mode -- printing only transsys \"%s\"\n", trlist->name);
    if (fontheight >= 0.0)
      style.fontheight = fontheight;
    set_default_postscript_style(1, trlist, &style);
    postscript_transsys_init(outfile, trlist->name, &style);
    postscript_transsys(outfile, trlist, &style);
    postscript_transsys_finish(outfile, &style);
  }
  else
  {
    set_default_postscript_style(0, NULL, &style);
    if (fontheight >= 0.0)
      style.fontheight = fontheight;
    if (transsys_name)
    {
      postscript_transsys_init(outfile, transsys_name, &style);
    }
    else
    {
      postscript_transsys_init(outfile, yyin_name, &style);
    }
    for (tr = trlist; tr; tr = tr->next)
    {
      if (transsys_name && strcmp(tr->name, transsys_name))
	continue;
      style.page_num++;
      postscript_transsys(outfile, tr, &style);
      if (transsys_name)
	break;
    }
    postscript_transsys_finish(outfile, &style);
  }
  free_transsys_list(parsed_transsys);
  if (yyin != stdin)
    fclose(yyin);
  if (outfile != stdout)
    fclose(outfile);
  return (EXIT_SUCCESS);
}

