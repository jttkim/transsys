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
 * Revision 1.1  2002/12/09 19:37:55  kim
 * added dot rendering of transsys network graphs (initial version)
 *
 */

#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>

#include "trconfig.h"
#include "trtypes.h"
#include "transsys.h"


char *prgname = "transdot";


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
  const char *transsys_name = NULL;

  while ((oc = getopt(argc, argv, "t:f:eh")) != -1)
  {
    switch(oc)
    {
    case 't':
      transsys_name = optarg;
      break;
    case 'h':
      printf("transdot -- render transsys in dot language\n");
      printf("usage: %s [options] <infile> <outfile>\n", argv[0]);
      printf("options:\n");
      printf("-t <name>: specify transsys by name\n");
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
      free_transsys_list(parsed_transsys);
      if (yyin != stdin)
	fclose(yyin);
      if (outfile != stdout)
	fclose(outfile);
      exit(EXIT_FAILURE);
    }
  }
  else
    trlist = parsed_transsys;
  for (tr = trlist; tr; tr = tr->next)
  {
    if (transsys_name && strcmp(tr->name, transsys_name))
      continue;
    dot_transsys(outfile, tr);
    if (transsys_name)
      break;
  }
  free_transsys_list(parsed_transsys);
  if (yyin != stdin)
    fclose(yyin);
  if (outfile != stdout)
    fclose(outfile);
  return (EXIT_SUCCESS);
}

