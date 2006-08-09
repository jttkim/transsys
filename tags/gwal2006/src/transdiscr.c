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
 * Revision 1.1  2002/12/09 12:59:20  kim
 * added (again) transdiscr.c
 *
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
  const char *outfile_name = NULL;
  FILE *outfile;
  int yyreturn;
  TRANSSYS *tr;
  LSYS *ls;

  while ((oc = getopt(argc, argv, "h")) != -1)
  {
    switch(oc)
    {
    case 'h':
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
  for (tr = parsed_transsys; tr; tr = tr->next)
  {
    fprint_transsys_as_discretenet(outfile, tr);
  }
  for (ls = parsed_lsys; ls; ls = ls->next)
  {
    fprintf(stderr, "lsys \"%s\"\n", ls->name);
  }
  for (ls = parsed_lsys; ls; ls = ls->next)
  {
    fprint_lsys(outfile, 0, ls);
  }
  free_lsys_list(parsed_lsys);
  free_transsys_list(parsed_transsys);
  if (yyin != stdin)
    fclose(yyin);
  if (outfile != stdout)
    fclose(outfile);
  return (EXIT_SUCCESS);
}

