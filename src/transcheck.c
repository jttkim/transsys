/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.1  2005/03/08 17:12:02  jtk
 * Initial revision
 *
 * Revision 1.5  2003/03/05 13:58:35  kim
 * added missing error messages upon failure exit
 *
 * Revision 1.4  2003/02/26 17:50:04  kim
 * fixed bug of returning success when parse errors occurred
 *
 * Revision 1.3  2003/02/04 23:48:12  kim
 * hack-fixed bug regarding ! by parentheses, no real fix!!
 *
 * Revision 1.2  2002/05/28 08:52:51  kim
 * added some discrete network stuff to sources
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
  const char *outfile_name = NULL, *discretenetfile_name = NULL;
  FILE *outfile, *discretenetfile = NULL;
  int yyreturn;
  TRANSSYS *tr;
  LSYS *ls;

  while ((oc = getopt(argc, argv, "d:yh")) != -1)
  {
    switch(oc)
    {
    case 'd':
      discretenetfile_name = optarg;
      break;
    case 'y':
      yydebug = 1;
      break;
    case 'h':
      printf("-d <filename>: specify name for discretenet topology output\n");
      printf("-y: generate parser debugging output by setting yydebug\n");
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
  if (discretenetfile_name && ((discretenetfile = fopen(discretenetfile_name, "w")) == NULL))
  {
    fprintf(stderr, "Failed to open \"%s\" for discrete network output -- exit\n", discretenetfile_name);
    if (yyin != stdin)
      fclose(yyin);
    if (outfile != stdout)
      fclose(outfile);
    exit(EXIT_FAILURE);
  }
  /* fprintf(stderr, "parsing from %s\n", yyin_name); */
  yyreturn = yyparse();
  if (yyreturn)
  {
    fprintf(stderr, "Failed to parse from \"%s\" -- exit\n", yyin_name);
    if (discretenetfile)
      fclose(discretenetfile);
    if (yyin != stdin)
      fclose(yyin);
    if (outfile != stdout)
      fclose(outfile);
    exit(EXIT_FAILURE);
  }
  for (tr = parsed_transsys; tr; tr = tr->next)
  {
    fprint_transsys(outfile, 0, tr);
    if (discretenetfile)
      fprint_transsys_as_discretenet(discretenetfile, tr);
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
  if (discretenetfile)
    fclose(discretenetfile);
  if (yyin != stdin)
    fclose(yyin);
  if (outfile != stdout)
    fclose(outfile);
  return (EXIT_SUCCESS);
}

