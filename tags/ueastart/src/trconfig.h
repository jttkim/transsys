/* Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de> */

/*
 * $Id$
 *
 * $Log$
 * Revision 1.1  2005/03/08 17:12:02  jtk
 * Initial revision
 *
 * Revision 1.1  2001/04/04 11:12:00  kim
 * Initial addition of files previously not CVS managed
 *
 */

#ifndef TRCONFIG_H
#define TRCONFIG_H

#define IDENTIFIER_MAX 1024

#ifdef FREE_DEADBEEF
#  define free(x) free_deadbeef(x)
#endif

#endif /* TRCONFIG_H */

