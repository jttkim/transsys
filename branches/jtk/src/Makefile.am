# Copyright (C) 2001 Jan T. Kim <kim@inb.mu-luebeck.de>
#
# This file is free software; as a special exception the author gives
# unlimited permission to copy and/or distribute it, with or without
# modifications, as long as this notice is preserved.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY, to the extent permitted by law; without even the
# implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

# $Id$
# $Log$
# Revision 1.1  2005/03/08 17:12:02  jtk
# Initial revision
#
# Revision 1.6  2004/04/13 17:36:37  kim
# some changes to adapt to new automake / autoconf under Debian
#
# Revision 1.5  2002/12/09 19:37:55  kim
# added dot rendering of transsys network graphs (initial version)
#
# Revision 1.4  2002/05/28 08:52:50  kim
# added some discrete network stuff to sources
#
# Revision 1.3  2002/01/25 03:35:03  kim
# Added gnuplot link functionality to transexpr, transscatter, improved
#     PostScript stuff
#
# Revision 1.2  2001/04/05 17:24:50  kim
# Added initialization of transsys_name = NULL
#
# Revision 1.1  2001/04/04 11:12:00  kim
# Initial addition of files previously not CVS managed
#

bin_PROGRAMS = transps transdot transcheck transexpr transscatter ltrcheck ltransps ltransgl transdiscr
lib_LIBRARIES = libtrans.a
include_HEADERS = transsys.h trconfig.h trtypes.h
noinst_HEADERS = trbison.h
libtrans_a_SOURCES = trlex.l trbison.y mem.c save.c urandom.c expr.c express.c pscript.c dot.c cell.c lsys.c transsys.h trconfig.h trtypes.h
transcheck_SOURCES = transcheck.c
transdiscr_SOURCES = transdiscr.c
transps_SOURCES = transps.c
transexpr_SOURCES = transexpr.c
transscatter_SOURCES = transscatter.c
ltrcheck_SOURCES = ltrcheck.c
ltransps_SOURCES = ltransps.c
ltransgl_SOURCES = ltransgl.c
LDADD = libtrans.a

CLEANFILES = *~
DISTCLEANFILES = trlex.c trbison.h trbison.c
EXTRA_DIST = trlex.c trbison.h trbison.c

trlex.c : trbison.c

