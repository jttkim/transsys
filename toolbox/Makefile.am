## Process this file with automake to produce Makefile.in

# Copyright (C) 2002 Jan T. Kim <kim@inb.mu-luebeck.de>
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
# Revision 1.5  2005/10/12 19:33:32  jtk
# finished adding unit test (transparse)
#
# Revision 1.4  2005/06/22 09:56:28  jtk
# first introduction of checking
#
# Revision 1.3  2005/04/06 21:02:08  jtk
# added transdisrupt.py to installation target
#
# Revision 1.2  2005/04/04 21:31:44  jtk
# added random_transsys_lsys to Makefile.am
#
# Revision 1.1.1.1  2005/03/08 17:12:02  jtk
# new cvs after loss at INB
#
# Revision 1.9  2004/04/13 17:36:37  kim
# some changes to adapt to new automake / autoconf under Debian
#
# Revision 1.8  2003/03/20 16:07:41  kim
# fixed error of reversed edge direction in create_disruption_network,
#     other small changes
#
# Revision 1.7  2003/02/27 10:36:48  kim
# added ArrayIntensityFunction class for simulating noisy arrays etc.
#
# Revision 1.6  2003/02/18 21:47:59  kim
# added choice of perturbers for simulating different environments,
#     ExpGauss and ExponentialReplacement perturbers are currently
#     available. *.rta file format was modified in the process.
#
# Revision 1.5  2003/02/17 21:44:15  kim
# included transarr R package in distribution
#
# Revision 1.4  2003/02/17 21:40:30  kim
# Transformed microarr.r into a R library called transarr. Added array_background
#     to regstruct system
#
# Revision 1.3  2003/02/14 16:31:27  kim
# regstruct project is nearing completion. regstruct_transsys for
#     generation of transsys programs and regstruct_transarray for
#     array simulation are (more or less) properly implemented.
#     Analysis with R has been assembled in microarr.r
#
# Revision 1.2  2003/02/04 23:41:08  kim
# added regstruct scripts and module
#
# Revision 1.1  2003/01/28 21:12:05  kim
# initial toolbox assembly
#

PYTHONMODULES	= arrtool.py regstruct.py transdisrupt.py transdecode.py
PYTHONSCRIPTS	= random_transsys random_transsys_lsys transarray regstruct_transsys regstruct_transarray mutant_transsys get_overlap
RPACKAGE	= transarr

pylibdir	= $(DESTDIR)/$(libdir)/python
rlibdir		= $(DESTDIR)/$(libdir)/R

bin_SCRIPTS	= $(PYTHONSCRIPTS)
pylib_DATA	= $(PYTHONMODULES)

EXTRA_DIST	= $(PYTHONSCRIPTS) $(PYTHONMODULES) $(RPACKAGE) htmldoc

CLEANFILES	= *~ *.pyc

htmldoc : $(PYTHONMODULES)
	epydoc --html -o htmldoc $(PYTHONMODULES)

install-data-local : transarr/R/transarr.r
	if test $(r_program) = yes ; then R_LIBS=$(rlibdir) R CMD INSTALL $(srcdir)/transsys ; fi
	if test $(r_program) = yes ; then R_LIBS=$(rlibdir) R CMD INSTALL $(srcdir)/transarr ; fi
	@if echo "$$PYTHONPATH" | grep "$(pylibdir)" ; then true ; else echo ; echo Please include $(pylibdir) in your PYTHONPATH environment variable ; echo ; fi

all-local : htmldoc

clean-local :
	rm -rf htmldoc

uninstall-hook :
	rm -rf $(rlibdir)/transarr

