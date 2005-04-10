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
# Revision 1.1  2005/03/08 17:12:01  jtk
# Initial revision
#
# Revision 1.7  2003/01/28 16:59:19  kim
# separated makefile for transsys_report from Makefile.am, simplified things
#
# Revision 1.6  2002/01/25 03:36:20  kim
# Wrote some paragraphs in transsys_report, hopefully removed necessity
#     of generating GL graphics after unpacking from archive
#
# Revision 1.4  2001/04/05 18:46:20  kim
# Added example files to EXTRA_DIST
#
# Revision 1.3  2001/04/05 17:24:27  kim
# Added examples (partially finished only)
#
# Revision 1.2  2001/04/05 15:41:39  kim
# moved transsys.tex to transsys_report.tex
#
# Revision 1.1  2001/04/04 11:11:55  kim
# Initial addition of files previously not CVS managed
#

TRANSSYSEXAMPLES = tr_ex00.tra tr_ex01.tra tr_ex02.tra tr_ex03.tra tr_ex04.tra tr_ex05.tra tr_ex06.tra \
trl_ex00.trl trl_ex01.trl trl_ex02.trl trl_ex03.trl trl_ex04.trl trl_ex05.trl trl_ex06.trl trl_ex07.trl trl_ex08.trl trl_ex09.trl trl_ex10.trl \
tr_ex02full.tra  tr_ex04full.tra  tr_ex05full.tra tr_ex06full.tra \
trl_ex04full.trl trl_ex05full.trl trl_ex06full.trl trl_ex09full.trl trl_ex10full.trl
EXTRA_DIST = transsys_report.tex transsys_report.ps transsys_report.pdf $(TRANSSYSEXAMPLES) report.makefile

CLEANFILES = transsys_report.aux transsys_report.log $(DATA) *~
MAINTAINERCLEANFILES = *.eps *.aux *.bbl *.blg *.ppm

transsys_report.ps transsys_report.pdf : report.makefile
	$(MAKE) -f report.makefile





