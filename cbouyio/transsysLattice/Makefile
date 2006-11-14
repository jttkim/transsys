# $Id: Makefile 223 2006-11-14 19:05:35Z cbouyio $
# $Date: 2006-11-14 19:05:35 +0000 (Tue, 14 Nov 2006) $
# $Author: cbouyio $

PYTHONMODULES	= translattice.py

html : $(PYTHONMODULES)
	rm -rf html/
	epydoc --html $(PYTHONMODULES)

install : $(PYTHONMODULES)
	python setup.py install --home=$(HOME)

clean :
	rm -rf html/ build/ *.pyc

.PHONY : install clean

