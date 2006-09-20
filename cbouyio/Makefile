# $Id$
# $Date$
# $Author$

PYTHONMODULES	= translattice.py

html : $(PYTHONMODULES)
	rm -rf html
	epydoc --html $(PYTHONMODULES)

install : $(PYTHONMODULES)
	python setup.py install --home=$(HOME)

clean :
	rm -rf html *.pyc build

.PHONY : install clean

