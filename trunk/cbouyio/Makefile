# $Id$
# $Date$
# $Author$

html : translattice.py
	rm -rf html
	epydoc --html translattice.py

install : translattice.py
	python setup.py install --home=$(HOME)

.PHONY : install

