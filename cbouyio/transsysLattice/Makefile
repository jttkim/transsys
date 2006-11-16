# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date$:  Date of last commit

PYTHONMODULES	= translattice.py

html : $(PYTHONMODULES)
	rm -rf html/
	epydoc --html $(PYTHONMODULES)

install : $(PYTHONMODULES)
	python setup.py install --home=$(HOME)

clean :
	rm -rf html/ build/ *.pyc

.PHONY : install clean

