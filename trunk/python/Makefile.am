# $Id$
# $Revision$
# $Author$
# $Date$

CLEANFILES	= transsys/*.pyc transsys/*.pyo
DISTCLEANFILES	= 
EXTRA_DIST	= setup.py transsys src

SUBDIRS		= scripts

# TESTS		= test_transsys


install-exec-local : pythonsetup-install

pythonsetup-install : setup.py
	cd $(<D) ; python setup.py install --$(python_scheme)=$(python_prefix)

pythonsetup-testinstall : setup.py
	rm -rf test
	mkdir test
	cd $(<D) ; python setup.py install --home=test

check-local : test_transsys pythonsetup-testinstall
	./test_transsys

clean-local : setup.py
	cd $(<D) ; python setup.py clean
	rm -rf build doc test

all-local : pythonsetup-build epydoc-html

pythonsetup-build : setup.py transsys
	cd $(<D) ; python setup.py build

# this is a weird kludge -- the problem is that the transsys
# subdirectory is not accepted...

epydoc-html : pythonsetup-build
	if test \! -d doc ; then mkdir doc ; fi
	if test \! -d doc/html ; then mkdir doc/html ; fi
	epydoc -v --html -o doc/html `find build -name transsys`

# checking after subversion update...
