# $Id$
# $Revision$
# $Author$
# $Date$

CLEANFILES	= transsys/*.pyc transsys/*.pyo
DISTCLEANFILES	= 
EXTRA_DIST	= setup.py test_transsys trsys2sbml trsystool

SUBDIRS		= scripts transsys src

# TESTS		= test_transsys
# kludge to avoid stalling on systems that don't have epydoc
if HAVE_EPYDOC
EPYDOC_HTML	= epydoc-html
endif

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

all-local : pythonsetup-build $(EPYDOC_HTML)

pythonsetup-build : setup.py transsys
	cd $(<D) ; python setup.py build

# this is a weird kludge -- the problem is that the transsys
# subdirectory is not accepted...

epydoc-html : pythonsetup-build
	if test \! -d doc ; then mkdir doc ; fi
	if test \! -d doc/html ; then mkdir doc/html ; fi
	epydoc -v --html -o doc/html `find build -name transsys`

