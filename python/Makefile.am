CLEANFILES = transsys/__init__.pyc
DISTCLEANFILES = 
EXTRA_DIST = setup.py transsys src

# TESTS		= test_transsys


install-exec-local : pythonsetup-install

pythonsetup-install :
	python setup.py install --$(python_scheme)=$(python_prefix)

pythonsetup-testinstall :
	rm -rf test
	mkdir test
	python setup.py install --home=test

check-local : test_transsys pythonsetup-testinstall
	./test_transsys

clean-local :
	python setup.py clean
	rm -rf build doc test

all-local : pythonsetup-build epydoc-html

pythonsetup-build : transsys
	python setup.py build

# this is a weird kludge -- the problem is that the transsys
# subdirectory is not accepted...

epydoc-html : pythonsetup-build
	if test \! -d doc ; then mkdir doc ; fi
	if test \! -d doc/html ; then mkdir doc/html ; fi
	epydoc --html -o doc/html `find build -name transsys`

