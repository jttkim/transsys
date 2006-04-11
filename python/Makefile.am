CLEANFILES = transsys/__init__pyc
DISTCLEANFILES = 
EXTRA_DIST = setup.py transsys src

install-exec-local :
	python setup.py install --$(python_scheme)=$(python_prefix)

clean-local :
	python setup.py clean
	rm -rf build

all-local :
	python setup.py build

