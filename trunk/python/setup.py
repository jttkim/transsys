#!/usr/bin/env python

import distutils.core
import os.path


# provisional determination of transsys installation prefix
# just assumes it's the user's login directory...
# FIXME: search for the include file and library and set parameters
# only according to what's found.

transsys_home = os.path.expanduser('~')
incdirs = ['../src']
# incdirs.append(os.path.join(transsys_home, 'include'))
libdirs = ['../src']
libdirs.append(os.path.join(transsys_home, 'lib'))
libs = []
libs.append('trans')


clib = distutils.core.Extension('transsys.clib',
                                sources = ['src/clib.c'],
                                include_dirs = incdirs,
                                library_dirs = libdirs,
                                libraries = libs,
                                extra_compile_args = ['-Wall', '-pedantic', '-Wno-long-long', '-fPIC'],
				extra_link_args = ['-fPIC']
                                )

distutils.core.setup(name = 'transsys',
                     version = '0.1',
                     description = 'transsys python module and utilities',
                     author = "Jan T. Kim",
                     author_email = "jtk@cmp.uea.ac.uk",
                     packages = ['transsys'],
		     scripts = ['trsys2sbml', 'trsystool'],
                     ext_modules = [clib])
