#!/usr/bin/env python

from distutils.core import setup

setup(name='translattice',
  version='0.1.0',
  description='A python module for lattice simulation of transsys programs',
  long_description='pending...',
  py_modules=['translattice'],
  scripts=['latticeSimulator', 'zeroTranssysDiffusibility'],
  author='Costas Bouyioukos',
  author_email='konsb@cmp.uea.ac.uk',
  )

