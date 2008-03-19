#!/usr/bin/env python

# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date$:  Date of last commit

from distutils.core import setup

setup(name='translattice',
  version='0.1.1',
  description='A python module for 2D structure simulation of transsys programs',
  #FIXME The discription should be changed after the introduction of the
  # optimisation modules.
  long_description='The setup apart from the python module translattice.py , installs several scripts, -in the <prefix>/bin directory. These scripts are essential for conducting simulation studies and they include: the latticeSimulator script (run the simulator according to the user specified control parameters) and the alterTranssysDiffusibility (produce the "control" transsys programs). All the experiments can be conducted automatically with the latticeExperinents.sh shell shell script (it is not part of the install).',
  py_modules=['translattice', 'randomLocalSearch', 'bimodalities'],
  scripts=['latticeSimulator', 'optimisation/runOptimiser', '../scripts/alterTranssysDiffusibility'],
  author='Costas Bouyioukos',
  author_email='konsb@cmp.uea.ac.uk',
  classifiers=[
   'Development Status :: Pending',
   'Environment :: Console',
   'Intended Audience :: Sci Software/Rearchers',
   'License :: GNU public Lisence v 2.0 or later',
   'Operating System :: POSIX :: Linux',
   'Programming Language :: Python :: Python2.5',
   'Topic :: Biological Simulators',
  ],
)

