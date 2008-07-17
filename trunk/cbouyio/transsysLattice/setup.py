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
  long_description='The setup apart from the python module translattice.py, installs anoe more script, -in the <prefix>/bin directory. These scripts are facilitating a range of simulation and control experiments, they include: The latticeSimulator python script which runs the simulator according to the user specified control parameters. Flags on this script run the two control experiments these are the well stirred reactor and the unstructured collection of cells. All the experiments can be conducted together with the same control parameters by the latticeExperinents.sh shell script which is not part of this installation. The installer also includes the modules required from the optimisation experiments. And installs at the <prefix>/bin the runOptimiser which conducts optimisation experiments.',
  py_modules=['translattice', 'randomLocalSearch', 'bimodalities'],
  scripts=['latticeSimulator', 'optimisation/runOptimiser'],
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

