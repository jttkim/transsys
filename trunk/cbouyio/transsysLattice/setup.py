#!/usr/bin/env python

# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date$:  Date of last commit

from distutils.core import setup

setup(name='translattice',
  version='0.1.1',
  description='A python module for 2D structure simulation of transsys programs',
  long_description='The setup apart from the python module translattice.py , installs several scripts, -in the <prefix>/bin directory. These scripts are essential for conducting simulation studies and they include: the latticeSimulator script (runs the simulator according to the user specified control parameters), the alterTranssysDiffusibility (produces the two "control" transsys programs) and the conductWholeExperiment.sh (wraps all the scripts together and produce an .R source file ready to load and analyse the data to R) the conductWholeExperiment.sh is obvious that is not a python but a shell script',
  py_modules=['translattice'],
  scripts=['latticeSimulator', '../scripts/alterTranssysDiffusibility', 'conductWholeExperiment.sh'],
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

