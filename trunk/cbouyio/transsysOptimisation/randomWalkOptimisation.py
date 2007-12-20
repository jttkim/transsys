#!/usr/bin/env python

# Subversion keywords.
# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date$:  Date of last commit

"""General description here.

@author: Costas Bouyioukos
@email: konsb@cmp.uea.ac.uk
@organisation: University of East Anglia
@since: <<DATE>>
@license: GNU General Public Lisence 2 or newer.
@contact: U{Costas Bouyioukos<mailto:konsb@cmp.uea.ac.uk>}
@version: $Id$"""


# Version Information.
__version__ = "$Id$"


import getopt
import sys

import transsys


# Name.
__name__ = sys.argv[0]

def usage() :
  print """Usage:"""
  print __name__, """[options] ARGUMENTS
Options:
  -h, --help : Print this screen and exit.
  -v, --version : Print the program revision information.
ARGUMENTS"""


## Command line parsing.
optlist, args = getopt.getopt(sys.argv[1:], 'hv', ['help', 'version'])

for opt, par in optlist :
  if opt in ('-h'):
    usage()
    sys.exit()
  if opt in ('--help'):
    print __name__
    print __doc__
    usage()
    sys.exit()
  if opt in ('-v', '--version'):
    print __version__
    print __doc__
    sys.exit()


# I/O implementation.
if len(args) > 0 :
  infile = open(args[0], 'r')
else :
  infile = sys.stdin
if len(args) > 1 :
  outfile = open(args[1], 'w')
else :
  outfile = sys.stdout


# Main program



# Close outfile and exit.
if outfile is not sys.stdout:
  outfile.close()

sys.exit

