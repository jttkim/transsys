#!/usr/bin/env python

# Subversion keywords.
# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date$:  Date of last commit

"""Python module containing all the functionalities to conduct the optimisation
experiments.

@author: Costas Bouyioukos
@email: konsb@cmp.uea.ac.uk
@organisation: University of East Anglia
@since: 20/12/2007
@license: GNU General Public Lisence 3 or newer.
@contact: U{Costas Bouyioukos<mailto:konsb@cmp.uea.ac.uk>}
@version: $Id$"""


# Version Information.
__version__ = "$Id$"


import copy

import transsyslattice
import bimodalities


# Random object already implemeted in the transsysLattice module.


# The perturbed transsys program object.

class PerturbTranssysProgram(transsys.TranssysProgram):
  """A class to implement the transsys program random walk

  Contains all the necessay methods to perform the transsys program
  randomisation.
  """

  def __init__(self, tp, cycle, rndSeed=1);
    """The constructor

    @param tp: the transsys program.
    @type tp: C{class 'transsys.TranssysProgram'}
    @param cycle: the optimisation cycle.
    @type cycle: C{int}
    @param rndSeed: the random seed.
    @type rndSeed: C{int}
    """

  def randomStep(self, tp, rndSeed):
    """A function wich perturbs a transsys program by one random step.

    @param tp: the transsys program.
    @type tp: C{class 'transsys.TranssysProgram'}
    @param rndSeed: the random seed.
    @type rndSeed: C{int}
    """


def runTranssysLatticeExperiment(transsysProgram, latticeSize, timesteps,
    experiment, samplingInterval=timesteps):
  """A function to run the transsyslattice simulator.

  The function deals with the running of the transsys lattice experiment and
  returns a transsyslattice.TranssysLatticeTimeseries object
  @param transsysProgram: A valid transsys program.
  @type transsysProgram: C{class 'transsys.TranssysProgram'}
  @param latticeSize: A list denoting the dimensions of the lattice.
  @type latticeSize: c{list} of C{int}
  @param timesteps: The number of timesteps for the lattice experiment.
  @type timesteps: C{int}
  @param experiment: The type of experiment that we need to conduct, lattice
  (with difussion) or the control.
  @type experiment: C{logical}
  @param samplingInterval: The interval between sampling. It is the same with
  the number of timesteps because we need to keep only the final state of the
  simulator.
  @type samplingInterval: C{int}
  @return: C{class 'transsyslattice.TranssysLatticeTimeseries'}
  """


def calculateBimodality(transsysLatticeTimeseries):
  """A function which wraps all the bimodality module functionalities and
  calculates the total bimodality of a collection of transsys instances.

  Returns a scalar.
  @param transsysLatticeTimeseries: A transsys lattice timeseries contains the
  transsys lattice mulator's output.
  @type transsysLatticeTimeseries: C{class 'transsyslattice.TranssysLatticeTimeseries'}
  """
