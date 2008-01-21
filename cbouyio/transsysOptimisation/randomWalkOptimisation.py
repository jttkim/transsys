#!/usr/bin/env python

# Subversion keywords.
# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date: 2007-12-20 20:55:15 +0000 (Thu, 20 Dec 2007) $:  Date of last commit

"""Python module containing all the functionalities to conduct the optimisation
experiments.

@author: Costas Bouyioukos
@email: konsb@cmp.uea.ac.uk
@organisation: University of East Anglia
@since: 20/12/2007
@license: GNU General Public Lisence 3 or newer.
@contact: U{Costas Bouyioukos<mailto:konsb@cmp.uea.ac.uk>}
@version: $Id: randomWalkOptimisation.py 275 2008-01-08 15:44:39Z cbouyio $"""

# Version Information.
__version__ = "$Id: randomWalkOptimisation.py 275 2008-01-08 15:44:39Z cbouyio $"

import copy
#import cPickle

import transsys
import translattice
import bimodalities

# Random objects already implemeted in the transsysLattice module.



def get_opt_list(transsysProgram):
  """A function to get the list of values subject to optimisation

  @param transsysProgram: A transsys program instance.
  @type transsysProgram: C{class 'transsys.TranssysProgram'}
  @ return: A list of size the gene number of lists with each list containing
  the constant values of each gene.
  @rtype: C{'list'} of C{'list'}s with the constant expression values.
  """
  optList = []
  for gene in transsysProgram.gene_list :
    optList.append(gene)
  optList = [g.getValueNodes() for g in optList]
  return optList



def random_step(transsysProgram, rndObj):
  """Perturb "on the fly" a transsys program by a random step.

  @param transsysProgram: A transsys program instance.
  @type transsysProgram: C{class 'transsys.TranssysProgram'}
  @param rndObj: A random object.
  @type rndObj: C{class 'random.Random'}
  """
  # Then strat the perturbations.
  for geneValues in get_opt_list(transsysProgram) :
    for exprValue in geneValues :
      exprValue.value = exprValue.value*rndObj.random_value(0, 2.001)



def zero_transsys_diffusibility(transsysProgram):
  """Set to zero the diffusibility expression of all factors.

  @param transsysProgram: An instance of a transsys program.
  """
  for diffusion in transsysProgram.getDiffusibilityValueNodes() :
    diffusion.value = 0



def expression_table(transsysLattice):
  """Return a factor concentration table (factors x instances) out of a transsys
  lattice.

  @param transsysLattice: A transsys instance lattice object.
  @type transsysLattice: C{class 'translattice.TranssysInstanceLattice'}
  @return: A factor x instances factor concentration table.
  @rtype; C{list} of {list}s of {float}s
  """
  exprTable = []
  for f in xrange(transsysLattice.transsysProgram.num_factors()):
    exprTable.append([])
  for ti in transsysLattice.transsys_instance_list() :
    for i in xrange(len(ti.factor_concentration)) :
      exprTable[i].append(ti.factor_concentration[i])
  return exprTable



def run_lattice(transsysProgram, latticeSize, timesteps, unifBorders, rndSeed):
  """Function to conduct the lattice experiment.

  It either runs the lattce or the control (diffusibility = 0), returns the
  factor expression table.
  """
  # Set the transsys program to a local variable tp) to protect it from
  # transforming the diffusion.
  transLat = translattice.TranssysInstanceLattice(transsysProgram, latticeSize, timesteps)
  transLat.initialise_lattice_concentrations(unifBorders[0], unifBorders[1], rndSeed)
  transsys.clib.srandom(rndSeed)
  transTS = translattice.TranssysLatticeTimeseries(transLat, timesteps, timesteps)
  lattice = transTS.latticeTimeseries.pop()
  expressionTable = expression_table(lattice)
  return expressionTable



def calculate_bimodality(expressionTable):
  """A function which wraps all the bimodalities module functions and
  calculate the total bimodality from an expression table.

  Returns a scalar.
  @param expressionTable: A table (as it is produced by the expression_table
  function) of size (factors x transsys instances) which contains the factor
  concentrations.
  @type expressionTable: C{'list'} od C{'list'}s of C{'float'}s
  @return: The total bimodality score.
  @rtype: C{'float'}
  """
  return bimodalities.total_bimodality(expressionTable)



def transsys_lattice_optimisation(transsysProgram, latticeSize, timesteps, unifBorders, rndSeed, optimisationCycles):
  """Conduct the optimisation experiment.

  The function implements all the necessary steps to run a complete optimisation
  experiment and returns an "optimised" transsys program.
  @param transsysProgram: A valid transsys program.
  @type transsysProgram: C{class 'transsys.TranssysProgram'}
  @param latticeSize: A 2n-tuple containing the dimensions of the lattice.
  @type latticeSize: c{tuple} of 2 c{int}s
  @param timesteps: The number of timesteps for the lattice experiment.
  @type timesteps: C{int}
  @param unifBorders: The borders out od which the uniform distribution will
  take random numbers.
  @type unifBorders: a 2n-C{'tuple'} of C{'float'}
  @param rndSeed: The random number generator seed.
  @type rndSeed: C{'int'}
  @param optimisationCycles: The number of optimisation cycles.
  @type optimisationCycles: C{'int'}
  @return: An "optimised" transsys program
  @rtype: C{class 'transsys.TranssysProgram'}
  """
  outputList = []
  # First set the random object.
  rndObj = translattice.UniformRNG(rndSeed)
  # Then calcullate the initial transsys program optimisation objective.
  # Run the lattice experiment.
  exprTableLattice = run_lattice(transsysProgram, latticeSize, timesteps, unifBorders, rndSeed)
  bmScoreLattice = calculate_bimodality(exprTableLattice)
  # Then the control.
  tpControl = copy.deepcopy(transsysProgram)
  zero_transsys_diffusibility(tpControl)
  exprTableControl = run_lattice(tpControl, latticeSize, timesteps, unifBorders, rndSeed)
  bmScoreControl = calculate_bimodality(exprTableControl)
  optimisationObjective = bmScoreLattice - bmScoreControl
  optimisationDict = {}
  walk = True
  for cycle in xrange(optimisationCycles) :
    # This condition implements the random local search instead of a pure random
    # walk.
    if walk :
      random_step(transsysProgram, rndObj)
      tpCopy = copy.deepcopy(transsysProgram)
    else :
      random_step(tpCopy, rndObj)
    # The lattice experiment.
    exprTableLattice = run_lattice(tpCopy, latticeSize, timesteps, unifBorders, rndSeed)
    bmScoreLattice = calculate_bimodality(exprTableLattice)
    # The control.
    tpControl = copy.deepcopy(tpCopy)
    zero_transsys_diffusibility(tpControl)
    exprTableControl = run_lattice(tpControl, latticeSize, timesteps, unifBorders, rndSeed)
    bmScoreControl = calculate_bimodality(exprTableControl)
    oo = bmScoreLattice - bmScoreControl
    # Actual optimisation procedure.
    if oo > optimisationObjective :
      optimisationObjective = oo
      walk = False
    else :
      walk = True
#    print cycle + 1, optimisationObjective, oo
    outputList.append((optimisationObjective, oo, tpCopy))
  return outputList

