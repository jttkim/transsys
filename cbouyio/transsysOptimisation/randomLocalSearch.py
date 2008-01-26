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

import transsys
import translattice
import bimodalities

# Random objects already implemeted in the transsysLattice module.


class TranssysProgramOprtimised(transsys.TranssysProgram):
  """Sybclass of transsys program to curry out optimisation needs
  """
#  def __init__(self) :
  pass


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



def perturb_transsys(transsysProgram, perturbObj):
  """Return a perturbed transsys program.

  Perturbation refers to multiply by arandom generated value the constant
  expresssion values of a transsys program.
  @param transsysProgram: A transsys program instance.
  @type transsysProgram: C{class 'transsys.TranssysProgram'}
  @param perturbObj: A random object.
  @type perturbObj: C{class 'random.Random'}
  """
  # First copy the transsys program.
  perturbedTP = copy.deepcopy(transsysProgram)
  for geneValues in get_opt_list(perturbedTP) :
    for exprValue in geneValues :
      exprValue.value = exprValue.value + exprValue.value*perturbObj.random_value()
  return perturbedTP



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
  if not isinstance(transsysLattice, translattice.TranssysInstanceLattice) :
    raise StandardError, 'Error  in expression_table(), called object is not a transsys lattice.'
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
  rndObj = translattice.UniformRNG(rndSeed, unifBorders[0], unifBorders[1])
  transLat.initialise_lattice(rndObj)
  transsys.clib.srandom(rndObj.rndSeed)
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


def calculate_objective(tp, latticeSize, timesteps, unifBorders, rndSeed):
  """A wrapper function which conducts the optimisation objective calculation.

  @param
  @type
  """
  # The lattice experiment.
  exprTableLattice = run_lattice(tp, latticeSize, timesteps, unifBorders, rndSeed)
  bmScoreLattice = calculate_bimodality(exprTableLattice)
  # The control.
  tpZero = copy.deepcopy(tp)
  zero_transsys_diffusibility(tpZero)
  exprTableControl = run_lattice(tpZero, latticeSize, timesteps, unifBorders, rndSeed)
  bmScoreControl = calculate_bimodality(exprTableControl)
  objective = bmScoreLattice - bmScoreControl
  return objective



def transsys_lattice_optimisation(transsysProgram, latticeSize, timesteps, unifBorders, rndSeed, perturbBorders, optimisationCycles):
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
  @param perturbBorders:
  @type perturbBorders:
  @param optimisationCycles: The number of optimisation cycles.
  @type optimisationCycles: C{'int'}
  @return: An "optimised" transsys program
  @rtype: C{class 'transsys.TranssysProgram'}
  """
  # Set the perturbation object.
  perturbObj = translattice.UniformRNG(rndSeed, perturbBorders[0], perturbBorders[1])
  # Calcullate the initial transsys program optimisation objective.
  optimisationObjective = calculate_objective(transsysProgram, latticeSize, timesteps, unifBorders, rndSeed)
  bestTP = transsysProgram
  outputList = [(0, optimisationObjective, optimisationObjective, bestTP, bestTP, [expr.value for expr in bestTP.getGeneValueNodes()], [expr.value for expr in bestTP.getGeneValueNodes()])]
  for cycle in xrange(optimisationCycles) :
    alternativeTP = perturb_transsys(bestTP, perturbObj)
    oo = calculate_objective(alternativeTP, latticeSize, timesteps, unifBorders, rndSeed)
    if oo > optimisationObjective :
      optimisationObjective = oo
      bestTP = alternativeTP
    outputList.append((cycle + 1, optimisationObjective, oo, bestTP, alternativeTP, [expr.value for expr in bestTP.getGeneValueNodes()], [expr.value for expr in alternativeTP.getGeneValueNodes()]))
  print bestTP
  return outputList

