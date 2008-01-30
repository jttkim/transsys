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
import math

# Random objects already implemeted in the transsysLattice module.


class TranssysProgramToy(transsys.TranssysProgram) :
  """Sybclass of transsys program to curry out optimisation needs

  """

  def __init__(self, name, factor_list = None, gene_list = None, resolve = True, comments = []) : #lowPoint = None, highPoint = None, circle = None) :
    self.name = name
    self.gene_list = gene_list
    self.factor_list = factor_list
    self.comments = comments


  def extract_coordinates(self) :
    """A function to retrieve the optimisation values from the toy model
    transsys program.

    @return: coordinates
    """
    decayList = [f.decay_expression.value for f in self.factor_list]
    coordinates = {}
    geneList = self.gene_list
    promoterList = [g.promoter[0] for g in geneList]
    expressionList = [p.expression for p in promoterList]
    nodesList = [e.getValueNodes() for e in expressionList]
    expressionValueList = [[]]*len(nodesList)
    for i in xrange(len(nodesList)) :
      expressionValueList[i] = [e.value for e in nodesList[i]]
    coordinates['lowPoint'] = [expressionValueList[0][0], expressionValueList[1][0]]
    coordinates['highPoint'] = [expressionValueList[0][3], expressionValueList[1][3]]
    coordinates['circle'] = [expressionValueList[0][5], expressionValueList[0][7], expressionValueList[0][9]]
    return coordinates


  def check_coordinates(self, coordinates) :
    """ A function to check weather a set of coordinates are fulfilling the
    designed criteria.

    """
    decayList = [f.decay_expression.value for f in self.factor_list]
    coordinates['lowPoint'] = [coordinates['lowPoint'][0] / float(decayList[0]), coordinates['lowPoint'][1] / float(decayList[1])]
    coordinates['highPoint'] = [coordinates['highPoint'][0] / float(decayList[0]), coordinates['highPoint'][1] / float(decayList[1])]
    if coordinates['circle'][0] - coordinates['circle'][2] < 0 or coordinates['circle'][1] - coordinates['circle'][2]  < 0 :
      return False
    if coordinates['lowPoint'][0] > coordinates['circle'][0] - coordinates['circle'][2] :
      return False
    if coordinates['lowPoint'][1] > coordinates['circle'][1] - coordinates['circle'][2] :
      return False
#    if coordinates['highPoint'][0] < coordinates['circle'][0] - coordinates['circle'][2] :
#      return False
#    if coordinates['highPocheck_coordinatesint'][1] < coordinates['circle'][1] - coordinates['circle'][2] :
#      return False
    if coordinates['highPoint'][0] > coordinates['circle'][0] + coordinates['circle'][2] and coordinates['highPoint'][1] > coordinates['circle'][1] + coordinates['circle'][2] :
      return False
    if coordinates['highPoint'][0] < coordinates['circle'][0] + coordinates['circle'][2] and coordinates['highPoint'][1] < coordinates['circle'][1] + coordinates['circle'][2] :
      return False
    else :
      return True


  def perturb_coordinates(self, coord, rndSeed, perturbRange) :
    """A function to perturb the coordinate set.

    """

    perturbObj = translattice.UniformRNG(rndSeed, perturbRange[0], perturbRange[1])
    for k, cValues in coord.iteritems() :
      value = []
      for v in cValues :
        value.append(v*math.exp(perturbObj.random_value()))
      coord[k] = value
#    if not self.check_coordinates(coord) :
#      self.perturb_coordinates(coord, rndSeed, perturbRange)



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
      # Multiply the current value with the exponential of the random value
      exprValue.value = exprValue.value*math.exp(perturbObj.random_value())
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
  @rtype: C{list} of {list}s of {float}s
  """
  if not isinstance(transsysLattice, translattice.TranssysInstanceLattice) :
    raise StandardError, 'Error  in expression_table(), called object is not a transsys lattice.'
  exprTable = [[]] * transsysLattice.transsysProgram.num_factors()
  for ti in transsysLattice.transsys_instance_list() :
    for i in xrange(len(ti.factor_concentration)) :
      exprTable[i].append(ti.factor_concentration[i])
  return exprTable



def run_lattice(transsysProgram, latticeSize, timesteps, rndObj):
  """Function to conduct the lattice experiment.

  It either runs the lattce or the control (diffusibility = 0), returns the
  factor expression table.
  """
  # Set the transsys program to a local variable tp) to protect it from
  # transforming the diffusion.
  transLat = translattice.TranssysInstanceLattice(transsysProgram, latticeSize, timesteps)
#  rndObj = translattice.UniformRNG(rndSeed, unifBorders[0], unifBorders[1])
  transLat.initialise_lattice(rndObj)
  # Set always the same random seed for the clib.
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


def calculate_objective(tp, latticeSize, timesteps, initialObj):
  """A wrapper function which conducts the optimisation objective calculation.

  @param
  @type
  """
  # Copy the random initalisation object to preserve the same initial conditions
  # for the lattice and the control.
  initialObj2 = copy.deepcopy(initialObj)
  # The lattice experiment.
  exprTableLattice = run_lattice(tp, latticeSize, timesteps, initialObj)
  bmScoreLattice = calculate_bimodality(exprTableLattice)
  # The control.
  tpZero = copy.deepcopy(tp)
  zero_transsys_diffusibility(tpZero)
  exprTableControl = run_lattice(tpZero, latticeSize, timesteps, initialObj2)
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
  initialObj = translattice.UniformRNG(rndSeed, unifBorders[0], unifBorders[1])
  # Calcullate the initial transsys program optimisation objective.
  optimisationObjective = calculate_objective(transsysProgram, latticeSize, timesteps, initialObj)
  bestTP = transsysProgram
  outputList = [(0, optimisationObjective, optimisationObjective, bestTP, bestTP, [expr.value for expr in bestTP.getGeneValueNodes()], [expr.value for expr in bestTP.getGeneValueNodes()])]
  for cycle in xrange(optimisationCycles) :
    alternativeTP = perturb_transsys(bestTP, perturbObj)
    oo = calculate_objective(alternativeTP, latticeSize, timesteps, initialObj)
    if oo > optimisationObjective :
      optimisationObjective = oo
      bestTP = alternativeTP
    outputList.append((cycle + 1, optimisationObjective, oo, bestTP, alternativeTP, [expr.value for expr in bestTP.getGeneValueNodes()], [expr.value for expr in alternativeTP.getGeneValueNodes()]))
    print optimisationObjective, oo
  return outputList

