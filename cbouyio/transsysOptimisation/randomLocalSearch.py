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
import sys

import transsys
import translattice
import bimodalities
import math

# Random objects already implemeted in the transsysLattice module.


class Point(translattice.TwoValueParameter) :
  """Class to represent a two dimensional point.
  """
  def __init__(self, x , y) :
    """The constructor of the class.

    @param x: The x-coordinate of the expression point.
    @type x: C{float}
    @param y: The y-coordinate of the expression point.
    @type y: C{float}
    """
    self.x = float(x)
    self.y = float(y)

  def get_x(self) :
    """Return the x- coordinate.
    """
    return self.x

  def get_y(self) :
    """Return the y- coordinate.
    """
    return self.y

  def set_values(self, x, y) :
    """Sets the parameter values.

    @param x: The x-coordinate of the point.
    @type x: C{float}
    @param y: The y-coordinate of the point.
    @type y: C{float}
    """
    self.x = copy.deepcopy(x)
    self.y = copy.deepcopy(y)


class Circle(translattice.TwoValueParameter) :
  """Class to represent a circle. Coordinates of the centre and radius
  """
  def __init__(self, centre, r) :
    """The constructor of the class.

    @param centre: An object of class Point representing the centre
    of the circle.
    @type centre: {'Point'} object.
    @param r: The radius of the circle.
    @type r: C{float'}
    """
    if not isinstance(centre, Point) :
      raise StandardError, 'Invalid circle centre. Centre should be a Point object.'
    self.centre = centre
    self.r = float(r)

  def get_centre(self) :
    """Return the centre of the circle.

    As a Point object.
    """
    return self.centre


  def get_r(self) :
    """Returns the radious of the circle.
    """
    return self.r

  def set_values(self, centre, r) :
    """Sets the parameter values.

    """
    if not isinstance(centre, Point) :
      raise StandardError, 'Invalide parameter setting for the circle. The centre should be of type "Point"'
    self.centre = copy.deepcopy(centre)
    self.r = copy.deepcopy(r)


class EngineeringParameters(translattice.ControlParameters) :
  """A class to aggregate the control parameters of the engineered transsys
  program.

  The transsys program has the structure:
  transsys engineered
  {
    factor A {decay: <decayA>; diffusibility: <diffusibilityA>;
  }
    factor B {decay: <decayB>; diffusibility: <diffusibilityB>;

  gene geneA
  {
    promoter
    {
      constitutive: (<lowPoint.x>/<decayA> + ((<highPoint.x>/<decayA> - <lowPoint.x>/<decayA>) * ((((A - <cirlce.x) * (A - <circle.x>)) + ((B - <circle.y>) * (B - <circle.y>))) <= (<circle.r> * <circle.r>))));
    }
    product
    {
      default: A;
    }
  }

  gene geneB
  {
    promoter
    {
      constitutive: (<lowPoint.y>/<decayB> + ((<highPoint.y>/<decayA> - <lowPoint.y>/<decayA>) * ((((A - <cirlce.x) * (A - <circle.x>)) + ((B - <circle.y>) * (B - <circle.y>))) <= (<circle.r> * <circle.r>))));
    }
    product
    {
      default: B;
    }
  }
  """

  def __init__(self, lowPoint, highPoint, circle, decayA, diffusibilityA, decayB, diffusibilityB) :
    """
    """
    self.lowPoint = lowPoint
    self.highPoint = highPoint
    self.circle = circle
    self.decayA = decayA
    self.decayB= decayB
    self.diffusibilityA = diffusibilityA
    self.diffusibilityB = diffusibilityB


  def get_lowPoint(self) :
    """Return the low point as a Point object.

    """
    return self.lowPoint


  def get_highPoint(self) :
    """Return the high point as a Point object.

    """
    return self.highPoint

  def get_circle(self) :
    """Return the circle as a Circle object.

    """
    return self.circle

  def get_decayA(self) :
    """Return the decay rate of factor A.

    """
    return self.decayA

  def get_decayB(self) :
    """Return the decay rate of factor B.

    """
    return self.decayB

  def get_diffusibilityA(self) :
    """Return the diffusibility rate of factor A.

    """
    return self.diffusibilityA

  def get_diffusibilityB(self) :
    """Return the diffusibility rate of factor B.

    """
    return self.diffusibilityB

  def set_perturbed_parameters(self, lowPoint, highPoint, circle) :
    """Sets the perturbed parameter values.

    """
    self.lowPoint = copy.deepcopy(lowPoint)
    self.highPoint = copy.deepcopy(highPoint)
    self.circle = copy.deepcopy(circle)



class TranssysProgramDummy(transsys.TranssysProgram) :
  """Subclass of transsys program to build "from scratch" a transsys program
  from the given control parameters set.

  """

  def __init__(self, name, engineeringParameters) :
    """The constructor of the class. Check the docstring of the class to see the
    structure of the transsys program.
    """
    self.name = name
    if isinstance(engineeringParameters, EngineeringParameters) :
      self.engineeringParameters = engineeringParameters
    else :
      raise StandardError, 'Error in enginnering parameters setting. The engineering parameters should be an EngineeringParameters object.'
    self.factor_list = [self.construct_factor('A', self.engineeringParameters.decayA, self.engineeringParameters.diffusibilityA), self.construct_factor('B', self.engineeringParameters.decayB, self.engineeringParameters.diffusibilityB)]
    self.gene_list = [self.construct_gene('geneA', 'A', self.engineeringParameters.lowPoint.x, self.engineeringParameters.highPoint.x, self.engineeringParameters.circle, self.engineeringParameters.decayA), self.construct_gene('geneB', 'B', self.engineeringParameters.lowPoint.y, self.engineeringParameters.highPoint.y, self.engineeringParameters.circle, self.engineeringParameters.decayB)]
    transsys.TranssysProgram.__init__(self, self.name, self.factor_list, self.gene_list)


  def construct_gene(self, name, product, lp, hp, circle, decay) :
    """Mehthod to construct the complex constitutive expression out of the
    system's points coordinates.

    Builds up all the complex expression nodes for each gene of the system, by
    constructing a ste af expression nodes insances.

    @return: A transsys gene, with the complex constitutive expression.
    @rtype: a C{'transsys.Gene'} object
    """
    exprH = transsys.ExpressionNodeValue(hp * float(decay))
    exprL = transsys.ExpressionNodeValue(lp * float(decay))
    circleX = transsys.ExpressionNodeValue(circle.centre.x)
    circleY = transsys.ExpressionNodeValue(circle.centre.y)
    circleR = transsys.ExpressionNodeValue(circle.r)
    factorA = transsys.ExpressionNodeIdentifier(self.factor_list[0])
    factorB = transsys.ExpressionNodeIdentifier(self.factor_list[1])
    expressionCircleX = transsys.ExpressionNodeMult(transsys.ExpressionNodeSubtract(factorA, circleX), transsys.ExpressionNodeSubtract(factorA, circleX))
    expressionCircleY = transsys.ExpressionNodeMult(transsys.ExpressionNodeSubtract(factorB, circleY), transsys.ExpressionNodeSubtract(factorB, circleY))
    expressionCircle = transsys.ExpressionNodeAdd(expressionCircleX, expressionCircleY)
    expressionRadious = transsys.ExpressionNodeMult(circleR, circleR)
    expressionBool = transsys.ExpressionNodeLowerEqual(expressionCircle, expressionRadious)
    expressionGeneExpr = transsys.ExpressionNodeMult(transsys.ExpressionNodeSubtract(exprH, exprL), expressionBool)
    constitutive = transsys.ExpressionNodeAdd(exprL, expressionGeneExpr)
    promoterConstitutive = transsys.PromoterElementConstitutive(constitutive)
    gene = transsys.Gene(name, product, [promoterConstitutive])
    return gene.unresolved_copy()


  def construct_factor(self, name, deacay, diffusibility) :
    """Method to construct the factor list of the transsys program.

    @return: The factor list of the transsys program.
    @rtype: a C{'transsys.Factor'} object
    """
    factor = transsys.Factor(name, transsys.ExpressionNodeValue(deacay), transsys.ExpressionNodeValue(diffusibility))
    return factor.unresolved_copy()



class OptimisationControlParameters(translattice.ControlParameters) :
  """Class to collect and represent the optimisation control parameters.

  @ivar cycles: The number of optimisation cycles.
  @type cycles: C{int}
  @ivar offset: The offset for the perturbation function. The perturbation
  function takes values from a uniform distribution of the interval [-<offset>,
  <offset>].
  @type offset: C{int} of C{float}
  """

  def __init__(self, cycles, offset) :
    """The constructor of the class.

    The two parameters comprise the instance variables as well.
    """
    self.cycles = cycles
    self.offset = offset



class OptimisationBookKeeping(object) :
  """Class to collect the bookkeeping information during the optimisation
  procedure.

  """

  def __init__(self) :
    self.result = []

  def append_record(self, optCycle, currentOO, optObjective, currentExprValues, bestExprValues, currentTP, bestTP) :
    """Method to aggregate the results from the optimiser.

    """
    self.result.append((optCycle, currentOO, optObjective, currentExprValues, bestExprValues, currentTP, bestTP,))


  def __str__(self) :
    """Prints the optimisation objective (current and best) and at the end
    append the best transsys program.

    """
    string = 'OptCycle:\tCurrentObjective:\tOptimisationObjective:\n'
    for i in xrange(len(self.result)) :
      string = string + '%i\t%f\t%f\n' % (self.result[i][0], self.result[i][1], self.result[i][2])
    string = string + str(self.result.pop()[6]) + '\n'
    return string


  def print_all(self) :
    """Prints all the contents of the OptimisationBookkeeping object.

    """
    string = ''
    for i in xrange(len(self.result)) :
      string = string + 'OptCycle: %i\nCurrentObjective: %f\nOptimisationObjective: %f\nCurrentExpressionValues: %s\nBestExpressionvalues: %s\nCurrentTranssysProgram:\n%s\nBestTranssysProgram:\n%s\n--------------------------------------------------------------------------------\n' % (self.result[i][0], self.result[i][1], self.result[i][2], self.result[i][3], self.result[i][4], self.result[i][5], self.result[i][6])
    return string




def get_perturbation_list(transsysProgram):
  """A function to get the list of values subject to random perturbations.

  Return a list of all ExpressionNodeValues of the genes of a transsys program.
  @param transsysProgram: A transsys program instance.
  @type transsysProgram: C{class 'transsys.TranssysProgram'}
  @ return: A list of all the expression node values of all genes.
  @rtype: C{list} of C{float}s.
  """
  return transsysProgram.getGeneValueNodes()


def perturb_value(value, perturbObj) :
  """Return a value perturbed by a random number.

  The value is perturbed according to the formula:
    - FC_perturbed = factor_concentration * exp(random_value)
    Where the random_value is get from a uniform [-z, z] distribution.
  @param value: A real value.
  @type value: C{float}
  @param perturbObj: An object of the L{translattice.RandomUniform} class
  @type perturbObj: C{class 'translattice.RandomUniform'}
  @return : the perturbed value.
  """
  perturbValue = value * math.exp(perturbObj.random_value())
  return perturbValue


def perturb_transsys_dummy(transsysProgram, perturbObj, cycle) :
  """Return a perturbed transsys program dummy. (the engineered transsys
  program.
  """
  engPar = transsysProgram.engineeringParameters
  perturbLowPoint = Point(perturb_value(engPar.get_lowPoint().get_x(), perturbObj), perturb_value(engPar.get_lowPoint().get_x(), perturbObj))
  perturbHighPoint = Point(perturb_value(engPar.get_highPoint().get_x(), perturbObj), perturb_value(engPar.get_highPoint().get_y(), perturbObj))
  perturbCentre = Point(perturb_value(engPar.get_circle().get_centre().get_x(), perturbObj), perturb_value(engPar.get_circle().get_centre().get_x(), perturbObj))
  perturbCircle = Circle(perturbCentre, perturb_value(engPar.get_circle().get_r(), perturbObj))
  engPar.set_perturbed_parameters(perturbLowPoint, perturbHighPoint, perturbCircle)
  return TranssysProgramDummy('engineered' + str(cycle), engPar)



def perturb_transsys(transsysProgram, perturbObj):
  """Return a perturbed transsys program.

  Perturbation refers to multiply by a random generated value the constant
  expression values of a transsys program.
  @param transsysProgram: A transsys program instance.
  @type transsysProgram: C{class 'transsys.TranssysProgram'}
  @param perturbObj: A random object.
  @type perturbObj: C{class 'random.Random'}
  """
  # First copy the transsys program.
  perturbedTP = copy.deepcopy(transsysProgram)
  for exprValue in get_perturbation_list(perturbedTP) :
    # Multiply the current value with the exponential of the random value
    exprValue.value = exprValue.value * math.exp(perturbObj.random_value())
  return perturbedTP



def zero_transsys_diffusibility(transsysProgram):
  """Set to zero the diffusibility expression of all factors.

  @param transsysProgram: An instance of a transsys program.
  """
  for diffusion in transsysProgram.getDiffusibilityValueNodes() :
    diffusion.value = 0



def expression_table(transsysLattice) :
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



def run_lattice(transsysProgram, latticeSize, timesteps, rndObj) :
  """Function to conduct the lattice experiment.

  It either runs the lattce or the control (diffusibility = 0), returns the
  factor expression table.
  """
  # Set the transsys program to a local variable tp) to protect it from
  # transforming the diffusion.
  transLat = translattice.TranssysInstanceLattice(transsysProgram, latticeSize, timesteps)
#  rndObj = translattice.UniformRNG(rndSeed, unifBorders[0], unifBorders[1])
  transLat.perturb_lattice(rndObj)
  # Set always the same random seed for the clib.
  transsys.clib.srandom(rndObj.rndSeed)
  transTS = translattice.TranssysLatticeTimeseries(transLat, timesteps, timesteps)
  lattice = transTS.latticeTimeseries.pop()
  return lattice


def bimodality(lattice) :
  """A function which wraps all the bimodalities module functions and
  calculate the total bimodality from an expression table.

  Returns a scalar.
  @param lattice: A TranssysInstanceLattice object.
  @type lattice: C{class translattice.TranssysInstanceLattice'}
  @return: The total bimodality score.
  @rtype: C{'float'}
  """
  expressionTable = expression_table(lattice)
  return bimodalities.total_bimodality(expressionTable)


def objective(tp, latticeSize, timesteps, initialObj):
  """A wrapper function which computes the optimisation objective.

  Calls the run_lattice function once for the lattice and once for the control
  experiment calculates the bimodalities and returns the optimisation
  objective.

  @return: The optimisation objective.
  @type: C{float}
  """
  # Copy the random initalisation object to preserve the same initial conditions
  # for the lattice and the control.
  initialObj2 = copy.deepcopy(initialObj)

  # The lattice experiment.
  exprTableLattice = run_lattice(tp, latticeSize, timesteps, initialObj)
  bmScoreLattice = bimodality(exprTableLattice)

  # The control.
  tpZero = copy.deepcopy(tp)
  zero_transsys_diffusibility(tpZero)
  exprTableControl = run_lattice(tpZero, latticeSize, timesteps, initialObj2)
  bmScoreControl = bimodality(exprTableControl)

  objective = bmScoreLattice - bmScoreControl
  return objective



def transsys_lattice_optimisation(transsysProgram, simulatorCP, optimiserCp, log = None) :
  """Conduct the optimisation experiment.

  The function implements all the necessary steps to run a complete optimisation
  experiment and returns an "optimised" transsys program.
  @param transsysProgram: A transsys program instance.
  @type transsysProgram: C{class 'transsys.TranssysProgram'}
  @param simulatorCP: An 'translattice.SimulatorControlParameters' object which
  contains the lattice simulator's control parameters. See docstring of
  L{translattice.SimulatorControlParameters}.
  @type simulatorCP: C{class 'translattice.SimulatorControlParameters'}
  @param optimiserCP: An 'OptimisationControlParameters' object which contains
  the optimisation procedure control parameters (i.e. the optimisation cycles
  and the perturbation offset). See L{OptimisatioControlParameters}.
  @type optimiserCP: C{class 'OptimisationControlParameters'}
  @param log: A L{OptimisationBookKeeping} instance. Holds the optimisation's
  results.
  @type log: C{class 'OptimisationBookKeeping'}
  @todo: Derive a log keeping class to store values during the optimisation
  procedure.
  """
  # Get the parameters.
  perturbOffset = optimiserCp.offset
  optimisationCycles = optimiserCp.cycles
  initialNoise = simulatorCP.initialisationvariables
  latticeSize = simulatorCP.latticeSize
  timesteps = simulatorCP.timesteps
  rndSeed = simulatorCP.randomSeed
  if not log :
    log = OptimisationBookKeeping()

  # Set the random objects.
  perturbObj = translattice.UniformRNG(rndSeed, -perturbOffset, perturbOffset)
  if isinstance(initialNoise, translattice.UniformParameters) :
    initialObj = translattice.UniformRNG(rndSeed, initialNoise.lower, initialNoise.upper)
  elif isinstance(initialNoise, translattice.GaussianParameters) :
    initialObj = translattice.GaussianRNG(rndSeed, initialNoise.mean, initialNoise.stddev)
  else :
    raise TypeError, 'The initial noise argument is not valid. Should be either an UniformParameters or a GaussianParameters instance.'

  # Calcullate the initial transsys program optimisation objective.
  optimisationObjective = objective(transsysProgram, latticeSize, timesteps, initialObj)
  bestTP = transsysProgram

  # BookKeeping.
  log.append_record(0, optimisationObjective, optimisationObjective, [expr.value for expr in bestTP.getGeneValueNodes()], [expr.value for expr in bestTP.getGeneValueNodes()], bestTP, bestTP,)

  # Basic optimisation loop.
  for cycle in xrange(optimisationCycles) :
    if isinstance(bestTP, TranssysProgramDummy) :
      alternativeTP = perturb_transsys_dummy(bestTP, perturbObj, cycle + 1)
    else :
      alternativeTP = perturb_transsys(bestTP, perturbObj)
    oo = objective(alternativeTP, latticeSize, timesteps, initialObj)
    if oo > optimisationObjective :
      optimisationObjective = oo
      bestTP = alternativeTP
    log.append_record(cycle + 1, oo, optimisationObjective, [expr.value for expr in bestTP.getGeneValueNodes()], [expr.value for expr in alternativeTP.getGeneValueNodes()], alternativeTP, bestTP,)

