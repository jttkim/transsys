#!/usr/bin/env python

# Subversion keywords.
# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date: 2007-12-20 20:55:15 +0000 (Thu, 20 Dec 2007) $:  Date of last commit

"""Python module implementing the random local search optimisation approach.

@author: Costas Bouyioukos
@organization: University of East Anglia
@since: 20/12/2007
@license: GNU General Public Licence 3 or newer.
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
import string

# Random objects already implemented in the transsysLattice module.


class Point(translattice.Parameter) :
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
    self.x = x
    self.y = y



class Circle(translattice.Parameter) :
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
    """Returns the radius of the circle.
    """
    return self.r


  def set_values(self, centre, r) :
    """Sets the parameter values.

    """
    if not isinstance(centre, Point) :
      raise StandardError, 'Invalid parameter setting for the circle. The centre should be of type "Point"'
    self.centre = copy.deepcopy(centre)
    self.r = r



class EngineeringParameters(translattice.ControlParameters) :
  """A class to aggregate the control parameters of the engineered transsys
  program.

  The transsys program has the structure:
  transsys engineered
  {factor A {decay: <decayA>; diffusibility: <diffusibilityA>;}
  {factor B {decay: <decayB>; diffusibility: <diffusibilityB>;}

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
      constitutive: (<lowPoint.y>/<decayB> + ((<highPoint.y>/<decayB> - <lowPoint.y>/<decayB>) * ((((A - <cirlce.x) * (A - <circle.x>)) + ((B - <circle.y>) * (B - <circle.y>))) <= (<circle.r> * <circle.r>))));
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
    self.decayB = decayB
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


  def perturb_eng_parameters(self, perturbObj) :
    """Perturb the parameters of a transsys program dummy

    Actualy the things that are get perturbed are the high and low expression
    points and the circle (centre and radius).
    @param perturbObj: The perturbation object.
    @type perturbObj: C{class 'translattice.RandomObject'}
    @param cycle: The optimisation cycle. Used only for naming purposes.
    @type cycle: C{int}
    @return: A perturbed copy of the engineered control parameters.
    @rtype: C{class 'EngineeringParameters'}
    """
    engP = copy.deepcopy(self)
    perturbLowPoint = Point(perturb_value(engP.get_lowPoint().get_x(), perturbObj), perturb_value(engP.get_lowPoint().get_x(), perturbObj))
    perturbHighPoint = Point(perturb_value(engP.get_highPoint().get_x(), perturbObj), perturb_value(engP.get_highPoint().get_y(), perturbObj))
    perturbCentre = Point(perturb_value(engP.get_circle().get_centre().get_x(), perturbObj), perturb_value(engP.get_circle().get_centre().get_x(), perturbObj))
    perturbCircle = Circle(perturbCentre, perturb_value(engP.get_circle().get_r(), perturbObj))
    engP.set_perturbed_parameters(perturbLowPoint, perturbHighPoint, perturbCircle)
    return engP



class TranssysProgramDummy(transsys.TranssysProgram) :
  """Subclass of transsys program to build "from scratch" a transsys program
  from the given control parameters set.

  """

  def __init__(self, name, engineeringParameters) :
    """The constructor of the class.

    Check the docstring of the EngineeringParameters class to see the
    structure of the transsys program.
    """
    self.name = name
    if isinstance(engineeringParameters, EngineeringParameters) :
      self.engineeringParameters = engineeringParameters
    else :
      raise StandardError, 'Error in engineering parameters setting. The engineering parameters should be an EngineeringParameters object.'
    self.factor_list = [self.construct_factor('A', self.engineeringParameters.decayA, self.engineeringParameters.diffusibilityA), self.construct_factor('B', self.engineeringParameters.decayB, self.engineeringParameters.diffusibilityB)]
    self.gene_list = [self.construct_gene('geneA', 'A', self.engineeringParameters.lowPoint.x, self.engineeringParameters.highPoint.x, self.engineeringParameters.circle, self.engineeringParameters.decayA), self.construct_gene('geneB', 'B', self.engineeringParameters.lowPoint.y, self.engineeringParameters.highPoint.y, self.engineeringParameters.circle, self.engineeringParameters.decayB)]
    transsys.TranssysProgram.__init__(self, self.name, self.factor_list, self.gene_list)


  def construct_gene(self, name, product, lp, hp, circle, decay) :
    """Method to construct the complex constitutive expression out of the
    system's points coordinates.

    Builds up all the complex expression nodes for each gene of the system, by
    constructing a set of expression nodes instances.

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


  def construct_factor(self, name, decay, diffusibility) :
    """Method to construct the factor list of the transsys program.

    @return: The factor list of the transsys program.
    @rtype: A C{'transsys.Factor'} object
    """
    factor = transsys.Factor(name, transsys.ExpressionNodeValue(decay), transsys.ExpressionNodeValue(diffusibility))
    return factor.unresolved_copy()


  def get_eng_parameters(self) :
    """Accessor to the engineered control parameters of a TranssysProgramDummy.

    @return : The control parameters of an engineered transsys program.
    @rtype : C{class 'EngineeringParameters'}
    """
    pass



class OptimisationControlParameters(translattice.ControlParameters) :
  """Class to collect and represent the optimisation control parameters.

  @ivar cycles: The number of optimisation cycles.
  @type cycles: C{int}
  @ivar offset: The offset for the perturbation function. The perturbation
  function takes values from a uniform distribution of the interval [-<offset>,
  <offset>].
  @type offset: C{int} of C{float}
  @ivar rndParam: the random number generator parameter for the optimisation
  procedure. The random seed is specified as a <rndSeed> = <optCycle> +
  <rndParam>
  @type rndParam: C{int} or C{float}
  """

  def __init__(self, cycles, offset, rndParam) :
    """The constructor of the class.

    The two parameters comprise the instance variables as well.
    """
    self.cycles = cycles
    self.offset = offset
    self.rndParam = rndParam



class OptimisationBookKeeping(object) :
  """Class to collect and print out information during the optimisation
  procedure.

  @ivar numericalFile: An output filename where all the optimisation's (current
  and best) numerical values will be printed out.
  @type numericalFile: C{str}
  @ivar transsysFile: An output filename where the transsys programs that
  have been generated by the optimisation process will be printed out.
  @type transsysFile: C{str}
  """

  def __init__(self, numericFilename, transsysFilename, scp, ocp) :
    self.numeric = numericFilename
    self.transsys = transsysFilename
    self.print_control_parameters(scp, ocp)
    self.print_header()


  def print_control_parameters(self, simulatorCP, optimiserCP) :
    """Print all the optimisation procedure control parameters.

    The control parameters are printed as comments.
    """
    nf = open(self.numeric, 'w')
    tf = open(self.transsys, 'w')
    nf.write(str(simulatorCP))
    nf.write(str(optimiserCP))
    tf.write(str(simulatorCP))
    tf.write(str(optimiserCP))
    nf.close()
    tf.close()


  def print_header(self) :
    """Print headers in the output files.

    The headers are printed in order to help parsing the results with other
    applications (i.e. R).
    """
    nf = open(self.numeric, 'a')
    nf.write('OptCycle\tCurrRNDSeed\tBestRNDSeed\tCurrObj\tCurrLatBM\tCurrCtrlBM\tBestObj\tBestLatBM\tBestCtrlBM\tCurrLowX\tCurrLowY\tCurrHighX\tCurrHighY\tCurrCircX\tCurrCircY\tCurrRadius\tBestLowX\tBestLowY\tBestHighX\tBestHighY\tBestCircX\tBestCircY\tBestRadius\n')
    nf.close()


  def write_log(self, optCycle, rndSeedCurrent, optStep, currentOO, optObjective, engParamCurrent, engParamBest, currentTP, bestTP) :
    """Curry out the printing of the optimisers results.

    Actually write_log is a wrapper method which calls print_objectives and
    print_parameters methods respectivelly.
    @param optCycle: The optimisation cycle.
    @type optCycle: C{int}
    @param rndSeedCurrent: The random number generator seed.
    @type rndSeed: C{int}
    @param optStep: A boolean designates whether the optimiser has actually made
    an optimisation step.
    @type currentOO: C{bool}
    @param currentOO: The optimisation objective score of the current
    alternative transsys program.
    @type currentOO: Numeric
    @param optObjective: The optimisation score of the current best transsys
    program.
    @type optObjective: Numeric
    @param engParamCurrent: The engineered control parameters of the current
    alternative transsys program.
    @type: C{class 'EngineeringParameters'}
    @param engBest: The engineered control parameters of the current best
    transsys program.
    @type engBest: C{class 'EngineeringParameters'}
    @param currentTP: The current alternative transsys program. That is the one
    generated in the particular optimisation cycle.
    @type currentTP: C{class 'transsys.TranssysProgram'}
    @param bestTP: The current best transsys program. That is the one with the
    highest score of the objective function.
    @type bestTP: C{class 'transsys.TranssysProgram'}
    """
    self.print_numerical(optCycle, rndSeedCurrent, optStep, currentOO, optObjective, engParamCurrent, engParamBest)
    self.print_transsys_programs(currentTP, bestTP)



  def print_numerical(self, optCycle, rndSeedCurrent, optStep, currentOO, optObjective, ePC, ePB) :
    """Method to print the optimisation objective values. Prints both the
    current objective value as well as the best objective value of the "best"
    transsys program found so far.
    @param optCycle: The optimisation cycle.
    @type optCycle: C{int}
    @param rndSeed: The random number generator seed.
    @type rndSeed: C{int}
    @param currentOO: The current (this particular optimisation cycle) value of
    the objective function.
    @type currentOO: Numeric
    @param optObjective: The overall optimisation objective value.
    @type optObjective: Numeric
    @param ePC: The current engineered control parameters.
    @type ePC: C{class 'EngineeringParameters'}
    @param ePB: The best engineered control parameters.
    @type ePB: C{class 'EngineeringParameters'}
    """
    nf = open(self.numeric, 'a')
    nf.write('%i\t%i\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' % (optCycle, rndSeedCurrent, str.upper(str(optStep)), currentOO[0], currentOO[1], currentOO[2], optObjective[0], optObjective[1], optObjective[2], ePC.get_lowPoint().x, ePC.get_lowPoint().y, ePC.get_highPoint().x, ePC.get_highPoint().y, ePC.get_circle().get_centre().x,  ePC.get_circle().get_centre().y, ePC.get_circle().r, ePB.get_lowPoint().x, ePB.get_lowPoint().y, ePB.get_highPoint().x, ePB.get_highPoint().y, ePB.get_circle().get_centre().x,  ePB.get_circle().get_centre().y, ePB.get_circle().r))
    nf.close()


  def print_transsys_programs(self, currentTP, bestTP) :
    """Print the transsys program parameters that are subject to perturbations
    (i.e. the ExpressionNodeValues of the genes). This method prints a list with
    all the perturbed values of both the current and the "best" perturbed
    program as well as the both (current and "best") transsys programs in every
    cycle of the optimiser.
    @param optCycle: The optimisation cycle.
    @type optCycle: C{int}
    @param rndSeed: The random number generator seed.
    @type rndSeed: C{int}
    @param currentExprValues: The current transsys program parameters that are
    subject to optimisation.
    @type currentExprValues: C{list} of {float}s
    @param bestExprValues: The parameters of the transsys program with the best
    optimisation score.
    @type bestExprValues: C{list} of {float}s
    @param currentTP: The transsys program that has generated in the particular
    optimisation cycle.
    @type currentTP: C{class 'transsys.TranssysProgram'}
    @param bestTP: The transsys program with the highest score of the objective
    function.
    @type bestTP: C{class 'transsys.TranssysProgram'}
    """
    tf = open(self.transsys, 'a')
    tf.write('%s\n%s\n' % (str(currentTP), str(bestTP)))
    tf.close()



def get_perturbation_list(transsysProgram):
  """A function to get the list of values subject to random perturbations.

  Return a list of all ExpressionNodeValues of the genes of a transsys program.
  @param transsysProgram: A transsys program instance.
  @type transsysProgram: C{class 'transsys.TranssysProgram'}
  @return: A list of all the expression node values of all genes.
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
  @param perturbObj: An object of an L{translattice.UniformRNG} class
  @type perturbObj: C{class 'translattice.UniformRNG'}
  @return: the perturbed value.
  @rtype: C{float}
  """
  perturbValue = value * math.exp(perturbObj.random_value())
  return perturbValue





def perturb_transsys(transsysProgram, perturbObj):
  """Return a perturbed transsys program.

  Perturbation refers to multiply by a random generated value the constant
  expression values of a transsys program (This function works only with
  TranssysProgram objects and not with any subclass of transsys programs).
  @param transsysProgram: A transsys program instance.
  @type transsysProgram: C{class 'transsys.TranssysProgram'}
  @param perturbObj: A random object.
  @type perturbObj: C{class 'random.Random'}
  @return: A transsys program with perturbed gene expression node values.
  @rtype: C{class 'transsys.TranssysProgram'}
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



def run_lattice(transsysProgram, latticeSize, timesteps, initialiseObj) :
  """Function to conduct the lattice experiment.

  It curries out the lattice simulation experiment with the specified
  parameters.
  @param transsysProgram: The transsys program.
  @type transsysProgram: C{class 'transsys.TranssyProgram'}
  @param latticeSize: The parameter holding the size of the lattice.
  @type latticeSize: C{class 'translattice.LatticeSize'}
  @param timesteps: The number of timesteps for the simulator.
  @type timesteps: C{int}
  @param initialiseObj: A randommnumber generator object specifying the initial
  condictions of the factor concentrations.
  @type initialiseObj: C{class 'translattice.RandomObject'}
  @return: A lattice after running a lattice timeseries for the specified number
  of timesteps.
  @rtype: C{class 'translattice.TranssysInstanceLattice'}
  """
  # Set the transsys program to a local variable tp) to protect it from
  # transforming the diffusion.
  transLat = translattice.TranssysInstanceLattice(transsysProgram, latticeSize, timesteps)
  transLat.perturb_lattice(initialiseObj)
  # Set always the same random seed for the clib.
  transsys.clib.srandom(initialiseObj.rndSeed)
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
  @rtype: C{float}
  """
  expressionTable = lattice.expression_table()
  return bimodalities.total_bimodality(expressionTable)


def objective(tp, latticeSize, timesteps, initialNoise, rndSeed):
  """A wrapper function which computes the optimisation objective.

  Calls the run_lattice function once for the lattice and once for the control
  experiment calculates the bimodalities and returns the optimisation
  objective.

  @return: The optimisation objective the lattice bimodality score and
  the control bimodality score.
  @rtype: C{tuple} of 3 C{float}s
  """
  # Construct the initial randomisation object for the lattice
  initialObj = initialNoise.getRNG(rndSeed)
  # Run the lattice experiment.
  lattice = run_lattice(tp, latticeSize, timesteps, initialObj)
  latticeBimodality = bimodality(lattice)

  # The control.
  # (re)Construct the initial randomisation object for the control. Note that
  # initialObj = initialObj2 to provide the same initial conditions for the
  # lattice and the control.
  initialObj2 = initialNoise.getRNG(rndSeed)
  # Copy the transsys program
  tpZero = copy.deepcopy(tp)
  # Zero the diffusibilities.
  zero_transsys_diffusibility(tpZero)
  # Run the control experiment.
  control = run_lattice(tpZero, latticeSize, timesteps, initialObj2)
  controlBimodality = bimodality(control)

  objective = latticeBimodality - controlBimodality
  return (objective, latticeBimodality, controlBimodality)



def optimisation(engineeredCP, simulatorCP, optimiserCP, logObj) :
  """Curry out the optimisation experiment.

  This function implements all the necessary steps to run a complete
  optimisation expreriment, it reports a bunch of optimisation results at the
  logging object.
  @param transsysProgram: A transsys program instance.
  @type transsysProgram: C{class 'transsys.TranssysProgram'}
  @param simulatorCP: An 'translattice.SimulatorControlParameters' object which
  contains the lattice simulator's control parameters. See docstring of
  L{translattice.SimulatorControlParameters}.
  @type simulatorCP: C{class 'translattice.SimulatorControlParameters'}
  @param optimiserCP: An 'OptimisationControlParameters' object which contains
  the optimisation procedure control parameters (i.e. the optimisation cycles
  and the perturbation offset). See L{OptimisationControlParameters}.
  @type optimiserCP: C{class 'OptimisationControlParameters'}
  @param logObj: An object to keep logging and print out the optimisation
  results.
  @type logObj: C{class 'OptimisationBookKeeping'}
  @rtype: C{None}
  """
  # Get the parameters.
  perturbOffset = optimiserCP.offset
  optimisationCycles = optimiserCP.cycles
  rndSeed = optimiserCP.rndParam
  initialNoise = simulatorCP.initialisationVariables
  latticeSize = simulatorCP.latticeSize
  timesteps = simulatorCP.timesteps
  transsysProgram = TranssysProgramDummy('engineered', engineeredCP)

  # Set the random objects.
  perturbObj = translattice.UniformRNG(rndSeed, -perturbOffset, perturbOffset)
  #Check for the inital noise parameter.
  if not isinstance(initialNoise, (translattice.UniformParameters, translattice.GaussianParameters,)) :
    raise TypeError, 'The initial noise argument is not valid. Should be either an UniformParameters or a GaussianParameters instance.'

  # Calculate the initial transsys program optimisation objective.
  objectiveInit = objective(transsysProgram, latticeSize, timesteps, initialNoise, rndSeed)
  engBest = engineeredCP
  tpBest = transsysProgram
  optStep = False

  # Print the optimisation's log for the initial run.
  logObj.write_log(0, rndSeed, optStep, objectiveInit, objectiveInit, engineeredCP, engBest, transsysProgram, tpBest)

  # Basic optimisation loop.
  for cycle in xrange(optimisationCycles) :
    # Boolean to check whether an optimisation step is succsefull or not.
    optStep = False
    # Set the current random seed.
    rndSeedCurrent = optimiserCP.rndParam + (cycle + 1)
    # Set the random objects.
    perturbObj = translattice.UniformRNG(rndSeedCurrent, -perturbOffset, perturbOffset)
    # Polymorphism will take care of the random object initalisation.
    initialObj = initialNoise.getRNG(rndSeedCurrent)

    # Perturb the engineered control parameters. (here is the actual
    # implementation of the random local search)
    engineeredCP = engBest.perturb_eng_parameters(perturbObj)
    # Generate the transsys programs.
    tpCurrent = TranssysProgramDummy('engineered' + str(cycle + 1), engineeredCP)
    tpBest = TranssysProgramDummy('engineeredBest' + str(rndSeedCurrent), engBest)

    # Calculate the objectives of both the best and the alternative transsys
    # programs, with the same initial conditions.
    objectiveBest = objective(tpBest, latticeSize, timesteps, initialNoise, rndSeedCurrent)
    objectiveCurrent = objective(tpCurrent, latticeSize, timesteps, initialNoise, rndSeedCurrent)

    # This is the optimisation condition.
    if objectiveCurrent[0] > objectiveBest[0] :
      engBest = copy.deepcopy(engineeredCP)
      optStep = True

    # Print the optimisation's log.
    logObj.write_log(cycle + 1, rndSeedCurrent, optStep, objectiveCurrent, objectiveBest, engineeredCP, engBest, tpCurrent, tpBest)

