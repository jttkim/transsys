#!/usr/bin/env python

# Subversion keywords.
# $Rev:: 305           $:  Revision of last commit
# $Author:: cbouyio    $:  Author of last commit
# $Date: 2007-12-20 20:55:15 +0000 (Thu, 20 Dec 2007) $:  Date of last commit


"""Python module implementing the random local search optimisation approach.

@author: Costas Bouyioukos
@organization: University of East Anglia
@since: 20/12/2007
@license: GNU General Public Licence 3 or newer.
@contact: U{Costas Bouyioukos<mailto:konsb@cmp.uea.ac.uk>}
@version: $Id$"""


# Version Information.
__version__ = "$Id$"


import copy
import sys
import math

import transsys
import translattice
import bimodalities
import transsys.optim


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



class EngineeringParameters(object) :
  """A class to aggregate the control parameters of the engineered transsys
  program.

  The transsys program has the following structure:

  transsys engineered
  {
    factor A
    {
      decay: <decayA>;
      diffusibility: <diffusibilityA>;
    }
    factor B
    {
      decay: <decayB>;
      diffusibility: <diffusibilityB>;
    }

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


  def __str__(self) :
    """
    """
    s = transsys.utils.dictionary_tablestring(self.__dict__)
    return s



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
    """Perturb the parameters of a transsys program dummy.

    The actual parameters that getting perturbed are the high and low
    expression points and the control circle (centre and radius).
    @param perturbObj: The perturbation object.
    @type perturbObj: C{class 'translattice.RandomObject'}
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


  def perturb_transsys(self, perturbObj) :
    """Method to perturb a TranssysProgramDummy.

    Works by calling the perturb_eng_parameters method.
    Polymorphic with the perturb_transsys function.
    """
    pass



  def get_eng_parameters(self) :
    """Accessor to the engineered control parameters of a TranssysProgramDummy.

    @return: The control parameters of an engineered transsys program.
    @rtype: C{class 'EngineeringParameters'}
    """
    pass



class OptimisationControlParameters(object) :
  """Class to collect and represent the optimisation control parameters.

  @ivar cycles: The number of optimisation cycles.
  @type cycles: C{int}
  @ivar perturbationOrder: The order for the perturbation function. The exponent
  of the perturbation function takes values from a random uiform distribution
  of the interval [-<perturbationOrder>, <perturbationOrder>].
  @type perturbationOrder: C{int} of C{float}
  @ivar rndParam: the random number generator parameter for the optimisation
  procedure. The random seed is specified as a <rndSeed> = <optCycle> +
  <rndParam>
  @type rndParam: C{int} or C{float}
  """

  def __init__(self, cycles, perturbationOrder, rndParam) :
    """The constructor of the class.

    The two parameters comprise the instance variables as well.
    """
    self.cycles = cycles
    self.perturbationOrder = perturbationOrder
    self.rndParam = rndParam


  def __str__(self) :
    """
    """
    s = transsys.utils.dictionary_tablestring(self.__dict__)
    return s



class OptimisationBookKeeping(object) :
  """Class to collect and print out information during the optimisation
  procedure.

  @ivar numericalFile: An output file where all the optimisation's (current
  and best) numerical values will be printed out.
  @type numericalFile: An open ready to write C{file}
  @ivar transsysFile: An output file where the transsys programs that
  have been generated by the optimisation process will be printed out.
  @type transsysFile: An open ready to write C{file}
  """

  def __init__(self, scp, ocp, tp, transsysFile = None, numericalFile = None) :
    if transsysFile :
      self.tf = transsysFile
    else :
      self.tf = sys.stdout
    if numericalFile :
      self.nf = numericalFile
    else :
      self.nf = None
    self.print_control_parameters(scp, ocp)
    self.print_header(tp)


  def print_control_parameters(self, simulatorCP, optimiserCP) :
    """Print all the optimisation procedure control parameters.

    The control parameters are printed as comments in the top of the output
    files.
    """
    if self.nf :
      self.nf.write(str(simulatorCP))
      self.nf.write(str(optimiserCP))
    self.tf.write(str(simulatorCP))
    self.tf.write(str(optimiserCP))


  def print_header(self, tp) :
    """Print headers in the output files.

    The headers are printed in order to help parsing the results with other
    applications (i.e. R).
    """
    if self.nf :
      if isinstance(tp, TranssysProgramDummy) :
        self.nf.write('OptCycle\tRNGSeed\tOptFlag\tAltObj\tAltLatBM\tAltCtrlBM\tBestObj\tBestLatBM\tBestCtrlBM\tAltLowX\tAltLowY\tAltHighX\tAltHighY\tAltCircX\tAltCircY\tAltRadius\tBestLowX\tBestLowY\tBestHighX\tBestHighY\tBestCircX\tBestCircY\tBestRadius\tLatBest_muA\tLatBest_stdevA\tLatBest_muB\tLatBest_stdevB\n')
      else :
        self.nf.write('OptCycle\tRNGSeed\tOptFlag\tAltObj\tAltLatBM\tAltCtrlBM\tBestObj\tBestLatBM\tBestCtrlBM')
        for fName in tp.factor_names() :
          self.nf.write('\tAverage_%s\tStddev_%s\tEntropy_%s\tMin_%s\tMax_%s' % (fName, fName, fName, fName, fName))
        self.nf.write('\n')


  def write_log(self, optCycle, rndSeedCurrent, optStep, altObj, bestObj, altTP, bestTP, altENGParam = None, bestENGParam = None) :
    """Curry out the printing of the optimisers results.

    Actually write_log is a wrapper method which calls print_objectives and
    print_parameters methods respectivelly.
    @param optCycle: The optimisation cycle.
    @type optCycle: C{int}
    @param rndSeedCurrent: The current run's random number generator seed.
    @type rndSeedCurrent: C{int}
    @param optStep: A boolean designates whether the optimiser has actually made
    an optimisation step.
    @type optStep: C{bool}
    @param altObj: The optimisation objective function score of the current
    alternative transsys program. (and some more things, a tuple)
    @type altObj: Some class (possibly objective class)
    @param bestObj: The optimisation objective function score of the current
    best transsys program. (and some more things, a tuple)
    @type bestObj: Some class (possibly objective class)
    @param altENGParam: The engineered control parameters of the current
    alternative transsys program.
    @type altENGParam: C{class 'EngineeringParameters'}
    @param bestENGParam: The engineered control parameters of the current best
    transsys program.
    @type bestENGParam: C{class 'EngineeringParameters'}
    @param altTP: The current alternative transsys program. That is the one
    generated in the particular optimisation cycle.
    @type altTP: C{class 'transsys.TranssysProgram'}
    @param bestTP: The current best transsys program. That is the one with the
    highest score of the objective function.
    @type bestTP: C{class 'transsys.TranssysProgram'}
    """
    self.print_transsys_programs(altTP, bestTP)
    if self.nf :
      self.print_numerical(optCycle, rndSeedCurrent, optStep, altObj, bestObj, altENGParam, bestENGParam)


  def print_numerical(self, optCycle, rndSeedCurrent, optStep, altObj, bestObj, altECP, bestECP) :
    """Method to print the optimisation objective values. Prints both the
    current objective value as well as the best objective value of the "best"
    transsys program found so far.
    @param optCycle: The optimisation cycle.
    @type optCycle: C{int}
    @param rndSeedCurrent: The random number generator seed.
    @type rndSeedCurrent: C{int}
    @param altObj: The current (this particular optimisation cycle) value of
    the objective function.
    @type altObj: Numeric
    @param bestObj: The overall optimisation objective value.
    @type bestObj: Numeric
    @param altECP: The current engineered control parameters.
    @type altECP: C{class 'EngineeringParameters'}
    @param bestECP: The best engineered control parameters.
    @type bestECP: C{class 'EngineeringParameters'}
    """
    if altECP and bestECP :
      statsBest = bestObj.lattice.statistics()
      self.nf.write('%i\t%i\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n' % (optCycle, rndSeedCurrent, str.upper(str(optStep)), altObj.objective, altObj.latticeBimodalities.totalBimodality, altObj.controlBimodalities.totalBimodality, bestObj.objective, bestObj.latticeBimodalities.totalBimodality, bestObj.controlBimodalities.totalBimodality, altECP.get_lowPoint().x, altECP.get_lowPoint().y, altECP.get_highPoint().x, altECP.get_highPoint().y, altECP.get_circle().get_centre().x, altECP.get_circle().get_centre().y, altECP.get_circle().r, bestECP.get_lowPoint().x, bestECP.get_lowPoint().y, bestECP.get_highPoint().x, bestECP.get_highPoint().y, bestECP.get_circle().get_centre().x, bestECP.get_circle().get_centre().y, bestECP.get_circle().r, statsBest.average[0], statsBest.standard_deviation[0], statsBest.average[1], statsBest.standard_deviation[1]))
    else :
      statsBest = bestObj.lattice.statistics()
      self.nf.write('%i\t%i\t%s\t%f\t%f\t%f\t%f\t%f\t%f' % (optCycle, rndSeedCurrent, str.upper(str(optStep)), altObj.objective, altObj.latticeBimodalities.totalBimodality, altObj.controlBimodalities.totalBimodality, bestObj.objective, bestObj.latticeBimodalities.totalBimodality, bestObj.controlBimodalities.totalBimodality))
      for i in xrange(statsBest.transsys_program.num_factors()) :
        self.nf.write('\t%f\t%f\t%f\t%f\t%f' % (statsBest.average[i], statsBest.standard_deviation[i], statsBest.shannon_entropy[i], statsBest.min_factor_concentration[i], statsBest.max_factor_concentration[i]))
      self.nf.write('\n')



  def print_transsys_programs(self, alternativeTP, currentBestTP) :
    """Print the transsys programs that are participating in to any optimisation
    cycle.

    This method writes in a transsys file (.tra) both the current alternative
    and the CURRENT best transsys programs of every cycle of the optimiser, the
    output file is fully compatible to the transsys program parser.
    @param alternativeTP: The transsys program that has generated in the
    particular optimisation cycle.
    @type alternativeTP: C{class 'transsys.TranssysProgram'}
    @param currentBestTP: The current best transsys program, that is the one
    that its objective is the one calculated in the current optimisation cycle.
    @type currentBestTP: C{class 'transsys.TranssysProgram'}
    """
    self.tf.write('%s\n%s\n' % (str(alternativeTP), str(currentBestTP)))



class OptimisationObjective(object) :
  """A class to wrap the calculations, represent the optimisation objective
  score and to keep other usefull ststistics (e.g. lattice and control
  bimodalities etc.).

  @ivar objective: The optimisation objective score.
  @type objective: Numeric
  @ivar latticeBimodalities: A bimodalities object for the lattice experiment.
  @type latticeBimodalities: C{class 'bimodalities.BimodalitiesCollection'}
  @ivar controlBimodalities: A bimodalities object for the control experiment.
  @type controlBimodalities: C{class 'bimodalities.BimodalitiesCollection'}
  @ivar lattice: The outcome of a transsys lattice simulation
  experiment.
  @type lattice: C{class 'translattice.TranssysInstanceLattice'}
  """

  def __init__(self, tp, latticeSize, timesteps, initialNoise, rndSeed) :
    """The constructor of the class.

    """
    # The lattice lattice simulation experiment.
    self.lattice = self.run_lattice(tp, latticeSize, timesteps, initialNoise, rndSeed)
    # Construct the lattice bimodalities object.
    self.latticeBimodalities = bimodalities.BimodalitiesCollection(self.lattice)
    # The control simulation experiment.
    zeroTP = set_transsys_diffusibilities(tp, 0)
    contolLattice = self.run_lattice(zeroTP, latticeSize, timesteps, initialNoise, rndSeed)
    # Construct the control experiment bimodalities object.
    self.controlBimodalities = bimodalities.BimodalitiesCollection(contolLattice)
    # Calculate the objective.
    self.objective = self.latticeBimodalities.totalBimodality - self.controlBimodalities.totalBimodality



  def run_lattice(self, tp, latticeSize, timesteps, initialNoise, rndSeed) :
    """Method to conduct the lattice experiment.

    The method curries out the lattice simulation experiment with the specified
    parameters.
    @param tp: The transsys program.
    @type tp: C{class 'transsys.TranssyProgram'}
    @param latticeSize: The parameter holding the size of the lattice.
    @type latticeSize: C{class 'translattice.LatticeSize'}
    @param timesteps: The number of timesteps for the simulator.
    @type timesteps: C{int}
    @param initialNoise: An UniformParameters object which holds the
    initialisation parameters for the lattice.
    @type initialNoise: C{class 'translattice.UniformParameters'}
    @param rndSeed: The random number generator seed.
    @type rndSeed: C{int}
    @return: A lattice after running the simulator for the specified number of
    timesteps.
    @rtype: C{class 'translattice.TranssysInstanceLattice'}
    """
    # Construct the initialisation object.
    initialiseObj = initialNoise.getRNG(rndSeed)
    # Instantiate a transsys lattice.
    transLat = translattice.TranssysInstanceLattice(tp, latticeSize, timesteps)
    transLat.perturb_lattice(initialiseObj)
    # Set always the same random seed for the clib.
    transsys.clib.srandom(initialiseObj.rndSeed)
    # Run the time-series.
    transTS = translattice.TranssysLatticeTimeseries(transLat, timesteps, timesteps)
    # Get the equilibration lattice (i.e. the last timestep)
    lattice = transTS.latticeTimeseries.pop()
    return lattice



def set_transsys_diffusibilities(transsysProgram, diffusionLimit) :
  """Assign a diffusibility expression of each factor of atranssys program.

  The diffusibility is calculated as a random number from the interval [0,
  diffusionLimit).
  If the diffusionLimit = 0 (the control experiment) then a zero (0)
  diffusibility expression is assigned to ALL the factors.
  Also notice that the seed for the random number generator is always 1, since
  we need to assigne diffusibilities only once in the begining of the
  optimisation.
  @param transsysProgram: An instance of a transsys program.
  @type transsysProgram: C{class 'transsys.TranssysProgram'}
  @param diffusionLimit: the upper limit for the diffusibility expression.
  Diffusibility is allowed to take values from the [0, 1] interval
  @type diffusionLimit: Numeric
  @return: A copy of the transsysProgram with zero diffusibilities.
  @rtype: C{class 'transsys.TranssysProgram'}
  """
  if diffusionLimit == 0 :
    zeroTP = copy.deepcopy(transsysProgram)
    for diffusion in zeroTP.getDiffusibilityValueNodes() :
      diffusion.value = 0
    return zeroTP
  else :
    diffList = []
    rng = translattice.UniformRNG(1, 0, diffusionLimit)
    for diffusion in transsysProgram.getDiffusibilityValueNodes() :
      diffusion.value = rng.random_value()
    return transsysProgram



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

  FC_perturbed = factor_concentration * exp(random_value)

  Where the random_value is got from a uniform [-perturbationOrder,
  perturbationOrder] distribution.
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



def optimisation(engineeredCP, simulatorCP, optimiserCP, logObj) :
  """Curry out the optimisation experiment.

  This is the wrapper function for the optimiser. this function integrates all
  the nessecary steps to run a random local search experiment on an engineered
  transsys program (a TranssysProgramDummy object, actually the function
  "optimises" the control parameters of the engineered transsys program). All
  the calculations are curried out in the relative methods/functions above here
  anly the evaluation is taking place. The function returns nothing but writes
  in the log object after each optimisation cycle.
  @param engineeredCP: The TranssyProgramDummy (engineered) control parameters.
  @type engineeredCP: C{class 'EngineeringParameters'}
  @param simulatorCP: An 'translattice.SimulatorControlParameters' object which
  contains the lattice simulator's control parameters. See docstring of
  L{translattice.SimulatorControlParameters}.
  @type simulatorCP: C{class 'translattice.SimulatorControlParameters'}
  @param optimiserCP: An 'OptimisationControlParameters' object which contains
  the optimisation procedure control parameters (i.e. the optimisation cycles
  and the perturbation perturbationOrder). See L{OptimisationControlParameters}.
  @type optimiserCP: C{class 'OptimisationControlParameters'}
  @param logObj: An object to keep logging and print out the optimisation
  results.
  @type logObj: C{class 'OptimisationBookKeeping'}
  @return: Nothing (But writes all the numerical output and the generated
  transsys programs in a numerical and a .tra output file respectively)
  @rtype: C{None}
  """
  # Get the parameters.
  perturbRange = optimiserCP.perturbationOrder
  optimisationCycles = optimiserCP.cycles
  rndSeed = optimiserCP.rndParam
  initialNoise = simulatorCP.initialisationVariables
  latticeSize = simulatorCP.latticeSize
  timesteps = simulatorCP.timesteps
  engBest = engineeredCP

  # Basic optimisation loop.
  for cycle in xrange(optimisationCycles) :
    # Boolean to check whether an optimisation step is succsefull or not.
    optStep = False
    # Set the current random seed.
    rndSeedCurrent = optimiserCP.rndParam + cycle
    # Set the random objects.
    perturbObj = translattice.UniformRNG(rndSeedCurrent, -perturbRange, perturbRange)
    # Polymorphism will take care of the random object initalisation.
    initialObj = initialNoise.getRNG(rndSeedCurrent)

    # Perturb the engineered control parameters. (here is the actual
    # implementation of the random local search)
    engPAlternative = engBest.perturb_eng_parameters(perturbObj)
    # Generate the transsys programs.
    tpAlternative = TranssysProgramDummy('tpAlternative_' + str(cycle + 1), engPAlternative)
    tpAlternative.comments = ['RNG Seed: %i' % rndSeedCurrent]
    tpBest = TranssysProgramDummy('currentBest_' + str(cycle + 1), engBest)
    tpBest.comments = ['RNG Seed: %i' % rndSeedCurrent]

    # Calculate the objectives of both the best and the alternative transsys
    # programs, with the same initial conditions.
    objectiveBest = OptimisationObjective(tpBest, latticeSize, timesteps, initialNoise, rndSeedCurrent)
    objectiveAlternative = OptimisationObjective(tpAlternative, latticeSize, timesteps, initialNoise, rndSeedCurrent)

    # If an optimisation step occurs.
    if objectiveAlternative.objective > objectiveBest.objective :
      optStep = True
    # Print the optimisation's log.
    logObj.write_log(cycle + 1, rndSeedCurrent, optStep, objectiveAlternative, objectiveBest, tpAlternative, tpBest, engPAlternative, engBest)
    # Change the control parameters (the random local search)
    if objectiveAlternative.objective > objectiveBest.objective :
      engBest = copy.deepcopy(engPAlternative)



def tp_optimisation(tp, sCP, oCP, logObj) :
  """Curry out the optimisation experiment.

  This is the wrapper function for the optimiser. this function integrates all
  the nessecary steps to run a random local search experiment on an engineered
  transsys program (a TranssysProgramDummy object, actually the function
  "optimises" the control parameters of the engineered transsys program). All
  the calculations are curried out in the relative methods/functions above here
  only the evaluation is taking place. The function returns nothing but writes
  in the log object after each optimisation cycle.
  @param tp: The transsys program that will be subject to optimisation.
  @type tp: C{class 'transsys.TranssysProgram'}
  @param sCP: An 'translattice.SimulatorControlParameters' object which
  contains the lattice simulator's control parameters. See docstring of
  L{translattice.SimulatorControlParameters}.
  @type sCP: C{class 'translattice.SimulatorControlParameters'}
  @param oCP: An 'OptimisationControlParameters' object which contains
  the optimisation procedure control parameters (i.e. the optimisation cycles
  and the perturbation perturbationOrder). See L{OptimisationControlParameters}.
  @type oCP: C{class 'OptimisationControlParameters'}
  @param logObj: An object to keep logging and print out the optimisation
  results.
  @type logObj: C{class 'OptimisationBookKeeping'}
  @return: Nothing (But writes all the numerical output and the generated
  transsys programs in a numerical and a .tra output file respectively)
  @rtype: C{None}
  """
  # Get the parameters.
  perturbRange = oCP.perturbationOrder
  optimisationCycles = oCP.cycles
  rndSeed = oCP.rndParam
  initialNoise = sCP.initialisationVariables
  latticeSize = sCP.latticeSize
  timesteps = sCP.timesteps
  tpBest = tp
  tpName = tpBest.name
  # Basic optimisation loop.
  for cycle in xrange(optimisationCycles) :
    # Boolean to check whether an optimisation step is succsefull or not.
    optStep = False
    # Set the current random seed.
    rndSeedCurrent = oCP.rndParam + cycle
    # Set the random object.
    perturbObj = translattice.UniformRNG(rndSeedCurrent, -perturbRange, perturbRange)
    # Perturbe the transsys program.
    tpAlternative = perturb_transsys(tpBest, perturbObj)
    tpAlternative.name = tpName + 'Alternative_' + str(cycle + 1)
    tpAlternative.comments = ['RNG Seed: %i' % rndSeedCurrent]
    tpBest.name = tpName + 'CurrentBest_' + str(cycle + 1)
    tpBest.comments = ['RNG Seed: %i' % rndSeedCurrent]
    # Calculate the objectives of both the best and the alternative transsys
    # programs, with the same initial conditions.
    objectiveBest = OptimisationObjective(tpBest, latticeSize, timesteps, initialNoise, rndSeedCurrent)
    objectiveAlternative = OptimisationObjective(tpAlternative, latticeSize, timesteps, initialNoise, rndSeedCurrent)
    # If an optimisation step occurs.
    if objectiveAlternative.objective > objectiveBest.objective :
      optStep = True
    # Print the optimisation's log.
    logObj.write_log(cycle + 1, rndSeedCurrent, optStep, objectiveAlternative, objectiveBest, tpAlternative, tpBest)
    # Alter the transsys programs (the random local search)
    if objectiveAlternative.objective > objectiveBest.objective :
      tpBest = copy.deepcopy(tpAlternative)

