#!/usr/bin/env python

"""
Here will be the translattice module description....
@author: Costas Bouyioukos
@organization: University of East Anglia
@since: Aprill 2006 (converted to a module by 08/09/2006)
@copyright: The program is coming as it is... You have the right to
redistribute, transform and change the source code presuming the apropriate
reference and the lisence is kept free.
@license: GNU GPL
@contact: U{Costas Bouyioukos<mailto:konsb@cmp.uea.ac.uk>}
"""

# latticeSimulator: A program that simulates a lattice of trannsys instances
# on a nxn lattice, each cell of the lattice cantains one transsys program.
# The program applies communication between neighbour cells and is an attempt
# to simulate the processes described in "Soile V. E. Keranen, 2004. Simulation
# study on effects of signaling network structure on the developmental
# increase in complexity."
# Start developing as a program by CBouyio on Aprill 2006 at UEA.
# Converted to a module (transsyslattice.py) at 08/09/2006 by CBouyio at UEA.


import random
import math

import transsys


class GaussianRNG :
  """
  Class of random objects out of a Gaussian distribution.

  By calling the random.gauss method, produces a "random" number out of a
  Gaussian distribution.
  @ivar rng: The random number object, generated by a call to the
  random.Random() method of Python's random module.
  @type rng: class 'random.Random'
  @ivar mu: The mean of the Gaussian distribution.
  @type mu: Float
  @ivar sigma: The standard deviation of the Gaussian distribution.
  @type sigma: Float
  """

  def __init__(self, rndseed, mu, sigma) :
    """
    Constrauctor of the class.

    Takes the random seed, the mean and the standard deviation.

    @param rndseed: The random seed of the random number generator.
    @type rndseed: Float
    @param mu: The mean of the Gaussian distribution.
    @type mu: Float
    @param sigma: The standard deviation of the Gaussian distribution.
    @type sigma: Float
    """
    self.rng   = random.Random(rndseed)
    self.mu    = mu
    self.sigma = sigma

  def random_value(self) :
    """
    Common method among all xxxxRNG classes.

    Returns a random number generated by the random.gauss Python function.
    @return: A random number from a Gaussian Distribution.
    @rtype: Float
    """
    return self.rng.gauss(self.mu, self.sigma)


class UniformRNG :
  """
  Class of random objects out of a Uniform distribution.
  By calling the random.uniform method, produces a "random" number out of a
  Uniform distribution.
  @ivar rng: The random number object, generated by a call to the
  random.Random() method of Python's random module.
  @type rng: Random Object
  """

  def __init__(self, rndseed) :
    """
    Constructor of the class.

    Takes only the random seed as parameter.
    @param rndseed: The random seed of the random number generator.
    @type rndseed: Float
    """
    self.rng = random.Random(rndseed)

  def random_value(self, a, b) :
    """
    Common method among all xxxxRNG classes.

    Returns a random number generated by the random.uniform Python function

    In this case the borders of the uniform distribution (a, b) are been
    provided from the function call parameters and not from the class's member
    variables.
    @param a: The left border of a Uniform distribution.
    @type a: Float
    @param b: The right border of a Uniform distribution.
    @type b: Float
    @return: A random number from a Gaussian Distribution.
    @rtype: Float
    """
    return self.rng.uniform(a, b)


class ConstantRNG :
  """
  Class of "random" objects with constant values!

  Basicaly for control reasons, the class actually generates constant numbers.
  @ivar constValue: A constant.
  @type constvalue: Float
  """
  # A "random" series of "constValue" will be produced!!!
  def __init__(self, constValue) :
    """
    Constructor of the class.

    The only instance variable is the constant value.
    """
    self.constValue = constValue

  def random_value(self) :
    """
    Returns the constant value.
    @return: The constant number.
    @rtype: Float
    """
    return self.constValue



class TranssysInstanceCoordinates(transsys.TranssysInstance) :
  """
  Extended TranssysInstance class.

  Provides the TranssysInstance objects with an additional functionality,
  theirs Cartesian coordinates on the lattice.
  @ivar transsys_program: A valid transsys program.
  @type transsys_program: class 'transsys.TranssysProgram'
  @ivar timestep: The number of timesteps that the simulator runs.
  @type timestep: Integer
  @ivar factor_concentration: A list with the factor concentrations.
  @type factor_concentration: List
  @ivar coordinates: A list denoting the coordinates on the lattice. (the
  length of the list equals the dimension of the lattice)
  @type coordinates: List
  """

  def __init__(self, tp, coords=None, timestep=None) :
    """
    Constructor of the class, overrides the TranssysInstance constructor.

    The coordinates of a transsys instance are stored in a list.
    """
#    transsys.TranssysInstance.__init__(self, tp, timestep)
    # No need to call the old constructor, fc_stddev and fc_entropy are
    # obsolate in this context. They have been implemented in the
    # transsys.statistics class.
    self.transsys_program = tp
    self.timestep = timestep
    self.factor_concentration = [0.0] * self.transsys_program.num_factors()
    self.coordinates = []

  def time_series(self, num_timesteps) :
    """
    Overrides the time_series method with a timeseries containing
    TranssysInstanceCoordinates instances instead.
    @param num_timesteps: The number of steps (runs) of the timeseries.
    @type num_timesteps: Integer
    @returns: A list of TranssysInstanceCoordinates instances (The length of
    list equals the number of timesteps).
    @rtype: List
    """
    timeseries = transsys.TranssysInstance.time_series(self, num_timesteps)
    # This list contains the timeseries with the coordinates.
    tsCoordinates = []
    for t in timeseries:
      a = TranssysInstanceCoordinates(t.transsys_program, t.timestep)
      a.transsys_program = t.transsys_program
      a.timestep = t.timestep
      a.factor_concentration = t.factor_concentration
     # No need to override these the new class doesn't need them, check
     # comments at the constructor of the class as well.
#      a.factor_concentration_stddev = t.factor_concentration_stddev
#      a.factor_concentration_entropy = t.factor_concentration_entropy
      a.coordinates = self.coordinates
      tsCoordinates.append(a)
    return tsCoordinates


class TranssysProgramLattice(transsys.TranssysInstanceCollection) :
  """
  The main class of the simulator, contains a lattice of transsys programs.

  Most of the methods for the simulation are contaned here. Is a subclass of
  the TranssysInstanceCollection super class of the transsys.py module.

  The lattice is forming a two dimensional toroidal form.
  @ivar name: The lattice name. Contains information of the type and the size.
  @type name: String
  @ivar lattice: A lattice containg a transsys programs in each cell.
  @type lattice: List of Lists
  @ivar transsysProgram: The transsys program.
  @type transsysProgram: class 'transsys.TranssysProgram'
  @ivar size: The size of the lattice.
  @type size: List
  """

  def __init__(self, tp, size) :
    """
    Constructor of the class.
    """
    # The name of the lattice.
    self.name = tp.name + '_on_a_toroidal_two_dimentional_' + str(size[0]) + 'x' + str(size[1]) + '_lattice'
    self.lattice = self.lattice_generator(tp, size)
    self.transsysProgram = tp
    self.size = size
    # Safeguard for acceptable diffusibility parameters.
    # Diffusibility ought to be a number between 0 and 1 [0, 1] to be accepted
    # by the update function.
    for factor in self.transsysProgram.factor_list :
      if 0.0 <= float(factor.diffusibility_expression.value) <= 1.0 :
        continue
      else :
        raise StandardError, 'Diffusibility expression error. Mind that latticeSimulator accepts diffusibility values within the range of 0.0 to 1.0 [0.0, 1.0]'


  def lattice_generator(self, tp, size) :
    """
    Method to generate the lattice of transsys instances.
    @param tp: The transsys program.
    @type tp: class 'transsys.TranssysProgram'
    @param size: The size of the lattice.
    @type size: List
    @return: A lattice populated with transsys program instances.
    @rtype: List of Lists
    """
    lattice = []
    for i in xrange(size[0]) :
      lattice.append([])
      for j in xrange(size[1]) :
        lattice[i].append(TranssysInstanceCoordinates(tp))
        lattice[i][j].coordinates.append(i + 1)
        lattice[i][j].coordinates.append(j + 1)
        lattice[i][j].timestep = 0
    return lattice


  def randomise_factor(self, fc, rng, rangeb, homogenise) :
    """
    Returns a factor concentration out of a uniform distribution.

    Checks for zero factor concentration and adapting the output.

    Also returns a constant value in case of the homogenisation switch is set.
    @param fc: The current factor concentration.
    @type fc: Float
    @param rng: The random number generator.
    @type rng: class 'xxxxxRNG'
    @param rangeb: The range notion for the Uniform distribution. (needs
    further explanation)
    @type rangeb: Float
    @param homogenise: The homogenisation switch.
    @type homogenise: Boolean
    @return: A randomize value for the factor concentration.
    @rtype: Float
  """
    if not homogenise :
      if fc != 0 :
        fc = rng.random_value((fc - fc * rangeb), (fc + fc * rangeb)) # This is a way to make sure that after a perturbation the factor concentration will remain positive.
      else :
        fc = rng.random_value(0, rangeb)
    else : # Populates the matrix with a "random" constant number
      fc = rng.random_value()
    return fc


  def randomise_lattice (self, rng, rangeb, homogenise) :
    """
    Method to randomise the factor concentrations of the lattice.

    Wraper of the previous method L{randomise_factor} for all the factors in a
    lattice.
    @param rng: The random number generator.
    @type rng: class 'xxxxxRNG'
    @param rangeb: The range notion for the Uniform distribution. (needs
    further explanation)
    @type rangeb: Float
    @param homogenise: The homogenisation switch.
    @type homogenise: Boolean
    """
#    k = self.transsysProgram.find_factor_index('Activator')
    for i in xrange(self.size[0]) :
      for j in xrange(self.size[1]) :
#        self.lattice[i][j].factor_concentration[k] = self.random_uniform_factor(self.lattice[i][j].factor_concentration[k], rng, rangeb)
        self.lattice[i][j].factor_concentration = map(lambda fc: self.randomise_factor(fc, rng, rangeb, homogenise), self.lattice[i][j].factor_concentration)


  def initialise_lattice_concentrations (self, homogenise, borderRange, rndSeed) :
    """
    Method to assign the initial factor concentrations on the lattice.
    """
    if not homogenise :
      rng = UniformRNG(rndSeed)
      self.randomise_lattice(rng, borderRange, homogenise)
    else :
      rng = ConstantRNG(borderRange * math.e)
      # Maybe is usefull (gives an homogenized matrix).
      self.randomise_lattice(rng, None, homogenise)


  def get_transsys_program(self) :
    """
    Overrides the get_transsys_program method of the superclass.

    Returns a transsys program.
    """
    return self.transsysProgram


  def transsys_instance_list(self) :
    """
    Overrides the transsys_instance_list method of the superclass.

    Returns a list of all the transsys program of a TranssysProgramCollection.
    """
    ti_list = []
    for i in xrange(len(self.lattice)) :
      for j in xrange(len(self.lattice[i])) :
        ti_list.append(self.lattice[i][j])
    return ti_list


  def write_table_header(self, f) :
    """
    Overrides the write_table_header method of the superclass.
    """
    f.write('# Table of coordinated (i,j) factor concentrations (Header)\n')
    f.write('timestep i j')
    for factor in self.transsysProgram.factor_list :
      f.write(' %s' % factor.name)
    f.write('\n')


  def write_table(self, f) :
    """
    Overrides the write_table method of the superclass.
    """
    for ti in self.transsys_instance_list() :
      if ti.timestep is None :
        f.write('NA')
      else :
        f.write('%i' % ti.timestep)
      f.write(' %i %i' % (ti.coordinates[0], ti.coordinates[1]))
      for fc in ti.factor_concentration :
        f.write(' %1.17e' % fc)
      f.write('\n')


  def update_factor_concentrations(self, currentState, size, rndseed) :
    """
    The update function of the simulator, calculates the instances of the next timestep.

    Careful usage of the factor concentrations...

    The method calculates the new factor concentrations after the calculation of all the diffusibility effects and returns the lattice with the "updated" factor concentrations.
    """
    # Calculate the updated factor concentrations.
    # All the instances in the lattice are interacting with the 4 neighbour
    # cells, the instances in the edges are interacting with the opposite
    # cells forming the toroidal.
    # First assign the dimensions of the lattice.
    x = size[0]
    y = size[1]
    latticeFactorConcentrations = []
    # Include some noise or not.
#    rng = GaussianRNG(rndseed, 0, 0.1)
    rng = ConstantRNG(0)
    for i in xrange(len(self.lattice)) :
      latticeFactorConcentrations.append([])
      for j in xrange(len(self.lattice[i])) :
        # The factor concentrations manipulation.
        factorConcentrationsUpdate = [0.0] * self.transsysProgram.num_factors()
        # Get the current lattice concentrations.
        factorConcentrationsOld = currentState[i][j].factor_concentration
        if len(factorConcentrationsOld) != len(factorConcentrationsUpdate) :
          raise StandardError, 'Factor number inconsistancy between the new and the old factor concentrations list.'
        for k in xrange(len(factorConcentrationsUpdate)) :
          # The main calculation of the new factor concentrations.
          # Introduce "per factor" diffusibility expression.
          p = self.transsysProgram.factor_list[k].diffusibility_expression.value
          # Check the toroidal transformations by using the modulo division!
          fc = (factorConcentrationsOld[k] + (rng.random_value() + (p * self.lattice[(x + i - 1) % x][j].factor_concentration[k])) + (rng.random_value() + (p * self.lattice[i][(y + j - 1) % y].factor_concentration[k])) + (rng.random_value() + (p * self.lattice[i][(j + 1) % y].factor_concentration[k])) + (rng.random_value() + (p * self.lattice[(i + 1) % x][j].factor_concentration[k]))) / (4 * p + 1)
          factorConcentrationsUpdate[k] = fc
        latticeFactorConcentrations[i].append(factorConcentrationsUpdate)
    return latticeFactorConcentrations


  def update_function(self, timesteps, rndseed) :
    """
    Calculates the new transsys instance.

    The method actually calls the transsys.time_series function for all the TranssysInstances on the lattice.
    """
    # Produce a copy of the current state of the simulator.
    currentState = self.lattice ## Possible Kkkkludge
    updateFactorConcentrations = self.update_factor_concentrations(currentState, self.size, rndseed)
    for i in xrange(len(self.lattice)) :
      for j in xrange(len(self.lattice[i])) :
        if len(self.lattice[i][j].factor_concentration) != len(updateFactorConcentrations[i][j]) :
          raise StandardError, 'update_function Error. Factor number inconsistancy'
        else :
          self.lattice[i][j].factor_concentration = updateFactorConcentrations[i][j]
          # The timeseries doesn't work with timestep 1, returns zero.
          # Thats why a timestep = 2 is used.
          # FIXME: The number of timesteps of each transsys instance on the
          # lattice timeseries is awlays one because of the usage of timestep
          # one all the time. Update it, to take the number of the current
          # timestep.
          self.lattice[i][j] = self.lattice[i][j].time_series(2)[1]
          self.lattice[i][j].timestep = (timesteps + 1)


  def signal_factor_concentration(self, xCenter, yCenter, signal, factor) :
    """
    Sets the factor concentration to signal.

    Method wich implements the introduce signal functionality.
    Sets the concentration of the factor to the given value.
    """
    for i in xCenter :
      for j in yCenter :
        if factor == None :
          self.lattice[i][j].factor_concentration = [signal for fc in self.lattice[i][j].factor_concentration]
        elif factor:
          k = self.transsysProgram.find_factor_index(factor)
          self.lattice[i][j].factor_concentration[k] = signal


  def introduce_signal(self, signalC, xFactor = None) :
    """
    Generic function for introduce signal.
    Has two functionalities depending on the call.

      1. Sets a specific factor concentration to the declared value in the begining (first timestep) of the simulator.
      If it is called according to the introduce signal factor functionality.

      2. Sets all the factor conscentrations the declared value on a declared timestep of the simulator.
      If it is called according to the introduce signal timestep functionallity.
    """
    # Check for the existance of the xFactor.
    if xFactor != None :
      if not xFactor in self.transsysProgram.factor_names() :
        raise StandardError, 'Factor %s is not a valid %s factor name. Please specify a valid factor name.' % (xFactor, self.transsysProgram.name)
    # Find the center of the lattice.
    # If abscissa is odd number.
    if divmod(self.size[0], 2)[1] == 1 :
      i = [divmod(self.size[0], 2)[0]]
    # If abscissa is even number.
    elif divmod(self.size[0], 2)[1] == 0 :
      q = divmod(self.size[0], 2)[0]
      i = range(q - 1, q + 1)
    # If ordinate is odd.
    if divmod(self.size[1], 2)[1] == 1 :
      j = [divmod(self.size[1], 2)[0]]
    # If ordinate is even.
    elif divmod(self.size[1], 2)[1] == 0 :
      q = divmod(self.size[1], 2)[0]
      j = range(q - 1, q + 1)
    # Call the signal_factor_concentration function.
    self.signal_factor_concentration(i, j, signalC, xFactor)
    #######
    # Alternatively, set a randomly choosen factor concentration.
#    k = random.randint(0, (self.ranssysProgram.num_factors() - 1))
#    self.lattice[i][j].factor_concentration[k] = random.gauss(signalC, 1)
#     # Randomise factor concentrations NOT USED anymore.
#     k = random.randint(0, (self.transsysProgram.num_factors() - 1))
      # Alternatively, set a randomly choosen factor concentration.
#      k = random.randint(0, (self.transsysProgram.num_factors() - 1))
#      self.lattice[i][j].factor_concentration[k] = random.gauss(signalC, 1)


def generate_pgm(fileObj, transsysLattice, maxCon) :
  """
  Produces the contents of the portable gray map (.pgm) image file.

  Suitable for view/animate/analyse with the Imagemagick graphics suite.
  """
  # Generate the apropriate .pgm image format.
  pgmMagic = 'P2' # .pgm file magic line, P2 for .pgm,
  pgmComment = '# ' + transsysLattice.name + ' .pgm image of ' + ' timestep'
  pgmSize = str(transsysLattice.size[0]) + ' ' + str(transsysLattice.size[1])
  pgmMaxval = 256
  pgmRaster = ''
  for i in xrange(transsysLattice.size[0]) :
    pgmRasterLine = ''
    for j in xrange(transsysLattice.size[1]) :
      #FIXME: this should simplified more and produce one image file for EACH
      # factor at more complicated transsys programs.
      p = pgmMaxval - int(round(sum(transsysLattice.lattice[i][j].factor_concentration) * pgmMaxval / maxCon))
      # Safeguard for negative color values.
      if p < 0 :
        raise StandardError, 'Netpbm raster image get a negative value. Please select a biger number for the maximum factor concentration (option -m)'
      pgmRasterLine = pgmRasterLine + str(p) + ' '
    pgmRaster = pgmRaster + pgmRasterLine + '\n'
  pgmText = pgmMagic + '\n' + pgmComment + '\n' + pgmSize + '\n' + str(pgmMaxval)  + '\n' + pgmRaster
  fileObj.write(pgmText)
  fileObj.close()

