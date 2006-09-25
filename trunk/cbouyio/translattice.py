#!/usr/bin/env python

"""
Here will be the translattice module description....

@author: Costas Bouyioukos
@organization: University of East Anglia
@since: Aprill 2006 (converted to a module on 08/09/2006)
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
import copy

import transsys


class GaussianRNG :
  """
  Class of random objects out of a Gaussian distribution.

  By calling the random.gauss method, produces a "random" number out of a
  Gaussian distribution.

  @ivar rng: The random number object, generated by a call to the
  random.Random() method of Python's random module.
  @type rng: C{class 'random.Random'}
  @ivar mu: The mean of the Gaussian distribution.
  @type mu: C{float} or C{int}
  @ivar sigma: The standard deviation of the Gaussian distribution.
  @type sigma: C{float} or C{int}
  """

  def __init__(self, rndseed, mu, sigma) :
    """
    Constrauctor of the class.

    Takes the random seed, the mean and the standard deviation.

    @param rndseed: The random seed of the random number generator.
    @type rndseed: C{int}
    @param mu: The mean of the Gaussian distribution.
    @type mu: C{float} or C{int}
    @param sigma: The standard deviation of the Gaussian distribution.
    @type sigma: C{float} or C{int}
    """
    if not isinstance(rndseed, int) :
      raise StandardError, 'The random number generator accepts only integers.'
    self.rng   = random.Random(rndseed)
    self.mu    = mu
    self.sigma = sigma

  def random_value(self) :
    """
    Common method among all xxxxRNG classes.

    Returns a random number generated by the random.gauss Python function.

    @return: A random number from a Gaussian Distribution.
    @rtype: C{float}
    """
    return self.rng.gauss(self.mu, self.sigma)


class UniformRNG :
  """
  Class of random objects out of a Uniform distribution.
  By calling the random.uniform method, produces a "random" number out of a
  Uniform distribution.

  @ivar rng: The random number object, generated by a call to the
  random.Random() method of Python's random module.
  @type rng: C{class 'random.Random'}
  """

  def __init__(self, rndseed) :
    """
    Constructor of the class.

    Takes only the random seed as parameter.

    @param rndseed: The random seed of the random number generator.
    @type rndseed: C{int}
    """
    if not isinstance(rndseed, int) :
      raise StandardError, 'The random number generator accepts only integers.'
    self.rng = random.Random(rndseed)

  def random_value(self, a, b) :
    """
    Common method among all xxxxRNG classes.

    Returns a random number generated by the random.uniform Python function

    In this case the borders of the uniform distribution (a, b) are been
    provided from the function call parameters and not from the class's member
    variables.

    @param a: The left border of a Uniform distribution.
    @type a: C{float} or C{int}
    @param b: The right border of a Uniform distribution.
    @type b: C{float} or C{int}
    @return: A random number from a Uniform Distribution.
    @rtype: C{float}
    """
    return self.rng.uniform(a, b)


class ConstantRNG :
  """
  Class of "random" objects with constant values.

  Basicaly for control reasons, the class actually generates constant numbers.

  @ivar constValue: A constant.
  @type constValue: C{float} or C{int}
  """
  # A "random" series of "constValue" will be produced!!!
  def __init__(self, constValue) :
    """
    Constructor of the class.

    The only instance variable is the constant value.

    @param constValue: A arbitarily chosen number.
    @type constValue: C{float} or C{int}
    """
    self.constValue = constValue

  def random_value(self) :
    """
    Returns the constant value.

    @return: The constant number.
    @rtype: C{float} or C{int}
    """
    return self.constValue



class TranssysInstanceCoordinated(transsys.TranssysInstance) :
  """
  Extended TranssysInstance class.

  Provides the TranssysInstance objects with an additional functionality,
  theirs Cartesian coordinates on the lattice.

  @ivar transsys_program: A valid transsys program.
  @type transsys_program: C{class 'transsys.TranssysProgram'}
  @ivar timestep: The number of timesteps that the simulator runs.
  @type timestep: C{int}
  @ivar factor_concentration: A list with the factor concentrations.
  @type factor_concentration: C{list} of floating point numbers
  @ivar coordinates: A list denoting the coordinates on the lattice. (the
  length of the list equals the dimension of the lattice)
  @type coordinates: C{list} of integers
  """

  def __init__(self, tp, coords=None, timestep=None) :
    """
    Constructor of the class, overrides the TranssysInstance constructor.

    The coordinates of a transsys instance are stored in a list.

    @param tp: A valid transsys program.
    @type tp: C{class 'transsys.TranssysProgram'}
    @param coords: A list denoting the coordinates on the lattice. (the
    length of the list equals the dimension of the lattice)
    @type coords: C{list} of integers
    @param timestep: The number of timesteps that the simulator runs.
    @type timestep: C{int}
    """
#    transsys.TranssysInstance.__init__(self, tp, timestep)
    # No need to call the old constructor, fc_stddev and fc_entropy are
    # obsolate in this context. They have been implemented in the
    # transsys.statistics class.
    self.transsys_program = tp
    self.timestep = timestep
    self.factor_concentration = [0.0] * self.transsys_program.num_factors()
    self.coordinates = []

  def time_series(self, num_timesteps, sampling_period = 1, lsys_lines = None, lsys_symbol = None ) :
    """
    Overrides the time_series method with a timeseries containing
    TranssysInstanceCoordinated instances instead.

    @param num_timesteps: The number of steps (runs) of the timeseries.
    @type num_timesteps: C{int}
    @returns: A list of TranssysInstanceCoordinated instances (The length of
    list equals the number of timesteps).
    @rtype: C{list} of C{TranssysInstanceCoordinated} objects
    """
    timeseries = transsys.TranssysInstance.time_series(self, num_timesteps)
    # This list contains the timeseries with the coordinates.
    tsCoordinated = []
    for t in timeseries:
      a = TranssysInstanceCoordinated(t.transsys_program, t.timestep)
      a.transsys_program = t.transsys_program
      a.timestep = t.timestep
      a.factor_concentration = t.factor_concentration
     # No need to override these the new class doesn't need them, check
     # comments at the constructor of the class as well.
#      a.factor_concentration_stddev = t.factor_concentration_stddev
#      a.factor_concentration_entropy = t.factor_concentration_entropy
      a.coordinates = self.coordinates
      tsCoordinated.append(a)
    return tsCoordinated



class TranssysInstanceLattice(transsys.TranssysInstanceCollection) :
  """
  The main class of the simulator, contains a lattice of transsys programs.

  Most of the methods for the simulation are contained here. Is a subclass of
  the TranssysInstanceCollection super class of the transsys module.

  The lattice is formated in a two dimensional toroidal form.

  @ivar name: The lattice name. Contains information of the type and the size
  of the structure.
  @type name: C{str}
  @ivar lattice: A lattice containg a transsys programs in each cell.
  @type lattice: C{list} of C{list} of C{transsys.TranssysProgram} objects
  @ivar transsysProgram: The transsys program.
  @type transsysProgram: C{class 'transsys.TranssysProgram'}
  @ivar size: The size of the lattice.
  @type size: C{list}
  @Ivar timestep: The timestep.
  @type timestep: C{int}
  """

  def __init__(self, tp, size, timestep=None) :
    """
    Constructor of the class.

    @param tp: The transsys program.
    @type tp: C{class 'transsys.TranssysProgram'}
    @param size: The size of the lattice.
    @type size: C{list}
    @param timestep: The timestep.
    @type timestep: C{int}
    """
    # The name of the lattice.
    self.name = tp.name + '_on_a_toroidal_two_dimentional_' + str(size[0]) + 'x' + str(size[1]) + '_lattice'
    self.lattice = self.lattice_generator(tp, size)
    self.transsysProgram = tp
    self.size = size
    self.timestep = timestep
    # Safeguard for acceptable diffusibility parameters.
    # Diffusibility ought to be a number between 0 and 1 [0, 1] to be accepted
    # by the update function.
    for factor in self.transsysProgram.factor_list :
      if 0.0 <= float(factor.diffusibility_expression.value) <= 1.0 :
        continue
      else :
        raise StandardError, 'Diffusibility expression error. The accepted diffusibility values resides within the range of 0.0 to 1.0 [0.0, 1.0]'


  def lattice_generator(self, tp, size) :
    """
    Method to generate the lattice of transsys instances.
    Populates the lattice with one transsys program instance to each cell.

    @param tp: The transsys program.
    @type tp: C{class 'transsys.TranssysProgram'}
    @param size: The size of the lattice.
    @type size: C{list} of integers
    @return: A lattice populated with transsys program instances.
    @rtype: C{list} of C{list} of C{transsys.TranssysProgram} objects
    """
    lattice = []
    for i in xrange(size[0]) :
      lattice.append([])
      for j in xrange(size[1]) :
        lattice[i].append(TranssysInstanceCoordinated(tp))
        lattice[i][j].coordinates.append(i + 1)
        lattice[i][j].coordinates.append(j + 1)
        lattice[i][j].timestep = 0
    return lattice


  def randomise_factor(self, fc, rng, rangeb) :
    """
    Returns a factor concentration out of a uniform distribution.

    Checks for zero factor concentration and adapting the output.

    Also returns a constant value in case of the homogenisation switch is set.

    @param fc: The current factor concentration.
    @type fc: C{float}
    @param rng: The random number generator.
    @type rng: C{class 'xxxxxRNG'}
    @param rangeb: The range notion for the Uniform distribution. (needs
    further explanation)
    @type rangeb: C{float}
    @return: A randomize value for the factor concentration.
    @rtype: C{float}
  """
    if fc != 0 :
      fc = rng.random_value((fc - fc * rangeb), (fc + fc * rangeb)) # This is a way to make sure that after a perturbation the factor concentration will remain positive.
    else :
      fc = rng.random_value(0, rangeb)
    return fc


  def randomise_lattice (self, rng, rangeb) :
    """
    Method to randomise the factor concentrations of the lattice.

    Wraper of the previous method L{randomise_factor} for all the factors in a
    lattice.

    @param rng: The random number generator.
    @type rng: C{class 'xxxxxRNG'}
    @param rangeb: The range notion for the Uniform distribution. (needs
    further explanation)
    @type rangeb: C{float}
    @rtype: C{None}
    """
    for i in xrange(self.size[0]) :
      for j in xrange(self.size[1]) :
        self.lattice[i][j].factor_concentration = map(lambda fc: self.randomise_factor(fc, rng, rangeb), self.lattice[i][j].factor_concentration)


  def initialise_lattice_concentrations (self, borderRange, rndSeed) :
    """
    Method to assign the initial factor concentrations on the lattice.

    Initialising the random seed and call the L{randomise_lattice}.

    @param borderRange: The upper limit for the Uniform distribution border.
    (Uniform returns random real number from the [0, borderRange) interval)
    @type borderRange: C{float}
    @param rndSeed: The seed for the random number genetaror.
    @type rndSeed: C{int}
    @rtype: C{None}
    @todo: This method might be merged with the L{randomise_lattice}
    """
    rng = UniformRNG(rndSeed)
    self.randomise_lattice(rng, borderRange)


  def get_transsys_program(self) :
    """
    Overrides the get_transsys_program method of the superclass.

    Returns a transsys program.

    @return: A C{transsys.TranssysProgram} object
    @rtype: C{class 'transsys.TranssysProgram'}
    """
    return self.transsysProgram


  def transsys_instance_list(self) :
    """
    Overrides the transsys_instance_list method of the superclass.

    Returns a list of all the transsys programs from a
    TranssysInstanceCollection.

    @return: A list of the transsys programs.
    @rtype: C{list} of C{transsys.TranssysProgram} objects
    """
    ti_list = []
    for i in xrange(len(self.lattice)) :
      for j in xrange(len(self.lattice[i])) :
        ti_list.append(self.lattice[i][j])
    return ti_list


  def write_table_header(self, f) :
    """
    Overrides the write_table_header method of the superclass.

    Writes to the output file the header of the factor's table.
    (compatible to R)

    @param f: An open for writing C{file} object.
    @type f: C{file}
    @rtype: C{None}
    """
    f.write('# Table of coordinated (i,j) factor concentrations (Header)\n')
    f.write('timestep x y')
    for factor in self.transsysProgram.factor_list :
      f.write(' %s' % factor.name)
    f.write('\n')


  def write_table(self, f) :
    """
    Overrides the write_table method of the superclass.

    Writes at the output file the factor table (factor concentrations) for each
    factor of each transsys program in each timestep of the simulator.

    @param f: An open for writing C{file} object.
    @type f: C{file}
    @rtype: C{None}
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
    The update function of the simulator, calculates the instances of the next
    timestep.

    Careful usage of the factor concentrations...

    The method calculates the new factor concentrations after the calculation
    of all the diffusibility effects and returns the lattice with the
    "updated" factor concentrations.

    @param currentState: A transsys program lattice.
    @type currentState: C{class 'TranssysInstanceLattice'}
    @param size: The size of the lattice.
    @type size: C{list} of C{int}
    @param rndseed: The random seed for the random number generator. (It is
    used only in cases that we want to incorporate some noise in the)
    @type rndseed: C{int}
    @return: A transsys program lattice containing the updated factor
    concentrations, after the diffusibility effects have been calculated.
    @rtype: C{class 'TranssysInstanceLattice'}
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
    if rndseed is None :
      rng = ConstantRNG(0)
    else :
      rng = GaussianRNG(rndseed, 0, 0.1)
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


  def update_function(self, timesteps, rndseed = None) :
    """
    Calculates the new transsys instance.

    The method first updates the factor concentrations by calling the
    L{update_factor_concentrations} function, and then calculates the new
    instance of the lattice (next timestep) by calling the
    L{TranssysInstanceCoordinated.time_series} for all the TranssysInstances
    in each cell of the lattice (it's a wrapper).

    @param timesteps: The number of timestep of the simulator proccedure.
    @type timesteps: C{int}
    @param rndseed: The random seed for the random number generator. (It is
    used only in cases that we want to incorporate some noise in the)
    @type rndseed: C{int}
    @rtype: C{None}
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
          self.lattice[i][j] = self.lattice[i][j].time_series(2)[1]
          self.lattice[i][j].timestep = timesteps


  def signal_factor_concentration(self, xCenter, yCenter, signal, factor = None) :
    """
    Set the factor's concentration to signal.

    Method wich implements the introduce signal functionality.
    Sets the concentration of the factor to the given value.

    @param xCenter: A list containing the coordiantes of the x-axis lattice
    center.
    @type xCenter: C{list}
    @param yCenter: A list containing the coordiantes of the y-axis lattice
    center.
    @type yCenter: C{list}
    @param signal: The signal factor concentration.
    @type signal: C{float}
    @param factor: The signal factor name.
    @type factor: C{str}
    @rtype: C{None}
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

      1. Sets a specific factor concentration to the declared value in the
      begining (first timestep) of the simulator.
      If it is called according to the introduce signal factor functionality.

      2. Sets all the factor conscentrations the declared value on a declared
      timestep of the simulator.
      If it is called according to the introduce signal timestep
      functionallity.

    The assignament of the signal actually takes place on the
    L{signal_factor_concentration}

    @param signalC: The signal factor concentration.
    @type signalC: C{float}
    @param xFactor: The signal factor name.
    @type xFactor: C{str}
    @rtype: C{None}
    """
    # Check for the existance of the xFactor.
    if xFactor != None :
      if not xFactor in self.transsysProgram.factor_names() :
        raise StandardError, 'Factor %s is not belonging to  %s transsys program. Please specify a valid factor name.' % (xFactor, self.transsysProgram.name)
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



class TranssysLatticeTimeseries :
  """
  Keep the whole timeseries of the simulator.

  Keep record of all the important things during the simulation process. The
  position of the transys instance on the lattice the timestep and the
  factor(s) concentration(s).

  Also calulates the maximum value of factor concentration that is observed
  during the simulation proccess.
  @ivar transsysLatticeName: The transsys lattice name. Contains information
  of the type and the size of the structure.
  @type transsysLatticeName: C{str}
  @ivar timesteps: The number of timesteps of the simulator.
  @type timesteps: C{int}
  @ivar latticeTimeseries: A list containing all the transsys lattice
  instances of each timestep of the simulator.
  @type latticeTimeseries: C{list} of C{TranssysInstanceLattice}
  objects
  @ivar maxFactorConcentration: The maximum factor concentration observed in
  the simulation.
  @type : C{float}
  """


  def __init__(self, transsysLattice, timesteps = None) :
    """
    Constructor of the class.
    @param transsysLattice: A transsys instance lattice object.
    L{TranssysInstanceLattice}
    @type transsysLattice: C{class 'TranssysInstanceLattice'}
    @param timesteps: The number of timesteps for the simulator.
    @type timesteps: C{int}
    """
#    self.factorTable = factor_table(fileObject)
    self.transsysLatticeName = transsysLattice.name
    self.timesteps = timesteps
    self.latticeTimeseries = self.transsys_lattice_timeseries(transsysLattice, timesteps)
    self.maxFactorConcentration = self.max_factor_concentration()


  def transsys_lattice_timeseries(self, transsysLattice, timestep = None) :
    """
    Keep all the L{TranssysInstanceLattice} instances in a list.

    @param transsysLattice: A transsys lattice.
    @type transsysLattice: C{class 'TranssysInstanceLattice'}
    @param timestep: The timesteps of the simulation.
    @type timestep: C{int}
    @return: A list with all the C{TranssysInstanceLattice} instances of the
    simulator.
    @rtype: C{list} of C{TranssysInstanceLattice} objects
    """
    # Put the first transsys lattice instance in the timeseries list.
    latticeTimeseries = [transsysLattice]
    # make a copy of it
    newLattice = copy.deepcopy(transsysLattice)
    if timestep :
      for i in xrange(timestep) :
        newLattice.update_function((i + 1))
        newLattice.timestep = (i + 1)
        latticeTimeseries.append(newLattice)
        # Copy itself for the next timestep.
        newLattice = copy.deepcopy(newLattice)
    return latticeTimeseries


  def write_factor_table(self, f) :
    """
    Writes the whole factor table in a file.

    Calls the L{TranssysInstanceLattice.write_table} for the whole transsys
    lattice timeseries.

    @param f: An open file object ready for writing.
    @type f: C{file}
    @rtype: C{None}
    """
    self.latticeTimeseries[0].write_table_header(f)
    for til in self.latticeTimeseries :
      til.write_table(f)


  def max_factor_concentration(self) :
    """
    Calculates the maximum factor concentration that is observed during the
    whole simulation proccess.

    @return: The maximal factor concentration of the whole simulation proccess.
    @rtype: C{float}
    """
    maxFC = 0.0
    for til in self.latticeTimeseries :
      for ti in til.transsys_instance_list() :
        for fc in ti.factor_concentration :
          if fc > maxFC :
            maxFC = fc
    return maxFC


  def factor_table(self, f) :
    """
    Returns the factor table out of a file that contain it in a string format.
    @param f: An open file object ready for reading.
    @type f: C{file}
    """
    self.write_factor_table(f)
    factorTable = f.read()
    return factorTable



def generate_pgm(fileObj, transsysLattice, maxCon) :
  """
  Produce the contents of the portable gray map (.pgm) image file.

  Write the to a .pgm file and then close it.

  Suitable for view/animate/analyse with the Imagemagick graphics suite.

  @param fileObj: An open raedy for writting file object.
  @type fileObj: C{file}
  @param transsysLattice: The whole transsys lattice.
  @type transsysLattice: C{class 'TranssysLattice'}
  @param maxCon: The maximum nuber of factor concentration that can be drawn
  by the imaging proccedure. (constant provided by the user)
  @type maxCon: C{int}
  @rtype: C{None}
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

