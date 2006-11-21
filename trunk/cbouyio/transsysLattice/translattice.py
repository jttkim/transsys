#!/usr/bin/env python

# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date$:  Date of last commit

"""
A module wich contains the essential tools for building structured 2D
transsys systems.

A 2D toroidal shaped grid (lattice) is the main structure that has been
implemented until now.

@author: Costas Bouyioukos
@organization: University of East Anglia
@since: Aprill 2006 (converted to a module on 08/09/2006)
@copyright: The program is coming as it is. You have the right to redistribute,
transform and change the source code presuming the apropriate reference and
the lisence is kept free.
@license: GNU GPL2 or newer.
@contact: U{Costas Bouyioukos<mailto:konsb@cmp.uea.ac.uk>}
@version: $Id$
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
import copy
import cPickle

import transsys


class GaussianRNG(object) :
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


class UniformRNG(object) :
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


  def random_value(self, lower, upper) :
    """
    Common method among all xxxxRNG classes.

    Returns a random number generated by the random.uniform Python function,
    the number belongs to the [lower, upper) interval.

    In this case the borders of the uniform distribution (lower, upper) are
    provided from the function call parameters and not from the class's member
    variables.

    @param lower: The left border of a Uniform distribution.
    @type lower: C{float} or C{int}
    @param upper: The right border of a Uniform distribution.
    @type upper: C{float} or C{int}
    @return: A random number from a Uniform Distribution.
    @rtype: C{float}
    """
    return self.rng.uniform(lower, upper)



class ConstantRNG(object) :
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
  @ivar factor_concentration: A list containing the factor concentrations.
  @type factor_concentration: C{list} of C{float}s
  @ivar coordinates: A list denoting the coordinates on the lattice. (the
  length of the list equals the dimension of the lattice)
  @type coordinates: C{list} of C{int}s
  """

  def __init__(self, tp, coords=None, timestep=None) :
    """
    Constructor of the class, overrides the TranssysInstance constructor.

    The coordinates of a transsys instance are stored in a list.

    @param tp: A valid transsys program.
    @type tp: C{class 'transsys.TranssysProgram'}
    @param coords: A list denoting the coordinates on the lattice. (the
    length of the list equals the dimension of the lattice)
    @type coords: C{list} of C{int}s
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
    self.coordinates = coords

  def time_series(self, num_timesteps, sampling_period=1, lsys_lines=None, lsys_symbol=None) :
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
    # This list will get the timeseries with the coordinates.
    tsCoordinated = []
    for ti in timeseries:
      coordTI = TranssysInstanceCoordinated(ti.transsys_program, ti.timestep)
      coordTI.transsys_program = ti.transsys_program
      coordTI.timestep = ti.timestep
      coordTI.factor_concentration = ti.factor_concentration
     # No need to override these the new class doesn't need them, check
     # comments at the constructor of the class as well.
#      a.factor_concentration_stddev = t.factor_concentration_stddev
#      a.factor_concentration_entropy = t.factor_concentration_entropy
      coordTI.coordinates = self.coordinates
      tsCoordinated.append(coordTI)
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
  @type size: C{tuple}
  @Ivar timestep: The timestep.
  @type timestep: C{int}
  """

  def __init__(self, tp, size, timestep=None) :
    """
    Constructor of the class.

    @param tp: The transsys program.
    @type tp: C{class 'transsys.TranssysProgram'}
    @param size: The size of the lattice.
    @type size: C{tuple}
    @param timestep: The timestep.
    @type timestep: C{int}
    """
    # Call the base class constructor.
    transsys.TranssysInstanceCollection.__init__(self)
    # The name of the lattice.
    self.name = tp.name + '_on_a_toroidal_two_dimentional_' + str(size[0]) + 'x' + str(size[1]) + '_lattice'
    self.size = size
    self.transsysProgram = tp
    self.timestep = timestep
    self.lattice = self.lattice_generator()
    # Safeguard for acceptable diffusibility parameters.
    # Diffusibility ought to be a number between 0 and 1 [0, 1] to be accepted
    # by the update function.
    for factor in self.transsysProgram.factor_list :
      if 0.0 <= float(factor.diffusibility_expression.value) <= 1.0 :
        continue
      else :
        raise StandardError, 'Diffusibility expression error. The accepted diffusibility values resides within the range of 0.0 to 1.0 [0.0, 1.0]'


  def lattice_generator(self) :
    """
    Method to generate the lattice of transsys instances.
    Populates the lattice with one transsys program instance to each cell.

    @return: A lattice populated with transsys program instances.
    @rtype: C{list} of C{list} of C{transsys.TranssysProgram} objects
    """
    lattice = []
    for i in xrange(self.size[0]) :
      lattice.append([])
      for j in xrange(self.size[1]) :
        coordinates = [i+1, j+1]
#        if self.timestep :
#          timestep = self.timestep
#        else :
#          timestep = 0
        lattice[i].append(TranssysInstanceCoordinated(self.transsysProgram, coordinates))
    return lattice


  def randomise_factor(self, fc, rng, lower, upper) :
    """
    Returns a factor concentration out of a uniform distribution.

    Checks for zero factor concentration and adapting the output.

    Also returns a constant value in case of the homogenisation switch is set.

    @param fc: The current factor concentration.
    @type fc: C{float}
    @param rng: The random number generator.
    @type rng: C{class 'xxxxxRNG'}
    @param lower: The lower limit for the Uniform distribution border.
    @type lower: C{float}
    @param upper: The upper limit for the Uniform distribution border.
    (Uniform returns random real number from the [a, b) interval)
    @type upper: C{float}
    @return: A randomize value for the factor concentration.
    @rtype: C{float}
    """
    rndValue = rng.random_value(lower, upper)
    if fc == 0 :
      fc = rndValue
    else :
      fc = fc + rndValue # This is a way to make sure that after a perturbation
                         # the factor concentration will remain positive.
    return fc


  def randomise_lattice(self, rng, lower, upper) :
    """
    Method to randomise the factor concentrations of the lattice.

    Wraper of the previous method L{randomise_factor} for all the factors in a
    lattice.

    @param rng: The random number generator.
    @type rng: C{class 'xxxxxRNG'}
    @param lower: The lower limit for the Uniform distribution border.
    @type lower: C{float}
    @param upper: The upper limit for the Uniform distribution border.
    (Uniform returns random real number from the [a, b) interval)
    @type upper: C{float}
    @rtype: C{None}
    """
    for i in xrange(self.size[0]) :
      for j in xrange(self.size[1]) :
        self.lattice[i][j].factor_concentration = map(lambda fc: self.randomise_factor(fc, rng, lower, upper), self.lattice[i][j].factor_concentration)


  def initialise_lattice_concentrations(self, lower, upper, rndSeed) :
    """
    Method to assign the initial factor concentrations on the lattice.

    Initialising the random seed and call the L{randomise_lattice}.

    @param lower: The lower limit for the Uniform distribution border.
    @type lower: C{float}
    @param upper: The upper limit for the Uniform distribution border.
    (Uniform returns random real number from the [a, b) interval)
    @type upper: C{float}
    @param rndSeed: The seed for the random number genetaror.
    @type rndSeed: C{int}
    @rtype: C{None}
    @todo: This method might be merged with the L{randomise_lattice}
    """
    rng = UniformRNG(rndSeed)
    self.randomise_lattice(rng, lower, upper)


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
    tiList = []
    for i in xrange(len(self.lattice)) :
      for j in xrange(len(self.lattice[i])) :
        tiList.append(self.lattice[i][j])
    return tiList


  def write_table_header(self, f, rndseed=None) :
    """
    Overrides the write_table_header method of the superclass.

    Writes to the output file the header of the factor's table.
    (compatible to R)

    @param f: An open for writing C{file} object.
    @type f: C{file}
    @param rndseed: The random seed of the simulator.
    @type rndseed: C{int}
    @rtype: C{None}
    @precondition: An open, ready for writting file object.
    """
    f.write('# Table of coordinated (x, y) factor concentrations (Header file)\n')
    f.write('# random seed: %i \n' % rndseed)
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
    @precondition: An open, ready for writting file object.
    """
    for ti in self.transsys_instance_list() :
      if self.timestep is None :
        f.write('NA')
      else :
        f.write('%i' % self.timestep)
      f.write(' %i %i' % (ti.coordinates[0], ti.coordinates[1]))
      for fc in ti.factor_concentration :
        f.write(' %1.17e' % fc)
      f.write('\n')


  def update_factor_concentrations(self, currentState, noiseSeed) :
    """
    The update function of the simulator, calculates the factor concentrations
    of the next timestep.

    Careful usage of factor concentrations diffusibilities and noise...!

    The method calculates the new factor concentrations after the calculation
    of all the diffusibility effects and returns the lattice list with the
    "updated" factor concentrations.

    @param currentState: A transsys program lattice.
    @type currentState: C{class 'TranssysInstanceLattice'}
    @param noiseSeed: The random seed for the random number generator. (It is
    used only in cases that we want to incorporate some noise in the)
    @type noiseSeed: C{int}
    @return: A list of lists of lists of factor concentrations, after the
    diffusibility effects have been calculated.
    @rtype: C{list} of C{list}s of C{list}s
    """
    # Calculate the updated factor concentrations.
    # All the instances in the lattice are interacting with the 4 neighbour
    # cells, the instances in the edges are interacting with the opposite
    # cells forming the toroidal.
    # First assign the dimensions of the lattice.
    x = self.size[0]
    y = self.size[1]
    latticeFactorConcentrations = []
    # Introduse some noise or not.
    if noiseSeed :
      rng = GaussianRNG(rndseed, 0, 0.01)
    else :
      rng = ConstantRNG(0)
    # Begin the calculations.
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
          diffRate = self.transsysProgram.factor_list[k].diffusibility_expression.value
          # Check the toroidal transformations by using the modulo division!
          fc = (factorConcentrationsOld[k] + (rng.random_value() + (diffRate * self.lattice[(x + i - 1) % x][j].factor_concentration[k])) + (rng.random_value() + (diffRate * self.lattice[i][(y + j - 1) % y].factor_concentration[k])) + (rng.random_value() + (diffRate * self.lattice[i][(j + 1) % y].factor_concentration[k])) + (rng.random_value() + (diffRate * self.lattice[(i + 1) % x][j].factor_concentration[k]))) / (4 * diffRate + 1)
          factorConcentrationsUpdate[k] = fc
        latticeFactorConcentrations[i].append(factorConcentrationsUpdate)
    return latticeFactorConcentrations


  def update_function(self, timesteps, noiseSeed=None) :
    """
    Calculates the new transsys instance.

    The method first updates the factor concentrations by calling the
    L{update_factor_concentrations} function, and then calculates the new
    instance of the lattice (next timestep) by calling the
    L{TranssysInstanceCoordinated.time_series} for all the TranssysInstances
    in each cell of the lattice (it's a wrapper).

    @param timesteps: The number of timestep that the simulator has reached.
    @type timesteps: C{int}
    @param noiseSeed: The random seed for the random number generator. (It is
    used only in cases that we want to incorporate some noise in the)
    @type noiseSeed: C{int}
    @rtype: C{None}
    """
    # Produce a copy of the current state of the simulator.
#    currentState = self.lattice ## Possible Kkkkludge
    updateFactorConcentrations = self.update_factor_concentrations(self.lattice, noiseSeed)
    for i in xrange(len(self.lattice)) :
      for j in xrange(len(self.lattice[i])) :
        if len(self.lattice[i][j].factor_concentration) != len(updateFactorConcentrations[i][j]) :
          raise StandardError, 'update_function Error. Factor number inconsistancy'
        else :
          self.lattice[i][j].factor_concentration = updateFactorConcentrations[i][j]
          # The timeseries doesn't work with timestep 1, returns zero.
          # Thats why a timestep = 2 is used.
          self.lattice[i][j] = self.lattice[i][j].time_series(2)[1]
          # Then set the timestep.
          self.lattice[i][j].timestep = timesteps


  def timestep_factor_concentration(self, signalC) :
    """
    Set all the factor concentrations on the lattice to signal concentration.

    @param signalC: The factor concentration of the signal.
    @type signalC: C{float}
    @rtype: C{None}
    """
    for i in xrange(self.size[0]) :
      for j in xrange(self.size[1]) :
        self.lattice[i][j].factor_concentration = [signalC for fc in self.lattice[i][j].factor_concentration]


  def signal_factor_concentration(self, sFactor, sConcentration, i=0, j=0) :
    """
    Function to set the factor concentration of a specified factor to the
    signal value.

    Function changes the factor concentration of the defiened factor before the
    begining of the simulator. Alters the factor concentration ONLY at the
    "first" (cell of the simulator (coordinates (1, 1)) wich acts as a signal.
    Thats why function's name includes signal.

    @param sFactor: The signal factor name.
    @type sFactor: C{str}
    @param sConcentration: The signal factor concentration.
    @type sConcentration: C{float}
    @param i: The abscissa of the coordinates.
    @type i: C{int}
    @param j: The ordinate of the coordinates.
    @type j: C{int}
    @rtype: C{None}
    @precondition: sFactor should be a valid factor of the corresponding
    transsys program, sConcentration should be a positive floating point
    number.
    """
    # Check for the existance of the signal factor.
    if not sFactor in self.transsysProgram.factor_names() :
      raise StandardError, 'Factor %s is not belonging to the %s transsys program. Please specify a valid factor name.' % (sFactor, self.transsysProgram.name)
    # Change the factor concentration.
    k = self.transsysProgram.find_factor_index(sFactor)
    self.lattice[i][j].factor_concentration[k] = sConcentration



class TranssysLatticeTimeseries(object) :
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
  @type maxFactorConcentration: C{float}
  """


  def __init__(self, transsysLattice, timesteps, timestepSignal=None) :
    """
    Constructor of the class.
    @param transsysLattice: A transsys instance lattice object.
    L{TranssysInstanceLattice}
    @type transsysLattice: C{class 'TranssysInstanceLattice'}
    @param timesteps: The number of timesteps for the simulator.
    @type timesteps: C{int}
    @param timestepSignal: The signal to alter factor concentrations on a
    specific timestep
    type timestepSignal: C{int}:C{float}
    """
#    self.factorTable = factor_table(fileObject)
    self.transsysLatticeName = transsysLattice.name
    self.timesteps = timesteps
    self.latticeTimeseries = self.transsys_lattice_timeseries(transsysLattice, timestepSignal)
    self.maxFactorConcentration = self.max_factor_concentration()


  def transsys_lattice_timeseries(self, transsysLattice, timestepSignal=None) :
    """
    Keep all the L{TranssysInstanceLattice} instances in a list.

    Returns a timeseries of C{TranssysInstanceLattice} objects for all the
    simulation in the form of a list.

    @param transsysLattice: A transsys lattice.
    @type transsysLattice: C{class 'TranssysInstanceLattice'}
    @param timestepSignal: The timestep signal parameter. Defines the factor
    concentration and the timestep that it will be introduced (to all factors
    untill now)
    @type timestepSignal: C{float}:C{int}
    @return: A list with all the C{TranssysInstanceLattice} instances of the
    simulator.
    @rtype: C{list} of C{TranssysInstanceLattice} objects
    """
    latticeTimeseries = []
    # Iterate for timestep + 1 to include the zero timestep. 
    for i in xrange(self.timesteps + 1) :
      transsysLattice.timestep = i
      if timestepSignal :
        if i == timestepSignal[1] :
          transsysLattice.timestep_factor_concentration(timestepSignal[0])
      # Append the deepcopy in each timestep.
      latticeTimeseries.append(copy.deepcopy(transsysLattice))
      transsysLattice.update_function(i)
#      transsysLattice.timestep = i
    return latticeTimeseries


  def write_factor_table(self, fileObj) :
    """
    Writes the whole factor table in a file.

    Calls the L{TranssysInstanceLattice.write_table} for the whole transsys
    lattice timeseries.

    @param fileObj: An open file object ready for writing.
    @type fileObj: C{file}
    @rtype: C{None}
    @precondition: An open, ready for writting file object.
    """
    # Then write the whole table.
    for til in self.latticeTimeseries :
      til.write_table(fileObj)


  def max_factor_concentration(self, factorName=None) :
    """
    The maximal observed factor concentration.

    Returns the maximum value of a factor concentration that is observed during
    the whole simulation proccess. If it is called without specifing a factor
    name calculate the maximum concentration value of all factors.

    @param factorName: The name of the factor that we maximum is wanted.
    @type factorName: C{str}
    @return: Maximum of factor concentration.
    @rtype: C{float}
    """
    maxFC = 0.0
    for til in self.latticeTimeseries :
      for ti in til.transsys_instance_list() :
        if not factorName :
          for fc in ti.factor_concentration :
            if fc > maxFC :
              maxFC = fc
        else :
          factorIndex = ti.transsys_program.find_factor_index(factorName)
          fc = ti.factor_concentration[factorIndex]
          if fc > maxFC :
            maxFC = fc
    return maxFC


  def factor_table(self, fileObj) :
    """
    Returns the factor table from a file that contains it.

    @todo: The return value might be usefull as an instance variable...
    @param fileObj: An open file object ready for reading.
    @type fileObj: C{file}
    @rtype: C{str}
    @precondition: An open, ready for reading file object.
    """
    factorTable = fileObj.read()
    return factorTable




def pickle(object) :
  """
  Return the pickled representation of the object as a string.
  @returns: A string representing the pickled object.
  @rtype: C{str}
  """
  return cPickle.dumps(object)


def unpickle(fileObj) :
  """
  Reconstructs and returns the original odject hierarchy from its pickled
  file.
  @returns: The pickled object.
  @rtype: C{class 'object'}
  @precondition: An open, ready for reading file object.
  """
  return cPickle.load(fileObj)


def generate_pgm(fileObj, transsysLattice, pgmFactor, maxCon) :
  """
  Produce the contents of the portable gray map (.pgm) image file.

  Write the to a .pgm file and then close it.

  Suitable for view/animate/analyse with the Imagemagick graphics suite.

  @param fileObj: An open raedy for writting file object.
  @type fileObj: C{file}
  @param transsysLattice: The whole transsys lattice.
  @type transsysLattice: C{class 'TranssysLattice'}
  @param pgmFactor: The name of the factor that will be imaged.
  @type pgmFactor: C{str}
  @param maxCon: The maximum nuber of factor concentration that can be drawn
  by the imaging proccedure. (constant provided by the user)
  @type maxCon: C{int}
  @rtype: C{None}
  @precondition: An open, ready for writting file object.
  """
  # Find the index (on  the transsys program level) of the apropriate factor.
  factorIndex = transsysLattice.transsysProgram.find_factor_index(pgmFactor)
  # Generate the apropriate .pgm image format.
  pgmMagic = 'P2' # .pgm file magic line, P2 for .pgm,
  pgmComment = '# Raster .pgm image representation of factor' + pgmFactor + ' on a' + transsysLattice.name
  pgmSize = str(transsysLattice.size[0]) + ' ' + str(transsysLattice.size[1])
  pgmMaxval = 256
  pgmRaster = ''
  for i in xrange(transsysLattice.size[0]) :
    pgmRasterLine = ''
    for j in xrange(transsysLattice.size[1]) :
      # Safeguard for zero maximum factor concentration.
      if maxCon == 0 :
        pixValue = pgmMaxval
      else :
        pixValue = int(pgmMaxval - round(transsysLattice.lattice[i][j].factor_concentration[factorIndex] * pgmMaxval / maxCon))
      # Safeguard for negative color values.
      if pixValue < 0 :
        raise StandardError, 'Netpbm raster image gets a negative value. Please select a biger number for the maximum factor concentration (option -m)'
      pgmRasterLine = pgmRasterLine + str(pixValue) + ' '
    pgmRaster = pgmRaster + pgmRasterLine + '\n'
  pgmText = pgmMagic + '\n' + pgmComment + '\n' + pgmSize + '\n' + str(pgmMaxval)  + '\n' + pgmRaster
  fileObj.write(pgmText)


def print_summary_statistics(tlt, fileObj) :
  """
  Function to calculate and print the collection statistics of the simulator.

  The calculation of the statistics is implemented in the
  L{transsys.CollectionStatistics} class.

  @param tlt: A transsys lattice timeseries.
  @type tlt: C{class 'TranssysLatticeTimeseries'}
  @param fileObj: An open ready for writting file object.
  @type fileObj: C{file}
  @precondition: An open, ready for writting file object.
  """
  statList = []
  for til in tlt :
    statList.append(til.statistics())
  fileObj.write('# Summary file for aggregate statistics (Header)\n')
  fileObj.write('timestep\tfactor\taverage\tstddev\tentropy\n')
  for i, stat in enumerate(statList) :
    for factor in stat.transsys_program.factor_list :
      fName = factor.name
      fileObj.write('%i\t%s\t%f\t%e\t%e\n' % (i, fName, stat.average[stat.transsys_program.find_factor_index(fName)], stat.standard_deviation[stat.transsys_program.find_factor_index(fName)], stat.shannon_entropy[stat.transsys_program.find_factor_index(fName)]))

