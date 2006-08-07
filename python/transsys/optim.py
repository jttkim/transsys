#!/usr/bin/env python

import sys
import copy
import math
import types
import os
import popen2
import random
import string
# import collections

import transsys
import utils


class Interval :
  """Represents an interval [a, b].
"""

  def __init__(self, a, b) :
    if a < b :
      self.lower = a
      self.upper = b
    else :
      self.lower = b
      self.upper = a


  def __str__(self) :
    return 'Interval(%f, %f)' % (self.lower, self.upper)


  def contains(self, other) :
    """Determine whether interval C{other} is a subset of this interval itself."""
    if type(other) is types.FloatType :
      return self.lower <= other and other <= self.upper
    if not isinstance(other, Interval) :
      raise StandardError, 'Interval.contains: illegal type'
    if self.lower <= other.lower and self.upper >= other.upper :
      return True
    else :
      return False


  def intersection(self, other) :
    """Compute the intersection of this interval itself with C{other}.
@return: An instance of C{Interval}, or C{None} if the intersection is
  empty.
"""
    if self.upper < other.lower or other.upper < self.lower :
      return None
    if self.contains(other) :
      return Interval(self.lower, self.upper)
    if other.contains(self) :
      return Interval(other.lower, other.upper)
    if self.lower < other.lower :
      return Interval(other.lower, self.upper)
    else :
      return Interval(self.lower, other.upper)


def active_and_inactive_set(rule, instance_series, factor_index) :
  active_set = []
  inactive_set = []
  for tsi in instance_series :
    if tsi.rule is rule :
      active_set.append(tsi.transsys_instance.factor_concentration[factor_index])
    else :
      inactive_set.append(tsi.transsys_instance.factor_concentration[factor_index])
  return active_set, inactive_set


def stddev_disparity(rule, instance_series, factor_index) :
  active_set, inactive_set = active_and_inactive_set(rule, instance_series, factor_index)
  if len(active_set) == 0 or len(inactive_set) == 0 :
    return 0.0
  m_active, s_active = utils.mean_and_stddev(active_set)
  m_inactive, s_inactive = utils.mean_and_stddev(inactive_set)
  m_diff = abs(m_active - m_inactive)
  if m_diff == 0.0 :
    return 1.0
  return 1.0 - 1.0 / (1.0 + (s_active + s_inactive) / m_diff)


def overlap_disparity(rule, instance_series, factor_index) :
  """Disparity of expression of factor in symbols activating rule and symbols
not activating rule."""
  active_set, inactive_set = active_and_inactive_set(rule, instance_series, factor_index)
  if len(active_set) == 0 or len(inactive_set) == 0 :
    return 1.0
  active_interval = Interval(min(active_set), max(active_set))
  inactive_interval = Interval(min(inactive_set), max(inactive_set))
  if active_interval.contains(inactive_interval) or inactive_interval.contains(active_interval) :
    return 1.0
  intersection = active_interval.intersection(inactive_interval)
  if intersection is None :
    return 0.0
  n = 0
  complete_set = active_set + inactive_set
  for x in complete_set :
    if intersection.contains(x) :
      n = n + 1
  return float(n) / float(len(complete_set))


class FitnessResult :
  """Base class for results of objective (or fitness) functions.

Calling objective functions should return instances of this class
or one of its subclasses. Optimisers use the C{fitness} instance
variable. Subclasses may introduce further instance variables to
describe the result achieved, break it down to individual components
of the transsys program etc.

@ivar fitness: the fitness value, i.e. the value computed by the
  objective function.
"""
  def __init__(self, fitness) :
    self.fitness = fitness


  def __str__(self) :
    return 'FitnessResult(%g)' % self.fitness


class LsysDisparityFitnessResult(FitnessResult) :

  def __init__(self, fitness, best_factor_list) :
    FitnessResult.__init__(self, fitness)
    self.best_factor_list = best_factor_list


def flat_symbol_instance_list(lstring_list, transsys_program = None) :

  def filter(si, transsys_program) :
    if transsys_program is None :
      return True
    if si.transsys_instance is None :
      return False
    if si.transsys_instance.transsys_program is transsys_program :
      return True
    return False
  
  flat_list = []
  for lstring in lstring_list :
    for si in lstring.symbol_list :
      if filter(si, transsys_program) :
        flat_list.append(si)
  return flat_list


def disparity_fitness(lsys_program, transsys_program, factor_names, num_timesteps, disparity) :
  if len(transsys_program.factor_list) == 0 :
    return 1.0
  fitness = 0.0
  best_factor_list = []
  factor_indices = []
  for factor_name in factor_names :
    factor_indices.append(transsys_program.find_factor_index(factor_name))
  instance_series = flat_symbol_instance_list(lsys_program.derivation_series(num_timesteps), transsys_program)
  for rule in lsys_program.rules :
    fmin = None
    for factor_index in factor_indices :
      f = disparity(rule, instance_series, factor_index)
      if fmin is None :
        fmin = f
        best_factors = [transsys_program.factor_list[factor_index].name]
      elif f < fmin :
        fmin = f
        best_factors = [transsys_program.factor_list[factor_index].name]
      elif f == fmin :
        best_factors.append(transsys_program.factor_list[factor_index].name)
    # print '  rule "%s": %f' % (rule.name, f)
    best_factor_list.append((rule.name, tuple(best_factors), ))
    fitness = fitness + fmin
  return LsysDisparityFitnessResult(fitness / float(len(lsys_program.rules)), best_factor_list)


class Mutator :

  def __init__(self, rng, replacement_rate, insertion_rate = 0.0, deletion_rate = 0.0, alphabet = 'acgt') :
    self.rng = rng
    self.replacement_rate = replacement_rate
    self.insertion_rate = insertion_rate
    self.deletion_rate = deletion_rate
    self.alphabet = alphabet


  def mutate(self, seq) :
    if self.insertion_rate != 0.0 or self.deletion_rate != 0.0 :
      raise StandardError, 'insertion / deletion not implemented yet'
    mutseq = ''
    for i in xrange(len(seq)) :
      if self.rng.random() < self.replacement_rate :
        mutseq = mutseq + self.alphabet[self.rng.randrange(len(self.alphabet))]
      else :
        mutseq = mutseq + seq[i]
    return mutseq


def hillclimb(objective_function, decoder, initial_dnaseq, num_generations, mutator) :
  dnaseq = initial_dnaseq
  fitness = objective_function(decoder.decode_transsys('implant', dnaseq)).fitness
  for g in xrange(num_generations) :
    print '%d: %f' % (g, fitness)
    mutant_dnaseq = mutator.mutate(dnaseq)
    mutant_fitness = objective_function(decoder.decode_transsys('implant', mutant_dnaseq)).fitness
    d = utils.hamming_distance(dnaseq, mutant_dnaseq)
    # print 'current: %f, mutant: %f, distance: %d (rel: %f)' % (fitness, mutant_fitness, d, float(d) / float(len(dnaseq)))
    if mutant_fitness <= fitness :
      dnaseq = mutant_dnaseq
      fitness = mutant_fitness
  return dnaseq


# moving to utils
def get_value_nodes(transsys_program, gene_name_list, factor_name_list) :
  value_expression_list = []
  if factor_name_list is None :
    value_expression_list.extend(transsys_program.getFactorValueNodes())
  else :
    for fn in factor_name_list :
      factor = transsys_program.find_factor(fn)
      value_expression_list.extend(factor.getValueNodes())
  if gene_name_list is None :
    value_expression_list.extend(transsys_program.getGeneValueNodes())
  else :
    for gn in gene_name_list :
      gene = transsys_program.find_gene(gn)
      value_expression_list.extend(gene.getValueNodes())
  return value_expression_list



# moving to utils
def randomise_transsys_values(transsys_program, random_function, gene_name_list = None, factor_name_list = None) :
  value_expression_list = get_value_nodes(transsys_program, gene_name_list, factor_name_list)
  for n in value_expression_list :
    n.value = random_function()


class ExponentialFunction(object) :
  """Exponential parameter transformation.

This is a callable class implementing the transformation function
y = minValue + exp(scaleFactor * x).
  """

  def __init__(self, scaleFactor = 1.0, minValue = 0.0) :
    """Constructor."""
    self.scaleFactor = scaleFactor
    self.minValue = minValue


  def __call__(self, x) :
    """Compute function value."""
    return self.minValue + math.exp(self.scaleFactor * x)


  def inverse(self, y) :
    """Compute inverse.

Notice that the inverse diverges as y approaches minValue.
Exactly at minValue, the inverse method attempts to return
the smallest value that can be transformed without a floating
point overflow.
"""
    if y < self.minValue :
      raise StandardError, 'argument %f < minimal value %f' % (y, self.minValue)
    if y == self.minValue :
      # 700 is 300 * log(10), and as floating point ranges typically
      # extend to around 1e300, this should be ok on normal platforms
      return 700.0 / self.scaleFactor
    return 1.0 / self.scaleFactor * math.log(y - self.minValue)


class SigmoidFunction(object) :
  """Sigmoid parameter transformation.

This is a callable class implementing the transformation function
y = minValue + (maxValue - minValue) / (1.0 + exp(-scaleFactor * x)).
"""
  def __init__(self, scaleFactor = 1.0, minValue = 0.0, maxValue = 1.0) :
    """Constructor."""
    self.scaleFactor = scaleFactor
    self.minValue = minValue
    self.valueRange = maxValue - minValue


  def __call__(self, x) :
    """Compute function value."""
    return self.minValue + self.valueRange / (1.0 + math.exp(-self.scaleFactor * x))


  def inverse(self, y) :
    """Compute inverse.

Notice that the inverse diverges towards the borders of the sigmoid
function's range. Exactly at the borders, the inverse method attempts
to return the closest value that can be transformed without a floating
point overflow.
"""
    maxValue = self.minValue + self.valueRange
    if y < self.minValue or y > maxValue :
      raise StandardError, 'argument %f out of range [%f, %f]' % (y, self.minValue, maxValue)
    # FIXME: proper DBL_MIN, DBL_MAX like values should be used here
    if y == self.minValue :
      return -700.0 / self.scaleFactor
    if y == maxValue :
      return 700.0 / self.scaleFactor
    return -1.0 / self.scaleFactor * math.log(self.valueRange / (y - self.minValue) - 1.0)


class ParameterTransformer(object) :
  """Base class for parameter transformers."""

  pass


class IdentityParameterTransformer(ParameterTransformer) :
  """No transformation, simply copy parameters from and to transsys program."""

  def __init__(self) :
    pass


  def getParameters(self, transsys_program) :
    nodes = transsys_program.getValueNodes()
    return map(lambda n : n.value, nodes)


  def setParameters(self, parameter_list, transsys_program) :
    nodes = transsys_program.getValueNodes()
    if len(nodes) != len(parameter_list) :
      raise StandardError,  'parameter list incompatible with transsys program'
    for i in xrange(len(parameter_list)) :
      nodes[i].value = parameter_list[i]


class ConstrainedParameterTransformer(ParameterTransformer) :
  """Transform unconstrained parameters to ranges that are sensible for
transsys programs.

Unconstrained parameters p in ]-Inf, Inf[ are transformed to
factor parameters (decay and diffusibility), constrained to ]0, 1[
by the sigmoid function y = 1 / (1 + exp(-x)), and to gene
parameters (constitutive, amax, aspec), constrained to ]0, +Inf[
by the exponential function y = exp(x).

Parameters are passed out and in using lists, and they are matched
to values in nodes contained within the transsys program by position.
Therefore, the transsys program must not be altered between getting
parameters and setting parameters.

Future versions of this class may use a proper indexable and iterable
class instead of a "flat", unstructured list.
"""

  # FIXME: should use proper class rather than unstructured list for parameters
  def __init__(self, scaleFactor = 1.0) :
    self.factorValueTransformation = SigmoidFunction(scaleFactor)
    self.geneValueTransformation = ExponentialFunction(scaleFactor)


  def getParameters(self, transsys_program) :
    """Get unconstrained parameters from a transsys program.

@param transsys_program: the program from which to get the parameters
@return: a list of the unconstrained parameters
    """
    parameter_list = []
    for p in map(lambda n : n.value, transsys_program.getFactorValueNodes()) :
      parameter_list.append(self.factorValueTransformation.inverse(p))
    for p in map(lambda n : n.value, transsys_program.getGeneValueNodes()) :
      parameter_list.append(self.geneValueTransformation.inverse(p))
    return parameter_list


  def setParameters(self, parameter_list, transsys_program) :
    """Set values in a transsys program based on a list of unconstrained parameters.

Notice that the parameter list must have the same length as that obtained
by C{getParameters}.
"""
    fn = transsys_program.getFactorValueNodes()
    gn = transsys_program.getGeneValueNodes()
    if len(parameter_list) != len(fn) + len(gn) :
      raise StandardError, 'parameter list incompatible with transsys program'
    i = 0
    for n in fn :
      n.value = self.factorValueTransformation(parameter_list[i])
      i = i + 1
    for n in gn :
      n.value = self.geneValueTransformation(parameter_list[i])
      i = i + 1


class TranssysTypedParameterTransformer(ParameterTransformer) :
  """
Parameter transformer differentiating decay, diffusibility, constitutive,
aspec, amax, rspec and rmax nodes.
"""
  def __init__(self, decayTransformation, diffusibilityTransformation, constitutiveTransformation, aspecTransformation, amaxTransformation, rspecTransformation, rmaxTransformation) :
    self.decayTransformation = decayTransformation
    self.diffusibilityTransformation = diffusibilityTransformation
    self.constitutiveTransformation = constitutiveTransformation
    self.aspecTransformation = aspecTransformation
    self.amaxTransformation = amaxTransformation
    self.rspecTransformation = rspecTransformation
    self.rmaxTransformation = rmaxTransformation
    self.transsys_program = None


  def getParameters(self, transsys_program) :
    """Get unconstrained parameters from a transsys program.

The C{transsys_program} is remembered and subsequent calls to
C{setParameters} are only accepted for the same C{TranssysProgram}
instance.

@param transsys_program: the program from which to get the parameters
@return: a list of the unconstrained parameters
    """
    parameter_list = []
    for p in map(lambda n : n.value, transsys_program.getDecayValueNodes()) :
      parameter_list.append(self.decayTransformation.inverse(p))
    for p in map(lambda n : n.value, transsys_program.getDiffusibilityValueNodes()) :
      parameter_list.append(self.diffusibilityTransformation.inverse(p))
    for p in map(lambda n : n.value, transsys_program.getConstitutiveValueNodes()) :
      parameter_list.append(self.constitutiveTransformation.inverse(p))
    for p in map(lambda n : n.value, transsys_program.getActivateSpecValueNodes()) :
      parameter_list.append(self.aspecTransformation.inverse(p))
    for p in map(lambda n : n.value, transsys_program.getActivateMaxValueNodes()) :
      parameter_list.append(self.amaxTransformation.inverse(p))
    for p in map(lambda n : n.value, transsys_program.getRepressSpecValueNodes()) :
      parameter_list.append(self.rspecTransformation.inverse(p))
    for p in map(lambda n : n.value, transsys_program.getRepressMaxValueNodes()) :
      parameter_list.append(self.rmaxTransformation.inverse(p))
    self.transsys_program = transsys_program
    return parameter_list


  def setParameters(self, parameter_list, transsys_program) :
    """Set values in a transsys program based on a list of unconstrained parameters.

Notice that the parameter list must have the same length as that obtained
by C{getParameters}.

@raise StandardError: when the C{transsys_program} is not the one used
  in the last C{getParameters} call or the parameter list is too long or
  too short
"""
    if transsys_program is not self.transsys_program :
      raise StandardError, 'parameters must be set on same program from which they were got'
    decay_nodes = transsys_program.getDecayValueNodes()
    diffusibility_nodes = transsys_program.getDiffusibilityValueNodes()
    constitutive_nodes = transsys_program.getConstitutiveValueNodes()
    aspec_nodes = transsys_program.getActivateSpecValueNodes()
    amax_nodes = transsys_program.getActivateMaxValueNodes()
    rspec_nodes = transsys_program.getRepressSpecValueNodes()
    rmax_nodes = transsys_program.getRepressMaxValueNodes()
    if len(parameter_list) != len(decay_nodes) + len(diffusibility_nodes) + len(constitutive_nodes) + len(aspec_nodes) + len(amax_nodes) + len(rspec_nodes) + len(rmax_nodes):
      raise StandardError, 'parameter list incompatible with transsys program'
    i = 0
    for n in decay_nodes :
      n.value = self.decayTransformation(parameter_list[i])
      i = i + 1
    for n in diffusibility_nodes :
      n.value = self.diffusibilityTransformation(parameter_list[i])
      i = i + 1
    for n in constitutive_nodes :
      n.value = self.constitutiveTransformation(parameter_list[i])
      i = i + 1
    for n in aspec_nodes :
      n.value = self.aspecTransformation(parameter_list[i])
      i = i + 1
    for n in amax_nodes :
      n.value = self.amaxTransformation(parameter_list[i])
      i = i + 1
    for n in rspec_nodes :
      n.value = self.rspecTransformation(parameter_list[i])
      i = i + 1
    for n in rmax_nodes :
      n.value = self.rmaxTransformation(parameter_list[i])
      i = i + 1


class OptimisationResult(object) :
  """Base class for returning results of an optimisation.

An optimisation result consists of the optimised version of the
transsys program, its objective function value, and, optionally,
an optimisation log.

The optimisation log is a list of optimisation records with a
structure that depends on the optimiser. Optimisation records
can be printed on a line, and they provide a table header. This
allows the optimisation log to be dumped in a format ready for
R's \C{read.table} function.

@ivar optimised_transsys_program: The transsys program resulting
  from the optimisation process.
@ivar objectiveOptimum: The result of evaluating the optimised transsys
  program using the objective function.
@ivar optimisation_log: A trace of the optimisation process, provided
  as a list of tuples.
"""

  def __init__(self, tp, objectiveOptimum, optimisation_log = None) :
    self.optimised_transsys_program = tp
    self.objectiveOptimum = objectiveOptimum
    self.optimisation_log = optimisation_log


  def write_log(self, f, column_prefix = '', write_header = True, header_prefix = '') :
    if self.optimisation_log is None :
      return
    if len(self.optimisation_log) > 0 :
      if header_prefix != '' :
        if not header_prefix[-1].isspace() :
          header_prefix = header_prefix + ' '
      if column_prefix != '' :
        if not column_prefix[-1].isspace() :
          column_prefix = column_prefix + ' '
      if write_header :
        f.write('%s%s\n' % (header_prefix, self.optimisation_log[0].table_header()))
      for l in self.optimisation_log :
        f.write('%s%s\n' % (column_prefix, str(l)))


class AbstractOptimisationRecord(object) :
  """Abstract base class for optimisation records"""

  def __init__(self) :
    pass
    

class SimulatedAnnealingRecord(AbstractOptimisationRecord) :

  def __init__(self, obj = None, obj_alt = None, temperature = None, stepsize = None, stepvector_entropy = None, stepvector_max = None, state_distance = None, p_accept = None, accept = None) :
    self.obj = obj
    self.obj_alt = obj_alt
    self.temperature = temperature
    self.stepsize = stepsize
    self.stepvector_entropy = stepvector_entropy
    self.stepvector_max = stepvector_max
    self.state_distance = state_distance
    self.p_accept = p_accept
    self.accept = accept


  def __str__(self) :
    return utils.table_row([self.obj, self.obj_alt, self.temperature, self.stepsize, self.stepvector_entropy, self.stepvector_max, self.state_distance, self.p_accept, self.accept])


  def table_header(self) :
    return 'obj obj_alt temperature stepsize stepvector_entropy stepvector_max state_distance p_accept accept'


class GradientOptimisationRecord(AbstractOptimisationRecord) :

  def __init__(self, obj, stepsize, numEvaluations, gradient) :
    self.obj = obj
    self.stepsize = stepsize
    self.numEvaluations = numEvaluations
    gradient_abs = map(abs, gradient)
    self.gradient_entropy = utils.shannon_entropy(gradient_abs)
    self.gradient_max = max(gradient_abs)


  def __str__(self) :
    return utils.table_row([self.obj, self.stepsize, self.numEvaluations, self.gradient_entropy, self.gradient_max])


  def table_header(self) :
    return 'obj stepsize numEvaluations gradient_entropy gradient_max'


class AbstractObjectiveFunction(object) :
  """Base class for objective functions.

This is an abstract class, methods just raise exceptions. Subclasses
should be callable and take a transsys program as a parameter. The
returned value should be an instance of C{FitnessResult} or a
subclass.
"""

  def __init__(self) :
    """Abstract constructor method, raises exception."""
    raise StandardError, 'cannot instantiate AbstractObjectiveFunction'


  def __call__(self, transsys_program) :
    """Abstract call method, raises exception.

@param transsys_program: the transsys program to be evaluated
    """
    raise StandardError, 'cannot call AbstractObjectiveFunction'


class ExpressionSeriesObjective(AbstractObjectiveFunction) :
  """Abstract base class for objective functions based on comparing observed
with desired expression profiles.

Inherits unimplemented call method from C{AbstractObjectiveFunction}.
"""

  def __init__(self, f = None) :
    """Construct an instance holding an empty series."""
    self.series = {}
    if f is not None :
      self.readProfiles(f)


  def __str__(self) :
    """Produce a string suitable for saving the profiles."""
    s = ''
    for factor_name in self.series.keys() :
      s = s + '%s:' % factor_name
      for xlevel in self.series[factor_name] :
        s = s + ' %f' % xlevel
      s = s + '\n'
    return s


  def series_length(self) :
    """Accessor to get the length of the expression profiles.

@return: Length of the expression profiles, or C{None} if no
  profiles contained in instance data.
"""
    if self.series is None :
      return None
    if len(self.series) == 0 :
      return None
    k = self.series.keys()[0]
    return len(self.series[k])


  def readProfiles(self, f) :
    """Read expression profiles from a file.

The format expected by this: One line per factor, each line consists
of whitespace separated words, the first word is the factor, followed
by the expression levels given as floating point numbers.

All profiles must have the same length.

@param f: file to read from
"""
    series_length = None
    s = f.readline()
    while s :
      w = s.split(':')
      factor_name = w[0]
      if factor_name == '' :
        raise StandardError, 'empty factor name'
      if factor_name in self.series.keys() :
        raise StandardError, 'multiple series for factor "%s"' % factor_name
      self.series[factor_name] = []
      xlist = w[1].split()
      if series_length is None :
        series_length = len(xlist)
      elif series_length != len(xlist) :
        raise StandardError, 'unequal series lengths (%d != %d)' % (series_length, len(xlist))
      self.series[factor_name] = map(float, xlist)
      s = f.readline()
    # print 'series_length:', series_length
    # continue here -- with work on trip to London...


class ExpressionSeriesSquareSumObjective(ExpressionSeriesObjective) :
  """The sum of squared Euclidean distances to a desired set of
expression profiles as an objective function."""

  def __init__(self, f = None) :
    super(ExpressionSeriesSquareSumObjective, self).__init__(f)


  def __call__(self, transsys_program) :
    l = self.series_length()
    if l is None :
      raise StandardError, 'call of uninitialised series objective function'
    if l == 0 :
      raise StandardError, 'call of objective function with series length 0'
    # this is a kludge -- will become much better once a time
    # series class is available
    for factor_name in self.series.keys() :
      if transsys_program.find_factor_index(factor_name) == -1 :
        raise StandardError, 'failed to find factor "%s" in transsys "%s"' % (factor_name, transsys_program.name)
    ti = transsys.TranssysInstance(transsys_program)
    tseries = ti.time_series(l)
    sq_sum = 0.0
    for factor_name in self.series.keys() :
      factor_index = transsys_program.find_factor_index(factor_name)
      # sys.stderr.write('factor "%s" --> index %d\n' % (factor_name, factor_index))
      v = []
      for i in xrange(l) :
        x = tseries[i].factor_concentration[factor_index]
        if utils.is_nan(x) :
          print str(tseries[i])
          raise StandardError, 'expression level of factor "%s" is NaN' % factor_name
        v.append(x)
      sq = transsys.utils.euclidean_distance_squared(self.series[factor_name], v)
      # sys.stderr.write('factor "%s": sq = %f = d2(%s, %s)\n' % (factor_name, sq, str(self.series[factor_name]), str(v)))
      sq_sum = sq_sum + sq
    return FitnessResult(sq_sum)


class ExpressionSeriesCorrelationObjective(ExpressionSeriesObjective) :
  """Correlation based objective function.

The objective value for each factor is 1 - r, where r is the
Pearson correlation coefficient of the desired and the observed
expression profile of the factor. The objective value of a
transsys program is the sum of the objective values of all
factors.

If all expression levels of a factor are identical in either
the desired or the observed profile, the objective value for
that factor is 2.
"""

  def __init__(self, f = None) :
    super(ExpressionSeriesCorrelationObjective, self).__init__(f)


  def __call__(self, transsys_program) :
    l = self.series_length()
    if l is None :
      raise StandardError, 'call of uninitialised series objective function'
    if l == 0 :
      raise StandardError, 'call of objective function with series length 0'
    # this is a kludge -- will become much better once a time
    # series class is available
    for factor_name in self.series.keys() :
      if transsys_program.find_factor_index(factor_name) == -1 :
        raise StandardError, 'failed to find factor "%s" in transsys "%s"' % (factor_name, transsys_program.name)
    ti = transsys.TranssysInstance(transsys_program)
    tseries = ti.time_series(l)
    cc_sum = 0.0
    for factor_name in self.series.keys() :
      factor_index = transsys_program.find_factor_index(factor_name)
      if min(self.series[factor_name]) < max(self.series[factor_name]) :
        # sys.stderr.write('factor "%s" --> index %d\n' % (factor_name, factor_index))
        v = []
        for i in xrange(l) :
          x = tseries[i].factor_concentration[factor_index]
          if utils.is_nan(x) :
            print str(tseries[i])
            raise StandardError, 'expression level of factor "%s" is NaN' % factor_name
          v.append(x)
        if min(v) < max(v) :
          cc = 1.0 - transsys.utils.correlation_coefficient(self.series[factor_name], v)
          cc_sum = cc_sum + cc
        else :
          cc_sum = cc_sum + 2.0
    return FitnessResult(cc_sum)


class AbstractOptimiser(object) :
  """Abstract base class for optimisers.

@ivar verbose: controls verbosity (of optimisation process)
  """

  def __init__(self, verbose = False) :
    """Partial constructor method.

This constructor just sets the C{verbose} instance variable
    """
    self.verbose = verbose


  def optimise(self, transsys_program_init, objective_function) :
    """Abstract optimise method.

Subclasses of C{AbstractOptimiser} must override this method and
conform to the parameter signature specified by this abstract method
in doing so. Implementations of the C{optimise} method may take
additional parameters, which must have suitable defaults so a call
with just the two parameters defined here remains valid. It is
generally preferable to provide such parameters through instance
variables.

@param transsys_program_init: The transsys program to be optimised.
  Optimisers should not modify this transsys program, modifications
  should be done to a copy.
@param objective_function: The objective function used to optimise the
  transsys program."""
    raise StandardError, 'call to unimplemented AbstractOptimiser optimise method'


  def check_savefile_magic(self, s) :
    raise StandardError, 'unimplemented abstract method'


  def parse_variables(self, f) :
    raise StandardError, 'unimplemented abstract method'


  def parse(self, f) :
    """Parse optimiser parameters from file C{f}."""
    l = f.readline()
    if l == '' :
      raise StandardError, 'expected magic but got EOF'
    l = l.strip()
    if not self.check_savefile_magic(l) :
      raise StandardError, 'bad magic "%s"' % l
    self.parse_variables(f)


class SimulatedAnnealer(AbstractOptimiser) :
  """Simple Simulated Annealing tool.

The parameterisation of stepping is based on the concepts outlined in
"Mathematical Optimization", with (rather major) modifications.

@cvar PERTURBATION_METHOD_UNIFORM: specifies local search by perturbation of
  current state by uniformly distributed random values
@cvar PERTURBATION_MOETHO_GAUSS: specifies local search by perturbation of
  current state by random values with Gaussian distribution
@ivar stepsize_learning_rate: learning rate for stepsize
@ivar target_improvement_ratio: target ratio of improvements in
  objective value among candidate solutions generated
@ivar stepsize_init: initial step size
@ivar stepvector_learning_rate: learning rate for updating stepvector
@ivar stepvector_learning_decay: parameter for decaying the "learned"
  stepvector direction by "learning" the unbiased stepvector
@ivar cooling_rate: current temperature is multiplied by this
  value in each iteration.
@ivar temperature_init: initial temperature. May be set to C{None},
  in this case, the initial value is the objective function
  value of the transsys program to be optimised.
@ivar termination_temperature: temperature threshold below which the annealing
  process is terminated.
@ivar termination_objective: objective function value at or below which
  process is terminated
@ivar termination_iteration: number of iteration after which process
  is terminated
@ivar transformer: parameter transformer.
@ivar rng: random number generator.
@ivar perturbation_method: perturbation method for local search
@ivar verbose: controls verbosity.
"""

  savefile_magic = 'SimulatedAnnealer-0.1'
  PERTURBATION_METHOD_UNIFORM = 1
  PERTURBATION_METHOD_GAUSS = 2

  def __init__(self, cooling_rate = 0.995, temperature_init = None, termination_temperature = None, termination_objective = None, termination_iteration = None, stepsize_learning_rate = 0.0, target_improvement_ratio = 0.2, stepsize_init = 1.0, stepvector_learning_rate = 0.0, stepvector_learning_decay = 0.0, transformer = None, rng = None, perturbation_method = PERTURBATION_METHOD_UNIFORM, verbose = False) :
    """Constructor."""
    super(SimulatedAnnealer, self).__init__(verbose)
    self.stepsize_learning_rate = stepsize_learning_rate
    self.target_improvement_ratio = target_improvement_ratio
    self.stepsize_init = stepsize_init
    self.stepvector_learning_rate = stepvector_learning_rate
    self.stepvector_learning_decay = stepvector_learning_decay
    self.cooling_rate = cooling_rate
    self.temperature_init = temperature_init
    self.termination_temperature = termination_temperature
    self.termination_objective = termination_objective
    self.termination_iteration = termination_iteration
    self.perturbation_method = perturbation_method
    if transformer is None :
      self.transformer = IdentityParameterTransformer()
    else :
      self.transformer = transformer
    if rng is None :
      self.rng = random.Random(1)
    else :
      self.rng = rng


  def setPerturbationMethod(self, pMethodString) :
    """Set the perturbation method, 'uniform' or 'gauss'."""
    if pMethodString == 'uniform' :
      self.perturbation_method = self.PERTURBATION_METHOD_UNIFORM
    elif pMethodString == 'gauss' :
      self.perturbationMethod = self.PERTURBATION_METHOD_GAUSS
    else :
      raise StandardError, 'unsupported perturbation method "%s"' % pMethodString


  def check_savefile_magic(self, s) :
    """Check whether C{s} is the "magic word" indicating the
beginning of a valid save file.
"""
    return s == self.savefile_magic


  def parse_variables(self, f) :
    """Get the values of instance variables from file C{f}."""
    self.stepsize_learning_rate = utils.parse_float(f, 'stepsize_learning_rate')
    self.target_improvement_ratio = utils.parse_float(f, 'target_improvement_ratio')
    self.stepsize_init = utils.parse_float(f, 'stepsize_init')
    self.stepvector_learning_rate = utils.parse_float(f, 'stepvector_learning_rate')
    self.stepvector_learning_decay = utils.parse_float(f, 'stepvector_learning_decay')
    self.cooling_rate = utils.parse_float(f, 'cooling_rate')
    self.temperature_init = utils.parse_float(f, 'temperature_init', allowNone = True)
    self.termination_temperature = utils.parse_float(f, 'termination_temperature', allowNone = True)
    self.termination_objective = utils.parse_float(f, 'termination_objective', allowNone = True)
    self.termination_iteration = utils.parse_int(f, 'termination_iteration', allowNone = True)
    self.setPerturbationMethod(utils.parse_string(f, 'perturbation_method').strip())
    rndseed = utils.parse_int(f, 'rndseed')
    self.rng = random.Random(rndseed)
    s = utils.parse_string(f, 'transformer').strip()
    if s == 'identity' :
      self.transformer = IdentityParameterTransformer()
    elif s == 'constrained' :
      self.transformer = ConstrainedParameterTransformer()
    else :
      raise StandardError, 'unknown parameter transformer "%s"' % s


  def __str__(self) :
    s = '%s\n' % self.savefile_magic
    s = s + '%s\n' % utils.name_value_pair(self.stepsize_learning_rate, 'stepsize_learning_rate')
    s = s + '%s\n' % utils.name_value_pair(self.target_improvement_ratio, 'target_improvement_ratio')
    s = s + '%s\n' % utils.name_value_pair(self.stepsize_init, 'stepsize_init')
    s = s + '%s\n' % utils.name_value_pair(self.stepvector_learning_rate, 'stepvector_learning_rate')
    s = s + '%s\n' % utils.name_value_pair(self.stepvector_learning_decay, 'stepvector_learning_decay')
    s = s + '%s\n' % utils.name_value_pair(self.cooling_rate, 'cooling_rate')
    s = s + '%s\n' % utils.name_value_pair(self.temperature_init, 'temperature_init')
    s = s + '%s\n' % utils.name_value_pair(self.termination_temperature, 'termination_temperature')
    s = s + '%s\n' % utils.name_value_pair(self.termination_objective, 'termination_objective')
    s = s + '%s\n' % utils.name_value_pair(self.termination_iteration, 'termination_iteration')
    if self.perturbation_method == self.PERTURBATION_METHOD_UNIFORM :
      s = s + '%s\n' % utils.name_value_pair('uniform', 'perturbation_method')
    elif self.perturbation_method == self.PERTURBATION_METHOD_GAUSS :
      s = s + '%s\n' % utils.name_value_pair('gauss', 'perturbation_method')
    else :
      raise StandardError, 'cannot encode perturbation method'
    s = s + '%s\n' % utils.name_value_pair('unknown', 'rndseed')
    if isinstance(self.transformer, IdentityParameterTransformer) :
      s = s + '%s\n' % utils.name_value_pair('identity', 'transformer')
    elif isinstance(self.transformer, ConstrainedParameterTransformer) :
      s = s + '%s\n' % utils.name_value_pair('constrained', 'transformer')
    else :
      raise StandardError, 'cannot encode parameter transformer'
    return s


  def verifyTermination(self) :
    """Verify that at least one termination condition is applicable.

If no termination condition is given, the annealer would run indefinitely.
Notice that setting C{termination_objective} only does constitute a
termination condition but there is no guarantee that it will ever be
satisfied. In contrast to this, setting C{termination_temperature}
implicitly limits the number of iterations and C{termination_iteration}
explicitly does so.

@raise StandardError: No termination condition applicable.
    """
    if self.termination_temperature is not None :
      return
    if self.termination_objective is not None :
      return
    if self.termination_iteration is not None :
      return
    raise StandardError, 'no termination condition set'


  def terminationCondition(self, temperature, objective, iteration) :
    """Check whether termination condition is satisfied.

@return: C{True} if termination condition is satisfied
@rtype: bool
    """
    if self.termination_temperature is not None :
      if temperature <= self.termination_temperature :
        if self.verbose :
          sys.stderr.write('SimulatedAnnealer: terminating because temperature %f <= %f\n' % (temperature, self.termination_temperature))
        return True
    if self.termination_objective is not None :
      if objective <= self.termination_objective :
        if self.verbose :
          sys.stderr.write('SimulatedAnnealer: terminating because objective %f <= %f\n' % (objective, self.termination_objective))
        return True
    if self.termination_iteration is not None :
      if iteration >= self.termination_iteration :
        if self.verbose :
          sys.stderr.write('SimulatedAnnealer: terminating because iteration %d >= %d\n' % (iteration, self.termination_iteration))
        return True
    return False


  def optimise(self, transsys_program_init, objective_function, stepvector_init = None, rndseed = None) :
    """Simulated annealing optimise method."""
    self.verifyTermination()
    if rndseed is not None :
      self.rng.seed(rndseed)
    transsys_program = copy.deepcopy(transsys_program_init)
    state = self.transformer.getParameters(transsys_program)
    if self.verbose :
      sys.stderr.write('SimulatedAnnealer: optimising %d values\n' % len(state))
    num_dimensions = len(state)
    unbiased_stepvector_component = 1.0 / math.sqrt(float(num_dimensions))
    if stepvector_init is None :
      stepvector = [unbiased_stepvector_component] * num_dimensions
    else :
      if len(stepvector_init) != num_dimensions :
        raise StandardError, 'stepvector incompatible with transsys program'
      stepvector = utils.normalised_vector(stepvector_init)
    stepsize = self.stepsize_init
    records = []
    obj = objective_function(transsys_program).fitness
    if self.temperature_init is None :
      temperature = obj
    else :
      temperature = self.temperature_init
    iteration = 0
    while not self.terminationCondition(temperature, obj, iteration) :
      if self.verbose :
        sys.stderr.write('starting: temp = %f, obj = %f\n' % (temperature, obj))
        # sys.stderr.write('  state: %s\n' % str(state))
      state_alt = []
      for i in xrange(num_dimensions) :
        if self.perturbation_method == self.PERTURBATION_METHOD_UNIFORM :
          # FIXME: this is marginally biased because the range of random
          # values includes -stepvector[i], but excludes +stepvector[i].
          state_alt.append(state[i] + (self.rng.random() - 0.5) * 2.0 * stepsize * stepvector[i])
        elif self.perturbation_method == self.PERTURBATION_METHOD_GAUSS :
          state_alt.append(self.rng.gauss(state[i], stepsize * stepvector[i]))
        else :
          raise StandardError, 'unimplemented perturbation method'
      delta_state = []
      for i in xrange(num_dimensions) :
        delta_state.append(state_alt[i] - state[i])
      state_distance = utils.euclidean_norm(delta_state)
      self.transformer.setParameters(state_alt, transsys_program)
      obj_alt = objective_function(transsys_program).fitness
      if self.verbose :
        # sys.stderr.write('  state_alt: %s, obj_alt = %f\n' % (str(state_alt), obj_alt))
        sys.stderr.write('  obj_alt = %f\n' % obj_alt)
        # sys.stderr.write(str(transsys_program))
      accept = obj_alt <= obj
      if accept :
        p_accept = 1.0
      else :
        p_accept = math.exp((obj - obj_alt) / temperature)
        if self.verbose :
          sys.stderr.write('delta_obj = %f, temperature = %f, p_accept = %f\n' % (obj_alt - obj, temperature, p_accept))
        accept = self.rng.random() < p_accept
      records.append(SimulatedAnnealingRecord(obj, obj_alt, temperature, stepsize, utils.shannon_entropy(stepvector), max(stepvector), state_distance, p_accept, accept))
      if obj_alt <= obj :
        stepsize = stepsize * (1.0 + self.stepsize_learning_rate / self.target_improvement_ratio)
        sdn = utils.normalised_vector(delta_state)
        v = []
        for i in xrange(num_dimensions) :
          v.append((1.0 - self.stepvector_learning_rate) * stepvector[i] + self.stepvector_learning_rate * abs(sdn[i]))
      else :
        stepsize = stepsize * (1.0 - self.stepsize_learning_rate)
        v = map(lambda x : (1.0 - self.stepvector_learning_decay) * x + self.stepvector_learning_decay * unbiased_stepvector_component, stepvector)
        stepvector = utils.normalised_vector(v)
      stepvector = utils.normalised_vector(v)
      if accept :
        if self.verbose :
          sys.stderr.write('  state_alt accepted (obj = %f, obj_alt = %f)\n' % (obj, obj_alt))
          # sys.stderr.write('  new stepvector: %s\n' % str(stepvector))
        state = state_alt
        obj = obj_alt
      temperature = temperature * self.cooling_rate
      iteration = iteration + 1
    self.transformer.setParameters(state, transsys_program)
    return OptimisationResult(transsys_program, obj, records)


class GradientOptimiser(AbstractOptimiser) :
  """Gradient descent optimiser.

This optimiser operates by iterating these steps until a termination
criterion is reached:

  1. Generate variants of the current state by adding or subtracting
      C{delta} to each component (unless C{eliminateFlatComponents} is
      set, see below) while keeping the remainding components
      fixed, and use these samples to estimate the local gradient.

  2. Construct a new state by adding C{stepsize} * the normalised
      gradient to the current state.

  3. If the new state improves the objective function, accept it. Then
      try whether performing a step of C{stepsize / stepsize_shrink}
      yields further improvement. If so, let
      C{stepsize = stepsize / stepsize_shrink} and iteratively continue
      doing so until no further improvement is attained or C{stepsize}
      exceeds C{stepsize_max}.

      If the new state does not improve the objective function, iteratively
      let C{stepsize = stepsize * stepsize_shrink} until improvement
      occurs or C{stepsize} falls below C{termination_stepsize}.

  4. Terminate if C{stepsize} has fallen below C{termination_stepsize} or if
      the improvement attained was below C{termination_improvement},
      otherwise start next cycle.

@ivar initial_stepsize: distance initially stepped in gradient
  direction. Subject to subsequent adaptation
@ivar delta: offset from current point in search space used
  for gradient computation (estimation)
@ivar stepsize_shrink: factor by which stepsize is multiplied / divided
  in stepsize adaptation
@ivar termination_stepsize: minimal stepsize, optimisation terminates if
  stepsize falls below this limit
@ivar termination_objective: optimisation terminates if objective function
  value falls to or below this
@ivar termination_iteration: optimisation terminates at this iteration
@ivar termination_numEvaluations: optimisation terminates if number of
  objective function evaluations reaches or exceeds this
@ivar termination_improvement: optimisation terminates if improvement
  drops below this threshold.
@ivar stepsize_max: maximal stepsize, adaptation will not make
  stepsize exceed this limit
@ivar transformer: parameter transformer
@ivar eliminateFlatComponents: whether components for which the
  partial derivative was 0 once should be excluded from subsequent
  optimisation iterations. This effectively freezes a component
  once its partial derivative becomes "flat", thus speeding up
  subsequent iterations but incurring the risk of freezing in a
  suboptimal state.
@ivar verbose: controls verbosity.
"""

  savefile_magic = 'GradientOptimiser-0.1'

  def __init__(self, initial_stepsize = 1.0, delta = 1.0e-6, stepsize_shrink = 0.5, termination_stepsize = 1.0e-30, termination_objective = None, termination_iteration = None, termination_numEvaluations = None, termination_improvement = 0.0, stepsize_max = 1.0e30, transformer = None, verbose = False) :
    """
@param initial_stepsize: distance initially stepped in gradient
  direction. Subject to subsequent adaptation
@param delta: offset from current point in search space used
  for gradient computation (estimation)
@param stepsize_shrink: factor by which stepsize is multiplied / divided
  in stepsize adaptation
@param termination_stepsize: minimal stepsize, optimisation terminates if
  stepsize falls below this limit
@param termination_objective: optimisation terminates if objective function
  value falls to or below this
@param termination_iteration: optimisation terminates at this iteration
@param termination_numEvaluations: optimisation terminates if number of
  objective function evaluations reaches or exceeds this
@param termination_improvement: optimisation terminates if improvement
  drops below this threshold.
@param stepsize_max: maximal stepsize, adaptation will not make
  stepsize exceed this limit
@param transformer: parameter transformer
"""
    super(GradientOptimiser, self).__init__(verbose)
    self.initial_stepsize = initial_stepsize
    self.delta = delta
    self.stepsize_shrink = stepsize_shrink
    self.termination_stepsize = termination_stepsize
    self.termination_objective = termination_objective
    self.termination_iteration = termination_iteration
    self.termination_numEvaluations = termination_numEvaluations
    self.termination_improvement = termination_improvement
    self.stepsize_max = stepsize_max
    if transformer is None :
      self.transformer = IdentityParameterTransformer()
    else :
      self.transformer = transformer
    self.eliminateFlatComponents = False


  def check_savefile_magic(self, s) :
    """Check whether C{s} is the "magic word" indicating the
beginning of a valid save file.
"""
    return s == self.savefile_magic


  def parse_variables(self, f) :
    """Get the values of instance variables from file C{f}."""
    self.initial_stepsize = utils.parse_float(f, 'initial_stepsize')
    self.delta = utils.parse_float(f, 'delta')
    self.stepsize_shrink = utils.parse_float(f, 'stepsize_shrink')
    self.termination_stepsize = utils.parse_float(f, 'termination_stepsize', allowNone = True)
    self.termination_objective = utils.parse_float(f, 'termination_objective', allowNone = True)
    self.termination_iteration = utils.parse_int(f, 'termination_iteration', allowNone = True)
    self.termination_numEvaluations = utils.parse_int(f, 'termination_numEvaluations', allowNone = True)
    self.termination_improvement = utils.parse_float(f, 'termination_improvement', allowNone = True)
    self.stepsize_max = utils.parse_float(f, 'stepsize_max')
    self.eliminateFlatComponents = utils.parse_boolean(f, 'eliminateFlatComponents')
    s = utils.parse_string(f, 'transformer').strip()
    if s == 'identity' :
      self.transformer = IdentityParameterTransformer()
    elif s == 'constrained' :
      self.transformer = ConstrainedParameterTransformer()
    else :
      raise StandardError, 'unknown parameter transformer "%s"' % s


  def __str__(self) :
    s = '%s\n' % self.savefile_magic
    s = s + '%s\n' % utils.name_value_pair(self.initial_stepsize, 'initial_stepsize')
    s = s + '%s\n' % utils.name_value_pair(self.delta, 'delta')
    s = s + '%s\n' % utils.name_value_pair(self.stepsize_shrink, 'stepsize_shrink')
    s = s + '%s\n' % utils.name_value_pair(self.termination_stepsize, 'termination_stepsize')
    s = s + '%s\n' % utils.name_value_pair(self.termination_objective, 'termination_objective')
    s = s + '%s\n' % utils.name_value_pair(self.termination_iteration, 'termination_iteration')
    s = s + '%s\n' % utils.name_value_pair(self.termination_numEvaluations, 'termination_numEvaluations')
    s = s + '%s\n' % utils.name_value_pair(self.termination_improvement, 'termination_improvement')
    s = s + '%s\n' % utils.name_value_pair(self.stepsize_max, 'stepsize_max')
    s = s + '%s\n' % utils.name_value_pair(self.eliminateFlatComponents, 'eliminateFlatComponents')
    if isinstance(self.transformer, IdentityParameterTransformer) :
      s = s + '%s\n' % utils.name_value_pair('identity', 'transformer')
    elif isinstance(self.transformer, ConstrainedParameterTransformer) :
      s = s + '%s\n' % utils.name_value_pair('constrained', 'transformer')
    else :
      raise StandardError, 'cannot encode parameter transformer'
    return s
    

  def terminationCondition(self, stepsize, objective, improvement, iteration, numEvaluations, gradient) :
    """Check whether termination condition is satisfied.

@return: C{True} if termination condition is satisfied
@rtype: bool
    """
    # TODO: print out satisfied condition when they are encountered
    if self.termination_stepsize is not None :
      if stepsize <= self.termination_stepsize :
        return True
    if self.termination_objective is not None :
      if objective <= self.termination_objective :
        return True
    if self.termination_improvement is not None and improvement is not None :
      if improvement <= self.termination_improvement :
        return True
    if self.termination_iteration is not None :
      if iteration >= self.termination_iteration :
        return True
    if self.termination_numEvaluations is not None :
      if numEvaluations >= self.termination_numEvaluations :
        return True
    if max(map(abs, gradient)) == 0.0 :
      return True
    

  def optimise(self, transsys_program, objective_function, gene_name_list = None, factor_name_list = None) :
    tp = copy.deepcopy(transsys_program)
    current_values = self.transformer.getParameters(tp)
    if self.verbose :
      sys.stderr.write('optimising %d values\n' % len(current_values))
    current_obj = objective_function(tp).fitness
    numEvaluations = 1
    stepsize = self.initial_stepsize
    old_gradient = [1.0 / math.sqrt(float(len(current_values)))] * len(current_values)
    iteration = 0
    optimisation_log = [GradientOptimisationRecord(current_obj, stepsize, numEvaluations, old_gradient)]
    improvement = None
    while not self.terminationCondition(stepsize, current_obj, improvement, iteration, numEvaluations, old_gradient) :
      gradient = []
      num_zeros = 0
      for i in xrange(len(current_values)) :
        if self.eliminateFlatComponents and old_gradient[i] == 0.0 :
          gradient.append(0.0)
          num_zeros = num_zeros + 1
        else :
          v = current_values[:]
          v[i] = current_values[i] - self.delta
          self.transformer.setParameters(v, tp)
          o_minus = objective_function(tp).fitness
          numEvaluations = numEvaluations + 1
          v[i] = current_values[i] + self.delta
          self.transformer.setParameters(v, tp)
          o_plus = objective_function(tp).fitness
          numEvaluations = numEvaluations + 1
          gradient.append(o_plus - o_minus)
          # restore current state
          self.transformer.setParameters(current_values, tp)
      if self.verbose and self.eliminateFlatComponents :
        sys.stderr.write('num_zeros: %d\n' % num_zeros)
      gradient_norm = utils.euclidean_norm(gradient)
      if gradient_norm == 0.0 :
        if self.verbose :
          sys.stderr.write('GradientOptimiser.optimise: flat gradient\n')
        break
      old_gradient = gradient
      gradient = map(lambda x : x / gradient_norm, gradient)
      # print gradient
      # if self.verbose :
      #   sys.stderr.write('obj = %g, values: %s\n' % (current_obj, str(current_values)))
      #   sys.stderr.write('  gradient = %s\n' % str(gradient))
      v = []
      for i in xrange(len(current_values)) :
        v.append(current_values[i] - gradient[i] * stepsize)
      self.transformer.setParameters(v, tp)
      new_obj = objective_function(tp).fitness
      numEvaluations = numEvaluations + 1
      if new_obj < current_obj :
        new_values = v[:]
        stepsize2 = stepsize / self.stepsize_shrink
        v = []
        for i in xrange(len(current_values)) :
          v.append(current_values[i] - gradient[i] * stepsize2)
        self.transformer.setParameters(v, tp)
        new_obj2 = objective_function(tp).fitness
        numEvaluations = numEvaluations + 1
        if self.verbose :
          sys.stderr.write('growth branch: tried stepsize %g, new_obj = %g, new_obj2 = %g\n' % (stepsize2, new_obj, new_obj2))
        while new_obj2 < new_obj and stepsize < self.stepsize_max :
          stepsize = stepsize2
          if self.verbose :
            sys.stderr.write('growing stepsize to %g\n' % stepsize)
          new_obj = new_obj2
          new_values = v[:]
          stepsize2 = stepsize / self.stepsize_shrink
          v = []
          for i in xrange(len(current_values)) :
            v.append(current_values[i] - gradient[i] * stepsize2)
          self.transformer.setParameters(v, tp)
          new_obj2 = objective_function(tp).fitness
          numEvaluations = numEvaluations + 1
          if self.verbose :
            sys.stderr.write('growth branch: tried stepsize %g, new_obj = %g, new_obj2 = %g\n' % (stepsize2, new_obj, new_obj2))
        if not new_obj2 < new_obj :
          self.transformer.setParameters(new_values, tp)
      else :
        while new_obj >= current_obj and stepsize > self.termination_stepsize :
          stepsize = stepsize * self.stepsize_shrink
          if self.verbose :
            sys.stderr.write('  %g: obj: current = %g, new = %g, d = %g\n' % (stepsize, current_obj, new_obj, current_obj - new_obj))
          new_values = []
          for i in xrange(len(current_values)) :
            new_values.append(current_values[i] - gradient[i] * stepsize)
          self.transformer.setParameters(new_values, tp)
          new_obj = objective_function(tp).fitness
          numEvaluations = numEvaluations + 1
      if stepsize > self.termination_stepsize :
        improvement = current_obj - new_obj
        if self.verbose :
          sys.stderr.write('  making step %g: improvement = %g, new_obj = %g\n' % (stepsize, improvement, new_obj))
        current_obj = new_obj
        current_values = new_values[:]
        optimisation_log.append(GradientOptimisationRecord(current_obj, stepsize, numEvaluations, old_gradient))
        iteration = iteration + 1
    return OptimisationResult(tp, current_obj, optimisation_log)


def parse_optimiser(f) :
  l = f.readline()
  if l == '' :
    return None
  l = l.strip()
  for o in [GradientOptimiser(), SimulatedAnnealer()] :
    if o.check_savefile_magic(l) :
      o.parse_variables(f)
      return o
  raise StandardError, 'unknown optimiser "%s"' % l


class LsysObjectiveFunction(AbstractObjectiveFunction) :

  def __init__(self, lsys, control_transsys, disparity_function, num_timesteps) :
    self.lsys = copy.deepcopy(lsys)
    self.control_transsys = copy.deepcopy(control_transsys)
    self.disparity_function = disparity_function
    self.num_timesteps = num_timesteps


  def __call__(self, tp) :
    transsys_program = copy.deepcopy(self.control_transsys)
    transsys_program.merge(tp)
    self.lsys.associate_transsys(transsys_program)
    # print self.lsys
    fitness = disparity_fitness(self.lsys, transsys_program, tp.factor_names(), self.num_timesteps, self.disparity_function)
    self.lsys.dissociate_transsys()
    return fitness


