#!/usr/bin/env python

import sys
import copy
import math
import types
import os
import random
import string
# import collections

import transsys
import transsys.optim
import transsys.utils


class Interval(object) :
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
  """Extract the set concentrations of a factor in instances that activated
a rule (the active set), and the complement set (the inactive set).

This method is rather specific to L-systems optimisation and should
perhaps not be in this module."""
  active_set = []
  inactive_set = []
  for tsi in instance_series :
    if tsi.rule is rule :
      active_set.append(tsi.transsys_instance.factor_concentration[factor_index])
    else :
      inactive_set.append(tsi.transsys_instance.factor_concentration[factor_index])
  return active_set, inactive_set


def stddev_disparity(rule, instance_series, factor_index) :
  """Standard deviation based disparity score.

This function applies some tricks for preventing limitations of
floating point precision from resulting in artificially low
scores. See comments in source for details.
"""
  active_set, inactive_set = active_and_inactive_set(rule, instance_series, factor_index)
  if len(active_set) == 0 or len(inactive_set) == 0 :
    return 0.0
  if len(active_set) > 1 :
    m_active, s_active = transsys.utils.mean_and_stddev(active_set)
  else :
    m_active, s_active = active_set[0], 0.0
  if len(inactive_set) > 1 :
    m_inactive, s_inactive = transsys.utils.mean_and_stddev(inactive_set)
  else :
    m_inactive, s_inactive = active_set[0], 0.0
  m_diff = abs(m_active - m_inactive)
  # if the difference of means falls below a threshold of 1e-150, the
  # squares calculated in the process of computing the standard deviation
  # may reach the order of magnitude of 1e-300, i.e. the limits of
  # floating point precision where they become indistinguishable from 0.
  # This condition prevents such imprecise values from being used.
  # (The standard deviations may be much larger than the difference
  # of the means and in that case not be affected by this problem, but
  # notice that the score will then be close to 1 anyway.)
  if m_diff < 1e-100 :
    sd_disp = 1.0
  else :
    # m_active and m_inactive may end up different even though all
    # values are identical, if the cardinalities of the active and
    # the inactive set are different (by orders of magnitude).
    # Therefore, we check the ratio of the mean and the mean difference
    # and treat the difference as 0 if it is too few orders of magnitude
    # above the numerical precision level
    # FIXME: currently, we use 1e-10 as a blanket threshold, this could
    # be made dependent on the ratio of set cardinalities.
    m_ratio = m_diff / (abs(m_active) + abs(m_inactive))
    if m_ratio < 1e-10 :
      sd_disp = 1.0
    else :
      sd_disp = 1.0 - 1.0 / (1.0 + (s_active + s_inactive) / m_diff)
    # print '%s, factor %d: active = %e +- %e, inactive = %e +- %e, m_ratio = %e, sd_disp = %e' % (rule.name, factor_index, m_active, s_active, m_inactive, s_inactive, m_ratio, sd_disp)
  return sd_disp


def difference_sign_disparity(rule, instance_series, factor_index) :
  """Disparity of sign values of differences between expression levels
in symbols activating C{rule} and other symbols.

Evaluates to 0 if all signs are the same (and hence activating and
inactivating instances have expression levels of the factor specified
by C{factor_index} in non-overlapping intervals and to 1 if sign ratio
is 50:50.
"""
  active_set, inactive_set = active_and_inactive_set(rule, instance_series, factor_index)
  if len(active_set) == 0 or len(inactive_set) == 0 :
    return 0.0
  n_minus = 0
  n_plus = 0
  for a in active_set :
    for i in inactive_set :
      if a < i :
        n_minus = n_minus + 1
      elif a > i :
        n_plus = n_plus + 1
  n_all = n_minus + n_plus
  if n_all == 0 :
    return 1.0
  return float(min(n_minus, n_plus)) / float(n_all) * 2.0



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


class FitnessResult(object) :
  """Base class for results of objective (or fitness) functions.

Calling objective functions (subclasses of
L{AbstractObjectiveFunction}) should return instances of this class or
one of its subclasses. Optimisers use the C{fitness} instance
variable. Subclasses may introduce further instance variables to
describe the result achieved, break it down to individual components
of the transsys program etc.

@ivar fitness: the fitness value, i.eI{.} the value computed by the
  objective function.
"""
  def __init__(self, fitness) :
    self.fitness = fitness


  def __str__(self) :
    return 'FitnessResult(%g)' % self.fitness


  def getTranssysProgramComments(self) :
    """Get a list of strings describing the objective value and
its computation.

The idea is to attach this list to the comments of the transsys
program that has been evaluated.

Notice that since the C{FitnessResult} does not keep track of the
transsys program it pertains to, so it's up to the user of this
method to make sure the comments are added to the correct transsys
program.
"""
    return ['objective: %g' % self.fitness]


class LsysDisparityFitnessResult(FitnessResult) :
  """Result of a disparity objective evaluation."""

  def __init__(self, fitness, best_factor_list, score_table) :
    FitnessResult.__init__(self, fitness)
    self.best_factor_list = best_factor_list
    self.score_table = score_table


  def getTranssysProgramComments(self) :
    comment_list = super(LsysDisparityFitnessResult, self).getTranssysProgramComments()
    comment_list.append('best factor list:')
    for rulename, best_factors in self.best_factor_list :
      s = '    %s' % rulename
      glue = ': '
      for fname in best_factors :
        s = s + '%s%s' % (glue, fname)
        glue = ', '
      comment_list.append(s)
    comment_list.append('objective by factor:')
    for st_entry in self.score_table :
      comment_list.append('    %s %s: %f' % st_entry)
    return comment_list


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
  """Objective function quantifying disparity of expression levels
between instances that activate and that do not activate rules.

This should be rewritten as a proper L{AbstractObjectiveFunction} subclass.
  """
  if len(transsys_program.factor_list) == 0 :
    return 1.0
  fitness = 0.0
  best_factor_list = []
  score_table = []
  factor_indices = []
  for factor_name in factor_names :
    factor_indices.append(transsys_program.find_factor_index(factor_name))
  instance_series = flat_symbol_instance_list(lsys_program.derivation_series(num_timesteps), transsys_program)
  for rule in lsys_program.rules :
    fmin = None
    for factor_index in factor_indices :
      f = disparity(rule, instance_series, factor_index)
      score_table.append((rule.name, transsys_program.factor_list[factor_index].name, f, ))
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
  return LsysDisparityFitnessResult(fitness / float(len(lsys_program.rules)), best_factor_list, score_table)


class Mutator(object) :

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
  """Experimental function for sequence encoding."""
  dnaseq = initial_dnaseq
  fitness = objective_function(decoder.decode_transsys('implant', dnaseq)).fitness
  for g in xrange(num_generations) :
    print '%d: %f' % (g, fitness)
    mutant_dnaseq = mutator.mutate(dnaseq)
    mutant_fitness = objective_function(decoder.decode_transsys('implant', mutant_dnaseq)).fitness
    d = transsys.utils.hamming_distance(dnaseq, mutant_dnaseq)
    # print 'current: %f, mutant: %f, distance: %d (rel: %f)' % (fitness, mutant_fitness, d, float(d) / float(len(dnaseq)))
    if mutant_fitness <= fitness :
      dnaseq = mutant_dnaseq
      fitness = mutant_fitness
  return dnaseq


# moving to utils
def get_value_nodes(transsys_program, gene_name_list, factor_name_list) :
  """Get value expression nodes from a transsys program.

@deprecated: use methods of L{TranssysProgram}.
  """
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
  """Randomise values in value expressions.

@deprecated: This functionality is now integrated into the optimising framework.
"""
  value_expression_list = get_value_nodes(transsys_program, gene_name_list, factor_name_list)
  for n in value_expression_list :
    n.value = random_function()


class TransformationFunction(object) :
  """Base class for transformation functions.

This base class implements the identity transformation.
"""

  functionName = 'TransformationFunction'


  def __call__(self, x) :
    """Compute the transformed value of C{x}."""
    return x


  def __str__(self) :
    return self.functionName


  def clip(self, y) :
    """Clip C{y} to the nearest admissible value."""
    return y


  def inverse(self, y) :
    """Compute the inverse value of C{y}.

This method should return the untransformed value C{x} corresponding to
C{y}, such that f(x) = y.
"""
    return y


  def check_savefile_magic(self, s) :
    return s.strip() == self.functionName


  def parse_variables(self, f) :
    pass


  def parse(self, f) :
    l = f.readline()
    if l == '' :
      raise StandardError, 'expected magic but got EOF'
    l = l.strip()
    if not self.check_savefile_magic(l) :
      raise StandardError, 'bad magic "%s" instead of "%s"' % (l, self.functionName)
    self.parse_variables(f)


class ExponentialFunction(TransformationFunction) :
  """Exponential parameter transformation.

This is a callable class implementing the transformation function
y = minValue + coefficient * exp(exponentMultiplier * x).
"""

  functionName = 'ExponentialFunction'

  def __init__(self, exponentMultiplier = 1.0, coefficient = 1.0, minValue = 0.0) :
    """Constructor."""
    self.exponentMultiplier = exponentMultiplier
    self.coefficient = coefficient
    self.minValue = minValue


  def __call__(self, x) :
    """Compute function value."""
    return self.minValue + self.coefficient * math.exp(self.exponentMultiplier * x)


  def __str__(self) :
    s = '%s\n' % self.functionName
    s = s + '%s\n' % transsys.utils.name_value_pair(self.exponentMultiplier, 'exponentMultiplier')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.coefficient, 'coefficient')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.minValue, 'minValue')
    return s[:-1]


  def parse_variables(self, f) :
    self.exponentMultiplier = transsys.utils.parse_float(f, 'exponentMultiplier')
    self.coefficient = transsys.utils.parse_float(f, 'coefficient')
    self.minValue = transsys.utils.parse_float(f, 'minValue')


  def clip(self, y) :
    return max(self.minValue, y)


  def inverse(self, y) :
    """Compute inverse.

Notice that the inverse diverges as y approaches minValue.
Exactly at minValue, the inverse method attempts to return
the smallest value that can be transformed without a floating
point overflow.


@return: inverse function value, or C{None} if C{y} is out of admissible range
@rtype: float or None
"""
    if self.coefficient == 0.0 :
      raise StandardError, 'exponential function with coefficient 0 cannot be inverted'
    if self.exponentMultiplier == 0.0 :
      raise StandardError, 'exponential function with exponentMultiplier 0 cannot be inverted'
    if y < self.minValue :
      return None
      # raise StandardError, 'argument %f < minimal value %f' % (y, self.minValue)
    xmin = -700.0 / self.exponentMultiplier
    v = (y - self.minValue) / self.coefficient
    # print '%1.17e %1.17e %s' % (y, v, str(v == 0.0))
    if v == 0.0 :
      # 700 is 300 * log(10), and as floating point ranges typically
      # extend to around 1e300, this should be ok on normal platforms
      return xmin
    return max(xmin, 1.0 / self.exponentMultiplier * math.log(v))


class ArctanFunction(TransformationFunction) :
  """Arcus tangens transformation function.

The arctan approaches its minimum and its maximum with O(1/x). From
a perspective of numeric computation, this makes it favourable over
the exponential sigmoid function, as x can take on values from the
greatest possible range without incurring floating point limitations.
"""

  functionName = 'ArctanFunction'

  def __init__(self, minValue = 0.0, maxValue = 1.0) :
    self.minValue = minValue
    self.valueRange = maxValue - minValue


  def __call__(self, x) :
    return self.minValue + (math.atan(x) / math.pi + 0.5) * self.valueRange


  def __str__(self) :
    s = '%s\n' % self.functionName
    s = s + '%s\n' % transsys.utils.name_value_pair(self.minValue, 'minValue')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.getMaxValue(), 'maxValue')
    return s[:-1]


  def parse_variables(self, f) :
    self.minValue = transsys.utils.parse_float(f, 'minValue')
    maxValue = transsys.utils.parse_float(f, 'maxValue')
    self.valueRange = maxValue - self.minValue


  def getMaxValue(self) :
    return self.minValue + self.valueRange


  def clip(self, y) :
    return min(self.getMaxValue(), max(self.minValue, y))


  def inverse(self, y) :
    maxValue = self.getMaxValue()
    if y < self.minValue or y > maxValue :
      return None
    if y == self.minValue :
      return -1e300
    if y == maxValue :
      return 1e300
    return math.tan(math.pi * ((y - self.minValue) / self.valueRange - 0.5))


class LowerBoundedFunction(TransformationFunction) :
  """A function asymptoting to a finite lower bound as x -> -Inf
and linearly diverging as x -> +Inf.

This implements the function 0.5 * x + sqrt(1.0 + 0.25 * x^2), which
approaches 1/x for small x and x for large x. Like arctan, this
function approaches its asymptote with O(1/x) and may therefore be
favourable over exponential functions. Linear divergence may also
be much less prone to extremes than exponential divergence.
"""

  functionName = 'LowerBoundedFunction'

  def __init__(self, minValue = 0.0) :
    self.minValue = minValue


  def __call__(self, x) :
    return self.minValue + 0.5 * x + math.sqrt(1.0 + 0.25 * x * x)


  def __str__(self) :
    s = '%s\n' % self.functionName
    s = s + '%s\n' % transsys.utils.name_value_pair(self.minValue, 'minValue')
    return s[:-1]


  def parse_variables(self, f) :
    self.minValue = transsys.utils.parse_float(f, 'minValue')


  def clip(self, y) :
    return max(self.minValue, y)


  def inverse(self, y) :
    if y < self.minValue :
      return None
    if y == self.minValue :
      return -1e300
    y1 = y - self.minValue
    return (y1 * y1 - 1.0) / y1


class SigmoidFunction(TransformationFunction) :
  """Sigmoid parameter transformation.

This is a callable class implementing the transformation function
y = minValue + (maxValue - minValue) / (1.0 + exp(-exponentMultiplier * x)).
"""

  # FIXME: this function may have a further parameter x0, a more general form
  # is minValue + valueRange / exp(-exponentMultiplier * x - x0)
  # Not clear whether this extension would be useful, could alternatively
  # use biased initialisation for x
  functionName = 'SigmoidFunction'

  def __init__(self, exponentMultiplier = 1.0, minValue = 0.0, maxValue = 1.0) :
    """Constructor."""
    self.exponentMultiplier = exponentMultiplier
    self.minValue = minValue
    self.valueRange = maxValue - minValue


  def __call__(self, x) :
    """Compute function value."""
    exp = -self.exponentMultiplier * x
    if exp < -700.0 :
      return self.getMaxValue()
    elif exp > 700.0 :
      return self.minValue
    return self.minValue + self.valueRange / (1.0 + math.exp(exp))


  def __str__(self) :
    s = '%s\n' % self.functionName
    s = s + '%s\n' % transsys.utils.name_value_pair(self.exponentMultiplier, 'exponentMultiplier')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.minValue, 'minValue')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.getMaxValue(), 'maxValue')
    return s[:-1]


  def parse_variables(self, f) :
    self.exponentMultiplier = transsys.utils.parse_float(f, 'exponentMultiplier')
    self.minValue = transsys.utils.parse_float(f, 'minValue')
    maxValue = transsys.utils.parse_float(f, 'maxValue')
    self.valueRange = maxValue - self.minValue


  def getMaxValue(self) :
    return self.minValue + self.valueRange


  def clip(self, y) :
    return min(self.getMaxValue(), max(self.minValue, y))


  def inverse(self, y) :
    """Compute inverse.

Notice that the inverse diverges towards the borders of the sigmoid
function's range. Exactly at the borders, the inverse method attempts
to return the closest value that can be transformed without a floating
point overflow.

@return: inverse function value, or C{None} if C{y} is out of admissible range
@rtype: float or None
"""
    maxValue = self.getMaxValue()
    if y < self.minValue or y > maxValue :
      return None
      # raise StandardError, 'argument %f out of range [%f, %f]' % (y, self.minValue, maxValue)
    # FIXME: proper DBL_MIN, DBL_MAX like values should be used here
    if y == self.minValue :
      return -700.0 / self.exponentMultiplier
    if y == maxValue :
      return 700.0 / self.exponentMultiplier
    return -1.0 / self.exponentMultiplier * math.log(self.valueRange / (y - self.minValue) - 1.0)


def parse_transformation_function(f) :
  l = f.readline()
  if l == '' :
    return None
  l = l.strip()
  for tf in [TransformationFunction(), ExponentialFunction(), SigmoidFunction(), LowerBoundedFunction(), ArctanFunction()] :
    if tf.check_savefile_magic(l) :
      tf.parse_variables(f)
      return tf
  raise StandardError, 'unknown transformation function "%s"' % l


class ParameterTransformer(object) :
  """Base class for parameter transformers.

Parameter transformation maps values from the domain ]-Inf, +Inf[ to
limited ranges of the real numbers, using L{TransformationFunction}s.

Furthermore, parameter transformers can restrict the set of factors
and genes for which optimisation is performed. Notionally, this is
equivalent to mapping all values from ]-Inf, +Inf[ to the constant
found in the original transsys program.

Notice that the factors and genes to be optimised must be specified
by their names, not by references.

This base class implements an identity transformer, i.e. the transformed
value of a parameter is just its actual value.

Technical note: Restriction of parameters that are subjected to
optimisation is not implemented by clamping them to a constant value but
by hiding the parameters of elements excluded from optimisation from the
optimiser altogether.  Effectively, this reduces the dimensionality of
the search space and therefore, it is much preferable from the notionally
equivalent clamping. From this perspective, the C{ParameterTransformer}
class currently combines the functions of transforming and selecting
parameters, and the optimisation framework may in the future be changed
to implement these separate function in separate classes.

This design change should not effect high level users of
C{ParameterTransformer} and its subclasses as the factor and gene
lists are normally received through the C{optimise} methods of
L{AbstractOptimiser} subclasses, and this interface will not be
affected by the design change.

@ivar transsys_program: the transsys program
@type transsys_program: C{TranssysProgram}
@ivar factor_name_list: list of names of factors for which parameters are to be
  optimised
@type factor_name_list: C{list} of C{string}s
@ivar gene_name_list: list of names of genes for which parameters are to be
  optimised
@type gene_name_list: C{list} of C{string}s
@ivar target_node_list: list of expression nodes holding the values subject to
  optimisation
@type target_node_list: C{list} of C{ExpressionNodeValue}
"""

  savefile_magic = 'ParameterTransformer'


  def __init__(self) :
    self.transsys_program = None
    self.target_node_list = None
    self.factor_name_list = None
    self.gene_name_list = None
    self.target_node_list = None


  def __str__(self) :
    return self.savefile_magic


  def check_savefile_magic(self, s) :
    return s.strip() == self.savefile_magic


  def parse_variables(self, f) :
    pass


  def parse(self, f) :
    l = f.readline()
    if l == '' :
      raise StandardError, 'expected magic but got EOF'
    l = l.strip()
    if not self.check_savefile_magic(l) :
      raise StandardError, 'bad magic "%s" instead of "%s"' % (l, self.savefile_magic)
    self.parse_variables(f)


  def setTranssysProgram(self, transsys_program, factor_name_list, gene_name_list) :
    """Set the transsys program on which the transformer is to operate.

Passing C{None} as factor name list selects all factor-related parameters for
optimisation, and likewise, passing C{None} as gene name list selects all
gene-related parameters for optimisation.

@param transsys_program: the transsys program
@type transsys_program: C{TranssysProgram}
@param factor_name_list: factors for which parameters are to be optimised
@type factor_name_list: C{list} of C{string}s, or C{None}
@param gene_name_list: genes for which parameters are to be optimised
@type gene_name_list: C{list} of C{string}s, or C{None}
"""
    self.transsys_program = transsys_program
    self.target_node_list = []
    if factor_name_list is None :
      self.target_node_list.extend(self.transsys_program.getFactorValueNodes())
    else :
      for factor_name in factor_name_list :
        factor = self.transsys_program.find_factor(factor_name)
        self.target_node_list.extend(factor.getValueNodes())
    if gene_name_list is None :
      self.target_node_list.extend(self.transsys_program.getGeneValueNodes())
    else :
      for gene_name in gene_name_list :
        gene = self.transsys_program.find_gene(gene_name)
        self.target_node_list.extend(gene.getValueNodes())


  def clipParameters(self) :
    """Clip the parameters in C{transsys_program} to the range
required by the respective transformer."""
    if self.transsys_program is None :
      raise StandardError, 'no transsys program, cannot clip'


  def randomiseParametersUniform(self, parameterRange, rng) :
    """Randomise the numeric values in C{transsys_program}.

Random values are drawn from a uniform distribution over
[-C{parameterRange}, C{parameterRange}[ and the numeric values
are set to the corresponding transformed values.
"""
    state = self.getParameters()
    for i in xrange(len(state)) :
      state[i] = (rng.random() * 2.0 - 1.0) * parameterRange
    self.setParameters(state)


  def getParameters(self) :
    """Get numeric parameters from the C{transsys_program}, converted to
the unconstrained optimiser space.
"""
    if self.transsys_program is None :
      raise StandardError, 'no transsys program, cannot get parameters'
    return map(lambda n : n.value, self.target_node_list)


  def setParameters(self, parameter_list) :
    """Set parameters in C{transsys_program} according to the optimiser space
values provided by C{parameter_list}.

This base class implementation does not do any transformation, so
it should be overridden by subclasses.
"""
    if self.transsys_program is None :
      raise StandardError, 'no transsys program, cannot set parameters'
    if len(self.target_node_list) != len(parameter_list) :
      raise StandardError,  'parameter list incompatible with transsys program'
    for i in xrange(len(parameter_list)) :
      self.target_node_list[i].value = parameter_list[i]


class TranssysTypedParameterTransformer(ParameterTransformer) :
  """
Parameter transformer differentiating decay, diffusibility, synthesis, constitutive,
aspec, amax, rspec and rmax nodes.

This class is intended for transsys programs with constant values
(i.e. L{ExpressionNodeValue} instances) for decay, diffusibility,
synthesis, constitutive, aspec, amax, rspec and rmax. Use with transsys
programs that do not have this standard, basic structure has undefined
effects and may be disabled in the future.

@ivar decayTransformation: transformation function for decay values.
@type decayTransformation: L{TransformationFunction}
@ivar diffusibilityTransformation: transformation function for diffusibility values.
@type diffusibilityTransformation: L{TransformationFunction}
@ivar synthesisTransformation: transformation function for synthesis values.
@type synthesisTransformation: L{TransformationFunction}
@ivar constitutiveTransformation: transformation function for constitutive values.
@type constitutiveTransformation: L{TransformationFunction}
@ivar aspecTransformation: transformation function for aspec values.
@type aspecTransformation: L{TransformationFunction}
@ivar amaxTransformation: transformation function for amax values.
@type amaxTransformation: L{TransformationFunction}
@ivar rspecTransformation: transformation function for rspec values.
@type rspecTransformation: L{TransformationFunction}
@ivar rmaxTransformation: transformation function for rmax values.
@type rmaxTransformation: L{TransformationFunction}
"""

  # version history:
  # TranssysTypedParameterTransformer: initial version (aka 1.0)
  # TranssysTypedParameterTransformer-1.1: added synthesistransformation
  savefile_magic = 'TranssysTypedParameterTransformer-1.1'

  def __init__(self, decayTransformation = None, diffusibilityTransformation = None, synthesisTransformation = None, constitutiveTransformation = None, aspecTransformation = None, amaxTransformation = None, rspecTransformation = None, rmaxTransformation = None) :
    super(TranssysTypedParameterTransformer, self).__init__()
    if decayTransformation is None :
      self.decayTransformation = TransformationFunction()
    else :
      self.decayTransformation = decayTransformation
    if diffusibilityTransformation is None :
      self.diffusibilityTransformation = TransformationFunction()
    else :
      self.diffusibilityTransformation = diffusibilityTransformation
    if synthesisTransformation is None :
      self.synthesisTransformation = TransformationFunction()
    else :
      self.synthesisTransformation = synthesisTransformation
    if constitutiveTransformation is None :
      self.constitutiveTransformation = TransformationFunction()
    else :
      self.constitutiveTransformation = constitutiveTransformation
    if aspecTransformation is None :
      self.aspecTransformation = TransformationFunction()
    else :
      self.aspecTransformation = aspecTransformation
    if amaxTransformation is None :
      self.amaxTransformation = TransformationFunction()
    else :
      self.amaxTransformation = amaxTransformation
    if rspecTransformation is None :
      self.rspecTransformation = TransformationFunction()
    else :
      self.rspecTransformation = rspecTransformation
    if rmaxTransformation is None :
      self.rmaxTransformation = TransformationFunction()
    else :
      self.rmaxTransformation = rmaxTransformation
    self.decay_nodes = None
    self.diffusibility_nodes = None
    self.synthesis_nodes = None
    self.constitutive_nodes = None
    self.aspec_nodes = None
    self.amax_nodes = None
    self.rspec_nodes = None
    self.rmax_nodes = None


  def __str__(self) :
    s = '%s\n' % self.savefile_magic
    s = s + 'decayTransformation\n'
    s = s + '%s\n' % str(self.decayTransformation)
    s = s + 'diffusibilityTransformation\n'
    s = s + '%s\n' % str(self.diffusibilityTransformation)
    s = s + 'synthesisTransformation\n'
    s = s + '%s\n' % str(self.synthesisTransformation)
    s = s + 'constitutiveTransformation\n'
    s = s + '%s\n' % str(self.constitutiveTransformation)
    s = s + 'aspecTransformation\n'
    s = s + '%s\n' % str(self.aspecTransformation)
    s = s + 'amaxTransformation\n'
    s = s + '%s\n' % str(self.amaxTransformation)
    s = s + 'rspecTransformation\n'
    s = s + '%s\n' % str(self.rspecTransformation)
    s = s + 'rmaxTransformation\n'
    s = s + '%s\n' % str(self.rmaxTransformation)
    return s[:-1]


  def setAllTransformations(self, transformation) :
    """Set all transformations to C{transformation}.

@param transformation: the transformation function to use
@type transformation: L{TransformationFunction}
"""
    self.decayTransformation = transformation
    self.diffusibilityTransformation = transformation
    self.synthesisTransformation = transformation
    self.constitutiveTransformation = transformation
    self.aspecTransformation = transformation
    self.amaxTransformation = transformation
    self.rspecTransformation = transformation
    self.rmaxTransformation = transformation
    self.transsys_program = None


  def setConstrainedTransformations(self) :
    """Configure transformer to be equivalent to the abandoned C{ConstrainedTransformer}.
@deprecated: assemble transformers manually, by assigning to instance variables.
"""
    self.decayTransformation = SigmoidFunction()
    self.diffusibilityTransformation = SigmoidFunction()
    self.synthesisTransformation = SigmoidFunction()
    self.constitutiveTransformation = ExponentialFunction()
    self.aspecTransformation = ExponentialFunction()
    self.amaxTransformation = ExponentialFunction()
    self.rspecTransformation = ExponentialFunction()
    self.rmaxTransformation = ExponentialFunction()
    self.transsys_program = None


  def parse_variables(self, f) :
    """Parse instance variables from file C{f}.

@param f: the file to read from
@type f: a C{file} object (file object emulators should generally work too)
"""
    l = f.readline().strip()
    if l != 'decayTransformation' :
      raise StandardError, 'expected "decauTransformation", got "%s"' % l
    self.decayTransformation = parse_transformation_function(f)
    l = f.readline().strip()
    if l != 'diffusibilityTransformation' :
      raise StandardError, 'expected "diffusibilityTransformation", got "%s"' % l
    self.diffusibilityTransformation = parse_transformation_function(f)
    l = f.readline().strip()
    if l != 'synthesisTransformation' :
      raise StandardError, 'expected "synthesisTransformation", got "%s"' % l
    self.synthesisTransformation = parse_transformation_function(f)
    l = f.readline().strip()
    if l != 'constitutiveTransformation' :
      raise StandardError, 'expected "constitutiveTransformation", got "%s"' % l
    self.constitutiveTransformation = parse_transformation_function(f)
    l = f.readline().strip()
    if l != 'aspecTransformation' :
      raise StandardError, 'expected "aspecTransformation", got "%s"' % l
    self.aspecTransformation = parse_transformation_function(f)
    l = f.readline().strip()
    if l != 'amaxTransformation' :
      raise StandardError, 'expected "amaxTransformation", got "%s"' % l
    self.amaxTransformation = parse_transformation_function(f)
    l = f.readline().strip()
    if l != 'rspecTransformation' :
      raise StandardError, 'expected "rspecTransformation", got "%s"' % l
    self.rspecTransformation = parse_transformation_function(f)
    l = f.readline().strip()
    if l != 'rmaxTransformation' :
      raise StandardError, 'expected "rmaxTransformation", got "%s"' % l
    self.rmaxTransformation = parse_transformation_function(f)
    self.transsys_program = None


  def setTranssysProgram(self, transsys_program, factor_name_list, gene_name_list) :
    self.transsys_program = transsys_program
    if factor_name_list is None :
      self.decay_nodes = transsys_program.getDecayValueNodes()
      self.diffusibility_nodes = transsys_program.getDiffusibilityValueNodes()
      self.synthesis_nodes = transsys_program.getSynthesisValueNodes()
    else :
      self.decay_nodes = []
      self.diffusibility_nodes = []
      self.synthesis_nodes = []
      for factor_name in factor_name_list :
        factor = self.transsys_program.find_factor(factor_name)
        self.decay_nodes.extend(factor.getDecayValueNodes())
        self.diffusibility_nodes.extend(factor.getDiffusibilityValueNodes())
        self.synthesis_nodes.extend(factor.getSynthesisValueNodes())
    if gene_name_list is None :
      self.constitutive_nodes = transsys_program.getConstitutiveValueNodes()
      self.aspec_nodes = transsys_program.getActivateSpecValueNodes()
      self.amax_nodes = transsys_program.getActivateMaxValueNodes()
      self.rspec_nodes = transsys_program.getRepressSpecValueNodes()
      self.rmax_nodes = transsys_program.getRepressMaxValueNodes()
    else :
      self.constitutive_nodes = []
      self.aspec_nodes = []
      self.amax_nodes = []
      self.rspec_nodes = []
      self.rmax_nodes = []
      for gene_name in gene_name_list :
        gene = self.transsys_program.find_gene(gene_name)
        self.constitutive_nodes.extend(gene.getConstitutiveValueNodes())
        self.aspec_nodes.extend(gene.getActivateSpecValueNodes())
        self.amax_nodes.extend(gene.getActivateMaxValueNodes())
        self.rspec_nodes.extend(gene.getRepressSpecValueNodes())
        self.rmax_nodes.extend(gene.getRepressMaxValueNodes())


  def clipParameters(self) :
    if self.transsys_program is None :
      raise StandardError, 'no transsys program, cannot clip'
    for n in self.decay_nodes :
      n.value = self.decayTransformation.clip(n.value)
    for n in self.diffusibility_nodes :
      n.value = self.diffusibilityTransformation.clip(n.value)
    for n in self.synthesis_nodes :
      n.value = self.synthesisTransformation.clip(n.value)
    for n in self.constitutive_nodes :
      n.value = self.constitutiveTransformation.clip(n.value)
    for n in self.aspec_nodes :
      n.value = self.aspecTransformation.clip(n.value)
    for n in self.amax_nodes :
      n.value = self.amaxTransformation.clip(n.value)
    for n in self.rspec_nodes :
      n.value = self.rspecTransformation.clip(n.value)
    for n in self.rmax_nodes :
      n.value = self.rmaxTransformation.clip(n.value)


  def getParameters(self) :
    """Get unconstrained parameters from a transsys program.

@return: a list of the unconstrained parameters
    """
    parameter_list = []
    for p in map(lambda n : n.value, self.decay_nodes) :
      parameter_list.append(self.decayTransformation.inverse(p))
    for p in map(lambda n : n.value, self.diffusibility_nodes) :
      parameter_list.append(self.diffusibilityTransformation.inverse(p))
    for p in map(lambda n : n.value, self.synthesis_nodes) :
      parameter_list.append(self.synthesisTransformation.inverse(p))
    for p in map(lambda n : n.value, self.constitutive_nodes) :
      parameter_list.append(self.constitutiveTransformation.inverse(p))
    for p in map(lambda n : n.value, self.aspec_nodes) :
      parameter_list.append(self.aspecTransformation.inverse(p))
    for p in map(lambda n : n.value, self.amax_nodes) :
      parameter_list.append(self.amaxTransformation.inverse(p))
    for p in map(lambda n : n.value, self.rspec_nodes) :
      parameter_list.append(self.rspecTransformation.inverse(p))
    for p in map(lambda n : n.value, self.rmax_nodes) :
      parameter_list.append(self.rmaxTransformation.inverse(p))
    return parameter_list


  def setParameters(self, parameter_list) :
    """Set values in a transsys program based on a list of unconstrained parameters.

Notice that the parameter list must have the same length as that obtained
by C{getParameters}.

@raise StandardError: if the parameter list is too long or too short
"""
    if len(parameter_list) != len(self.decay_nodes) + len(self.diffusibility_nodes) + len(self.synthesis_nodes) + len(self.constitutive_nodes) + len(self.aspec_nodes) + len(self.amax_nodes) + len(self.rspec_nodes) + len(self.rmax_nodes):
      raise StandardError, 'parameter list incompatible with transsys program'
    i = 0
    for n in self.decay_nodes :
      n.value = self.decayTransformation(parameter_list[i])
      i = i + 1
    for n in self.diffusibility_nodes :
      n.value = self.diffusibilityTransformation(parameter_list[i])
      i = i + 1
    for n in self.synthesis_nodes :
      n.value = self.synthesisTransformation(parameter_list[i])
      i = i + 1
    for n in self.constitutive_nodes :
      n.value = self.constitutiveTransformation(parameter_list[i])
      i = i + 1
    for n in self.aspec_nodes :
      n.value = self.aspecTransformation(parameter_list[i])
      i = i + 1
    for n in self.amax_nodes :
      n.value = self.amaxTransformation(parameter_list[i])
      i = i + 1
    for n in self.rspec_nodes :
      n.value = self.rspecTransformation(parameter_list[i])
      i = i + 1
    for n in self.rmax_nodes :
      n.value = self.rmaxTransformation(parameter_list[i])
      i = i + 1


def parse_parameter_transformer(f) :
  """Read a parameter transformer from a file.

This function uses the "magic" first line to determine what subtype
of L{ParameterTransformer} is specified.

@return: the parameter transformer
@rtype: L{ParameterTransformer}
"""
  l = f.readline()
  if l == '' :
    return None
  l = l.strip()
  for tf in [ParameterTransformer(), TranssysTypedParameterTransformer()] :
    if tf.check_savefile_magic(l) :
      tf.parse_variables(f)
      return tf
  raise StandardError, 'unknown parameter transformer "%s"' % l


class OptimisationRecord(object) :
  """Base class for optimisation records.

@cvar table_header: string providing a suitable header for a table
  containing instance representations generated by the C{__str__} method
"""

  table_header = ''


  def __init__(self) :
    pass


  def __str__(self) :
    return ''


class SimulatedAnnealingRecord(OptimisationRecord) :
  """Record containing data from a Simulated Annealing update.

@cvar table_header: string providing a suitable header for a table
  containing instance representations generated by the C{__str__} method
@ivar obj: current objective function value
@ivar obj_alt: objective function value of alternative solution
@ivar temperature: temperature
@ivar stepsize: stepsize
@ivar stepvector_entropy: entropy of entries of the step vector
@ivar stepvector_max: largest step vector component
@ivar state_distance: distance between current and new state
@ivar p_accept: accaptance probability
@ivar accept: boolean, whether alternative solution was accepted
"""

  table_header = 'obj obj_alt temperature stepsize stepvector_entropy stepvector_max state_distance p_accept accept'


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


  # FIXME: using __str__ to generate R table rows is bending the rules a bit...
  def __str__(self) :
    return transsys.utils.table_row([self.obj, self.obj_alt, self.temperature, self.stepsize, self.stepvector_entropy, self.stepvector_max, self.state_distance, self.p_accept, self.accept])


class GradientOptimisationRecord(OptimisationRecord) :
  """Record containing data from a step of gradient optimisation.

@cvar table_header: string providing a suitable header for a table
  containing instance representations generated by the C{__str__} method
@ivar obj: current objective function value
@ivar numEvaluations: number of evaluations (cumulative)
@ivar gradient_entropy: entropy of the gradient vector
@ivar gradient_max: maximum of the gradient vector
"""


  table_header = 'obj stepsize numEvaluations gradient_entropy gradient_max'


  def __init__(self, obj, stepsize, numEvaluations, gradient) :
    """Constructor.

@param obj: current objective function value
@param numEvaluations: cumulative number of evaluations
@param gradient: gradient vector as a list, normalised
"""
    self.obj = obj
    self.stepsize = stepsize
    self.numEvaluations = numEvaluations
    gradient_abs = map(abs, gradient)
    # print transsys.utils.euclidean_norm(gradient)
    # print gradient_abs
    self.gradient_entropy = transsys.utils.shannon_entropy(gradient_abs)
    self.gradient_max = max(gradient_abs)


  def __str__(self) :
    return transsys.utils.table_row([self.obj, self.stepsize, self.numEvaluations, self.gradient_entropy, self.gradient_max])


class OptimisationResult(object) :
  """Base class for returning results of an optimisation.

An optimisation result consists of the optimised version of the
transsys program, its objective function value, and, optionally,
an optimisation log.

The optimisation log is a list of optimisation records with a
structure that depends on the optimiser. Optimisation records
can be printed on a line, and they provide a table header. This
allows the optimisation log to be dumped in a format ready for
R's C{read.table} function.

@cvar recordClass: class of the records in the log list
@ivar optimised_transsys_program: The transsys program resulting
  from the optimisation process.
@ivar objectiveOptimum: The result of evaluating the optimised transsys
  program using the objective function.
@ivar optimisation_log: A trace of the optimisation process, provided
  as a list of tuples.
"""

  recordClass = OptimisationRecord


  def __init__(self, tp, objectiveOptimum, optimisation_log = None) :
    self.optimised_transsys_program = tp
    self.objectiveOptimum = objectiveOptimum
    self.optimisation_log = optimisation_log


  def write_log(self, f, column_prefix = '', write_header = True, header_prefix = '') :
    if self.optimisation_log is None :
      return
    if len(self.optimisation_log) > 0 :
      if write_header :
        if header_prefix != '' :
          if not header_prefix[-1].isspace() :
            header_prefix = header_prefix + ' '
        f.write('%s%s\n' % (header_prefix, self.recordClass.table_header))
      if column_prefix != '' :
        if not column_prefix[-1].isspace() :
          column_prefix = column_prefix + ' '
      for l in self.optimisation_log :
        f.write('%s%s\n' % (column_prefix, str(l)))


class SimulatedAnnealingResult(OptimisationResult) :
  """Result of optimisation by Simulated Annealing."""


  recordClass = SimulatedAnnealingRecord


  def __init__(self, tp, objectiveOptimum, optimisation_log = None) :
    # print 'annealing', objectiveOptimum
    super(SimulatedAnnealingResult, self).__init__(tp, objectiveOptimum, optimisation_log)


class GradientOptimisationResult(OptimisationResult) :
  """Result of optimisation by Simulated Annealing."""


  recordClass = GradientOptimisationRecord


  def __init__(self, tp, objectiveOptimum, optimisation_log = None) :
    # print 'gradient', objectiveOptimum
    super(GradientOptimisationResult, self).__init__(tp, objectiveOptimum, optimisation_log)


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
        if transsys.utils.is_nan(x) :
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
          if transsys.utils.is_nan(x) :
            print str(tseries[i])
            raise StandardError, 'expression level of factor "%s" is NaN' % factor_name
          v.append(x)
        if min(v) < max(v) :
          cc = 1.0 - transsys.utils.uncentered_correlation(self.series[factor_name], v)
          cc_sum = cc_sum + cc
        else :
          cc_sum = cc_sum + 2.0
    return FitnessResult(cc_sum)


class AbstractOptimiser(object) :
  """Abstract base class for optimisers.

@ivar rng: random number generator, used to initialise optimisation, but
  may subsequently be used by stochastic optimisers.
@ivar transformer: parameter transformer
@ivar randomInitRange: if not C{None}, parameters in transformed space
  will be initialised with uniform random values from
  [-C{randomInitRange}, +C{randomInitRange}[
@ivar verbose: controls verbosity (of optimisation process)
  """

  def __init__(self, rng, transformer = None, randomInitRange = None, verbose = 0) :
    """Partial constructor method.

This constructor is to be invoked from subclass constructors. Instances
of C{AbstractOptimiser} do not have any use.

If the C{rng} parameter is C{None} (the default), the constructor initialises
the C{rng} instance variable to C{random.Random(1)}.

@param rng: random number generator (see C{rng} instance variable)
@type rng: C{random.Random}
@param transformer: parameter transformer (see C{transformer} instance variable)
@type transformer: C{ParameterTransformer}
@param randomInitRange: range for initialising parameters (see C{randomInitRange} instance variable)
@type randomInitRange: C{float}
@param verbose: verbosity level
@type verbose: C{int}
"""
    if rng is None :
      self.rng = random.Random(1)
    else :
      self.rng = rng
    if transformer is None :
      self.transformer = ParameterTransformer()
    else :
      self.transformer = transformer
    self.randomInitRange = randomInitRange
    self.verbose = verbose


  def setRandomNumberGenerator(self, rng) :
    """Set the random number generator.

@param rng: the random number generator to be used.
@type rng: C{random.Random}
"""
    self.rng = rng


  def setRandomNumberSeed(self, rndseed) :
    """Set a new random number generator with the specified seed.

This is a convenience method equivalent to calling
C{setRandomNumberGenerator(random.Random(1))}.

@param rndseed: the random seed
@type rndseed: C{int}
"""
    self.rng = random.Random(rndseed)


  def initialiseParameterTransformer(self, transsys_program, factor_name_list, gene_name_list) :
    """Initialise the parameter transformer and randomise values if requested.

This method should only be called from the optimise method, the
startup of each optimise method should contain something like
C{self.initialiseParameterTransformer(transsys_program,
factor_name_list, gene_name_list)}
"""
    self.transformer.setTranssysProgram(transsys_program, factor_name_list, gene_name_list)
    if self.randomInitRange is None :
      self.transformer.clipParameters()
    else :
      self.transformer.randomiseParametersUniform(self.randomInitRange, self.rng)


  def optimise(self, transsys_program_init, objective_function, factor_name_list = None, gene_name_list = None) :
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
@type transsys_program_init: C{TranssysProgram}
@param objective_function: The objective function used to optimise the
  transsys program.
@type objective_function: L{AbstractObjectiveFunction} (or subclass)
@param factor_name_list: factors for which parameters are to be optimised
@type factor_name_list: C{list} of C{string}s, or C{None}
@param gene_name_list: genes for which parameters are to be optimised
@type gene_name_list: C{list} of C{string}s, or C{None}
"""
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
@cvar PERTURBATION_METHOD_GAUSS: specifies local search by perturbation of
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
@ivar perturbation_method: perturbation method for local search
"""

  savefile_magic = 'SimulatedAnnealer-0.1'
  PERTURBATION_METHOD_UNIFORM = 1
  PERTURBATION_METHOD_GAUSS = 2

  def __init__(self, cooling_rate = 0.995, temperature_init = None, termination_temperature = None, termination_objective = None, termination_iteration = None, stepsize_learning_rate = 0.0, target_improvement_ratio = 0.2, stepsize_init = 1.0, stepvector_learning_rate = 0.0, stepvector_learning_decay = 0.0, perturbation_method = PERTURBATION_METHOD_UNIFORM, rng = None, transformer = None, randomInitRange = None, verbose = 0) :
    """Constructor."""
    super(SimulatedAnnealer, self).__init__(rng, transformer, randomInitRange, verbose)
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
    self.stepsize_learning_rate = transsys.utils.parse_float(f, 'stepsize_learning_rate')
    self.target_improvement_ratio = transsys.utils.parse_float(f, 'target_improvement_ratio')
    self.stepsize_init = transsys.utils.parse_float(f, 'stepsize_init')
    self.stepvector_learning_rate = transsys.utils.parse_float(f, 'stepvector_learning_rate')
    self.stepvector_learning_decay = transsys.utils.parse_float(f, 'stepvector_learning_decay')
    self.cooling_rate = transsys.utils.parse_float(f, 'cooling_rate')
    self.temperature_init = transsys.utils.parse_float(f, 'temperature_init', allowNone = True)
    self.termination_temperature = transsys.utils.parse_float(f, 'termination_temperature', allowNone = True)
    self.termination_objective = transsys.utils.parse_float(f, 'termination_objective', allowNone = True)
    self.termination_iteration = transsys.utils.parse_int(f, 'termination_iteration', allowNone = True)
    self.setPerturbationMethod(transsys.utils.parse_string(f, 'perturbation_method').strip())
    rndseed = transsys.utils.parse_int(f, 'rndseed')
    self.rng = random.Random(rndseed)
    self.transformer = parse_parameter_transformer(f)


  def __str__(self) :
    s = '%s\n' % self.savefile_magic
    s = s + '%s\n' % transsys.utils.name_value_pair(self.stepsize_learning_rate, 'stepsize_learning_rate')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.target_improvement_ratio, 'target_improvement_ratio')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.stepsize_init, 'stepsize_init')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.stepvector_learning_rate, 'stepvector_learning_rate')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.stepvector_learning_decay, 'stepvector_learning_decay')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.cooling_rate, 'cooling_rate')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.temperature_init, 'temperature_init')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.termination_temperature, 'termination_temperature')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.termination_objective, 'termination_objective')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.termination_iteration, 'termination_iteration')
    if self.perturbation_method == self.PERTURBATION_METHOD_UNIFORM :
      s = s + '%s\n' % transsys.utils.name_value_pair('uniform', 'perturbation_method')
    elif self.perturbation_method == self.PERTURBATION_METHOD_GAUSS :
      s = s + '%s\n' % transsys.utils.name_value_pair('gauss', 'perturbation_method')
    else :
      raise StandardError, 'cannot encode perturbation method'
    s = s + '%s\n' % transsys.utils.name_value_pair('unknown', 'rndseed')
    s = s + '%s\n' % str(self.transformer)
    return s[:-1]


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


  def optimise(self, transsys_program_init, objective_function, factor_name_list = None, gene_name_list = None, stepvector_init = None) :
    """Simulated annealing optimise method."""
    self.verifyTermination()
    transsys_program = copy.deepcopy(transsys_program_init)
    self.initialiseParameterTransformer(transsys_program, factor_name_list, gene_name_list)
    # slightly silly to set the parameters in initialising and then getting them next thing
    # may help to maintain consistency, though.
    state = self.transformer.getParameters()
    if self.verbose :
      sys.stderr.write('SimulatedAnnealer: optimising %d values\n' % len(state))
      # sys.stderr.write('state: %s\n' % str(state))
      sys.stderr.write('%s\n' % str(transsys_program))
    num_dimensions = len(state)
    unbiased_stepvector_component = 1.0 / math.sqrt(float(num_dimensions))
    if stepvector_init is None :
      stepvector = [unbiased_stepvector_component] * num_dimensions
    else :
      if len(stepvector_init) != num_dimensions :
        raise StandardError, 'stepvector incompatible with transsys program'
      stepvector = transsys.utils.normalised_vector(stepvector_init)
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
      state_distance = transsys.utils.euclidean_norm(delta_state)
      self.transformer.setParameters(state_alt)
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
      records.append(SimulatedAnnealingRecord(obj, obj_alt, temperature, stepsize, transsys.utils.shannon_entropy(stepvector), max(stepvector), state_distance, p_accept, accept))
      if obj_alt <= obj :
        stepsize = stepsize * (1.0 + self.stepsize_learning_rate / self.target_improvement_ratio)
        sdn = transsys.utils.normalised_vector(delta_state)
        v = []
        for i in xrange(num_dimensions) :
          v.append((1.0 - self.stepvector_learning_rate) * stepvector[i] + self.stepvector_learning_rate * abs(sdn[i]))
      else :
        stepsize = stepsize * (1.0 - self.stepsize_learning_rate)
        v = map(lambda x : (1.0 - self.stepvector_learning_decay) * x + self.stepvector_learning_decay * unbiased_stepvector_component, stepvector)
      stepvector = transsys.utils.normalised_vector(v)
      if accept :
        if self.verbose :
          sys.stderr.write('  state_alt accepted (obj = %f, obj_alt = %f)\n' % (obj, obj_alt))
          # sys.stderr.write('  new stepvector: %s\n' % str(stepvector))
        state = state_alt
        obj = obj_alt
      temperature = temperature * self.cooling_rate
      iteration = iteration + 1
    self.transformer.setParameters(state)
    # FIXME: this one last evaluation of the objective function can be saved
    return SimulatedAnnealingResult(transsys_program, objective_function(transsys_program), records)


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

  4. Terminate if C{stepsize} has fallen below
      C{termination_stepsize}, if the improvement attained was below
      C{termination_improvement} or the relative improvement fell
      below the C{termination_relative_improvement} threshold, or if
      the gradient is flat (i.e. all components are 0), otherwise
      start next cycle.

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
@ivar termination_relative_improvement: optimisation terminates if
  relative improvement drops below this threshold.
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

  savefile_magic = 'GradientOptimiser-0.1.1'

  def __init__(self, initial_stepsize = 1.0, delta = 1.0e-6, stepsize_shrink = 0.5, termination_stepsize = 1.0e-30, termination_objective = None, termination_iteration = None, termination_numEvaluations = None, termination_improvement = None, termination_relative_improvement = None, stepsize_max = 1.0e30, rng = None, transformer = None, randomInitRange = None, verbose = 0) :
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
@param termination_relative_improvement: optimisation terminates if relative
  improvement drops below this threshold.
@param stepsize_max: maximal stepsize, adaptation will not make
  stepsize exceed this limit
@param transformer: parameter transformer
"""
    super(GradientOptimiser, self).__init__(rng, transformer, randomInitRange, verbose)
    self.initial_stepsize = initial_stepsize
    self.delta = delta
    self.stepsize_shrink = stepsize_shrink
    self.termination_stepsize = termination_stepsize
    self.termination_objective = termination_objective
    self.termination_iteration = termination_iteration
    self.termination_numEvaluations = termination_numEvaluations
    self.termination_improvement = termination_improvement
    self.termination_relative_improvement = termination_relative_improvement
    self.stepsize_max = stepsize_max
    self.eliminateFlatComponents = False


  def check_savefile_magic(self, s) :
    """Check whether C{s} is the "magic word" indicating the
beginning of a valid save file.
"""
    return s == self.savefile_magic


  def parse_variables(self, f) :
    """Get the values of instance variables from file C{f}."""
    self.initial_stepsize = transsys.utils.parse_float(f, 'initial_stepsize')
    self.delta = transsys.utils.parse_float(f, 'delta')
    self.stepsize_shrink = transsys.utils.parse_float(f, 'stepsize_shrink')
    self.termination_stepsize = transsys.utils.parse_float(f, 'termination_stepsize', allowNone = True)
    self.termination_objective = transsys.utils.parse_float(f, 'termination_objective', allowNone = True)
    self.termination_iteration = transsys.utils.parse_int(f, 'termination_iteration', allowNone = True)
    self.termination_numEvaluations = transsys.utils.parse_int(f, 'termination_numEvaluations', allowNone = True)
    self.termination_improvement = transsys.utils.parse_float(f, 'termination_improvement', allowNone = True)
    self.termination_relative_improvement = transsys.utils.parse_float(f, 'termination_relative_improvement', allowNone = True)
    self.stepsize_max = transsys.utils.parse_float(f, 'stepsize_max')
    self.eliminateFlatComponents = transsys.utils.parse_boolean(f, 'eliminateFlatComponents')
    self.transformer = parse_parameter_transformer(f)


  def __str__(self) :
    s = '%s\n' % self.savefile_magic
    s = s + '%s\n' % transsys.utils.name_value_pair(self.initial_stepsize, 'initial_stepsize')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.delta, 'delta')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.stepsize_shrink, 'stepsize_shrink')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.termination_stepsize, 'termination_stepsize')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.termination_objective, 'termination_objective')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.termination_iteration, 'termination_iteration')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.termination_numEvaluations, 'termination_numEvaluations')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.termination_improvement, 'termination_improvement')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.termination_relative_improvement, 'termination_relative_improvement')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.stepsize_max, 'stepsize_max')
    s = s + '%s\n' % transsys.utils.name_value_pair(self.eliminateFlatComponents, 'eliminateFlatComponents')
    s = s + '%s\n' % str(self.transformer)
    return s[:-1]


  def terminationCondition(self, stepsize, objective, improvement, iteration, numEvaluations, gradient) :
    """Check whether termination condition is satisfied.

@return: C{True} if termination condition is satisfied
@rtype: bool
    """
    # TODO: print out satisfied condition when they are encountered
    if self.termination_stepsize is not None :
      if stepsize <= self.termination_stepsize :
        if self.verbose :
          sys.stderr.write('GradientOptimiser: terminating because stepsize %f <= %f\n' % (stepsize, self.termination_stepsize))
        return True
    else :
      if stepsize == 0.0 :
        if self.verbose :
          sys.stderr.write('GradientOptimiser: terminating because stepsize == 0.0 (implicit termination criterion)\n')
        return True
    if self.termination_objective is not None :
      if objective <= self.termination_objective :
        if self.verbose :
          sys.stderr.write('GradientOptimiser: terminating because objective %f <= %f\n' % (objective, self.termination_objective))
        return True
    if self.termination_improvement is not None and improvement is not None :
      if improvement <= self.termination_improvement :
        if self.verbose :
          sys.stderr.write('GradientOptimiser: terminating because improvement %f <= %f\n' % (improvement, self.termination_improvement))
        return True
    if self.termination_relative_improvement is not None and improvement is not None :
      # FIXME: just a patch, not a solution, relative improvement will be numerically unstable around objective = 0.0 anyway
      if objective == 0.0 :
        raise StandardError, 'relative improvement diverges because objective is 0.0 (consider setting termination_objective)'
      relative_improvement = improvement / objective
      if relative_improvement <= self.termination_relative_improvement :
        if self.verbose :
          sys.stderr.write('GradientOptimiser: terminating because relative improvement %f <= %f\n' % (relative_improvement, self.termination_relative_improvement))
        return True
    if self.termination_iteration is not None :
      if iteration >= self.termination_iteration :
        if self.verbose :
          sys.stderr.write('GradientOptimiser: terminating because iteration %d <= %d\n' % (iteration, self.termination_iteration))
        return True
    if self.termination_numEvaluations is not None :
      if numEvaluations >= self.termination_numEvaluations :
        if self.verbose :
          sys.stderr.write('GradientOptimiser: terminating because numEvaluations %d <= %d\n' % (numEvaluations, self.termination_numEvaluations))
        return True
    if max(map(abs, gradient)) == 0.0 :
      return True


  def optimise(self, transsys_program, objective_function, factor_name_list = None, gene_name_list = None) :
    tp = copy.deepcopy(transsys_program)
    self.initialiseParameterTransformer(tp, factor_name_list, gene_name_list)
    current_values = self.transformer.getParameters()
    if self.verbose :
      sys.stderr.write('optimising %d values\n' % len(current_values))
    current_obj = objective_function(tp).fitness
    numEvaluations = 1
    stepsize = self.initial_stepsize
    # FIXME: this is some kludge, should use a proper NA representation
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
          self.transformer.setParameters(v)
          o_minus = objective_function(tp).fitness
          numEvaluations = numEvaluations + 1
          v[i] = current_values[i] + self.delta
          self.transformer.setParameters(v)
          o_plus = objective_function(tp).fitness
          numEvaluations = numEvaluations + 1
          gradient.append(o_plus - o_minus)
          # restore current state
          self.transformer.setParameters(current_values)
      if self.verbose and self.eliminateFlatComponents :
        sys.stderr.write('num_zeros: %d\n' % num_zeros)
      gradient_norm = transsys.utils.euclidean_norm(gradient)
      if gradient_norm == 0.0 :
        if self.verbose :
          sys.stderr.write('GradientOptimiser.optimise: flat gradient\n')
        break
      gradient = map(lambda x : x / gradient_norm, gradient)
      old_gradient = gradient
      # print gradient
      # if self.verbose :
      #   sys.stderr.write('obj = %g, values: %s\n' % (current_obj, str(current_values)))
      #   sys.stderr.write('  gradient = %s\n' % str(gradient))
      v = []
      for i in xrange(len(current_values)) :
        v.append(current_values[i] - gradient[i] * stepsize)
      self.transformer.setParameters(v)
      new_obj = objective_function(tp).fitness
      numEvaluations = numEvaluations + 1
      # 2010-09-23: moved line below out of if new_obj < current_obj block
      # new_values can otherwise end up unassigned in else branch (or even worse, may carry over values from previous iteration...)
      # FIXME: code is rather unmaintainable and thus prone to such problems, needs modularisation
      new_values = v[:]
      if new_obj < current_obj :
        stepsize2 = stepsize / self.stepsize_shrink
        v = []
        for i in xrange(len(current_values)) :
          v.append(current_values[i] - gradient[i] * stepsize2)
        self.transformer.setParameters(v)
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
          self.transformer.setParameters(v)
          new_obj2 = objective_function(tp).fitness
          numEvaluations = numEvaluations + 1
          if self.verbose :
            sys.stderr.write('growth branch: tried stepsize %g, new_obj = %g, new_obj2 = %g\n' % (stepsize2, new_obj, new_obj2))
        if not new_obj2 < new_obj :
          self.transformer.setParameters(new_values)
      else :
        while new_obj > current_obj :
          # FIXME: this ought to be done in the terminationCondition method -- this is just a quick fix to prevent infinite looping
          if (self.termination_stepsize is not None and stepsize < self.termination_stepsize ) or stepsize == 0.0 :
            break
          stepsize = stepsize * self.stepsize_shrink
          if self.verbose :
            sys.stderr.write('  %g: obj: current = %g, new = %g, d = %g\n' % (stepsize, current_obj, new_obj, current_obj - new_obj))
          new_values = []
          for i in xrange(len(current_values)) :
            new_values.append(current_values[i] - gradient[i] * stepsize)
          self.transformer.setParameters(new_values)
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
    return GradientOptimisationResult(tp, objective_function(tp), optimisation_log)


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


