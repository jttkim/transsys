#!/usr/bin/env python

import sys
import copy
import math
import types
import os
import popen2
import random
import string

import transsys
import transdecode


class Interval :

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
    if type(other) is types.FloatType :
      return self.lower <= other and other <= self.upper
    if not isinstance(other, Interval) :
      raise StandardError, 'Interval.contains: illegal type'
    if self.lower <= other.lower and self.upper >= other.upper :
      return True
    else :
      return False


  def intersection(self, other) :
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


def mean_and_stddev(l) :
  """compute mean and standard deviation of a list of floating point values."""
  if len(l) == 1 :
    return l[0], 0.0
  m = sum(l) / float(len(l))
  d = map(lambda x : x - m, l)
  d2 = map(lambda x : x * x, d)
  v = sum(d2) / float(len(l) - 1)
  # print l
  # print v
  sd = math.sqrt(v)
  return m, sd


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
  m_active, s_active = mean_and_stddev(active_set)
  m_inactive, s_inactive = mean_and_stddev(inactive_set)
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


def hamming_distance(s1, s2) :
  if len(s1) != len(s2) :
    raise StandardError, 'length mismatch'
  n = 0
  for i in xrange(len(s1)) :
    if s1[i] != s2[i] :
      n = n + 1
  return n


def hillclimb(objective_function, decoder, initial_dnaseq, num_generations, mutator) :
  dnaseq = initial_dnaseq
  fitness = objective_function(decoder.decode_transsys('implant', dnaseq)).fitness
  for g in xrange(num_generations) :
    print '%d: %f' % (g, fitness)
    mutant_dnaseq = mutator.mutate(dnaseq)
    mutant_fitness = objective_function(decoder.decode_transsys('implant', mutant_dnaseq)).fitness
    d = hamming_distance(dnaseq, mutant_dnaseq)
    # print 'current: %f, mutant: %f, distance: %d (rel: %f)' % (fitness, mutant_fitness, d, float(d) / float(len(dnaseq)))
    if mutant_fitness <= fitness :
      dnaseq = mutant_dnaseq
      fitness = mutant_fitness
  return dnaseq


class ObjectiveFunction :

  def __init__(self) :
    pass


  def __call__(self, tp) :
    return None


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


class OptimisationResult :
  """Class for returning results of an optimisation.

The members are:
- optimised_transsys_program: The transsys program resulting
  from the optimisation process
- objectiveOptimum: The result of evaluating the optimised transsys
  program using the objective function
- optimisation_log: A trace of the optimisation process, provided
  as a list of tuples
The detailed structure of optimisation_log member may vary with the
nature of the particular optimiser. As an emerging standard, the
first element of each tuple (i.e. the leftmost column) should
be the objective value of the current solution. (For population
based optimisers, this might be replaced with the average objective
value.)"""

  def __init__(self, tp, objectiveOptimum, optimisation_log) :
    self.optimised_transsys_program = tp
    self.objectiveOptimum = objectiveOptimum
    self.optimisation_log = optimisation_log


class GradientOptimiser :

  def __init__(self, initial_stepsize, delta = 1.0e-6, stepsize_shrink = 0.5, stepsize_min = 1.0e-30, stepsize_max = 1.0e30, improvement_threshold = 0.0) :
    """
@param initial_stepsize: distance initially stepped in gradient
  direction. Subject to subsequent adaptation
@param delta: offset from current point in search space used
  for gradient computation (estimation)
@param stepsize_shrink: factor by which stepsize is multiplied / divided
  in stepsize adaptation
@param stepsize_min: minimal stepsize, optimisation terminates if
  stepsize falls below this limit
@param stepsize_max: maximal stepsize, adaptation will not make
  stepsize exceed this limit
@param improvement_threshold: optimisation terminates if improvement
  drops below this threshold.
"""
    self.initial_stepsize = initial_stepsize
    self.delta = delta
    self.stepsize_shrink = stepsize_shrink
    self.stepsize_min = stepsize_min
    self.stepsize_max = stepsize_max
    self.improvement_threshold = improvement_threshold
    self.eliminateFlatComponents = False
    self.verbose = False


  def optimise(self, transsys_program, objective_function, gene_name_list = None, factor_name_list = None) :
    optimisation_log = []
    tp = copy.deepcopy(transsys_program)
    value_expression_list = get_value_nodes(tp, gene_name_list, factor_name_list)
    if self.verbose :
      sys.stderr.write('optimising %d values\n' % len(value_expression_list))
    initial_values = map(lambda n : n.value, value_expression_list)
    current_values = initial_values[:]
    current_obj = objective_function(tp).fitness
    optimisation_log = [current_obj]
    stepsize = self.initial_stepsize
    old_gradient = [1.0] * len(current_values)
    while stepsize > self.stepsize_min :
      gradient = []
      num_zeros = 0
      for i in xrange(len(current_values)) :
        if self.eliminateFlatComponents and old_gradient[i] == 0.0 :
          gradient.append(0.0)
          num_zeros = num_zeros + 1
        else :
          node = value_expression_list[i]
          node.value = current_values[i] - self.delta
          o_minus = objective_function(tp).fitness
          node.value = current_values[i] + self.delta
          o_plus = objective_function(tp).fitness
          gradient.append(o_plus - o_minus)
          node.value = current_values[i]
      # print 'num_zeros: %d' % num_zeros
      gradient_norm = math.sqrt(reduce(lambda x, y : x + y * y, gradient, 0.0))
      old_gradient = gradient
      if gradient_norm == 0.0 :
        if self.verbose :
          sys.stderr.write('GradientOptimiser.optimise: flat gradient\n')
        break
      gradient = map(lambda x : x / gradient_norm, gradient)
      # print gradient
      if self.verbose :
        sys.stderr.write('obj = %g, values: %s\n' % (current_obj, str(current_values)))
        sys.stderr.write('  gradient = %s\n' % str(gradient))
      for i in xrange(len(current_values)) :
        value_expression_list[i].value = current_values[i] - gradient[i] * stepsize
      new_obj = objective_function(tp).fitness
      if new_obj < current_obj :
        v = map(lambda n : n.value, value_expression_list)
        stepsize2 = stepsize / self.stepsize_shrink
        for i in xrange(len(current_values)) :
          value_expression_list[i].value = current_values[i] - gradient[i] * stepsize2
        new_obj2 = objective_function(tp).fitness
        if self.verbose :
          sys.stderr.write('growth branch: tried stepsize %g, new_obj = %g, new_obj2 = %g\n' % (stepsize2, new_obj, new_obj2))
        while new_obj2 < new_obj and stepsize < self.stepsize_max :
          stepsize = stepsize2
          if self.verbose :
            sys.stderr.write('growing stepsize to %g\n' % stepsize)
          new_obj = new_obj2
          v = map(lambda n : n.value, value_expression_list)
          stepsize2 = stepsize / self.stepsize_shrink
          for i in xrange(len(current_values)) :
            value_expression_list[i].value = current_values[i] - gradient[i] * stepsize2
          new_obj2 = objective_function(tp).fitness
          if self.verbose :
            sys.stderr.write('growth branch: tried stepsize %g, new_obj = %g, new_obj2 = %g\n' % (stepsize2, new_obj, new_obj2))
        if not new_obj2 < new_obj :
          for i in xrange(len(value_expression_list)) :
            value_expression_list[i].value = v[i]
      else :
        while new_obj >= current_obj and stepsize > self.stepsize_min :
          stepsize = stepsize * self.stepsize_shrink
          if self.verbose :
            sys.stderr.write('  %g: obj: current = %g, new = %g, d = %g\n' % (stepsize, current_obj, new_obj, current_obj - new_obj))
          for i in xrange(len(current_values)) :
            value_expression_list[i].value = current_values[i] - gradient[i] * stepsize
          new_obj = objective_function(tp).fitness
      if stepsize > self.stepsize_min :
        d_obj = current_obj - new_obj
        if self.verbose :
          sys.stderr.write('  making step %g: d_obj = %g, new_obj = %g\n' % (stepsize, d_obj, new_obj))
        current_obj = new_obj
        optimisation_log.append(current_obj)
        current_values = map(lambda n : n.value, value_expression_list)
        if d_obj < self.improvement_threshold :
          if self.verbose :
            sys.stderr.write('  terminating because d_obj = %f < threshold %f\n' % (d_obj, self.improvement_threshold))
          break
    return OptimisationResult(tp, current_obj, optimisation_log)


class LsysObjectiveFunction(ObjectiveFunction) :

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


