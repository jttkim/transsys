#!/usr/bin/env python

# $Id$

# $Log$
# Revision 1.1  2005/03/08 17:12:02  jtk
# Initial revision
#
# Revision 1.2  2003/02/04 23:46:21  kim
# added RouletteWheel
#
# Revision 1.1  2003/01/28 21:12:05  kim
# initial toolbox assembly
#


"""Miscellaneous utilities.
"""

import types
import math
import random


def is_nan(x) :
  """Determine whether x is NaN.

This function is a temporary solution, to be replaced if and when
python provides this functionality. The current implementation
returns x != x; this works with current versions of python and Linux,
but cannot be guaranteed to work generally.

@param x: the floating point object to be tested
@return: C{True} if C{x} is NaN
@rtype: boolean
"""
  return x != x


def tablecell(x) :
  """Render C{x} as a string suitable as an R table element.

Currently supported types are:
  - C{None}, rendered as C{NA}
  - C{int}
  - C{float}, formatted using C{'%1.17e'} as a format string
  - C{boolean}, rendered as C{TRUE} or C{FALSE}.

@param x: the value to be converted to a table element
@return: a string representing C{x}
@rtype: C{String}
"""
  if x is None :
    return 'NA'
  elif type(x) is types.IntType :
    return '%d' % x
  elif type(x) is types.FloatType :
    return '%1.17e' % x
  elif type(x) is types.BooleanType :
    if x :
      return 'TRUE'
    else :
      return 'FALSE'
  raise StandardError, 'unsupported type %s' % str(type(x))


def hamming_distance(s1, s2) :
  """Compute the Hamming distance (number of different elements)
between two strings (or other sequence types).

Notice that this function compares objects for identity, not for
equality.

@param s1: first sequence
@param s2: second sequence
@return: Hamming distance between C{s1} and C{s2}
@rtype: int
"""
  if len(s1) != len(s2) :
    raise StandardError, 'length mismatch'
  n = 0
  for i in xrange(len(s1)) :
    if s1[i] != s2[i] :
      n = n + 1
  return n


def inner_product(v1, v2) :
  """Compute the inner product of two vectors.

@param v1: the first vector (any sequence consisting of numeric elements)
@param v2: the second vector (a sequence of the same length)
@return: the inner product of C{v1} and C{v2}
@rtype: float
"""
  if len(v1) != len(v2) :
    raise StandardError, 'unequal vector lengths'
  s = 0.0
  for i in xrange(len(v1)) :
    s = s + v1[i] * v2[i]
  return s


def euclidean_distance_squared(v1, v2) :
  """Compute the square of the Euclidean distance between C{v1} and C{v2}."""
  if len(v1) != len(v2) :
    raise StandardError, 'length mismatch'
  d2 = 0.0
  for i in xrange(len(v1)) :
    d = v1[i] - v2[i]
    d2 = d2 + d * d
  return d2


def euclidean_distance(v1, v2) :
  """Compute the Euclidean distance between C{v1} and C{v2}."""
  return math.sqrt(euclidean_distance_squared(v1, v2))


def euclidean_norm_squared(v) :
  """Comute the square of the euclidean norm of C{v}."""
  return sum(map(lambda x: x * x, v))


def euclidean_norm(v) :
  """Comute the euclidean norm of C{v}."""
  return math.sqrt(euclidean_norm_squared(v))


def mean_and_stddev(l) :
  """Compute mean and standard deviation of a list of floating point values.

@return: the mean and the standard deviation
@rtype: tuple
"""
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


def correlation_coefficient(x, y) :
  """Compute the Pearson correlation coefficient of C{x} and C{y}."""
  if len(x) != len(y) :
    raise StandardError, 'x and y have unequal length'
  nx = euclidean_norm(x)
  if nx == 0.0 :
    raise StandardError, 'standard deviation of x is zero'
  ny = euclidean_norm(y)
  if ny == 0.0 :
    raise StandardError, 'standard deviation of y is zero'
  r = inner_product(x, y) / euclidean_norm(x) / euclidean_norm(y)
  return r


class UniformRNG :
  """Callable class producing random floating point values from
a uniform distribution.

Based on the standard random module.
"""

  def __init__(self, rndseed, min_value = 0.0, max_value = 1.0) :
    """
@param rndseed: random seed
@param min_value: minimal value (inclusive)
@param max_value: maximal value (exclusive)
@return: a random value r from [min_value, max_value[
@rtype: float
"""    
    self.rng = random.Random(rndseed)
    self.min_value = min_value
    self.max_value = max_value


  def __call__(self) :
    return self.rng.uniform(self.min_value, self.max_value)


def randomise_transsys_values(transsys_program, random_function, gene_name_list = None, factor_name_list = None) :
  """Randomise numerical values in value nodes by replacing them
with values obtained from the random_function.

gene_name_list and factor_name_list identify the genes and the
factors to have their numerical values randomised. The default value
of None indicates that all genes / factors should be randomised.
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
  for n in value_expression_list :
    n.value = random_function()


class transrnd :
  """Random number generator class.

This is a legacy random number generator which was introduced before
the standard Python random package appeared. It is a Python implementation
of the random number generator implemented by the C{random} function
in the GNU C standard library.

This class should not be used anymore in new projects.
"""

  def __init__(self, seed = 1, state = None, fptr = None, rptr = None) :
    self.rand_type = 4
    self.rand_deg = 63
    self.rand_sep = 1
    self.RAND_MAX = 0x7fffffff
    self.gset = None
    if state is not None and fptr is not None and rptr is not None :
      raise StandardError, 'transrnd::__init__: restarting from saved state info not implemented'
      self.state = copy.deepcopy(state)
      self.fptr = fptr
      self.rptr = rptr
    else :
      self.fptr = 2
      self.rptr = 1
      self.seed = seed
      self.srandom(seed)


  def srandom(self, seed) :
    self.state = [0L] * self.rand_deg
    self.state[0] = long(seed)
    self.gset = None
    for i in range(1, self.rand_deg) :
      self.state[i] = 1103515245L * self.state[i - 1] + 12345L;
      # self.state[i] = 1103515145L * self.state[i - 1] + 12345L;  # linux libc version
    self.fptr = self.rand_sep
    self.rptr = 0
    for i in range(10 * self.rand_deg) :
      self.random()


  def random(self) :
    self.state[self.fptr] = (self.state[self.fptr] + self.state[self.rptr]) & 0xffffffffL
    i = (self.state[self.fptr] >> 1) & 0x7fffffff
    self.fptr = self.fptr + 1
    if self.fptr >= self.rand_deg :
      self.fptr = 0
      self.rptr = self.rptr + 1
    else :
      self.rptr = self.rptr + 1
      if self.rptr >= self.rand_deg :
        self.rptr = 0
    return int(i)


  def rnd(self) :
    return float(self.random()) / (float(self.RAND_MAX) + 1)


  def random_range(self, rceil) :
    m = long(self.RAND_MAX) + 1L
    m = m - m % rceil
    r = self.random()
    while r >= m :
      r = self.random()
    return r % rceil


  def gauss(self) :
    """adapted from Numerical Recipes in C"""
    if self.gset is None :
      v1 = 2.0 * self.random() / 0x7fffffff - 1.0
      v2 = 2.0 * self.random() / 0x7fffffff - 1.0
      rsq = v1 * v1 + v2 * v2
      while rsq >= 1.0 or rsq == 0.0 :
        v1 = 2.0 * self.random() / 0x7fffffff - 1.0
        v2 = 2.0 * self.random() / 0x7fffffff - 1.0
        rsq = v1 * v1 + v2 * v2
      fac = math.sqrt(-2.0 * math.log(rsq) / rsq)
      self.gset = v1 * fac
      return v2 * fac
    else :
      gset = self.gset
      self.gset = None
      return gset
    

class RouletteWheel :
  """Roulette wheel for evolutionary algorithms and similar purposes.
"""
  def __init__(self, rng, d = None) :
    if d is None :
      d = [1.0]
    self.set_wheel(d)
    self.rng = rng


  def gnuplot(self, fname) :
    f = open(fname, 'w')
    for w in self.wheel :
      f.write('%f\n' % w)


  def set_wheel(self, d) :
    x = 0.0
    for v in d :
      x = x + v
    self.wheel = [0.0]
    for v in d :
      self.wheel.append(self.wheel[-1] + v / x)


  def pocket(self) :
    x = self.rng.rnd()
    i = 1
    j = len(self.wheel)
    while i < j :
      k = (i + j) / 2
      # print x, self.wheel[k - 1], self.wheel[k], i, j, k
      if self.wheel[k - 1] <= x and x < self.wheel[k] :
        return k - 1
      elif x >= self.wheel[k] :
        i = k
      else :
        j = k
    raise ValueError, 'RouletteWheel::slot: bad value ' % x


