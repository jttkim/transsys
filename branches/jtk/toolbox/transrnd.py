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


import math


class transrnd :

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
