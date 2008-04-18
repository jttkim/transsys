#!/usr/bin/env python

# Subversion keywords.
# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date: 2008-01-08 15:44:39 +0000 (Tue, 08 Jan 2008) $:  Date of last commit

"""bimodalities.py

A pure python module for calculating bimodalirties.

The module contains all the high level functions to conduct all the necessary
bimodalities calculations. It calculates the bimodality score of a
TranssysInstanceCollection object. Returns also (as instance variables) the
thresholds and the partial bimodality scores.

@author: Costas Bouyioukos
@organization: University of East Anglia
@since: 20/12/2007
@license: GNU General Public Lisence 3 or newer.
@contact: U{Costas Bouyioukos<mailto:konsb@cmp.uea.ac.uk>}
@version: $Id$
"""


# Version Information.
__version__ = "$Id$"

import sys
import math
import statlib.stats


class BimodalitiesCollection(object) :
  """A class wich gathers the bimodality calculations.

  @ivar totalBimodality: The total bimodality score of an extire lattice.
  vector.
  @type totalBimodality: Numeric
  @ivar thresholds: A dictionary of factors (keys) and threshold values (values)
  by wich the bimodality scores has been calculated.
  @type thresholds: A C{dict} of <factor><thresholds> key/value pairs.
  @ivar bimodalities: A dictionary of factors (keys) and bimodality scores
  (values).
  @type bimodalities:  A C{dict} of <factor><bimodalities> key/value pairs.
  """


  def __init__(self, transsysLattice) :
    """The constructor

    @param transsysLattice: A transsys instance lattice object. (this class can
    be constructed for any other TranssysInstanceCollection object).
    @type transsysLattice: C{class 'transsys.Transsyslattice'}
    """
    self.bimodalities = {}
    self.thresholds = {}
    self.totalBimodality = self.total_bimodality(transsysLattice)


  def total_bimodality(self, lattice) :
    """Calculates the bimodality score for the relevant thresold. Stores the
    bimodalities and the thresholds into the relative dictionaries and returns
    the total bimodality score of the lattice.

    The function raises an exception when the length of the expression profile
    is less than 4 (a bimodality score cannot be defined for lists smaller than
    four.
    @return: The bimodality score of the lattice (TranssysInstanceCollection).
    @rtype: C{float}
    """
    totalBM = 0
    for factorName in lattice.transsysProgram.factor_names() :
      bmScore = 0
      optThres = 0
      factorExpression = lattice.get_factor_expression_list(factorName)
      if len(factorExpression) < 4 :
        raise StandardError, 'Bimodality calculation does not make sense for such a small data set.'
      # iterate over set... set is a non reduntant set of the factorExpression
      # list.
      for thres in set(factorExpression) :
        # Actual calculation of the bimodality score.
        bm = bimodality(factorExpression, thres)
        if bm > bmScore :
          bmScore = bm
          optThres = thres
      self.bimodalities[factorName] = bmScore
      self.thresholds[factorName] = optThres
    for bm in self.bimodalities.itervalues() :
      totalBM = totalBM + bm
    return totalBM



def bimodality(data, threshold):
  """Function to calculate the bimodality of an one dimensional data set for a
  given threshold.

  First the data set is partitioned in two subsets S1 and S2 (for threshold t
  then S1 = {x_i, if x_i <= t} and S2 = {x_i, if x_i > t}).
  Then, bimodality is calculated as the quotient of the distance between the
  means of each subset over the sum of the standard deviations.
  Returns zero if all the elements are equal and also zero if any subset
  is smaller than one.
  @param data: A one dimentional data set.
  @type data: C{list}
  @param threshold: A scalar variable.
  @type threshold: Numeric
  @return: The bimodality.
  @rtype: C{float}
  """
  subset1 = []
  subset2 = []
  # Partition the data.
  for e in data :
    if e <= threshold :
      subset1.append(e)
    else :
      subset2.append(e)
  # Check the length of the subsets.
  if len(subset1) <= 1 or len(subset2) <= 1 :
    return 0
  mu1 = statlib.stats.lmean(subset1)
  mu2 = statlib.stats.lmean(subset2)
  sigma1 = statlib.stats.lstdev(subset1)
  sigma2 = statlib.stats.lstdev(subset2)
  bimodality = float(abs(mu1 - mu2)) / float(sigma1 + sigma2)
  return bimodality

