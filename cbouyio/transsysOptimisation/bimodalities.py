#!/usr/bin/env python

# Subversion keywords.
# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date: 2008-01-08 15:44:39 +0000 (Tue, 08 Jan 2008) $:  Date of last commit

"""bimodalities.py

A pure python module for calculating bimodalirties.

The module contains all the high level functions to conduct all the necessary
bimodalities calculations. It calculates the bimodality score of a data set for
a given threshold, the optimal threshold wich maximises bimodality score and the
sum of bimodalities for a data table/frame.

@author: Costas Bouyioukos
@email: konsb@cmp.uea.ac.uk
@organisation: University of East Anglia
@since: 20/12/2007
@license: GNU General Public Lisence 3 or newer.
@contact: U{Costas Bouyioukos<mailto:konsb@cmp.uea.ac.uk>}
@version: $Id$"""


# Version Information.
__version__ = "$Id$"

import sys
import math
from statlib import stats


def bimodality(data, threshold):
  """Function to calculate bimodality for a given threshold.

  Bimodality is calculated as the distance between two means over the sum of
  their standard deviations.
  @param data: A one dimentional data set.
  @type data: C{list}
  @param threshold: A scalar variable.
  @type threshold: C{float}
  @precondition: threshold should be an element of the list containing the data.
  @return: The bimodality (C{float})
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
  mu1 = stats.lmean(subset1)
  mu2 = stats.lmean(subset2)
  sigma1 = stats.lstdev(subset1)
  sigma2 = stats.lstdev(subset2)
#  # Treat the special case of the sum of stdevs to be zero.
#  if sigma1 + sigma2 == 0 :
#    return 1e3000
  bimodality = float(abs(mu1 - mu2)) / float(sigma1 + sigma2)
  return bimodality


def bimodality_score(data):
  """Function to calculate the bimodality score (maximum bimodality) of an one
  dimensional data set.

  Bimodality score is defined as the maximum bimodality that a partition of the
  data set can have.
  @param data: A one dimentional data set.
  @type data: C{list}
  @return: The bimodality (C{float})
  """
  if len(data) < 4 :
    raise StandardError, 'Bimodality calculation does not make sense for such a small data set.'
  bmScore = 0
  for thres in set(data) :
    # Actual calculation of the bimodality score.
    bm = bimodality(data, thres)
    if bm > bmScore :
      bmScore = bm
  return bmScore


def total_bimodality(dataFrame):
  """Function to calculate the total bimodality score of a multi dimensional
  data set/frame.

  The total bimodality score is the sum of the bimodality scores for all the
  data vectors of the data frame.
  @param dataFrame: A mutli-dimenstional data frame.
  @type dataFrame: C{list} of C{list}s
  @return: The total bimodality score of a data set/frame. (C{float})
  """
  totalBM = 0
  for data in dataFrame :
    bm = bimodality_score(data)
    totalBM = totalBM + bm
  return totalBM

