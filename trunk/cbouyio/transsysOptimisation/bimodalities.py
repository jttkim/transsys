#!/usr/bin/env python

# Subversion keywords.
# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date$:  Date of last commit

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
  data.sort()
  n = len(data)
  # Check the existance of a threshold in the data list.
  if threshold not in data :
    raise StandardError, 'Specified threshold "%s" does not belong to the data. Cannot calculate bimodality' % str(threshold)
  # Safeguard for the minimum and the maximum value!
  if threshold == data[0] or threshold == data[n - 1] or threshold == data[n - 2] :
    raise StandardError, 'Specified threshold is not suiteble for bimodality calculations. It is either the smallest, the bigest or the second bigest element of the data list.'
  index = data.index(threshold) + 1
  subset1 = data[index:]
  subset2 = data[:index]
  mu1 = stats.mean(subset1)
  mu2 = stats.mean(subset2)
  sigma1 = stats.stdev(subset1)
  sigma2 = stats.stdev(subset2)
  bimodality = float(abs(mu1 - mu2))/float(sigma1 + sigma2)
  return bimodality


def bimodalityScore(data):
  """Function to calculate the bimodality score (maximum bimodality) of a one
  dimensional data set.

  Bimodality score is defined as the maximum bimodality that a partition of the
  data set can have.
  @param data: A one dimentional data set.
  @type data: C{list}
  @return: The bimodality (C{float})
  """


def bimodalityTotal(dataFrame):
  """Function to calculate the total bimodality score of a multi dimensional
  data set/frame.

  The total bimodality score is the sum of the bimodality scores for all the
  data vectors of the data frame.
  @param dataFrame: A mutli-dimenstional data frame.
  @type dataFrame: C{list} of C{list}s
  @return: The total bimodality score of a data set/frame. (C{float})
  """

