g#!/usr/bin/env python

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


import statlib


def bimodality(data, threshold):
  """Function to calculate bimodality for a given threshold.

  Bimodality is calculated as the distance between two means over the sum of
  their standard deviations.
  @param data: A one dimentional data set.
  @type data: C{list}
  @param threshold: A scalar variable.
  @type threshold: C{float}
  @precondition: threshold should be an element of the data set.
  @return: The bimodality (C{float})
  """


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

