#!/usr/bin/env python

# Subversion keywords.
# $Rev::               $:  Revision of last commit
# $Author::            $:  Author of last commit
# $Date$:  Date of last commit

"""General description here.

@author: Costas Bouyioukos
@email: konsb@cmp.uea.ac.uk
@organisation: University of East Anglia
@since: <<DATE>>
@license: GNU General Public Lisence 2 or newer.
@contact: U{Costas Bouyioukos<mailto:konsb@cmp.uea.ac.uk>}
@version: $Id$"""


# Version Information.
__version__ = "$Id$"


import transsys
import igraph


# Parameter Initialization.
decay         = transsys.ExpressionNodeValue(0.2)
diffusibility = transsys.ExpressionNodeValue(0.0)
v_max         = transsys.ExpressionNodeValue(2.0)
k_m           = transsys.ExpressionNodeValue(1.0)
constitutive  = transsys.ExpressionNodeValue(0.1)



def igraphTOtranssys (graph, tpName='transsysFROMigraph') :
  """Function to convert an igraph object to a transsys program.

  @param graph: An igraph.Graph object.
  @type graph: C{class 'igraph.Graph'}
  @return: A transsys program
  @rtype: C{class 'transsys.TranssysProgram'}
  """

  transsysGenes = []
  transsysFactors = []

  for vertex in xrange(graph.vcount()) :
    transsysFactors.append(transsys.Factor('factor_%04d' % vertex, decay, diffusibility))

  for edge in graph.es : # This is the edge set of the graph.
    pass


def transsysTOigraph (tp) :
  """Function to generate an igraph object from a transsys program.

  @param tp: A transsys program object.
  @type tp: C{class 'transsys.TranssysProgram'}
  @return : An igraph object representing the network.
  @rtype : C{class 'igraph.Graph'}
  """

  geneInteractionList = [] # a list of 2n-tuples representing the edges.
  for gene in tp.gene_list :
    geneIndex = tp.find_gene_index(gene.name)
    factor = gene.product
    regulatedGenesIndex = []
    for rgene in tp.regulated_genes(factor) :
      regulatedGenesIndex.append(tp.find_gene_index(rgene.name))
    for rIndex in regulatedGenesIndex :
      geneInteractionList.append((geneIndex, rIndex))


  g = igraph.Graph(n=tp.num_genes(), edges=geneInteractionList, directed=True)
  return g

