#!/usr/bin/env python

import string
import copy
import math
import transrnd
import transsys
import types
import sys

class ErrorLogger :

  def __init__(self, fname = None, vlevel = 0) :
    """might still need some modification...
    The higher the level the more gets printed.
    """
    self.errorfilename = fname
    if fname is None :
      self.errfile = sys.stderr
    else :
      self.errfile = open(fname, 'w')
    self.vlevel = vlevel
    self.errors=0


  def message(self, msg, level = 0) :
    """print errormessage to screen or file fname, if it has lower or equall verbositylevel than vlevel.
    """
    if level < self.vlevel:
      self.errors=self.errors + 1
    if self.vlevel >= level :
      self.errfile.write('%s\n' % msg)
      self.errfile.flush()


def uniform_cmp(element1, element2) :
  if not isinstance(element1, int) :
    raise StandardError, 'uniform_cmp:: type of element1 is not INT'
  if not isinstance(element2, int) :
    raise StandardError, 'uniform_cmp:: type of element2 is not INT'
  if element1 == element2 :
    if element1 == 1 :
      return 1
  return 0


class ExpressionSimilarityMatrix :
  """The ExpressionSimilarityMatrix is basically an implementation of a matrix, that is suited for our purposes. This means that entries can easily be extracted by factornames instead of indices. 
  """

  def __init__(self, name, gf_dict, default_element = 0.0) :
    """Constructor. Initialize ExpressionSimilarityMatrix with all factors given as values in the dictionary gf_dict. When an error is raised the network has factors that are encoded by multiple genes. this functionality is not (yet) implemented.
    """
    self.name = name
    self.matrix = []
    self.gf_dict = gf_dict
    self.factor_list = []
    for fn in self.gf_dict.values() :
      if fn in self.factor_list :
        raise StandardError, 'multiple genes encode same factor'
      self.factor_list.append(fn)
    self.factor_list.sort()
    self.fill_all(default_element)


  def __str__(self, factor_fmt = '%10s', element_fmt = ' %7.3f', none_symbol = ' None  ') :
    """Prints the matrix.
    """
    s = ''
    flist = self.get_fnamelist()
    for i in xrange(len(flist)) :
      s = s + (factor_fmt % flist[i])
      glue = ': '
      for j in xrange(len(flist)) :
        s = s + glue
        if self.matrix[i][j] is None:
          s = s + none_symbol
        else :
          s = s + (element_fmt % float(self.matrix[i][j]))
        glue = ' '
      s = s + '\n'
    return s


  def gnu_comparison(self, other, cmp_function) :
    """this function seems to be half-implemented -- it does not do anything
reasonable or useful"""
    raise StandardError, 'ExpressionSimilarityMatrix::gnu_comparison: unimplemented'
    if self.factor_list != other.factor_list :
      raise StandardError, 'ExpressionSimilarityMatrix.gnu_comparison:: matrices have not got the same factor_list'
    for factor1 in self.factor_list :
      for factor2 in self.factor_list :
        element1 = self.get_element(factor1, factor2)
        element2 = other.get_element(factor1, factor2)
        if cmp_function(element1, element2) == 1 :
          pass
          

  def gene_factor_dictionary(self):
    """Returns the gene-factor dictionary.
    """
    return self.gf_dict


  def factor_gene_dictionary(self) :
    fg_dict = {}
    for k in self.gf_dict.keys() :
      v = self.gf_dict[k]
      if v in fg_dict.keys() :
        raise StandardError, 'ExpressionSimilarityMatrix.factor_gene_dictionary: multiple genes encode factor "%s"' % v
      fg_dict[v] = k
    return fg_dict


  def return_as_list(self) :
    """Returns the matrix as a list where rows follow each other.
    """
    matrix2list = []
    flist = self.get_fnamelist()
    for i in flist :
      for j in flist :
        matrix2list.append(self.get_element(i,j))
    return matrix2list


  def fill_all(self, v) :
    """Writes v in each entry of the matrix.
    """
    flist = self.get_fnamelist()
    self.matrix = []
    for i in xrange(len(flist)) :
      row = [v] * len(flist)
      self.matrix.append(row)


  def get_index(self, fname) :
    """Returns the index of the factor fname in the membervariable factor_list owned by the matrix.
    """
    if self.factor_list.count(fname) != 1 :
      raise StandardError, 'factor "%s" not in matrix' % fname
    return self.factor_list.index(fname)


  def get_fname(self, gname) :
    """Returns the product of the gene gname.
    """
    return self.gf_dict[gname]


  def get_element(self, f1, f2) :
    """Returns the entry of the matrix associated with the factors f1 and f2.
    """
    i = self.get_index(f1)
    j = self.get_index(f2)
    return self.matrix[i][j]


  def set_element(self, f1, f2, element) :
    """Fill entry of the matrix associated with the factors f1, f2 with element.
    """
    i = self.get_index(f1)
    j = self.get_index(f2)
    self.matrix[i][j] = element


  def get_gnamelist(self) :
    """Returns the gene names in a list.
    """
    return self.gf_dict.keys()


  def get_fnamelist(self) :
    """Returns the factor names in a list.
    """
    return self.factor_list


  def replace_none(self, element) :
    flist = len(self.get_fnamelist())
    for i in xrange(flist) :
      for j in xrange(flist) :
        if self.matrix[i][j] is None :
          self.matrix[i][j] = element


class KnockoutTranssysProgram(transsys.TranssysProgram) :
  """KnockoutTranssysProgram
  """

  def __init__(self, tp, nonfunc_name = "nonfunc") :
    """docu missing
    """
    self.basename = tp.name
    self.name = tp.name
    self.factor_list = tp.factor_list[:]
    self.gene_list = tp.gene_list[:]
    self.comments = tp.comments
    self.nonfunc_name = nonfunc_name
    self.add_nonfunctional_factor()
    self.functional_factor = None
    self.knockout_index = None


  def add_nonfunctional_factor(self) :
    """docu missing
    """
    # FIXME: should check for existence of nonfunc factor somewhere
    self.nonfunc = transsys.Factor(self.nonfunc_name)
    self.factor_list.append(self.nonfunc)


  def do_knockout(self, knockout_index) :
    """docu missing
    """
    if self.knockout_index is not None :
      raise StandardError, 'multiple knockout backtrace not (yet) implemented'
    self.name = '%s_ko_%d' % (self.basename, knockout_index)
    self.knockout_index = knockout_index
    self.functional_factor = self.gene_list[knockout_index].product
    self.gene_list[knockout_index].product = self.nonfunc


  def undo_knockout(self) :
    """docu missing
    """
    if self.knockout_index is None :
      raise StandardError, 'no knockout done which could be undone'
    self.name = self.basename
    self.gene_list[self.knockout_index].product = self.functional_factor
    self.knockout_index = None
    self.functional_factor = None


def link_scorefunc(promoter) :
  if isinstance(promoter, transsys.PromoterElementLink) :
    return 1
  else :
    # scoreMatrix uses NONE to identify set matrix elements.
    # Errors 
    return None


def score_matrix(network, promoter_element_scorefunc) :
  """compute a distance adjacency matrix in which entries connecting
a regulator to a cognate regulatee contain a score value computed by
the scorefunc on the promoter element. All other entries are None.
  """
  gf_dict = gene_factor_dictionary(network)
  adj_matrix = ExpressionSimilarityMatrix('adj_matrix', gf_dict)
  fg_dict = adj_matrix.factor_gene_dictionary()
  adj_matrix.fill_all(None)
  for gene_i in gf_dict.keys():
    fi = gf_dict[gene_i]
    gene = network.find_gene(gene_i)
    for p in gene.promoter :
      if isinstance(p, transsys.PromoterElementLink) :
        for f in p.factor_list :
          fj = f.name
          #
          # FIXME ??? what's wrong here?
          if adj_matrix.get_element(fj, fi) is not None :
            errlogger.message('score_matrix_error: multiple links %s --> %s' % (fj, fi), 0)
            raise StandardError, 'score_matrix: multiple links %s --> %s' % (fj, fi)
          #
          #         
          adj_matrix.set_element(fj, fi, promoter_element_scorefunc(p))
  return adj_matrix


def pathscore_matrix(tp, scorefunc, no_connection_default_value) :
  """get cost matrix of a transsys network. Calculation rules are implemented in a generic score_function().
  """
  gf_dict = gene_factor_dictionary(tp)
  scoreMatrix = score_matrix(tp, scorefunc)
  scoreMatrix.name = 'scoreMatrix'
  # obviously distance from gene i to itself is zero
  for i in gf_dict.keys() :
    scoreMatrix.set_element(gf_dict[i], gf_dict[i], 0)
  for n in gf_dict.keys() :
    for i in gf_dict.keys() :
      for j in gf_dict.keys() :
        s1 = scoreMatrix.get_element(gf_dict[i], gf_dict[j])
        s2 = scoreMatrix.get_element(gf_dict[i], gf_dict[n])
        s3 = scoreMatrix.get_element(gf_dict[n], gf_dict[j])
        if s1 is not None :
          if s2 is not None :
            if s3 is not None :
              if s1 > s2 + s3 :
                scoreMatrix.set_element(gf_dict[i], gf_dict[j], s2 + s3)
            else :
              scoreMatrix.set_element(gf_dict[i], gf_dict[j], s1)
          else :
            scoreMatrix.set_element(gf_dict[i], gf_dict[j], s1)
        else : 
          if s2 is not None :
            if s3 is not None :
              scoreMatrix.set_element(gf_dict[i], gf_dict[j], s2 + s3)
  scoreMatrix.replace_none(no_connection_default_value)
  return scoreMatrix


# -1 := activating
# +1 := repressing
# must be consistent with matrix D !!!
def adjacency_matrix_scorefunc(promoter) :
  if isinstance(promoter, transsys.PromoterElementLink) :
    if isinstance(promoter, transsys.PromoterElementActivate) :
      return -1
    elif isinstance(promoter, transsys.PromoterElementRepress) :
      return 1
    else :
      errlogger.message('adjacency_error: unknown PromoterElementLink subclass "%s"' % promoter.__class__.__name__ , 0)
      raise StandardError, 'adjacency_matrix_scorefunc: unknown PromoterElementLink subclass "%s"' % promoter.__class__.__name__
  else :
    errlogger.message('adjacency_error: unknown PromoterElement class "%s"' % promoter.__class__.__name__ , 0)
    raise StandardError, 'adjacency_matrix_scorefunc: unknown PromoterElement class "%s"' % promoter.__class__.__name__


def adjacency_matrix (network) :
  adj_matrix = score_matrix(network, adjacency_matrix_scorefunc)
  adj_matrix.replace_none(0)
  return adj_matrix


def gene_factor_dictionary(tp) :
  """make a dictionary, where gene is the key and the product is the value, for easy and convenient access.  It also provides a more understandable code and prevents errors, which could occur when using indices.
  """
  d = {}
  for gene in tp.gene_list :
    d[gene.name] = gene.product.name
  return d


def identity_perturber(c, fn) :
  """docu missing
  """
  return c


#FIXME: a_spec, a_max is both used for activation and repression
def create_disruption_network(name, matrix_D, gene_order, decay_rate, a_spec, a_max, constit) :
  """docu missing
  """
  factor_list = []
  gene_list = []
  for gname in gene_order :
    fname = matrix_D.get_fname(gname)
    factor_list.append(transsys.Factor(fname, transsys.ExpressionNodeValue(decay_rate)))
  for gname in gene_order :
    fname = matrix_D.get_fname(gname)
    temp_promoter_list = []
    for reg_gname in gene_order :
      reg_fname = matrix_D.get_fname(reg_gname)
      if matrix_D.get_element(reg_fname, fname) == -1 :
        temp_promoter_list.append(transsys.PromoterElementActivate(transsys.ExpressionNodeValue(a_spec), transsys.ExpressionNodeValue(a_max), [reg_fname]))
      if matrix_D.get_element(reg_fname, fname) ==  1 :
        temp_promoter_list.append(transsys.PromoterElementRepress(transsys.ExpressionNodeValue(a_spec), transsys.ExpressionNodeValue(a_max), [reg_fname]))
    temp_promoter_list.append(transsys.PromoterElementConstitutive(transsys.ExpressionNodeValue(constit)))
    gene_list.append(transsys.Gene(gname, fname, temp_promoter_list))
  disruption_network = transsys.TranssysProgram(name, factor_list, gene_list)
  return disruption_network


# matrix_D:
# ---------
# dim(matrix_D) = n x n
#
# value_k  e {-1, 0, 1}
# list_k   = [value_0, .., value_n]
# matrix_D = [[list_0], .., [list_n]]
#
# column_i of matrix_D = list_i
# row_j of matrix_D    = list_1[j], .., list_n[j]

def reconstructed_adjacency_matrix(matrix_E, gamma, sigma) :
  """docu missing
  """
  # initialize matrix D
  if sigma is None :
    sigma = {}
    for fname in matrix_E.get_fnamelist() :
      sigma[fname] = 1.0
  if len(sigma) != len(matrix_E.get_fnamelist()) :
    errlogger.message('reconstructed_adjacency_matrix: sigma_error: invalid size ',errorlogger_verbositylevel)
    raise StandardError, 'matrix_D: invalid size of vector sigma'
  matrix_D = copy.deepcopy(matrix_E)
  matrix_D.name = 'matrix_D'
  # build matrix D from matrix E
  for fname in matrix_D.get_fnamelist() :
    for reg_fname in matrix_D.get_fnamelist() :
      # "normalize"
      x = matrix_E.get_element(fname, reg_fname) / sigma[fname]
      # d_ij = ...
      # ignore selfregulation, meaning edges with in/out gene is the same
      if fname == reg_fname :
        matrix_D.set_element(fname, reg_fname, 0)
      # gene is activating
      elif x <= -gamma :
        matrix_D.set_element(fname, reg_fname, -1)
      # gene is repressing
      elif x >= gamma :
        matrix_D.set_element(fname, reg_fname, 1)
      # gene is not regulating
      else :
        matrix_D.set_element(fname, reg_fname, 0)
  return matrix_D    


#FIXME: i want to be an array not a matrix anymore...!
####
def referencestate_concentration_matrix(network, timesteps) :
  """docu missing
  """
  gfd = gene_factor_dictionary(network)  
  ti = transsys.TranssysInstance(network)
  tseries = ti.time_series(timesteps, 1)
  reference_state = tseries[max(tseries.keys())]
  ref_conc_matrix = ExpressionSimilarityMatrix('ref_conc_matrix', gfd )
  for fi in network.factor_names() :
    for fj in network.factor_names() :
      ref_conc_matrix.set_element(fj, fi, reference_state.factor_concentration[network.find_factor_index(fi)])
  return ref_conc_matrix


####
#### FIXME: see referencestate_concentration_matrix
def knockout_concentration_matrix(knockout_network, ref_conc_matrix, timesteps) :
  """compute a matrix of expression values obtained from single gene
knockout networks. The matrix element m[i][j] contains the expression
level of factor j observed in a network in which factor i is knocked out
(by replacement with a nonfunc product in the encoding gene). Expression
starts with the reference state, taken from ref_conc_matrix, and is
carried out for timesteps timesteps.
  """
  gfd_wildtype = gene_factor_dictionary(knockout_network)
  conc_matrix = ExpressionSimilarityMatrix('conc_matrix', gfd_wildtype)
  ref_state = transsys.TranssysInstance(knockout_network)
  wildtype_factor_names = knockout_network.factor_names()[:]
  wildtype_factor_names.remove(knockout_network.nonfunc_name)
  for factor in wildtype_factor_names :
    ref_state.factor_concentration[knockout_network.find_factor_index(factor)] = ref_conc_matrix.get_element(factor, factor)
  # make 'experiments'  
  for gi in knockout_network.gene_names() :
    # mutate gene i
    knockout_network.do_knockout(knockout_network.find_gene_index(gi))
    initial_state = ref_state.perturbed_copy(identity_perturber)
    initial_state.transsys_program = knockout_network
    m_tseries = initial_state.time_series(timesteps, 1)
    this_mutant = m_tseries[max(m_tseries.keys())]
    for fi in wildtype_factor_names :
      # expression of fi is affected by knockout of gi as quantified by concentration of fi
      conc_matrix.set_element(gfd_wildtype[gi], fi, this_mutant.factor_concentration[knockout_network.find_factor_index(fi)])
    # demutate gene i
    knockout_network.undo_knockout()
  return conc_matrix


def knockout_concentration_matrix_lsys(knockout_network, aux_network, lsys_lines, timesteps) :
  """docu missing
  """
  # gfd_wildtype contains only the wildtype genes because
  # KnockoutTranssysProgram has only one factor added, genes are
  # unchanged (unless in knockout state, which we don't expect at the start).
  gfd_wildtype = gene_factor_dictionary(knockout_network)
  conc_matrix = ExpressionSimilarityMatrix('conc_matrix', gfd_wildtype)
  wildtype_factor_names = knockout_network.factor_names()[:]
  wildtype_factor_names.remove(knockout_network.nonfunc_name)
  # make 'experiments'  
  for gi in knockout_network.gene_names() :
    # mutate gene i
    # sys.stderr.write('knocking out gene "%s"\n' % gi)
    knockout_network.do_knockout(knockout_network.find_gene_index(gi))
    if aux_network is None :
      merged_network = knockout_network
    else :
      merged_network = copy.deepcopy(aux_network)
      merged_network.merge(knockout_network)
    # print merged_network
    initial_state = transsys.TranssysInstance(merged_network)
    m_tseries = initial_state.time_series(timesteps, 1, lsys_lines)
    # FIXME: all the intervening entries in the time series are not needed...
    this_mutant = m_tseries[max(m_tseries.keys())]
    for fi in wildtype_factor_names :
      # expression of fi is affected by knockout of gi as quantified by concentration of fi
      i = merged_network.find_factor_index(fi)
      c = this_mutant.factor_concentration[i]
      # sys.stderr.write('  [%s] (index %d) is now %g\n' % (fi, i, c))
      conc_matrix.set_element(gfd_wildtype[gi], fi, c)
    # demutate gene i
    knockout_network.undo_knockout()
  return conc_matrix


####
####
def logratio_matrix(ref_conc_matrix, conc_matrix) :
  """docu missing
  """
  overflows = 0
  log2 = math.log(2)
  # initializing logratio matrix (matrix E in paper)
  matrix_E = ExpressionSimilarityMatrix('matrix_E', conc_matrix.gene_factor_dictionary())
  for fi in conc_matrix.get_fnamelist() :
    for fj in conc_matrix.get_fnamelist() :
      # FIXME: no selfregulation yet
      if fi == fj :
        matrix_E.set_element(fi, fj, 0.0)
      elif ((ref_conc_matrix.get_element(fi, fj) == 0.0) and (conc_matrix.get_element(fi, fj) == 0.0)) :
        matrix_E.set_element(fi, fj, 0.0)
      elif (ref_conc_matrix.get_element(fi, fj) == 0.0) :
        matrix_E.set_element(fi, fj, float(sys.maxint))
      elif conc_matrix.get_element(fi, fj) == (0.0) :
        matrix_E.set_element(fi, fj, float(-sys.maxint))
      # calculate logratio (basis 2)
      else :
        try :
          # !!! caution with log(0)
          # positive overflow should not be possible, since the arguments of the logarithm should by then
          # cause an overflow. so if an OverflowError is caught it surely is going to be a negative one.
          # to bypass this problem we chose to set matrix_E[i][j] to some arbitrary floating point number,
          # that makes sure that the according edge is going to be activating.
          # underflow seems to be caught by python internally.
          matrix_E.set_element(fi, fj, (math.log(conc_matrix.get_element(fi, fj) / ref_conc_matrix.get_element(fi, fj)))/log2 )
        except OverflowError :
          overflows = overflows + 1
          matrix_E.set_element(fi, fj, float(-sys.maxint))
  if overflows > 0 :
    errlogger.message('number of overflows, when calculating the logarithm: '+str(overflows), errorlogger_verbositylevel)
  return matrix_E


def overlaps (net_a, net_b) :
  """docu missing
  """
  promoter_only_net_a = []
  promoter_only_net_b = []
  promoter_net_a_and_net_b = []
  #collect equivalent (correct) edges, count original edges and reconstructed edges
  genelist = net_a.gene_names()
  for gene in genelist :
    gene_a = net_a.find_gene(gene)
    gene_b = net_b.find_gene(gene)
    for promoter_element_net_a in gene_a.promoter :
      if not isinstance(promoter_element_net_a, transsys.PromoterElementConstitutive):
        promoter_only_net_a.append(promoter_element_net_a)
    for promoter_element_net_b in gene_b.promoter :
      if not isinstance(promoter_element_net_b, transsys.PromoterElementConstitutive):
        promoter_only_net_b.append(promoter_element_net_b)
    for promoter_element_net_a in gene_a.promoter :
      if isinstance(promoter_element_net_a, transsys.PromoterElementLink) :
        factor_list_net_a = []
        for factor in promoter_element_net_a.factor_list :
          factor_list_net_a.append(factor.name)
        factor_list_net_a.sort()
        for promoter_element_net_b in gene_b.promoter :
          if isinstance(promoter_element_net_b, transsys.PromoterElementLink) :
            factor_list_net_b = []
            for factor in promoter_element_net_b.factor_list :
              factor_list_net_b.append(factor.name)
            factor_list_net_b.sort()
            # Activate
            if isinstance(promoter_element_net_a, transsys.PromoterElementActivate) and isinstance(promoter_element_net_b, transsys.PromoterElementActivate) :
              if factor_list_net_a == factor_list_net_b :
                promoter_net_a_and_net_b.append((promoter_element_net_a, promoter_element_net_b))
                promoter_element_net_b.dot_attributes["color"]="green"
                promoter_element_net_a.dot_attributes["color"]="green"                
                promoter_only_net_a.remove(promoter_element_net_a)
                promoter_only_net_b.remove(promoter_element_net_b)
            # Repress
            elif isinstance(promoter_element_net_a, transsys.PromoterElementRepress) and isinstance(promoter_element_net_b, transsys.PromoterElementRepress) :
              if factor_list_net_a == factor_list_net_b :
                promoter_net_a_and_net_b.append((promoter_element_net_a, promoter_element_net_b))
                promoter_element_net_b.dot_attributes["color"]="green"
                promoter_element_net_a.dot_attributes["color"]="green"                
                promoter_only_net_a.remove(promoter_element_net_a)
                promoter_only_net_b.remove(promoter_element_net_b)
  return [ promoter_net_a_and_net_b, promoter_only_net_a, promoter_only_net_b ]

