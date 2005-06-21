#!/usr/bin/env python

# $Id$

# $Log$
# Revision 1.15  2005/06/21 16:38:04  jtk
# initial success with lsys parser (no serious debugging yet)
#
# Revision 1.14  2005/06/21 09:34:04  jtk
# lsys parser accepts syntax (of one file), semantics still to be implemented
#
# Revision 1.13  2005/06/20 21:10:03  jtk
# in process of implementing lsys parsing capability
#
# Revision 1.12  2005/06/14 18:33:37  jtk
# added comments to factors and genes
#
# Revision 1.11  2005/06/14 09:57:17  jtk
# added some safeguards against bad types in __init__ functions
#
# Revision 1.10  2005/04/14 18:50:17  jtk
# added entropy computation, introduced magic line for transexpr
#
# Revision 1.9  2005/04/10 19:38:00  jtk
# fixed powerlaw_linkist off by one bug
#
# Revision 1.8  2005/04/08 21:07:58  jtk
# fixed random transsys generation to prevent generation of multiple links
#
# Revision 1.7  2005/04/08 20:05:18  jtk
# added indegree and outdegree methods to TranssysProgram
#
# Revision 1.6  2005/04/06 21:04:54  jtk
# minor fixes
#
# Revision 1.5  2005/04/06 09:54:38  jtk
# debugged transdisrupt.py
#
# Revision 1.4  2005/04/05 20:27:37  jtk
# allowed snippets in CyclicSequence, adapted time_series, lsys kludge
#
# Revision 1.3  2005/04/05 10:13:19  jtk
# fixed unresolved_copy for promoter element
#
# Revision 1.2  2005/04/04 21:26:50  jtk
# added merging and unresolved_copy functions to TranssysProgram
#
# Revision 1.1.1.1  2005/03/08 17:12:02  jtk
# new cvs after loss at INB
#
# Revision 1.16  2003/04/03 13:31:30  kim
# removed superfluous print statement in PromoterElementLink.__init__
#
# Revision 1.15  2003/03/20 21:56:54  kim
# added regulated_genes() and regulating_factors()
#
# Revision 1.14  2003/03/12 19:22:12  kim
# fixed wait() for both transexpr child and python child
#     in TranssysInstance.time_series()
#
# Revision 1.13  2003/03/05 13:59:16  kim
# improved error reporting, other minor changes
#
# Revision 1.12  2003/02/28 00:46:14  kim
# improved intensity computation
#
# Revision 1.11  2003/02/27 10:36:48  kim
# added ArrayIntensityFunction class for simulating noisy arrays etc.
#
# Revision 1.10  2003/02/26 19:19:45  schwarzd
# fixed syntax by eliminating tabs...
#
# Revision 1.9  2003/02/26 19:17:24  schwarzd
# fixed problems with dot_attribute members set to None
#
# Revision 1.7  2003/02/26 18:22:53  schwarzd
# systematically implemented dot_attributes output
#
# Revision 1.6  2003/02/17 21:40:31  kim
# Transformed microarr.r into a R library called transarr. Added array_background
#     to regstruct system
#
# Revision 1.5  2003/02/14 16:31:27  kim
# regstruct project is nearing completion. regstruct_transsys for
#     generation of transsys programs and regstruct_transarray for
#     array simulation are (more or less) properly implemented.
#     Analysis with R has been assembled in microarr.r
#
# Revision 1.4  2003/02/05 14:58:05  kim
# minor updates
#
# Revision 1.3  2003/02/04 23:46:58  kim
# added parser, major extension and reworking of RandomTranssysParameters
#
# Revision 1.2  2003/01/30 12:14:33  kim
# initial proper implementation of regstruct_transsys
#
# Revision 1.1  2003/01/28 21:12:05  kim
# initial toolbox assembly
#

# old log from python/jtkmods
# Revision 1.1  2002/03/22 17:23:50  kim
# added transrnd.py and transsys.py to CVS repository
#

import copy
import types
import os
import sys
import popen2
import string
import StringIO
import math
import re

import transrnd


def dot_attribute_string(d) :
  s = ''
  glue = ''
  for k in d.keys() :
    s = s + '%s%s="%s"' % (glue, k, d[k])
    glue = ' '
  return s


def comment_lines(clines, prefix = '') :
  s = ''
  for l in clines :
    s = s + prefix + '// ' + l + '\n'
  return s


# note on expression evaluation: this is not yet properly designed.
# Implementation fragments should be taken to be no more than random notices.
# the transsys_instance argument may be replaced with some more adequate
# evaluation_context thingy...

class ExpressionNode :

  def resolve(self, tp) :
    pass


  def unresolved_copy(self) :
    return copy.deepcopy(self)


class ExpressionNodeValue(ExpressionNode) :

  def __init__(self, v) :
    self.value = float(v)


  def __str__(self) :
    return str(self.value)


  def evaluate(self, transsys_instance_list) :
    return self.value


class ExpressionNodeIdentifier(ExpressionNode) :

  def __init__(self, factor, transsys_label = None) :
    self.factor = factor
    self.transsys_label = transsys_label


  def __str__(self) :
    if self.transsys_label is None :
      return self.factor_name()
    else :
      return '%s.%s' % (self.transsys_label, self.factor_name())


  def factor_name(self) :
    if type(self.factor) is types.StringType :
      return self.factor
    else :
      return self.factor.name


  def resolve(self, tp) :
    self.factor = tp.find_factor(self.factor)


  def unresolved_copy(self) :
    return self.__class__(self.factor_name(), self.transsys_label)
    

  def evaluate(self, transsys_instance) :
    """transsys_instance may either be a single instance or a list of
instances (in a L-transsys context)"""
    if self.transsys_label is None :
      return transsys_instance.factor_concentration[self.factor.index]
    else :
      return transsys_instance[self.transsys_label].factor_concentration[self.factor.index]


class BinaryExpressionNode(ExpressionNode) :

  def __init__(self, op1, op2, operator_sym = None) :
    if not isinstance(op1, ExpressionNode) :
      raise StandardError, 'BinaryExpressionNode::__init__: bad type of op1'
    if not isinstance(op2, ExpressionNode) :
      raise StandardError, 'BinaryExpressionNode::__init__: bad type of op2'
    self.operand1 = op1
    self.operand2 = op2
    self.operator_sym = operator_sym


  def __str__(self) :
    return '(%s %s %s)' % (str(self.operand1), self.operator_sym, str(self.operand2))


  def resolve(self, tp) :
    self.operand1.resolve(tp)
    self.operand2.resolve(tp)


  def unresolved_copy(self) :
    op1 = self.operand1.unresolved_copy()
    op2 = self.operand1.unresolved_copy()
    return self.__class__(op1, op2, self.operator_sym)


class ExpressionNodeAdd(BinaryExpressionNode) :

  def __init__(self, op1, op2) :
    BinaryExpressionNode.__init__(self, op1, op2, '+')


  def evaluate(self, transsys_instance) :
    return self.operand1.evaluate(transsys_instance) + self.operand2.evaluate(transsys_instance)


class ExpressionNodeSubtract(BinaryExpressionNode) :

  def __init__(self, op1, op2) :
    BinaryExpressionNode.__init__(self, op1, op2, '-')


  def evaluate(self, transsys_instance) :
    return self.operand1.evaluate(transsys_instance) - self.operand2.evaluate(transsys_instance)


class ExpressionNodeMult(BinaryExpressionNode) :

  def __init__(self, op1, op2) :
    BinaryExpressionNode.__init__(self, op1, op2, '*')


  def evaluate(self, transsys_instance) :
    return self.operand1.evaluate(transsys_instance) * self.operand2.evaluate(transsys_instance)


class ExpressionNodeDiv(BinaryExpressionNode) :

  def __init__(self, op1, op2) :
    BinaryExpressionNode.__init__(self, op1, op2, '/')


  def evaluate(self, transsys_instance) :
    return self.operand1.evaluate(transsys_instance) / self.operand2.evaluate(transsys_instance)


class ExpressionNodeLower(BinaryExpressionNode) :

  def __init__(self, op1, op2) :
    BinaryExpressionNode.__init__(self, op1, op2, '<')


  def evaluate(self, transsys_instance) :
    return self.operand1.evaluate(transsys_instance) < self.operand2.evaluate(transsys_instance)


class ExpressionNodeLowerEqual(BinaryExpressionNode) :

  def __init__(self, op1, op2) :
    BinaryExpressionNode.__init__(self, op1, op2, '<=')


  def evaluate(self, transsys_instance) :
    return self.operand1.evaluate(transsys_instance) <= self.operand2.evaluate(transsys_instance)


class ExpressionNodeGreater(BinaryExpressionNode) :

  def __init__(self, op1, op2) :
    BinaryExpressionNode.__init__(self, op1, op2, '>')


  def evaluate(self, transsys_instance) :
    return self.operand1.evaluate(transsys_instance) > self.operand2.evaluate(transsys_instance)


class ExpressionNodeGreaterEqual(BinaryExpressionNode) :

  def __init__(self, op1, op2) :
    BinaryExpressionNode.__init__(self, op1, op2, '>=')


  def evaluate(self, transsys_instance) :
    return self.operand1.evaluate(transsys_instance) >= self.operand2.evaluate(transsys_instance)


class ExpressionNodeEqual(BinaryExpressionNode) :

  def __init__(self, op1, op2) :
    BinaryExpressionNode.__init__(self, op1, op2, '==')


  def evaluate(self, transsys_instance) :
    return self.operand1.evaluate(transsys_instance) == self.operand2.evaluate(transsys_instance)


class ExpressionNodeUnequal(BinaryExpressionNode) :

  def __init__(self, op1, op2) :
    BinaryExpressionNode.__init__(self, op1, op2, '!=')


  def evaluate(self, transsys_instance) :
    return self.operand1.evaluate(transsys_instance) != self.operand2.evaluate(transsys_instance)


class ExpressionNodeNot(ExpressionNode) :

  def __init__(self, operand) :
    if not isinstance(operand, ExpressionNode) :
      raise StandardError, 'ExpressionNodeNot::__init__: bad type of operand'
    self.operand = operand


  def __str__(self) :
    return '(!(%s))' % str(self.operand)


  def evaluate(self, transsys_instance) :
    return not self.operand.evaluate(transsys_instance)


  def resolve(self, tp) :
    self.operand.resolve(tp)


  def unresolved_copy(self) :
    op = self.op.unresolved.copy()
    return self.__class__(op)


class ExpressionNodeAnd(BinaryExpressionNode) :

  def __init__(self, op1, op2) :
    BinaryExpressionNode.__init__(self, op1, op2, '&&')


  def evaluate(self, transsys_instance) :
    return self.operand1.evaluate(transsys_instance) and self.operand2.evaluate(transsys_instance)


class ExpressionNodeOr(BinaryExpressionNode) :

  def __init__(self, op1, op2) :
    BinaryExpressionNode.__init__(self, op1, op2, '||')


  def evaluate(self, transsys_instance) :
    return self.operand1.evaluate(transsys_instance) or self.operand2.evaluate(transsys_instance)


class FunctionExpressionNode(ExpressionNode) :

  def __init__(self, operand, funcname = 'undefined_function') :
    for op in operand :
      if not isinstance(op, ExpressionNode) :
        raise StandardError, 'FunctionExpressionNode::__init__: bad type of operand'
    self.operand = operand
    self.funcname = funcname


  def __str__(self) :
    s = self.funcname
    glue = '('
    for op in self.operand :
      s = s + glue + str(op)
      glue = ', '
    s = s + ')'
    return s


  def evaluate(self, transsys_instance) :
    raise StandardError, 'cannot evaluate undefined function'


  def resolve(self, tp) :
    for op in self.operand :
      op.resolve(tp)


  def unresolved_copy(self) :
    u = []
    for op in self.operand :
      u.append(op.unresolved_copy())
    return self.__class__(u)


class ExpressionNodeUniformRandom(FunctionExpressionNode) :

  def __init__(self, op1, op2) :
    FunctionExpressionNode.__init__(self, [op1, op2], 'random')


  def evaluate(self, transsys_instance) :
    arg1 = self.operand[0].evaluate()
    arg2 = self.operand[1].evaluate()
    if arg1 < arg2 :
      return arg1 + (arg2 - arg1) * transsys_instance[0].rng.rnd()
    else :
      return arg2 + (arg1 - arg2) * transsys_instance[0].rng.rnd()


class ExpressionNodeGaussianRandom(FunctionExpressionNode) :

  def __init__(self, op1, op2) :
    FunctionExpressionNode.__init__(self, [op1, op2], 'gauss')


  def evaluate(self, transsys_instance) :
    arg1 = self.operand[0].evaluate()
    arg2 = self.operand[1].evaluate()
    return arg1 + arg2 * transsys_instance[0].rng.gauss()


class PromoterElement :
  pass


class PromoterElementConstitutive(PromoterElement) :

  def __init__(self, expr) :
    if not isinstance(expr, ExpressionNode) :
      raise StandardError, 'ExpressionNodeNot::__init__: bad type of expr'
    self.expression = expr


  def __str__(self) :
    return 'constitutive: %s' % str(self.expression)


  def resolve(self, tp) :
    self.expression.resolve(tp)


  def unresolved_copy(self) :
    return self.__class__(self.expression.unresolved_copy())
    

  def write_dot_edge(self, f, target_name, dot_parameters, transsys) :
    pass # constitutive expression of genes should probably be indicated
         # by node shape / border or the like...
  

class PromoterElementLink(PromoterElement) :

# expr1, expr2, ... should perhaps become a list at some stage...
  def __init__(self, expr1, expr2, factor_list, dot_attrs = None) :
    # print factor_list
    if not isinstance(expr1, ExpressionNode) :
      raise StandardError, 'ExpressionNodeNot::__init__: bad type of expr1'
    if not isinstance(expr2, ExpressionNode) :
      raise StandardError, 'ExpressionNodeNot::__init__: bad type of expr2'
    self.factor_list = factor_list[:]
    if dot_attrs is None :
      self.dot_attributes = {}
    else :
      self.dot_attributes = copy.deepcopy(dot_attrs)
    self.expression1 = expr1
    self.expression2 = expr2


  def __str__(self) :
    s = ''
    glue = ''
    for f in self.factor_list :
      if type(f) is types.StringType :
        fn = f
      else :
        fn = f.name
      s = s + glue + fn
      glue = ' + '
      return s


  def resolve(self, tp) :
    fi_list = []
    for f_str in self.factor_list :
      fi_list.append(tp.find_factor(f_str))
    self.factor_list = fi_list
    self.expression1.resolve(tp)
    self.expression2.resolve(tp)


  def unresolved_copy(self) :
    flist = []
    for f in self.factor_list :
      if type(f) is types.StringType :
        flist.append(f)
      else :
        flist.append(f.name)
    return self.__class__(self.expression1.unresolved_copy(), self.expression2.unresolved_copy(), flist, self.dot_attributes)


  def write_dot_edge(self, f, target_name, display_factors, arrowhead, transsys) :

    def factor_name(f) :
      if type(f) is types.StringType :
        return f
      else :
        return f.name

    if display_factors :
      for factor in self.factor_list :
        f.write('  %s -> %s' % (factor_name(factor), target_name))
        f.write(' [arrowhead="%s"' % arrowhead)
        f.write(' %s];\n' % dot_attribute_string(self.dot_attributes))
    else :
      for factor in self.factor_list :
        gene_list = transsys.encoding_gene_list(factor_name(factor))
        for gene in gene_list :
          f.write('  %s -> %s' % (gene.name, target_name))
          f.write(' [arrowhead="%s"' % arrowhead)
          f.write(' %s];\n' % dot_attribute_string(self.dot_attributes))


class PromoterElementActivate(PromoterElementLink) :

  def __init__(self, expr1, expr2, factor_list, dot_attrs = None) :
    PromoterElementLink.__init__(self, expr1, expr2, factor_list, dot_attrs)


  def __str__(self) :
    return '%s: activate(%s, %s)' % (PromoterElementLink.__str__(self), str(self.expression1), str(self.expression2))


  def write_dot_edge(self, f, target_name, dot_parameters, transsys) :
    PromoterElementLink.write_dot_edge(self, f, target_name, dot_parameters.display_factors, dot_parameters.activate_arrowhead, transsys)


class PromoterElementRepress(PromoterElementLink) :

  def __init__(self, expr1, expr2, factor_list, dot_attrs = None) :
    PromoterElementLink.__init__(self, expr1, expr2, factor_list, dot_attrs)


  def __str__(self) :
    return '%s: repress(%s, %s)' % (PromoterElementLink.__str__(self), str(self.expression1), str(self.expression2))


  def write_dot_edge(self, f, target_name, dot_parameters, transsys) :
    PromoterElementLink.write_dot_edge(self, f, target_name, dot_parameters.display_factors, dot_parameters.repress_arrowhead, transsys)


class Factor :

  def __init__(self, name, decay_expr = None, diffusibility_expr = None, dot_attrs = None) :
    if not isinstance(decay_expr, ExpressionNode) :
      raise StandardError, 'Factor::__init__: bad decay_expr type'
    if not isinstance(diffusibility_expr, ExpressionNode) :
      raise StandardError, 'Factor::__init__: bad diffusibility_expr type'
    self.name = name
    self.decay_expression = decay_expr
    if decay_expr is None :
      self.decay_expression = ExpressionNodeValue(1.0)
    self.diffusibility_expression = diffusibility_expr
    if diffusibility_expr is None :
      self.diffusibility_expression = ExpressionNodeValue(1.0)
    self.comments = []
    if dot_attrs is None :
      self.dot_attributes = {}
    else :
      self.dot_attributes = copy.deepcopy(dot_attrs)


  def __str__(self) :
    return """  factor %s
  {
%s    decay: %s;
    diffusibility: %s;
  }
""" % (self.name, comment_lines(self.comments, '    '), str(self.decay_expression), str(self.diffusibility_expression))


  def resolve(self, tp) :
    self.decay_expression.resolve(tp)
    self.diffusibility_expression.resolve(tp)


  def unresolved_copy(self) :
    return self.__class__(self.name, self.decay_expression.unresolved_copy(), self.diffusibility_expression.unresolved_copy(), self.dot_attributes)


  def write_dot_node(self, f, dot_parameters, transsys) :
    f.write('  %s' % self.name)
    if len(self.dot_attributes) > 0 :
      f.write(' [%s]' % dot_attribute_string(self.dot_attributes))
    f.write(';\n')


class Gene :

  def __init__(self, name, product, promoter = None, dot_attrs = None) :
    self.name = name
    self.product = product
    if promoter is None :
      self.promoter = []
    else :
      for p in promoter :
        if not isinstance(p, PromoterElement) :
          raise StandardError, 'Gene::__init__: bad element in promoter'
      self.promoter = copy.deepcopy(promoter)
    self.comments = []
    if dot_attrs is None :
      self.dot_attributes = {}
    else :
      self.dot_attributes = copy.deepcopy(dot_attrs)


  def __str__(self) :
    s = '  gene %s\n' % self.name
    s = s + '  {\n'
    s = s + comment_lines(self.comments, '    ')
    s = s + '    promoter\n'
    s = s + '    {\n'
    for p in self.promoter :
      s = s + '      %s;\n' % str(p)
    s = s + '    }\n'
    s = s + '    product\n'
    s = s + '    {\n'
    s = s + '      default: %s;\n' % self.product_name()
    s = s + '    }\n'
    s = s + '  }\n'
    return s


  def product_name(self) :
    if type(self.product) is types.StringType :
      return self.product
    else :
      return self.product.name


  def resolve(self, tp) :
    for pe in self.promoter :
      pe.resolve(tp)
    self.product = tp.find_factor(self.product)


  def unresolved_copy(self) :
    unresolved_promoter = []
    for pe in self.promoter :
      unresolved_promoter.append(pe.unresolved_copy())
    return self.__class__(self.name, self.product_name(), unresolved_promoter, self.dot_attributes)


  def write_dot_node(self, f, dot_parameters, transsys) :
    f.write('  %s' % self.name)
    if len(self.dot_attributes) > 0 :
      f.write(' [%s]' % dot_attribute_string(self.dot_attributes))
    f.write(';\n')


  def write_dot_edges(self, f, dot_parameters, transsys) :
    """write edges directed towards this gene"""
    for p in self.promoter :
      p.write_dot_edge(f, self.name, dot_parameters, transsys)


  def regulating_factors(self) :
    """return a list of factors that regulate expression of this gene"""
    rfactors = []
    for p in self.promoter :
      if isinstance(p, PromoterElementLink) :
        for f in p.factor_list :
          if f not in rfactors :
            rfactors.append(f)
    return rfactors


class TranssysProgram :
  """Note: unresolved_copy returns a TranssysProgram instance, also for
classes derived from TranssysProgram. If unresolved_copy for such classes
is required to return the derived class, these classes will have to
re-implement this method.
  """

  def __init__(self, name, factor_list = None, gene_list = None, resolve = True) :
    self.name = name
    if factor_list is None :
      self.factor_list = []
    else :
      for f in factor_list :
        if not isinstance(f, Factor) :
          raise StandardError, 'TranssysProgram::__init__: unsuitable element in factor_list'
      # FIXME: I don't think this deepcopy works usefully for a list of resolved factors
      self.factor_list = copy.deepcopy(factor_list)
    if gene_list is None :
      self.gene_list = []
    else :
      for g in gene_list :
        if not isinstance(g, Gene) :
          raise StandardError, 'TranssysProgram::__init__: unsuitable element in gene_list'
      # FIXME: same deepcopy problem as above
      self.gene_list = copy.deepcopy(gene_list)
    self.comments = []
    if resolve :
      self.resolve()


  def __str__(self) :
    s = 'transsys %s\n' % self.name
    s = s + comment_lines(self.comments)
    s = s + '{\n'
    glue = ''
    for f in self.factor_list :
      s = s + glue + str(f)
      glue = '\n'
    for g in self.gene_list :
      s = s + glue + str(g)
      glue = '\n'
    s = s + '}\n'
    return s


  def num_factors(self) :
    return len(self.factor_list)


  def num_genes(self) :
    return len(self.gene_list)


  def factor_names(self) :
    fnames = []
    for f in self.factor_list :
      fnames.append(f.name)
    return fnames


  def gene_names(self) :
    gnames = []
    for g in self.gene_list :
      gnames.append(g.name)
    return gnames


  def find_factor_index(self, f_name) :
    for i in xrange(len(self.factor_list)) :
      if self.factor_list[i].name == f_name :
        return i
    return -1


  def find_factor(self, f_name) :
    if type(f_name) is not types.StringType :
      f_name = f_name.name
    i = self.find_factor_index(f_name)
    if i < 0 :
      raise StandardError, 'TranssysProgram::find_factor: no factor "%s" in transsys %s' % (f_name, self.name)
    return self.factor_list[i]


  def find_gene_index(self, g_name) :
    for i in xrange(len(self.gene_list)) :
      if self.gene_list[i].name == g_name :
        return i
    return -1


  def find_gene(self, g_name) :
    if type(g_name) is not types.StringType :
      g_name = g_name.name
    i = self.find_gene_index(g_name)
    if i < 0 :
      raise StandardError, 'TranssysProgram::find_gene: no gene "%s" in transsys %s' % (g_name, self.name)
    return self.gene_list[i]


  def resolve(self) :
    """The resolve() method replaces all gene and factor specifications
by strings with references to the gene object or factor object, respectively.
Strings which do not properly resolve trigger an exception. resolve()
should be called whenever alterations are done to a TranssysProgram instance in
which strings are used to specify genes or factors. Only after resolving,
gene names and factor names can be altered without changing the network."""
    for factor in self.factor_list :
      factor.resolve(self)
    for gene in self.gene_list :
      gene.resolve(self)


  def unresolved_copy(self) :
    factor_list = []
    for f in self.factor_list :
      factor_list.append(f.unresolved_copy())
    gene_list = []
    for g in self.gene_list :
      gene_list.append(g.unresolved_copy())
    # since we don't know what might be derived from TranssysProgram and what
    # constructors they may have, we don't use the self.__class__() constructor
    # call here.
    return TranssysProgram(self.name, factor_list, gene_list, False)


  def encoding_gene_list(self, factor_name) :
    """return a list of all genes which encode factor named factor_name"""
    # print 'encoding_gene_list("%s")' % factor_name
    glist = []
    for gene in self.gene_list :
      if gene.product_name() == factor_name :
        glist.append(gene)
    return glist


  def write_dot(self, f, dot_parameters) :
    """write transsys program as graph in the dot language"""
    f.write('digraph %s\n' % self.name)
    f.write('{\n')
    for gene in self.gene_list :
      gene.write_dot_node(f, dot_parameters, self)
    if dot_parameters.display_factors :
      for factor in self.factor_list :
        factor.write_dot_node(f, dot_parameters, self)
    for gene in self.gene_list :
      # print 'write_dot: edges for gene "%s"' % gene.name
      gene.write_dot_edges(f, dot_parameters, self)
    f.write('}\n')


  def dot_positions_circle(self, x0, y0, radius) :
    for i in xrange(len(self.gene_list)) :
      angle = 2.0 * math.pi * i / len(self.gene_list)
      self.gene_list[i].dot_attributes['pos'] = '%f,%f!' % (x0 + radius * math.cos(angle), y0 + radius * math.sin(angle))


  def regulated_genes(self, factor) :
    if type(factor) is types.StringType :
      i = self.find_factor_index(factor)
      if i == -1 :
        raise StandardError, 'TranssysProgram.regulated_genes: transsys program "%s" does not contain factor "%s"' % (self.name, factor)
      factor = self.factor_list[i]
    if not isinstance(factor, Factor) :
      raise TypeError, 'TranssysProgram.regulated_genes: factor argument not a Factor instance'
    rgenes = []
    for gene in self.gene_list :
      if factor in gene.regulating_factors() :
        rgenes.append(gene)
    return rgenes


  def indegree_list(self) :
    dlist = []
    for gene in self.gene_list :
      dlist.append(len(gene.regulating_factors()))
    return dlist


  def outdegree_list(self) :
    dlist = []
    for gene in self.gene_list :
      factor = gene.product
      dlist.append(len(self.regulated_genes(factor)))
    return dlist


  def merge(self, other) :
    other_unresolved = other.unresolved_copy()
    for f in other_unresolved.factor_list :
      if self.find_factor_index(f.name) == -1 :
        self.factor_list.append(f)
      else :
        print 'factor %s already present -- skipping' % f.name
    for g in other_unresolved.gene_list :
      if self.find_gene_index(g.name) == -1 :
        self.gene_list.append(g)
      else :
        print 'gene %s already present -- skipping' % g.name
    self.comments.append('merged with transsys %s' % other.name)


class GraphicsPrimitive :

  def __init__(self) :
    pass


class GraphicsPrimitiveMove(GraphicsPrimitive) :

  def __init__(self, expression) :
    self.expression = expression


  def __str__(self) :
    return 'move(%s);' % str(self.expression)


class GraphicsPrimitivePush(GraphicsPrimitive) :

  def __init__(self) :
    pass


  def __str__(self) :
    return 'push();'


class GraphicsPrimitivePop(GraphicsPrimitive) :

  def __init__(self) :
    pass


  def __str__(self) :
    return 'pop();'


class GraphicsPrimitiveTurn(GraphicsPrimitive) :

  def __init__(self, expression) :
    self.expression = expression


  def __str__(self) :
    return 'turn(%s);' % str(self.expression)
    

class GraphicsPrimitiveRoll(GraphicsPrimitive) :

  def __init__(self, expression) :
    self.expression = expression


  def __str__(self) :
    return 'roll(%s);' % str(self.expression)
    

class GraphicsPrimitiveBank(GraphicsPrimitive) :

  def __init__(self, expression) :
    self.expression = expression


  def __str__(self) :
    return 'bank(%s);' % str(self.expression)


class GraphicsPrimitiveSphere(GraphicsPrimitive) :

  def __init__(self, expression) :
    self.expression = expression


  def __str__(self) :
    return 'sphere(%s);' % str(self.expression)


class GraphicsPrimitiveCylinder(GraphicsPrimitive) :

  def __init__(self, diameterExpression, lengthExpression) :
    self.diameterExpression = diameterExpression
    self.lengthExpression = lengthExpression


  def __str__(self) :
    return 'cylinder(%s, %s);' % (str(self.diameterExpression), str(self.lengthExpression))


class GraphicsPrimitiveBox(GraphicsPrimitive) :

  def __init__(self, xExpression, yExpression, zExpression) :
    self.xExpression = xExpression
    self.yExpression = yExpression
    self.zExpression = zExpression


  def __str__(self) :
    return 'box(%s, %s, %s);' % (str(self.xExpression), str(self.yExpression), str(self.zExpression))


class GraphicsPrimitiveColor(GraphicsPrimitive) :

  def __init__(self, redExpression, greenExpression, blueExpression) :
    self.redExpression = redExpression
    self.greenExpression = greenExpression
    self.blueExpression = blueExpression


  def __str__(self) :
    return 'color(%s, %s, %s);' % (str(self.redExpression), str(self.greenExpression), str(self.blueExpression))


class Symbol :

  def __init__(self, name, transsys, gpl = None) :
    self.name = name
    self.transsys = transsys
    self.graphicsPrimitiveList = gpl


  def __str__(self) :
    if self.transsys is None :
      return 'symbol %s' % self.name
    else :
      if isinstance(self.transsys, TranssysProgram) :
        return 'symbol %s(%s)' % (self.name, self.transsys.name)
      elif type(self.transsys) is types.StringType :
        return 'symbol %s(%s)' % (self.name, self.transsys)


  def graphics_string(self) :
    s = ''
    for g in self.graphicsPrimitiveList :
      s = s + '      %s\n' % str(g)
    return s
      


class Assignment :

  def __init__(self, transsys, factor, value) :
    self.transsys = transsys
    self.factor = factor
    self.value = value


  def __str__(self) :
    return '%s = %s' % (self.factor, str(self.value))


class LhsSymbol :

  def __init__(self, symbol, transsys_name) :
    self.symbol = symbol
    self.transsys_name = transsys_name


  def __str__(self) :
    if self.transsys_name :
      return '%s(%s)' % (self.symbol.name, self.transsys_name)
    else :
      return self.symbol_name


class ProductionElement :

  def __init__(self, symbol, template_label, assignments) :
    self.symbol = symbol
    self.template_label = template_label
    self.assignments = assignments


  def __str__(self) :
    tl = ''
    if self.template_label is not None :
      tl = 'transsys %s: ' % self.template_label
    as = ''
    glue = ''
    for a in self.assignments :
      as = as + glue + str(a)
      glue = ', '
    if tl != '' or as != '' :
      return '%s(%s%s)' % (self.symbol.name, tl, as)
    else :
      return self.symbol.name


class Rule :

  def __init__(self, name, lhs, condition, rhs) :
    self.name = name
    self.lhs = lhs
    self.condition = condition
    self.rhs = rhs


  def __str__(self) :
    s = '  rule %s\n' % self.name
    s = s + '  {\n'
    s = s + '    '
    for l in self.lhs :
      s = s + ' %s' % str(l)
    if self.condition is not None :
      s = s + ':\n'
      s = s + '    %s' % str(self.condition)
    s = s + '\n'
    s = s + '    -->\n'
    for r in self.rhs :
      s = s + '    %s\n' % str(r)
    s = s + '  }\n'
    return s

class LsysProgram :

  def __init__(self, name, symbols, axiom, diffusionrange, rules) :
    self.name = name
    self.symbols = symbols
    self.axiom = axiom
    self.diffusionrange = diffusionrange
    self.rules = rules


  def __str__(self) :
    s = 'lsys %s\n' % self.name
    s = s + '{\n'
    s = s + '  diffusionrange: %d;\n' % self.diffusionrange
    for sym in self.symbols :
      s = s + '  %s;\n' % str(sym)
    s = s + '  axiom'
    for a in self.axiom :
      s = s + ' %s;\n' % str(a)
    for rule in self.rules :
      s = s + str(rule)
    s = s + '  graphics\n'
    s = s + '  {\n'
    for sym in self.symbols :
      s = s + sym.graphics_string()
    s = s + '  }\n'
    s = s + '}\n'
    return s

class TranssysInstance :

  # the 'magic first line' of transexpr output. If this is not found
  # an error is triggered, as a safeguard against nonsense due to
  # incompatiblilities between TranssysInstance and the transexpr
  # program invoked by the time_series method.

  magic = '# transsys expression records 1.0'

  def __init__(self, tp) :
    if not isinstance(tp, TranssysProgram) :
      raise StandardError, 'TranssysInstance::__init__: unsuitable type of tp (no TranssysProgram)'
    self.transsys_program = tp
    self.factor_concentration = [0.0] * self.transsys_program.num_factors()
    # FIXME (conceptual issue):
    # the stddev and entropy members are optional and they do not pertain
    # to an individual instance (obviously). As long as TranssysInstance
    # instances are only used for getting levels from transexpr, this is
    # not too much of a problem, but eventually, instances and statistical
    # results should be handled by separated classes.
    self.factor_concentration_stddev = [0.0] * self.transsys_program.num_factors()
    self.factor_concentration_entropy = [0.0] * self.transsys_program.num_factors()


  def __str__(self) :
    s = 'transsys instance of %s\n' % self.transsys_program.name
    glue = ''
    for f in xrange(len(self.factor_concentration)) :
      s = s + '  %s: %f\n' % (self.transsys_program.factor_list[f].name, self.factor_concentration[f])
    return s


  def set_uniform_random(self, rng, c_min, c_max) :
    self.factor_concentration = []
    c_range = c_max - c_min
    for i in xrange(len(self.transsys_program.factor_list)) :
      self.factor_concentration.append(c_min + rng.rnd() * c_range)


  def perturb(self, perturb_func) :
    for i in xrange(len(self.transsys_program.factor_list)) :
      self.factor_concentration[i] = perturb_func(self.factor_concentration[i], self.transsys_program.factor_list[i].name)


  def perturbed_copy(self, perturb_func) :
    pc = TranssysInstance(self.transsys_program)
    pc.factor_concentration = copy.deepcopy(self.factor_concentration)
    pc.perturb(perturb_func)
    return pc


  def write_gnuplot(self, f, time_step) :
    f.write('%d' % time_step)
    for c in self.factor_concentration :
      f.write(' %f' % c)
    f.write('\n')


  def time_series(self, num_timesteps, sampling_period = 1, lsys_lines = None, lsys_symbol = None) :
    # FIXME: lsys_lines is a kludge that should be removed once lsys parsing
    # is available.
    cmd = 'transexpr -n %d -d %d' % (num_timesteps, sampling_period)
    glue = ''
    s = ''
    for x in self.factor_concentration :
      s = s + '%s%f' % (glue, x)
      glue = ' '
    cmd = cmd + ' -F \'%s\'' % s
    if lsys_lines is not None :
      cmd = cmd + ' -l -t \'%s\'' % self.transsys_program.name
    if lsys_symbol is not None :
      cmd = cmd + ' -y \'%s\'' % lsys_symbol
    # sys.stderr.write('%s\n' % cmd)
    p = popen2.Popen3(cmd, 1)
    sys.stdout.flush()
    sys.stderr.flush()
    pid = os.fork()
    if pid == 0 :
      p.fromchild.close()
      p.tochild.write(str(self.transsys_program))
      if lsys_lines is not None :
        for l in lsys_lines :
          p.tochild.write('%s\n' % l)
      #sys.stdout.write(str(self.transsys_program))
      p.tochild.close()
      sys.exit()
    p.tochild.close()
    tseries_dict = {}
    line = p.fromchild.readline()
    magic = line.strip()
    if magic != self.magic :
      sys.stderr.write('bad "magic" first line: got "%s", expected "%s"\n' % (magic, self.magic))
      raise StandardError, 'TranssysInstance.time_series: bad header (incompatible transexpr version?)'
    while line :
      line = string.strip(line)
      if line :
        if line[0] != '#' :
          l = string.split(line)
          if len(l) != 3 * self.transsys_program.num_factors() + 2 :
            raise StandardError, 'TranssysInstance.time_series: bad line format: expected %d words, got %d' % (3 * self.transsys_program.num_factors() + 2, len(l))
          t = string.atoi(l[0])
          ti = TranssysInstance(self.transsys_program)
          for i in xrange(len(self.transsys_program.factor_list)) :
            ti.factor_concentration[i] = float(l[3 * i + 2])
            ti.factor_concentration_stddev[i] = float(l[3 * i + 3])
            ti.factor_concentration_stddev[i] = float(l[3 * i + 4])
          tseries_dict[t] = ti
      line = p.fromchild.readline()
    p.fromchild.close()
    status = p.wait()
    if status != 0 :
      errmsg = p.childerr.readline()
      while errmsg :
        sys.stderr.write(errmsg)
        errmsg = p.childerr.readline()
      raise StandardError, 'TranssysInstance::time_series: transexpr exit status %d ("%s")' % (status, errmsg.strip())
    os.wait()
    return tseries_dict


class DotParameters :

  def __init__(self) :
    self.display_genes = 1
    self.display_factors = 1
    self.activate_arrowhead = 'normal'
    self.repress_arrowhead = 'tee'


class CyclicSequence :
  """A convenience class for use with RandomTranssysParameters. Implements an
object from which numerical values can be pulled out indefinitely using the
nextval() method. Values are taken sequentially from a list, repeating from
the start when the end is reached.

One could subclass from this class in order to design other value
generators. But then, this might be of limited use, considering that
all methods will have to be overridden anyway."""

  def __init__(self, l, start = None) :
    if len(l) == 0 :
      raise StandardError, 'CyclicSequence::__init__: cannot init with empty list'
    self.l = copy.deepcopy(l)
    if start is None :
      self.i = len(self.l) - 1
    else :
      self.i = start % len(self.l)


  def __str__(self) :
    if len(self.l) == 1 :
      return str(self.l[0])
    else :
      return str(self.l)


  def from_string(self, s) :
    # fixme (?) this implementation allows a trailing comma
    # ... but then, python does this too in several contexts...
    quoted_string_re = re.compile(' *\'([^\']+)\'( *, *)?')
    string_re = re.compile(' *([^,]*)(, *)?')
    s = s.strip()
    l = []
    while s :
      m = quoted_string_re.match(s)
      if m is None :
        m = string_re.match(s)
        if m is None :
          raise StandardError, 'CyclicSequence::from_string: malformed string "%s"' % s
      l.append(m.group(1))
      s = s[len(m.group(0)):]
    self.i = len(l) - 1
    self.l = l


  def currval(self) :
    return self.l[self.i]


  def nextval(self) :
    self.i = (self.i + 1) % len(self.l)
    return self.l[self.i]


  def reset(self) :
    self.i = len(self.l) - 1


class RandomTranssysParameters :
  """A factory for generating "random" transsys programs.
The name "RandomTranssysParameters" really is a historical accident,
this class started out as a dumb class with just a bunch of member
variables for conveniently storing an ever growing bunch of parameters
for random program generation. But finally, the implementation of
classes for the transsys program components has permitted writing a
reasonable generate_program() method.

Usage of this class: Instantiate with a random seed.
Set the member variables by direct assignment as you
see fit. Finally, call generate_program() to pull out endless
variations of transsys programs..."""

  savefile_magic = 'RandomTranssysParameters-1.1'

  def __init__(self, rndseed = 1) :
    self.set_seed(rndseed)
    self.serial = 0
    self.topology = None
    self.constitutive = CyclicSequence([0.0])
    self.km_activation = CyclicSequence([0.0])
    self.km_repression = CyclicSequence([0.0])
    self.vmax_activation = CyclicSequence([0.0])
    self.vmax_repression = CyclicSequence([0.0])
    self.decay = CyclicSequence([0.0])
    self.diffusibility = CyclicSequence([0.0])


  def __str__(self) :
    """display transsys program generation parameters in a string, in a
format that can later on be loaded into a RandomTranssysParameter instance
using the parse() method.

Important note: The state of the RandomTranssysParameters instance is
not completely represented; the RNG state as well as the current value positions
in the CyclicSequence instances are lost. Reproducible transsys generation
after parsing a saved file should work if __str__() is used before generating
any programs, and if no fancy tricks with CyclicSequence instances are
played."""
    s = 'topology: %s\n' % self.topology
    if self.topology == 'random_nk' :
      s = s + 'n: %d\n' % self.n
      s = s + 'k: %d\n' % self.k
    elif self.topology == 'random_uniform' :
      s = s + 'n: %d\n' % self.n
      s = s + 'num_edges: %d\n' % self.num_edges
    elif self.topology == 'random_powerlaw' :
      s = s + 'n: %d\n' % self.n
      s = s + 'num_edges: %d\n' % self.num_edges
      s = s + 'power_exp: %g\n' % self.power_exp
      s = s + 'power_base: %g\n' % self.power_base
    elif self.topology == 'linklist' :
      for in_list in self.linklist :
        s = s + 'in_list: %s\n' % str(in_list)
    else :
      raise StandardError, 'RandomTranssysParameters::write: unknown topology type "%s"' % self.topology
    s = s + 'topology: end\n'
    s = s + 'constitutive: %s\n' % str(self.constitutive)
    s = s + 'km_activation: %s\n' % str(self.km_activation)
    s = s + 'km_repression: %s\n' % str(self.km_repression)
    s = s + 'vmax_activation: %s\n' % str(self.vmax_activation)
    s = s + 'vmax_repression: %s\n' % str(self.vmax_repression)
    s = s + 'decay: %s\n' % str(self.decay)
    s = s + 'diffusibility: %s\n' % str(self.diffusibility)
    s = s + 'rndseed: %d\n' % self.rng.seed
    return s


  def write(self, f) :
    """write the transsys program generation parameters in file f, in a
format that can later on be loaded into a RandomTranssysParameter instance
using the parse() method.

Important note: The state of the RandomTranssysParameters instance is
not completely saved; the RNG state as well as the current value positions
in the CyclicSequence instances are lost. Reproducible transsys generation
after parsing a saved file should work if write() is used before generating
any programs, and if no fancy tricks with CyclicSequence instances are
played."""
    f.write('%s\n' % self.savefile_magic)
    f.write(str(self))
    f.write('\n')


  def parse(self, f) :
    """set parameter members by parsing from file f. This should work with files
generated by the write() method, but of course, files can also be manually written.

Notes on the format: The format is of the simple <identifier>: <value> type.
All parameters must be specified in the order in which they are written by the
write() method. Permuting this order was deliberately not permitted because
that would make it more easy to omit parameters, and more difficult to check
against this. Incomplete specifications are considered dangerous because
unspecified parameters may end up in unclear defaults or even undefined states,
so it's best to require explicit specification of all parameters."""

    def parse_int(f, label) :
      r = '%s\\s*:\\s*([0-9]+)' % label
      line = f.readline()
      m = re.match(r, line.strip())
      if m is None :
        raise StandardError, 'RandomTranssysParameters::parse: failed to obtain int "%s" in "%s"' % (label, line.strip())
      return int(m.group(1))

    def parse_float(f, label) :
      r = '%s\\s*:\\s*([+-]?([0-9]+(\\.[0-9]+)?)|(\\.[0-9]+)([Ee][+-]?[0-9]+)?)' % label
      line = f.readline()
      m = re.match(r, line.strip())
      if m is None :
        raise StandardError, 'RandomTranssysParameters::parse: failed to obtain float "%s" in "%s"' % (label, line.strip())
      return float(m.group(1))

    def parse_cyclicseq(f, label) :
      r = '%s\\s*:\\s*(.*)' % label
      line = f.readline()
      m = re.match(r, line.strip())
      if m is None :
        raise StandardError, 'RandomTranssysParameters::parse: failed to obtain CyclicSequence "%s"' % (label, line.strip())
      c = CyclicSequence([0.0])
      c.from_string(m.group(1))
      return c


    line = f.readline()
    if string.strip(line) != self.savefile_magic :
      raise StandardError, 'RandomTranssysParameters::parse: bad magic "%s"' % string.strip(line)

    line = f.readline()
    m = re.match('topology\\s*:\\s*([A-Za-z_]+)', line)
    if m is None :
      raise StandardError, 'RandomTranssysParameters::parse: failed to obtain topology'
    # print 'topology type: "%s"' % m.group(1)
    self.topology = m.group(1)
    if self.topology == 'random_nk' :
      self.n = parse_int(f, 'n')
      self.k = parse_int(f, 'k')
    elif self.topology == 'random_uniform' :
      self.n = parse_int(f, 'n')
      self.num_edges = parse_int(f, 'num_edges')
    elif self.topology == 'random_powerlaw' :
      self.n = parse_int(f, 'n')
      self.num_edges = parse_int(f, 'num_edges')
      self.power_exp = parse_float(f, 'power_exp')
      self.power_base = parse_float(f, 'power_base')
    elif self.topology == 'linklist' :
      self.topology = 'linklist'
      in_list_re = re.compile('in_list\\s*:\\s*\\[([-0-9,\\s]*)\\]')
      csvsplit_re = re.compile('\\s*,\\s*')
      linklist = []
      line = f.readline()
      m = in_list_re.match(line)
      while m :
        l = csvsplit_re.split(m.group(1))
        linklist.append([])
        for x in l :
          if len(x) > 0 :
            linklist[-1].append(int(x))
        line = f.readline()
        m = in_list_re.match(line)
      print 'parsed topology linklist:', str(linklist)
      self.topology_linklist = linklist
    else :
      raise StandardError, 'RandomTranssysParameters::parse: unknown topology type %s' % self.topology
    line = f.readline()
    m = re.match('topology\\s*:\\s*end', line.strip())
    if m is None :
      raise StandardError, 'RandomTranssysParameters::parse: topology improperly terminated by "%s"' % line.strip()
    self.constitutive = parse_cyclicseq(f, 'constitutive')
    self.km_activation = parse_cyclicseq(f, 'km_activation')
    self.km_repression = parse_cyclicseq(f, 'km_repression')
    self.vmax_activation = parse_cyclicseq(f, 'vmax_activation')
    self.vmax_repression = parse_cyclicseq(f, 'vmax_repression')
    self.decay = parse_cyclicseq(f, 'decay')
    self.diffusibility = parse_cyclicseq(f, 'diffusibility')
    rndseed = parse_int(f, 'rndseed')
    self.set_seed(rndseed)
    line = f.readline()
    if line.strip() != '' :
      raise StandardError, 'RandomTranssysParameters::parse: trailing garbage line "%s"' % line.strip()


  def set_seed(self, rndseed = 1) :
    self.rng = transrnd.transrnd(rndseed)


  def set_constitutive(self, v) :
    if type(v) is types.ListType :
      self.constitutive = CyclicSequence(v)
    else :
      self.constitutive = CyclicSequence([v])


  def set_km_activation(self, v) :
    if type(v) is types.ListType :
      self.km_activation = CyclicSequence(v)
    else :
      self.km_activation = CyclicSequence([v])


  def set_km_repression(self, v) :
    if type(v) is types.ListType :
      self.km_repression = CyclicSequence(v)
    else :
      self.km_repression = CyclicSequence([v])


  def set_vmax_activation(self, v) :
    if type(v) is types.ListType :
      self.vmax_activation = CyclicSequence(v)
    else :
      self.vmax_activation = CyclicSequence([v])


  def set_vmax_repression(self, v) :
    if type(v) is types.ListType :
      self.vmax_repression = CyclicSequence(v)
    else :
      self.vmax_repression = CyclicSequence([v])


  def set_decay(self, v) :
    if type(v) is types.ListType :
      self.decay = CyclicSequence(v)
    else :
      self.decay = CyclicSequence([v])


  def set_diffusibility(self, v) :
    if type(v) is types.ListType :
      self.diffusibility = CyclicSequence(v)
    else :
      self.diffusibility = CyclicSequence([v])


  def check_linklist(self) :
    for n in self.linklist :
      for i in n :
        if (i < 0 and -i - 1 >= len(t)) or i >= len(t) :
          raise StandardError, 'RandomTranssysParameters::check_linklist: illegal index %d' % i


  def num_genes(self) :
    if self.topology == 'random_nk' or self.topology == 'random_uniform' or self.topology == 'random_powerlaw' :
      return self.n
    elif self.topology == 'linklist' :
      return len(self.linklist)
    else :
      raise StandardError, 'RandomTranssysParameters::num_genes: unknown topology type "%s"' % self.topology


  def num_links(self) :
    if self.topology == 'random_nk' :
      return self.n * self.k
    elif self.topology == 'random_uniform' or self.topology == 'random_powerlaw' :
      return self.num_edges
    elif self.topology == 'linklist' :
      nl = 0
      for x in self.linklist :
        nl = nl + len(x)
      return nl
    else :
      raise StandardError, 'RandomTranssysParameters::num_links: unknown topology type "%s"' % self.topology


  def write_parameter_comments(self, f) :
    for s in str(self).split('\n') :
      f.write('// %s\n' % s)


  def generate_transsys(self, name) :
    """Generate a random transsys program. A large variety of parameters controlling
the "random" generation of the transsys program can be set in the
RandomTranssysParameters instance rtp passed as the third argument to this
function."""

    def factor_name(i) :
      return 'f%04d' % i

    def gene_name(i) :
      return 'g%04d' % i


    def random_nk_linklist(n, k, rng) :
      if rng is None :
        raise StandardError, 'RandomTranssysParameters::generate_transsys: cannot generate random NK linklist without RNG'
      linklist = []
      for i in xrange(n) :
        linklist.append([])
        in_list = range(n)
        for j in xrange(k) :
          x = rng.random_range(len(in_list))
          t = in_list[x]
          del in_list[x]
          if rng.random_range(2) :
            linklist[i].append(-t -1)
          else :
            linklist[i].append(t)
      return linklist

    def random_powerlaw_linklist(num_genes, num_edges, power_exp, power_base, rng) :
      linklist = []
      rw = []
      for i in xrange(num_genes) :
        linklist.append([])
        try :
          x = (power_base * i)**power_exp
        except OverflowError :
          if a < 0.0 :
            x = 0.0
          else :
            raise StandardError, 'RandomTranssysParameters::generate_transsys: overflow in power distribution'
        rw.append(x)
      rw_jumbled = []
      rwi = range(num_genes)
      for i in xrange(num_genes) :
        j = rng.random_range(len(rwi))
        r = rwi[j]
        del rwi[j]
        rw_jumbled.append(rw[r])
      wheel = transrnd.RouletteWheel(rng, rw_jumbled)
      for i in xrange(num_edges) :
        g0 = wheel.pocket()
        g1 = wheel.pocket()
        # FIXME: this loop may run very often or even infinitely often when
        # generating dense graphs
        while g1 in linklist[g0] or -g1 - 1 in linklist[g0] :
          g0 = wheel.pocket()
          g1 = wheel.pocket()
        if rng.random_range(2) :
          linklist[g0].append(-g1 - 1)
        else :
          linklist[g0].append(g1)
      return linklist

    def random_uniform_linklist(num_genes, num_edges, rng) :
      all_links = []
      for i in xrange(num_genes) :
        for j in xrange(num_genes) :
          if i != j :
            all_links.append((i, j))
      linklist = []
      for i in xrange(num_genes) :
        linklist.append([])
      for i in xrange(num_edges) :
        r = rng.random_range(len(all_links))
        g0, g1 = all_links[r]
        del all_links[r]
        if rng.random_range(2) :
          linklist[g0].append(-g1 - 1)
        else :
          linklist[g0].append(g1)
      return linklist

    def next_expression(s) :
      p = TranssysProgramParser(StringIO.StringIO(s.nextval()))
      return p.parse_expr()

    if self.topology == 'random_nk' :
      linklist = random_nk_linklist(self.n, self.k, self.rng)
    elif self.topology == 'random_uniform' :
      linklist = random_uniform_linklist(self.n, self.num_edges, self.rng)
    elif self.topology == 'random_powerlaw' :
      linklist = random_powerlaw_linklist(self.n, self.num_edges, self.power_exp, self.power_base, self.rng)
    elif self.topology == 'linklist' :
      linklist = self.linklist
    flist = []
    glist = []
    for i in xrange(self.num_genes()) :
      flist.append(Factor(factor_name(i), next_expression(self.decay), next_expression(self.diffusibility)))
      promoter = []
      promoter.append(PromoterElementConstitutive(next_expression(self.constitutive)))
      for t in linklist[i] :
        if t < 0 :
          promoter.append(PromoterElementRepress(next_expression(self.km_repression), next_expression(self.vmax_repression), [factor_name(abs(t) - 1)]))
        else :
          promoter.append(PromoterElementActivate(next_expression(self.km_activation), next_expression(self.vmax_activation), [factor_name(t)]))
      glist.append(Gene(gene_name(i), factor_name(i), promoter))
    tp = TranssysProgram(name, flist, glist)
    self.serial = self.serial + 1
    tp.comments.extend(str(self).split('\n'))
    tp.comments.append('serial #%d' % self.serial)
    return tp


class TranssysProgramScanner :

  def __init__(self, f) :
    self.infile = f
    self.buffer = ''
    self.lineno = 0
    self.keywords = ['factor', 'gene', 'promoter', 'product', 'constitutive', 'activate', 'repress', 'default', 'gauss', 'random', 'transsys', 'decay', 'diffusibility', 'lsys', 'symbol', 'axiom', 'rule', 'diffusionrange', '-->', 'graphics', 'move', 'sphere', 'cylinder', 'box', 'turn', 'roll', 'bank', 'color', 'push', 'pop', '<=', '>=', '==', '!=', '&&', '||']
    self.identifier_re = re.compile('[A-Za-z_][A-Za-z0-9_]*')
    self.realvalue_re = re.compile('[+-]?(([0-9]+(\\.[0-9]*)?)|(\\.[0-9]+))([Ee][+-]?[0-9]+)?')
    self.next_token = self.get_token()


  def lookahead(self) :
    return self.next_token[0]


  def token(self) :
    return_token = self.next_token
    self.next_token = self.get_token()
    return return_token


  def get_token(self) :
    # FIXME: lineno reflects line of lookahead token
    if len(self.buffer) > 0 :
      if self.buffer[0] == '#' or self.buffer[0:2] == '//' :
        self.buffer = ''
    while self.buffer == '' :
      self.buffer = self.infile.readline()
      self.lineno = self.lineno + 1
      if self.buffer == '' :
        return None, None
      self.buffer = string.strip(self.buffer)
      if len(self.buffer) > 0 :
        if self.buffer[0] == '#' or self.buffer[0:2] == '//' :
          self.buffer = ''
    for kw in self.keywords :
      if self.buffer[:len(kw)] == kw :
        self.buffer = string.strip(self.buffer[len(kw):])
        return kw, None
    m = self.identifier_re.match(self.buffer)
    if m :
      s = m.group()
      self.buffer = string.strip(self.buffer[len(s):])
      return ('identifier', s)
    m = self.realvalue_re.match(self.buffer)
    if m :
      s = m.group()
      v = float(s)
      self.buffer = string.strip(self.buffer[len(s):])
      return ('realvalue', v)
    c = self.buffer[0]
    self.buffer = string.strip(self.buffer[1:])
    return c, None
    raise StandardError, 'line %d: scanner stalled at "%s"' % (self.lineno, self.buffer)


  def get_lines(self) :
    if self.next_token[0] == 'identifier' :
      l = self.next_token[1]
    elif self.next_token[0] == 'realvalue' :
      l = str(self.next_token[1])
    elif self.next_token[0] is None :
      l = ''
    else :
      l = self.next_token[0] + ' '
    l = l + self.buffer
    lines = []
    if l :
      lines.append(l)
    l = self.infile.readline()
    while l :
      lines.append(l[:-1])
      l = self.infile.readline()
    return lines


class TranssysProgramParser :

  def __init__(self, infile) :
    if type(infile) is types.StringType :
      f = open(infile, 'r')
    else :
      f = infile
    self.scanner = TranssysProgramScanner(f)


  def expect_token(self, expected_token) :
    t, v = self.scanner.token()
    if t != expected_token :
      raise StandardError, 'line %d: expected token "%s" but got "%s"' % (self.scanner.lineno, expected_token, t)
    return v


  def consume_if_found(self, expected_token) :
    if self.scanner.lookahead() == expected_token :
      self.expect_token(expected_token)
      return True
    else :
      return False


  def parse_realvalue_expr(self) :
    v = self.expect_token('realvalue')
    return ExpressionNodeValue(v)


  def parse_identifier_expr(self, transsys_label_list) :
    id = self.expect_token('identifier')
    if self.consume_if_found('.') :
      if id not in transsys_label_list :
        raise StandardError, 'unknown transsys label "%s"' % id
      transsys_id = id
      id = self.expect_token('identifier')
    else :
      transsys_id = None
    return ExpressionNodeIdentifier(id, transsys_id)


  def parse_argument_list(self, transsys_label_list) :
    arglist = []
    self.expect_token('(')
    if self.consume_if_found(')') :
      return arglist
    sep = ','
    while sep == ',' :
      arglist.append(self.parse_expr(transsys_label_list))
      sep, v = self.scanner.token()
      if sep != ',' and sep != ')' :
        raise StandardError, 'line %d: expected token "," or ")" in arglist but got "%s"' % (self.scanner.lineno, sep)
    return arglist


  def parse_gauss_expr(self, transsys_label_list) :
    self.expect_token('gauss')
    arglist = self.parse_argument_list(transsys_label_list)
    if len(arglist) != 2 :
      raise StandardError, 'line %d: gauss() takes 2 parameters, but got %d' % (self.scanner.lineno, len(arglist))
    return ExpressionNodeGaussianRandom(arglist[0], arglist[1])


  def parse_random_expr(self, transsys_label_list) :
    self.expect_token('random')
    arglist = self.parse_argument_list(transsys_label_list)
    if len(arglist) != 2 :
      raise StandardError, 'line %d: gauss() takes 2 parameters, but got %d' % (self.scanner.lineno, len(arglist))
    return ExpressionNodeUniformRandom(arglist[0], arglist[1])


  def parse_value_expr(self, transsys_label_list) :
    l = self.scanner.lookahead()
    if l == 'realvalue' :
      return self.parse_realvalue_expr()
    elif l == 'identifier' :
      return self.parse_identifier_expr(transsys_label_list)
    elif l == 'random' :
      return self.parse_random_expr(transsys_label_list)
    elif l == 'gauss' :
      return self.parse_gauss_expr(transsys_label_list)
    elif l == '(' :
      self.expect_token('(')
      expr = self.parse_expr(transsys_label_list)
      self.expect_token(')')
      return expr
    else :
      raise StandardError, 'line %d: got token "%s", but value must be either <realvalue>, <identifier>, "random" or "gauss"' % (self.scanner.lineno, l)


  def parse_term_expr(self, transsys_label_list) :
    expr1 = self.parse_value_expr(transsys_label_list)
    while 1:
      l = self.scanner.lookahead()
      if l == '*' :
        self.expect_token('*')
        expr2 = self.parse_value_expr(transsys_label_list)
        expr1 = ExpressionNodeMult(expr1, expr2)
      elif l == '/' :
        self.expect_token('/')
        expr2 = self.parse_value_expr(transsys_label_list)
        expr1 = ExpressionNodeDiv(expr1, expr2)
      else :
        break
    return expr1


  def parse_arithmetic_expr(self, transsys_label_list) :
    expr1 = self.parse_term_expr(transsys_label_list)
    while 1:
      l = self.scanner.lookahead()
      if l == '+' :
        self.expect_token('+')
        expr2 = self.parse_term_expr(transsys_label_list)
        expr1 = ExpressionNodeAdd(expr1, expr2)
      elif l == '-' :
        self.expect_token('-')
        expr2 = self.parse_term_expr(transsys_label_list)
        expr1 = ExpressionNodeSubtract(expr1, expr2)
      else :
        break
    return expr1


  def parse_cmp_expr(self, transsys_label_list) :
    expr1 = self.parse_arithmetic_expr(transsys_label_list)
    while 1 :
      l = self.scanner.lookahead()
      if l == '<' :
        self.expect_token('<')
        expr2 = self.parse_arithmetic_expr(transsys_label_list)
        expr1 = ExpressionNodeLower(expr1, expr2)
      elif l == '<=' :
        self.expect_token('<=')
        expr2 = self.parse_arithmetic_expr(transsys_label_list)
        expr1 = ExpressionNodeLowerEqual(expr1, expr2)
      elif l == '>' :
        self.expect_token('>')
        expr2 = self.parse_arithmetic_expr(transsys_label_list)
        expr1 = ExpressionNodeGreater(expr1, expr2)
      elif l == '>=' :
        self.expect_token('?=')
        expr2 = self.parse_arithmetic_expr(transsys_label_list)
        expr1 = ExpressionNodeGreaterEqual(expr1, expr2)
      elif l == '==' :
        self.expect_token('==')
        expr2 = self.parse_arithmetic_expr(transsys_label_list)
        expr1 = ExpressionNodeEqual(expr1, expr2)
      elif l == '!=' :
        self.expect_token('==')
        expr2 = self.parse_arithmetic_expr(transsys_label_list)
        expr1 = ExpressionNodeUnequal(expr1, expr2)
      else :
        break
    return expr1


  def parse_not_expr(self, transsys_label_list) :
    if self.consume_if_found('!') :
      expr = self.parse_cmp_expr(transsys_label_list)
      return ExpressionNodeNot(expr)
    else :
      return self.parse_cmp_expr(transsys_label_list)


  def parse_and_expr(self, transsys_label_list) :
    expr1 = self.parse_not_expr(transsys_label_list)
    while 1 :
      if self.consume_if_found('&&') :
        expr2 = self.parse_not_expr(transsys_label_list)
        expr1 = ExpressionNodeAnd(expr1, expr2)
      else :
        break
    return expr1


  def parse_expr(self, transsys_label_list) :
    expr1 = self.parse_and_expr(transsys_label_list)
    while 1 :
      if self.consume_if_found('||') :
        expr2 = self.parse_and_expr(transsys_label_list)
        expr1 = ExpressionNodeOr(expr1, expr2)
      else :
        break
    return expr1


  def parse_product_statement(self) :
    self.expect_token('default')
    self.expect_token(':')
    product_id = self.expect_token('identifier')
    self.expect_token(';')
    return product_id


  def parse_product_component(self) :
    self.expect_token('product')
    self.expect_token('{')
    p = self.parse_product_statement()
    self.expect_token('}')
    return p


  def parse_factor_combination(self) :
    fc = [self.expect_token('identifier')]
    while self.consume_if_found('+') :
      fc.append(self.expect_token('identifier'))
    return fc


  def parse_promoter_statement(self) :
    if self.consume_if_found('constitutive') :
      self.expect_token(':')
      expr = self.parse_expr([])
      self.expect_token(';')
      return PromoterElementConstitutive(expr)
    fc = self.parse_factor_combination()
    self.expect_token(':')
    if self.consume_if_found('activate') :
      arglist = self.parse_argument_list([])
      if len(arglist) != 2 :
        raise StandardError, 'line %d: activate() takes exactly 2 parameters, got %d' % (self.scanner.lineno, len(arglist))
      self.expect_token(';')
      return PromoterElementActivate(arglist[0], arglist[1], fc)
    elif self.consume_if_found('repress') :
      arglist = self.parse_argument_list([])
      if len(arglist) != 2 :
        raise StandardError, 'line %d: repress() takes exactly 2 parameters, got %d' % (self.scanner.lineno, len(arglist))
      self.expect_token(';')
      return PromoterElementRepress(arglist[0], arglist[1], fc)
    raise StandardError, 'line %d: bad token "%s" in promoter statement' % (self.scanner.lineno, self.scanner.lookahead())


  def parse_promoter_component(self) :
    self.expect_token('promoter')
    self.expect_token('{')
    promoter_elements = []
    while self.scanner.lookahead() != '}' :
      promoter_elements.append(self.parse_promoter_statement())
    self.expect_token('}')
    return promoter_elements


  def parse_factor_definition(self) :
    self.expect_token('factor')
    factor_name = self.expect_token('identifier')
    self.expect_token('{')
    decay_expr = ExpressionNodeValue(1.0)
    diffusibility_expr = ExpressionNodeValue(1.0)
    while 1 :
      l = self.scanner.lookahead()
      if l == 'decay' :
        self.expect_token('decay')
        self.expect_token(':')
        decay_expr = self.parse_expr([])
      elif l == 'diffusibility' :
        self.expect_token('diffusibility')
        self.expect_token(':')
        diffusibility_expr = self.parse_expr([])
      if self.scanner.lookahead() == ';' :
        self.expect_token(';')
      else :
        break
    self.expect_token('}')
    return Factor(factor_name, decay_expr, diffusibility_expr)


  def parse_gene_definition(self) :
    self.expect_token('gene')
    gene_name = self.expect_token('identifier')
    self.expect_token('{')
    promoter = self.parse_promoter_component()
    product = self.parse_product_component()
    self.expect_token('}')
    return Gene(gene_name, product, promoter)


  def parse_transsys_elements(self) :
    factors = []
    genes = []
    l = self.scanner.lookahead()
    while l == 'factor' or l == 'gene' :
      #print 'parsing "%s", lookahead is "%s"' % (l, self.scanner.lookahead())
      if l == 'factor' :
        factors.append(self.parse_factor_definition())
      elif l == 'gene' :
        genes.append(self.parse_gene_definition())
      l = self.scanner.lookahead()
    return factors, genes


  def parse_transsys(self) :
    t, v = self.scanner.token()
    if t is None :
      return None
    transsys_name = self.expect_token('identifier')
    self.expect_token('{')
    factors, genes = self.parse_transsys_elements()
    self.expect_token('}')
    return TranssysProgram(transsys_name, factors, genes)


  def parse_diffusionrange(self) :
    self.expect_token('diffusionrange')
    self.expect_token(':')
    d = self.expect_token('realvalue')
    self.expect_token(';')
    return int(d)


  def parse_symbol(self) :
    transsys_name = None
    self.expect_token('symbol')
    symbol_name = self.expect_token('identifier')
    if self.scanner.lookahead() == '(' :
      self.expect_token('(')
      transsys_name = self.expect_token('identifier')
      self.expect_token(')')
    self.expect_token(';')
    return Symbol(symbol_name, transsys_name)


  def parse_assignment(self, transsys, transsys_label_list) :
    factor = self.expect_token('identifier')
    self.expect_token('=')
    value = self.parse_expr(transsys_label_list)
    a = Assignment(transsys, factor, value)
    print str(a)
    return a


  def parse_assignment_list(self, transsys, transsys_label_list) :
    alist = []
    while self.scanner.lookahead() == 'identifier' :
      alist.append(self.parse_assignment(transsys, transsys_label_list))
    return alist


  def parse_lhs_symbol(self, symbol_dict) :
    symbol_name = self.expect_token('identifier')
    if symbol_name not in symbol_dict.keys() :
      raise StandardError, 'unknown symbol "%s"' % symbol_name
    symbol = symbol_dict[symbol_name]
    if self.scanner.lookahead() == '(' :
      self.expect_token('(')
      transsys_name = self.expect_token('identifier')
      self.expect_token(')')
    else :
      transsys_name = None
    return LhsSymbol(symbol, transsys_name)


  def parse_lhs(self, symbol_dict) :
    lhs = [self.parse_lhs_symbol(symbol_dict)]
    l = self.scanner.lookahead()
    while l != ':' and l != '-->' :
      lhs.append(self.parse_lhs_symbol(symbol_dict))
      l = self.scanner.lookahead()
    return lhs


  def parse_transsys_initializer(self, transsys, transsys_label_list) :
    if self.consume_if_found('transsys') :
      template_label = self.expect_token('identifier')
      self.expect_token(':')
    else :
      template_label = None
    assignments = self.parse_assignment_list(transsys, transsys_label_list)
    return template_label, assignments
      

  def parse_production_element(self, symbol_dict, transsys_label_list) :
    symbol_name = self.expect_token('identifier')
    if symbol_name not in symbol_dict.keys() :
      raise StandardError, 'unknown symbol "%s"' % symbol_name
    symbol = symbol_dict[symbol_name]
    if self.consume_if_found('(') :
      template_label, assignments = self.parse_transsys_initializer(symbol.transsys, transsys_label_list)
      self.expect_token(')')
    else :
      template_label = None
      assignments = []
    return ProductionElement(symbol, template_label, assignments)


  def parse_rhs(self, symbol_dict, transsys_label_list) :
    rhs = []
    while self.scanner.lookahead() == 'identifier' :
      rhs.append(self.parse_production_element(symbol_dict, transsys_label_list))
    return rhs


  def parse_rule(self, symbol_dict) :
    self.expect_token('rule')
    rule_name = self.expect_token('identifier')
    self.expect_token('{')
    lhs = self.parse_lhs(symbol_dict)
    transsys_label_list = map(lambda x: x.transsys_name, lhs)
    if self.consume_if_found(':') :
      condition = self.parse_expr(transsys_label_list)
    else :
      condition = None
    self.expect_token('-->')
    rhs = self.parse_rhs(symbol_dict, transsys_label_list)
    self.expect_token('}')
    return Rule(rule_name, lhs, condition, rhs)


  def parse_axiom(self, symbol_dict) :
    self.expect_token('axiom')
    axiom = self.parse_rhs(symbol_dict, [])
    self.expect_token(';')
    return axiom


  def parse_move(self, symbol) :
    self.expect_token('move')
    arglist = self.parse_argument_list([])
    if len(arglist) != 1 :
      raise StandardError, 'move takes 1 argument, got %d' % len(arglist)
    return GraphicsPrimitiveMove(arglist[0])


  def parse_push(self, symbol) :
    self.expect_token('push')
    arglist = self.parse_argument_list([])
    if len(arglist) != 0 :
      raise StandardError, 'push does not take arguments, got %d' % len(arglist)
    return GraphicsPrimitivePush()


  def parse_pop(self, symbol) :
    self.expect_token('pop')
    arglist = self.parse_argument_list([])
    if len(arglist) != 0 :
      raise StandardError, 'pop does not take arguments, got %d' % len(arglist)
    return GraphicsPrimitivePop()


  def parse_turn(self, symbol) :
    self.expect_token('turn')
    arglist = self.parse_argument_list([])
    if len(arglist) != 1 :
      raise StandardError, 'turn takes 1 argument, got %d' % len(arglist)
    return GraphicsPrimitiveTurn(arglist[0])


  def parse_bank(self, symbol) :
    self.expect_token('bank')
    arglist = self.parse_argument_list([])
    if len(arglist) != 1 :
      raise StandardError, 'bank takes 1 argument, got %d' % len(arglist)
    return GraphicsPrimitiveBank(arglist[0])


  def parse_roll(self, symbol) :
    self.expect_token('roll')
    arglist = self.parse_argument_list([])
    if len(arglist) != 1 :
      raise StandardError, 'roll takes 1 argument, got %d' % len(arglist)
    return GraphicsPrimitiveRoll(arglist[0])


  def parse_sphere(self, symbol) :
    self.expect_token('sphere')
    arglist = self.parse_argument_list([])
    if len(arglist) != 1 :
      raise StandardError, 'sphere takes 1 argument, got %d' % len(arglist)
    return GraphicsPrimitiveSphere(arglist[0])


  def parse_cylinder(self, symbol) :
    self.expect_token('cylinder')
    arglist = self.parse_argument_list([])
    if len(arglist) != 2 :
      raise StandardError, 'cylinder takes 2 argument, got %d' % len(arglist)
    return GraphicsPrimitiveCylinder(arglist[0], arglist[1])


  def parse_box(self, symbol) :
    self.expect_token('box')
    arglist = self.parse_argument_list([])
    if len(arglist) != 3 :
      raise StandardError, 'box takes 3 argument, got %d' % len(arglist)
    return GraphicsPrimitiveBox(arglist[0], arglist[1], arglist[2])


  def parse_color(self, symbol) :
    self.expect_token('color')
    arglist = self.parse_argument_list([])
    if len(arglist) != 3 :
      raise StandardError, 'color takes 3 argument, got %d' % len(arglist)
    return GraphicsPrimitiveColor(arglist[0], arglist[1], arglist[2])


  def parse_symbol_graphics(self, symbol_dict) :
    graphics_primitives = ['move', 'push', 'pop', 'turn', 'bank', 'roll', 'sphere', 'cylinder', 'box', 'color']
    symbol_name = self.expect_token('identifier')
    if symbol_name not in symbol_dict.keys() :
      raise StandardError, 'unknown symbol "%s"' % symbol_name
    symbol = symbol_dict[symbol_name]
    self.expect_token('{')
    plist = []
    l = self.scanner.lookahead()
    while l in graphics_primitives :
      if l == 'move' :
        p = self.parse_move(symbol)
      elif l == 'push' :
        p = self.parse_push(symbol)
      elif l == 'pop' :
        p = self.parse_pop(symbol)
      elif l == 'turn' :
        p = self.parse_turn(symbol)
      elif l == 'bank' :
        p = self.parse_bank(symbol)
      elif l == 'roll' :
        p = self.parse_roll(symbol)
      elif l == 'sphere' :
        p = self.parse_sphere(symbol)
      elif l == 'cylinder' :
        p = self.parse_cylinder(symbol)
      elif l == 'box' :
        p = self.parse_box(symbol)
      elif l == 'color' :
        p = self.parse_color(symbol)
      self.expect_token(';')
      plist.append(p)
      l = self.scanner.lookahead()
    self.expect_token('}')
    if symbol.graphicsPrimitiveList is not None :
      sys.stderr.write('symbol "%s": superseding graphics\n' % symbol_name)
    symbol.graphicsPrimitiveList = plist


  def parse_graphics(self, symbol_dict) :
    self.expect_token('graphics')
    self.expect_token('{')
    while self.scanner.lookahead() == 'identifier' :
      self.parse_symbol_graphics(symbol_dict)
    self.expect_token('}')


  def parse_lsys_elements(self) :
    lsys_elements = ['symbol', 'axiom', 'diffusionrange', 'rule', 'graphics']
    symbols = []
    symbol_dict = {}
    axiom = None
    diffusionrange = None
    rules = []
    l = self.scanner.lookahead()
    while l in lsys_elements :
      if l == 'symbol' :
        symbol = self.parse_symbol()
        symbols.append(symbol)
        symbol_dict[symbol.name] = symbol
      elif l == 'axiom' :
        if axiom is not None :
          raise StandardError, 'multiple axioms'
        axiom = self.parse_axiom(symbol_dict)
      elif l == 'diffusionrange' :
        if diffusionrange is not None :
          raise StandardError, 'multiple diffusionrange specs'
        diffusionrange = self.parse_diffusionrange()
      elif l == 'rule' :
        rules.append(self.parse_rule(symbol_dict))
      elif l == 'graphics' :
        self.parse_graphics(symbol_dict)
      l = self.scanner.lookahead()
    return symbols, axiom, diffusionrange, rules


  def parse_lsys(self) :
    self.expect_token('lsys')
    lsys_name = self.expect_token('identifier')
    self.expect_token('{')
    symbols, axiom, diffusionrange, rules = self.parse_lsys_elements()
    self.expect_token('}')
    return LsysProgram(lsys_name, symbols, axiom, diffusionrange, rules)


  def parse(self) :
    t = self.scanner.lookahead()
    if t is None :
      return None
    elif t == 'transsys' :
      return self.parse_transsys()
    elif t == 'lsys' :
      return self.parse_lsys()
    else :
      raise StandardError, 'line %d: expected transsys or lsys, got %s' % (self.scanner.lineno, t)

  # kludges for handling lsys code in lists of lines

  def lsys_lines(self) :
    """read lsys code lines (no parsing)"""
    return self.scanner.get_lines()


def write_lsys_code(f, lines) :
  for l in lines :
    f.write('%s\n' % l)


def write_random_transsys(f, name, rtp) :
  raise StandardError, 'write_random_transsys is not implemented anymore. Use RandomTranssysParameters::generate_transsys() method'


def expression_level_dict(transsys_filename, num_timesteps, sampling_period, initial_levels = None) :
  cmd = 'transexpr -n %d -d %d ' % (num_timesteps, sampling_period)
  if type(initial_levels) is types.FloatType :
    cmd = cmd + '-f %f ' % initial_levels
  elif type(initial_levels) is types.ListType :
    glue = ''
    s = ''
    for x in initial_levels :
      s = s + '%s%f' % (glue, x)
      glue = ' '
    cmd = cmd + '-F \'%s\' ' % s
  cmd = cmd + transsys_filename
  # print cmd
  p = os.popen(cmd, 'r')
  level_dict = {}
  factor_list = None
  line = p.readline()
  while line :
    line = string.strip(line)
    l = string.split(line)
    if line[0] == '#' :
      # this is a hack to figure out the factor names
      if l[1] == 'time' :
        factor_list = l[2:]
      else :
        raise StandardError, 'expression_level_dict: encountered unknown comment "%s"' % line
    else :
      t = string.atoi(l[0])
      levels = []
      for s in l[1:] :
        levels.append(string.atof(s))
      level_dict[t] = levels
    line = p.readline()
  status = p.close()
  if status :
    raise StandardError, 'expression_level_dict: transexpr exit status %d' % status
  return level_dict, factor_list


class ArrayIntensityFunction :
  """A basic intensity function which computes the intensity by
taking the expression level and adding background and noise
"""

  def __init__(self, offset, dispersion, rndseed) :
    self.offset = offset
    self.dispersion = dispersion
    self.rng = transrnd.transrnd(rndseed)


  def __call__(self, factor_name, expression_level) :
    intensity = expression_level + self.offset * math.exp(self.rng.gauss() * self.dispersion)
    return intensity


class ArraySeries :

  def __init__(self, transsys_filename, initial_state, num_timesteps_init, exp_threshold = 0.0) :
    self.transsys_program = initial_state.transsys_program
    self.exp_threshold = exp_threshold
    self.array_description = []
    self.series_list = []
    self.series_directory = {}
    self.initial_state = TranssysInstance(self.transsys_program)
    self.initial_state.factor_concentration = copy.deepcopy(initial_state.factor_concentration)
    self.num_timesteps_init = num_timesteps_init
    self.array_data = {}
    ts = self.initial_state.time_series(self.num_timesteps_init, self.num_timesteps_init - 1)
    self.reference_state = ts[self.num_timesteps_init - 1]
    for f in self.transsys_program.factor_list :
      self.array_data[f.name] = []


  def simulate_timeseries(self, series_name, initial_state, num_timesteps, sampling_period, reference_intensity_func, array_intensity_func) :

    def log2(x) :
      return math.log(x) / math.log(2)

    if initial_state.transsys_program is not self.transsys_program :
      raise StandardError, 'ArraySeries::simulate_timeseries: transsys programs not identical'
    ts = initial_state.time_series(num_timesteps, sampling_period)
    time_list = ts.keys()
    time_list.sort()
    self.series_list.append(series_name)
    self.series_directory[series_name] = []
    for t in time_list :
      arr_name = '%s_t%d' % (series_name, t)
      self.series_directory[series_name].append(arr_name)
      self.array_description.append(arr_name)
      for i in xrange(len(self.transsys_program.factor_list)) :
        ref_intensity = reference_intensity_func(self.transsys_program.factor_list[i].name, self.reference_state.factor_concentration[i])
        intensity = array_intensity_func(self.transsys_program.factor_list[i].name, ts[t].factor_concentration[i])
        if ref_intensity <= self.exp_threshold or intensity <= self.exp_threshold :
          ratio = None
        else :
          ratio = log2(intensity / ref_intensity)
        self.array_data[self.transsys_program.factor_list[i].name].append(ratio)


  def write_experimentspecs(self, f) :
    for n in self.series_list :
      f.write('%s\t%s\n' % (n, string.join(self.series_directory[n], '\t')))


  def write_data(self, f, factor_list = None, series_name = None, complete_only = 0, fieldsep = '\t') :
    if factor_list is None :
      factor_list = self.transsys_program.factor_names()
      # factor_list.sort()
    if series_name is None :
      array_index_list = range(len(self.array_description))
    else :
      if series_name not in self.series_directory.keys() :
        raise StandardError, 'ArraySeries::write_data: unknown array series "%s"' % series_name
      array_index_list = self.series_directory[series_name]
    s = 'Factor'
    for i in array_index_list :
      s = s + fieldsep + self.array_description[i]
    f.write('%s\n' % s)
    for factor in factor_list :
      s = factor
      if complete_only :
        if not self.factor_complete(factor) :
          f.write('# %s: incomplete data\n' % factor)
          continue
      for i in array_index_list :
        v = self.array_data[factor][i]
        s = s + fieldsep
        if v is not None :
          s = s + '%f' % v
      f.write('%s\n' % s)


  def array_complete(self, array_index) :
    for f in self.array_data.keys() :
      if self.array_data[f][array_index] is None :
        return 0
    return 1


  def factor_complete(self, factor_name) :
    for v in self.array_data[factor_name] :
      if v is None :
        return 0
    return 1


  def num_arrays(self) :
    f = self.array_data.keys()[0]
    return len(self.array_data[f])


  def num_genes(self) :
    return self.array_data.num_genes()


  def remove_array(self, array_index) :
    del self.array_description[array_index]
    for f in self.array_data.keys() :
      del self.array_data[f][array_index]


  def remove_incomplete_arrays(self) :
    i = 0
    while i < self.num_arrays() :
      if self.array_complete(i) :
        i = i + 1
      else :
        self.remove_array(i)


  def remove_incomplete_factors(self) :
    for f in self.array_data.keys() :
      if not self.factor_complete(f) :
        del self.array_data[f]


if __name__ == '__main__' :
  import sys
  import transrnd

  rtp = RandomTranssysParameters()
  rtp.rng = transrnd.transrnd(1)
  rtp.n = 3
  rtp.k = 2
  write_random_transsys(sys.stdout, 'rndtrans0', rtp)

  rtp.rng = transrnd.transrnd(1)
  rtp.set_constitutive(1.0)
  write_random_transsys(sys.stdout, 'rndtrans_c1', rtp)

  rtp.rng = transrnd.transrnd(1)
  rtp.set_km_activation([0.0, 1.0])
  rtp.set_km_repression([2.0, 3.0])
  rtp.set_vmax_activation([47.11, 8.15])
  rtp.set_vmax_repression([1.2, 3.4, 5.6, 7.8, 9.0, 10.11])
  rtp.set_decay([0.1, 0.2])
  rtp.set_diffusibility([0.88, 0.99])
  write_random_transsys(sys.stdout, 'rndtrans_cycles', rtp)

  rtp.rng = transrnd.transrnd(1)
  rtp.set_topology([[0, 1], [], [0, 4], [0, -1, 1], [-1, -2, -3, -4, -5]])
  write_random_transsys(sys.stdout, 'rndtrans_topo', rtp)
