#!/usr/bin/env python

# $Id$

# $Log$
# Revision 1.1  2005/03/08 17:12:02  jtk
# Initial revision
#
# Revision 1.8  2003/03/14 01:20:18  kim
# added stuff for Pearson-like distance accroding to Eisen et.al.
#
# Revision 1.7  2003/02/28 00:46:14  kim
# improved intensity computation
#
# Revision 1.6  2003/02/27 10:36:48  kim
# added ArrayIntensityFunction class for simulating noisy arrays etc.
#
# Revision 1.5  2003/02/26 11:46:07  kim
# added uniform perturber
#


import sys
import math
import getopt
import copy
import string
import StringIO
import re
import os

import transsys
import transrnd


def get_namevalue(f) :
  line = f.readline()
  while line != '' :
    line = line.strip()
    if line != '' :
      if line[0] != '#' :
        break
    line = f.readline()
  if line == '' :
    return None, None
  nv = line.split(':', 1)
  if len(nv) == 1 :
    return nv[0].strip(), None
  else :
    return nv[0].strip(), nv[1].strip()


def parse_int_csv(s) :
  l = []
  for v in s.split(',') :
    l.append(int(v.strip()))
  return l


def parse_string_csv(s) :
  l = []
  for v in s.split(',') :
    l.append(v.strip())
  return l


class ClusterNode :

  def __init__(self, name = '???', ancestor = None, distance = 0.0, descendant_list = None) :
    self.name = name
    self.ancestor = ancestor
    if ancestor is not None :
      ancestor.descendant_list.append(self)
    self.distance = distance
    if descendant_list is None :
      descendant_list = []
    self.descendant_list = copy.copy(descendant_list)
    for d in self.descendant_list :
      d.ancestor = self


  def set_ancestor(self, ancestor, distance = None) :
    if self.ancestor is not None :
      i = self.ancestor.descendant_list.index(self)
      del self.ancestor.descendant_list[i]
    self.ancestor = ancestor
    if distance is not None :
      self.distance = distance
    ancestor.descendant_list.append(self)


  def __str__(self) :
    if len(self.descendant_list) == 0 :
      return '%s:%f' % (self.name, self.distance)
    else :
      s = '('
      glue = ''
      for d in self.descendant_list :
        s = s + glue + str(d)
        glue = ', '
      s = s + ')' 
      return s + '%s:%f' % (self.name, self.distance)


  def cluster_size(self) :
    s = 1
    for d in self.descendant_list :
      s = s + d.cluster_size()
    return s


  def set_ultrametric_distance(self, d) :
    self.distance = d
    dl = self.descendant_list
    while len(dl) > 0 :
      self.distance = self.distance - dl[0].distance
      dl = dl[0].descendant_list


def parse_cluster(f) :

  def get_cluster(s, cluster_dict) :
    m = re.match('[0-9]+', s)
    if m :
      s = int(s)
    else :
      cluster_dict[s] = ClusterNode(s)
    return cluster_dict[s]

  int_re = re.compile('[0-9]+')
  cluster_dict = {}
  cluster_index = 1
  line = f.readline()
  a = None
  while line :
    if line.strip() == '' :
      break
    l = line.split()
    d0 = get_cluster(l[0], cluster_dict)
    d1 = get_cluster(l[1], cluster_dict)
    distance = float(l[2])
    a = ClusterNode('', descendant_list = [d0, d1])
    d0.set_ultrametric_distance(distance)
    d1.set_ultrametric_distance(distance)
    cluster_dict[cluster_index] = a
    cluster_index = cluster_index + 1
    line = f.readline()
  if a.cluster_size() != len(cluster_dict) :
    raise StandardError, 'parse_cluster: unconnected components'
  return a


def regulon_name(r) :
  return 'regulon_%02d' % r


def struct_protein_name(r, i) :
  return 'strp%04d_%02d' % (r, i)


def struct_gene_name(r, i) :
  return 'sgene%04d_%02d' % (r, i)


def regcontrol_gene_name(r) :
  return 'gctrl%04d' % r


def regcontrol_factor_name(r) :
  return 'fctrl%04d' % r


class RegulonDescriptor :

  def __init__(self, name, controlling_factors = None, coregulated_factors = None) :
    if controlling_factors is None :
      controlling_factors = []
    if coregulated_factors is None :
      coregulated_factors = []
    self.name = name
    self.controlling_factors = controlling_factors[:]
    self.coregulated_factors = coregulated_factors[:]


  def cluster(self) :
    r = ClusterNode(string.join(self.controlling_factors, '/'))
    for f in self.controlling_factors :
      ClusterNode(f, r, 2.0)
    c = ClusterNode('', r, 1.0)
    for f in self.coregulated_factors :
      ClusterNode(f, c, 1.0)
    return r


  def __str__(self) :
    s = 'RegulonDescriptor %s\n' % self.name
    s = s + 'controlling_factors: %s\n' % string.join(self.controlling_factors, ', ')
    s = s + 'coregulated_factors: %s\n' % string.join(self.coregulated_factors, ', ')
    return s


  def parse(self, f) :

    def parse_factor_list(f, label) :
      l = '%s:' % label
      line = f.readline()
      if line[:len(l)] != l :
        raise StandardError, 'RegulonDescriptor::parse: expected label "%s" but found "%s"' % (label, line.strip())
      return parse_string_csv(line[len(l):])

    line = f.readline()
    m = re.match('RegulonDescriptor (.*)', line)
    if m is None :
      raise StandardError, 'RegulonDescriptor::parse: bad header line "%s"' % line.strip()
    self.name = m.group(1).strip()
    self.controlling_factors = parse_factor_list(f, 'controlling_factors')
    self.coregulated_factors = parse_factor_list(f, 'coregulated_factors')


  def genegroup_line(self) :
    s = self.name
    for f in self.controlling_factors :
      s = s + '\t%s' % f
    for f in self.coregulated_factors :
      s = s + '\t%s' % f
    return s


class RegStructTranssys :

  def __init__(self, tp_name, reg_rtp, struct_rtp, regulon_size_list, radius) :
    x0 = radius
    y0 = radius
    r_radius = 0.6 * radius
    self.tp = reg_rtp.generate_transsys(tp_name)
    self.noncontrolling_factors = self.tp.factor_names()
    self.controlling_factors = []
    self.perturbed_factors = []
    self.structural_factors = []
    self.regulon_list = []
    self.regulons = {}
    self.tp.dot_positions_circle(x0, y0, r_radius)
    a_step = 2.0 * math.pi / len(self.tp.gene_list)
    a_range = 0.6 * a_step
    rlist = range(len(self.tp.gene_list))
    for num_structural_genes in regulon_size_list :
      r = struct_rtp.rng.random_range(len(rlist))
      g = rlist[r]
      del rlist[r]
      angle = g * a_step
      self.add_regulon(g, num_structural_genes, struct_rtp, angle, a_range, radius)
    self.tp.resolve()


  def add_regulon(self, g, num_structural_genes, struct_rtp, angle, a_range, radius) :
    """add a regulon, consisting of num_structural_genes controlled by
gene number g, to transsys program self.tp. Parameters for the promoter and
the factor are generated using RandomTranssysParameters instance struct_rtp,
but not all facilities are used."""
    x0 = radius
    y0 = radius
    a_start = angle - 0.5 * a_range
    self.tp.gene_list[g].name = regcontrol_gene_name(g)
    del self.noncontrolling_factors[self.noncontrolling_factors.index(self.tp.gene_list[g].product.name)]
    reg_name = regulon_name(g)
    ctrl_name = regcontrol_factor_name(g)
    self.controlling_factors.append(ctrl_name)
    self.tp.gene_list[g].product.name = ctrl_name
    self.regulons[reg_name] = [ctrl_name]
    # c_radius = 0.7 * radius
    # self.tp.gene_list[g].dot_attributes['pos'] = '%f,%f!' % (x0 + c_radius * math.cos(angle), y0 + c_radius * math.sin(angle))
    cfactor = self.tp.gene_list[g].product_name()
    angle = a_start
    coregulated_factors = []
    for i in xrange(num_structural_genes) :
      struct_name = struct_protein_name(g, i)
      self.regulons[reg_name].append(struct_name)
      self.tp.factor_list.append(transsys.Factor(struct_name, transsys.ExpressionNodeValue(struct_rtp.decay.nextval()), transsys.ExpressionNodeValue(struct_rtp.diffusibility.nextval())))
      coregulated_factors.append(struct_name)
      p = []
      p.append(transsys.PromoterElementConstitutive(transsys.ExpressionNodeValue(struct_rtp.constitutive.nextval())))
      p.append(transsys.PromoterElementActivate(transsys.ExpressionNodeValue(struct_rtp.km_activation.nextval()), transsys.ExpressionNodeValue(struct_rtp.vmax_activation.nextval()), [cfactor]))
      sgene = transsys.Gene(struct_gene_name(g, i), struct_protein_name(g, i), p)
      sgene.dot_attributes['pos'] = '%f,%f!' % (x0 + radius * math.cos(angle), y0 + radius * math.sin(angle))
      if num_structural_genes > 1 :
        angle = angle + a_range / (num_structural_genes - 1)
      self.tp.gene_list.append(sgene)
    self.regulon_list.append(RegulonDescriptor(reg_name, [cfactor], coregulated_factors))
    self.structural_factors.extend(coregulated_factors)


  def regnet_factors(self) :
    return self.noncontrolling_factors + self.controlling_factors

  def choose_perturbed_factors(self, num_perturbed_factors, rng) :
    if num_perturbed_factors is None :
      self.perturbed_factors = self.noncontrolling_factors[:]
      return
    if num_perturbed_factors > len(self.noncontrolling_factors) :
      raise StandardError, 'RegStructTranssys::choose_perturbed_factors: cannot choose %d factors out of %d' % (num_perturbed_factors, len(self.noncontrolling_factors))
    l = range(len(self.noncontrolling_factors))
    self.perturbed_factors = []
    for i in xrange(num_perturbed_factors) :
      r = rng.random_range(len(l))
      j = l[r]
      del l[r]
      self.perturbed_factors.append(self.noncontrolling_factors[j])


  def write_genegroups(self, f) :
    f.write('perturbed')
    for factor in self.perturbed_factors :
      f.write('\t%s' % factor)
    f.write('\n')
    f.write('nonperturbed')
    for factor in self.nonperturbed_factor_list() :
      f.write('\t%s' % factor)
    f.write('\n')
    for r in self.regulons.keys() :
      f.write(r)
      for factor in self.regulons[r] :
        f.write('\t%s' % factor)
      f.write('\n')


class ExpGaussPerturber :

  def __init__(self, skew, dispersion, rndseed, perturbed_factor_list = None) :
    self.perturbed_factor_list = copy.deepcopy(perturbed_factor_list)
    self.skew = skew
    self.dispersion = dispersion
    self.rng = transrnd.transrnd(rndseed)


  def __call__(self, c, factor_name) :

    def perturb(c, factor_name, skew, dispersion, rng) :
      x = rng.gauss() * dispersion
      c_new = (c + skew) * 2.0**(x)
      # print 'perturbing "%s" by 2**%f: %f --> %f' % (factor_name, x, c, c_new)
      return c_new
      
    if self.perturbed_factor_list is None :
      return perturb(c, factor_name, self.skew, self.dispersion, self.rng)
    if factor_name in self.perturbed_factor_list :
      return perturb(c, factor_name, self.skew, self.dispersion, self.rng)
    else :
      return c


class ExponentialReplacementPerturber :

  def __init__(self, halfmax, rndseed, perturbed_factor_list = None) :
    self.perturbed_factor_list = copy.deepcopy(perturbed_factor_list)
    self.halfmax = halfmax
    self.rng = transrnd.transrnd(rndseed)


  def __call__(self, c, factor_name) :

    def perturb(c, factor_name, halfmax, rng) :
      x = rng.rnd()
      c_new = halfmax * math.log(x) / math.log(0.5)
      # print 'perturbing "%s" by 2**%f: %f --> %f' % (factor_name, x, c, c_new)
      return c_new
      
    if self.perturbed_factor_list is None :
      return perturb(c, factor_name, self.halfmax, self.rng)
    if factor_name in self.perturbed_factor_list :
      return perturb(c, factor_name, self.halfmax, self.rng)
    else :
      return c


class UniformReplacementPerturber :

  def __init__(self, c_min, c_max, rndseed, perturbed_factor_list = None) :
    self.perturbed_factor_list = copy.deepcopy(perturbed_factor_list)
    self.c_min = c_min
    self.c_max = c_max
    self.c_range = c_max - c_min
    self.rng = transrnd.transrnd(rndseed)


  def __call__(self, c, factor_name) :
    if self.perturbed_factor_list is None :
      return self.rng.rnd() * self.c_range + self.c_min
    if factor_name in self.perturbed_factor_list :
      return self.rng.rnd() * self.c_range + self.c_min
    else :
      return c


class RTAControlParameters :

  def __init__(self) :
    self.transsys_program = None
    self.regnet_factors = None
    self.structural_factors = None
    self.perturbed_factors = None
    self.regulon_list = []
    self.min_initial_concentration = None
    self.max_initial_concentration = None
    self.num_environments = None
    self.perturbation_type = None
    self.delta_t = 1
    self.samples_per_environment = None
    self.array_offset = 0.0
    self.array_dispersion = 1.0
    self.num_timesteps_init = 0
    self.rndseed = 1
    self.histo_max = 1.0
    self.histo_nbins = 0
    self.basename = None
    self.exp_threshold = 1e-8


  def environment_perturber(self) :
    if self.perturbed_factors is None :
      raise StandardError, 'RTAControlParameters::environment_perturber: perturbed_factors list missing'
    if self.perturbation_type == 'ExpGauss' :
      return ExpGaussPerturber(self.perturbation_skew, self.perturbation_dispersion, self.rndseed, self.perturbed_factors)
    elif self.perturbation_type == 'ExponentialReplacement' :
      return ExponentialReplacementPerturber(self.perturbation_halfmax, self.rndseed, self.perturbed_factors)
    elif self.perturbation_type == 'UniformReplacement' :
      return UniformReplacementPerturber(self.perturbation_min, self.perturbation_max, self.rndseed, self.perturbed_factors)
    else :
      raise StandardError, 'RTAControlParameters::environment_perturber: unknown type "%s"' % self.perturbation_type


  def nonperturbed_factor_list(self) :
    l = []
    for f in self.transsys_program.factor_names() :
      if f not in self.perturbed_factors :
        l.append(f)
    return l


  def nonperturbed_regnet_factor_list(self) :
    l = []
    for f in self.regnet_factors :
      if f not in self.perturbed_factors :
        inlist = 1
        for r in self.regulon_list :
          inlist = inlist and (f not in r.controlling_factors)
          inlist = inlist and (f not in r.coregulated_factors)
          if not inlist :
            break
        if inlist :
          l.append(f)
    return l


  def nonperturbed_cluster(self) :
    c = ClusterNode('root')
    for r in self.regulon_list :
      rc = r.cluster()
      rc.set_ancestor(c, 1.0)
    for f in self.nonperturbed_regnet_factor_list() :
      ClusterNode(f, c, 3.0)
    return c


  def complete_arrtool_basename(self) :
    return '%s_complete' % self.basename


  def nonperturbed_arrtool_basename(self) :
    return '%s_nonperturbed' % self.basename


  def transsys_fname(self) :
    return '%s.tra' % self.basename


  def dot_fname(self) :
    return '%s.dot' % self.basename


  def gnuplot_command_fname(self) :
    return '%s.gpc' % self.basename


  def complete_array_fname(self, env = None) :
    if env is None :
      return '%s.dat' % self.complete_arrtool_basename()
    else :
      return '%s_%03d.dat' % (self.complete_arrtool_basename(), env)


  def nonperturbed_array_fname(self, env = None) :
    if env is None :
      return '%s.dat' % self.nonperturbed_arrtool_basename()
    else :
      return '%s_%03d.dat' % (self.nonperturbed_arrtool_basename(), env)


  def nonperturbed_clustertree_fname(self) :
    return '%s_nonperturbed.tre' % self.basename


  def nonperturbed_clustertree_basename(self) :
    return '%s_nonperturbed' % self.basename


  def treelist_fname(self) :
    return '%s_clusterlist.txt' % self.basename


  def nonperturbed_euclidean_average_fname(self) :
    return '%s_nonperturbed_euclidean_avg.tre' % self.basename


  def nonperturbed_euclidean_single_fname(self) :
    return '%s_nonperturbed_euclidean_single.tre' % self.basename


  def nonperturbed_euclidean_complete_fname(self) :
    return '%s_nonperturbed_euclidean_complete.tre' % self.basename


  def experimentspecs_fname(self) :
    return '%s_expspec.dat' % self.basename


  def genegroups_fname(self) :
    return '%s_ggroups.dat' % self.basename


  def r_fname(self) :
    return '%s.r' % self.basename


  def write_nonperturbed_cluster(self) :
    f = open(self.nonperturbed_clustertree_fname(), 'w')
    f.write(str(self.nonperturbed_cluster()))


  def write_genegroups(self) :
    f = open(self.genegroups_fname(), 'w')
    f.write('perturbed\t%s\n' % string.join(self.perturbed_factors, '\t'))
    f.write('nonperturbed\t%s\n' % string.join(self.nonperturbed_factor_list(), '\t'))
    for r in self.regulon_list :
      f.write('%s\n' % r.genegroup_line())


  def write_r_commands(self) :
    f = open(self.r_fname(), 'w')
    f.write('library(transarr)\n')
    f.write('distfunc <- list();\n')
    f.write('distfunc[["euclidean"]] <- function(m) { dist(m, method = "euclidean", upper = TRUE); };\n')
    f.write('distfunc[["eisen"]] <- function(m) { dist.eisen(m); };\n');
    f.write('clustfunc <- list();\n');
    f.write('clustfunc[["average"]] <- function(dm) { hclust(dm, method = "average"); };\n')
    f.write('clustfunc[["single"]] <- function(dm) { hclust(dm, method = "single"); };\n')
    f.write('clustfunc[["complete"]] <- function(dm) { hclust(dm, method = "complete"); };\n')
    f.write('%sarr <- read.transarr("%s", espfile = "%s", ggfile = "%s");\n' % (self.basename, self.complete_array_fname(), self.experimentspecs_fname(), self.genegroups_fname()))
    f.write('expnames <- setdiff(names(attr(%sarr, "expspec")), "ref_init");\n'% self.basename)
    f.write('allenv <- c();\n')
    f.write('for (env in expnames)\n')
    f.write('{\n')
    f.write('  allenv <- c(allenv, attr(%sarr, "expspec")[[env]]);\n' % self.basename)
    f.write('}\n')
    f.write('attr(%sarr, "expspec")[["allenv"]] <- allenv;\n' % self.basename)
    f.write('ggnames <- "nonperturbed";\n')
    f.write('for (df in names(distfunc))\n')
    f.write('{\n')
    f.write('  for (cl in names(clustfunc))\n')
    f.write('  {\n')
    f.write('    fname <- sprintf("%s_%%s_%%s.tre", df, cl);\n' % self.nonperturbed_clustertree_basename())
    f.write('    write(fname, file = "");\n')
    f.write('    clustset <- cluster.transarr(%sarr, distfunc[[df]], clustfunc[[cl]], expnames = expnames, ggnames = ggnames);\n' % self.basename)
    f.write('    write(clusterstring(clustset), fname);\n')
    f.write('    fname <- sprintf("%s_%%s_%%s_allenv.tre", df, cl);\n' % self.nonperturbed_clustertree_basename())
    f.write('    allclust <- cluster.transarr(%sarr, distfunc[[df]], clustfunc[[cl]], expnames = "allenv", ggnames = ggnames);\n' % self.basename)
    f.write('    write(clusterstring(allclust), fname);\n')
    f.write('  }\n')
    f.write('}\n')



  def run_r(self) :
    self.write_r_commands()
    return os.system('R --no-save --slave < %s > %s' % (self.r_fname(), self.treelist_fname()))


  def run_plottree(self) :
    f = open(self.treelist_fname(), 'r')
    line = f.readline()
    while line :
      fname = line.strip()
      if fname != '' :
        if fname[-4:] != '.tre' :
          raise StandardError, 'RTAControlParameters::run_plottree: malformed file name "%s"' % fname
        psname = '%s.ps' % fname[:-4]
        cmd = 'plottree -f 10 -u -l r -B -i %s -b %s -o %s' % (self.nonperturbed_clustertree_fname(), fname, psname)
        print cmd
        os.system(cmd)
      line = f.readline()


  def adjust_histogram_parameters(self, num_factors) :
    if self.histo_nbins == 0 :
      self.histo_nbins = num_factors / 20
    if self.histo_nbins < 5 :
      self.histo_nbins = 5


  def complete(self) :
    if self.transsys_program is None :
      return 0
    if self.regnet_factors is None :
      return 0
    if self.structural_factors is None :
      return 0
    if self.perturbed_factors is None :
      return 0
    if self.perturbation_type is None :
      return 0
    if self.regulon_list is None :
      return 0
    if self.min_initial_concentration is None :
      return 0
    if self.max_initial_concentration is None :
      return 0
    if self.num_environments is None :
      return 0
    if self.samples_per_environment is None :
      return 0
    if self.basename is None :
      return 0
    return 1


  def write(self, f) :

    def write_perturbation_parameters(f, rta) :
      f.write('perturbation\n')
      f.write('perturbation_type: %s\n' % rta.perturbation_type)
      if rta.perturbation_type == 'ExpGauss' :
        f.write('perturbation_skew: %g\n' % rta.perturbation_skew)
        f.write('perturbation_dispersion: %g\n' % rta.perturbation_dispersion)
      elif rta.perturbation_type == 'ExponentialReplacement' :
        f.write('perturbation_halfmax: %g\n' % rta.perturbation_halfmax)
      elif rta.perturbation_type == 'UniformReplacement' :
        f.write('perturbation_min: %g\n' % rta.perturbation_min)
        f.write('perturbation_max: %g\n' % rta.perturbation_max)
      else :
        raise StandardError, 'RTAControlParameters::write: unknown perturbation type "%s"' % rta.perturbation_type        

    if not self.complete() :
      sys.stderr.write('RTAControlParameters::write: writing incomplete instance\n')
    f.write('min_initial_concentration: %g\n' % self.min_initial_concentration)
    f.write('max_initial_concentration: %g\n' % self.max_initial_concentration)
    f.write('num_environments: %d\n' % self.num_environments)
    write_perturbation_parameters(f, self)
    f.write('delta_t: %d\n' % self.delta_t)
    f.write('samples_per_environment: %d\n' % self.samples_per_environment)
    f.write('array_offset: %g\n' % self.array_offset)
    f.write('array_dispersion: %g\n' % self.array_dispersion)
    f.write('exp_threshold: %g\n' % self.exp_threshold)
    f.write('num_timesteps_init: %d\n' % self.num_timesteps_init)
    f.write('rndseed: %d\n' % self.rndseed)
    f.write('histo_max: %g\n' % self.histo_max)
    f.write('histo_nbins: %d\n' % self.histo_nbins)
    f.write('basename: %s\n' % self.basename)
    f.write('regnet_factors: %s\n' % string.join(self.regnet_factors, ', '))
    f.write('structural_factors: %s\n' % string.join(self.structural_factors, ', '))
    f.write('perturbed_factors: %s\n' % string.join(self.perturbed_factors, ', '))
    for r in self.regulon_list :
      f.write('regulon\n')
      f.write(str(r))
    f.write('transsys_program\n')
    f.write(str(self.transsys_program))
    f.write('transsys_program: end\n')


  def parse(self, f) :

    def parse_transsys_program(f) :
      line = f.readline()
      tp_string = ''
      while line :
        if line.strip() == 'transsys_program: end' :
          break
        tp_string = tp_string + line
        line = f.readline()
      if line == '' :
        raise StandardError, 'no end tag for transsys program -- premature end of file'
      sf = StringIO.StringIO(tp_string)
      p = transsys.TranssysProgramParser(sf)
      return p.parse_transsys()

    def parse_perturbation_spec(f, rta) :
      n, v = get_namevalue(f)
      if n != 'perturbation_type' :
        raise StandardError, 'RTAControlParameters::parse: expected perturbation_type but got "%s"' % n
      rta.perturbation_type = v
      if rta.perturbation_type == 'ExpGauss' :
        n, v = get_namevalue(f)
        if n != 'perturbation_skew' :
          raise StandardError, 'RTAControlParameters::parse: expected perturbation_skew but got "%s"' % n
        rta.perturbation_skew = float(v)
        n, v = get_namevalue(f)
        if n != 'perturbation_dispersion' :
          raise StandardError, 'RTAControlParameters::parse: expected perturbation_dispersion but got "%s"' % n
        rta.perturbation_dispersion = float(v)
      elif rta.perturbation_type == 'ExponentialReplacement' :
        n, v = get_namevalue(f)
        if n != 'perturbation_halfmax' :
          raise StandardError, 'RTAControlParameters::parse: expected perturbation_halfmax but got "%s"' % n
        rta.perturbation_halfmax = float(v)
      elif rta.perturbation_type == 'UniformReplacement' :
        n, v = get_namevalue(f)
        if n != 'perturbation_min' :
          raise StandardError, 'RTAControlParameters::parse: expected perturbation_min but got "%s"' % n
        rta.perturbation_min = float(v)
        n, v = get_namevalue(f)
        if n != 'perturbation_max' :
          raise StandardError, 'RTAControlParameters::parse: expected perturbation_max but got "%s"' % n
        rta.perturbation_max = float(v)
      else :
        raise StandardError, 'RTAControlParameters::parse: unknown perturbation type "%s"' % n

    float_re = re.compile('[+-]?([0-9]+(\\.[0-9]+)?)|(\\.[0-9]+)([Ee][+-]?[0-9]+)?')
    int_re = re.compile('[0-9]+')
    n, v = get_namevalue(f)
    while n is not None :
      if n == 'transsys_program' :
        self.transsys_program = parse_transsys_program(f)
      elif n == 'regnet_factors' :
        self.regnet_factors = parse_string_csv(v)
      elif n == 'structural_factors' :
        self.structural_factors = parse_string_csv(v)
      elif n == 'perturbed_factors' :
        self.perturbed_factors = parse_string_csv(v)
      elif n == 'regulon' :
        r = RegulonDescriptor('')
        r.parse(f)
        self.regulon_list.append(r)
      elif n == 'min_initial_concentration' :
        self.min_initial_concentration = float(v)
      elif n == 'max_initial_concentration' :
        self.max_initial_concentration = float(v)
      elif n == 'num_environments' :
        self.num_environments = int(v)
      elif n == 'perturbation' :
        parse_perturbation_spec(f, self)
      elif n == 'delta_t' :
        self.delta_t = int(v)
      elif n == 'samples_per_environment' :
        self.samples_per_environment = int(v)
      elif n == 'array_offset' :
        self.array_offset = float(v)
      elif n == 'array_dispersion' :
        self.array_dispersion = float(v)
      elif n == 'exp_threshold' :
        self.exp_threshold = float(v)
      elif n == 'num_timesteps_init' :
        self.num_timesteps_init = int(v)
      elif n == 'rndseed' :
        self.rndseed = int(v)
      elif n == 'histo_max' :
        self.histo_max = float(v)
      elif n == 'histo_nbins' :
        self.histo_nbins = int(v)
      elif n == 'basename' :
        self.basename = v
      else :
        raise StandardError, 'RTAControlParameters::parse: unknown attribute "%s"' % n
      n, v = get_namevalue(f)

