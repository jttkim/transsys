#!/usr/bin/env python

# $Id$

# $Log$
# Revision 1.1  2005/03/08 17:12:02  jtk
# Initial revision
#
# Revision 1.3  2003/02/14 16:31:27  kim
# regstruct project is nearing completion. regstruct_transsys for
#     generation of transsys programs and regstruct_transarray for
#     array simulation are (more or less) properly implemented.
#     Analysis with R has been assembled in microarr.r
#
# Revision 1.2  2003/02/04 23:42:17  kim
# made experiment_plots function from transarray a method of MicroarrayData
#
# Revision 1.1  2003/01/28 21:12:05  kim
# initial toolbox assembly
#


import sys
import string
import types
import re
import copy


def cmp_by_score(x1, x2) :
  if x1[1] < x2[1] :
    return -1
  elif x1[1] > x2[1] :
    return 1
  return 0


class ExperimentSpec :

  def __init__(self, col_index, time_list, title = '(untitled)') :
    if len(col_index) != len(time_list) :
      raise StandardError, 'ExperimentSpec::__init__: lengths of column_index (%d) and time_list (%d) differ' % (len(col_index), len(time_list))
    self.column_index = copy.deepcopy(col_index)
    self.time_list = copy.deepcopy(time_list)
    self.title = title


  def num_columns(self) :
    return len(self.time_list)


  def gene_plot_data(self, arrdata, row_index) :
    row = arrdata.data[row_index]
    pd = []
    for i in range(self.num_columns()) :
      ci = self.column_index[i]
      if type(row[ci]) is types.FloatType :
        pd.append([self.time_list[i], row[ci]])
    return pd


  def write_gnuplot_file(self, f, arrdata) :
    for row_index in range(arrdata.num_rows()) :
      pd = self.gene_plot_data(arrdata, row_index)
      f.write('# %s\n' % arrdata.data[row_index][0])
      if len(pd) > 0 :
        for p in pd :
          f.write('%f %f\n' % (p[0], p[1]))
      f.write('\n')
          

class MicroarrayData :

  def __init__(self) :
    self.column_label = []
    self.data = []
    self.gene_dict = {}


  def setup_gene_dict(self) :
    self.gene_dict = {}
    for i  in range(len(self.data)) :
      gene_id = self.data[i][0]
      if type(gene_id) is not types.StringType :
        raise StandardError, 'MicroarrayData::setup_gene_dict: non-string gene ID'
      if gene_id == '' :
        raise StandardError, 'MicroarrayData::setup_gene_dict: empty gene ID'
      if gene_id in self.gene_dict.keys() :
        raise StandardError, 'MicroarrayData::setup_gene_dict: duplicate gene ID %s' % gene_id
      self.gene_dict[gene_id] = i


  def num_rows(self) :
    return len(self.data)


  def num_columns(self) :
    return len(self.column_label)


  def read_file(self, f) :

    def get_line(f) :
      line = f.readline()
      while line != '' :
        if line.strip()[0] == '#' :
          line = f.readline()
        else :
          break
      return line

    int_re = re.compile('[+-]?[0-9]+')
    float_re = re.compile('[+-]?([0-9]+(\\.[0-9]+)?)|(([0-9]+)?\\.[0-9]+)([Ee][+-]?[0-9]+)?')
    line = get_line(f)
    if line == '' :
      raise StandardError, 'MicroarrayData::read_file: no header line'
    if line[-1] == '\n' :
      line = line[:-1]
    self.column_label = line.split('\t')
    line = get_line(f)
    while line != '' :
      if line[-1] == '\n' :
        line = line[:-1]
      raw_vlist = line.split('\t')
      vlist = []
      for raw_v in raw_vlist :
        m = float_re.match(raw_v)
        if m :
          vlist.append(string.atof(raw_v))
        else :
          vlist.append(raw_v)
      self.data.append(vlist)
      line = get_line(f)
    self.setup_gene_dict()


  def column_type(self, col) :
    for row in self.data :
      if len(row) < col :
        raise StandardError, 'MicroarrayData::column_type: column index %d out of range' % col
      if type(row[col]) is not types.FloatType :
        return types.StringType
    return types.FloatType


  def write_gnuplot_file(self, f) :
    floatrows = []
    for i in xrange(len(self.column_label)) :
      if self.column_type is types.FloatType :
        f.write('# %d: %s\n' % (i, self.column_label[i]))
        floatrows.append(i)
    for row in self.data :
      for i in floatrows :
        f.write('%1.12g ' % row[i])
      f.write('\n')


  def regulation_index(self, row_index) :
    n = 0
    s = 0.0
    for x in self.data[row_index] :
      if type(x) is types.FloatType :
        n = n + 1
        s = s + abs(x)
    if n == 0 :
      return None
    return s / n


  def regulation_index_list(self) :
    l = []
    for row_index in range(self.num_rows()) :
      r = self.regulation_index(row_index)
      if r is not None :
        l.append((self.data[row_index][0], r))
    l.sort(cmp_by_score)
    return l


  def experiment_plots(self, gpcfile, espec_list, histo_max, histo_nbins, basename, cmp_basename = None) :
    """write a gnuplot file containing the expression profiles of all genes, and
also a histogram showing a profile of regulatory strengths"""
    for e in espec_list :
      fname = '%s_%s.plt' % (basename, e.title)
      f = open(fname, 'w')
      e.write_gnuplot_file(f, self)
      f.close()
      gpcfile.write('plot \'%s\' with linespoints' % fname)
      if cmp_basename :
        cmp_fname = '%s_%s.plt' % (cmp_basename, e.title)
        gpcfile.write(', \'%s\' with linespoints\n' % cmp_fname)
      gpcfile.write('\n')
      gpcfile.write('pause -1 \'Hit return\'\n')
    fname = '%s_strength.plt' % basename
    f = open(fname, 'w')
    f.write('# regulation index values\n')
    histogram = histo_nbins * [0]
    for r in self.regulation_index_list() :
      hi = int((r[1] * histo_nbins) / histo_max)
      if hi < histo_nbins :
        histogram[hi] = histogram[hi] + 1
      else :
        sys.stderr.write('arraydata_plots: regulation strength %f out of histogram range %f: %d >= %d\n' % (r[1], histo_max, hi, histo_nbins))
      f.write('# %s\n' % r[0])
      f.write('%f\n' % r[1])
    f.close()
    gpcfile.write('plot \'%s\' with boxes' % fname)
    #    if cmp_basename :
    #      cmp_fname = '%s_strength.plt' % cmp_basename
    #      gpcfile.write(', \'%s\' with impulses\n' % cmp_fname)
    gpcfile.write('\n')
    gpcfile.write('pause -1 \'Hit return\'\n')
    fname = '%s_shist.plt' % basename
    f = open(fname, 'w')
    f.write('# regulation strength histogram\n')
    for i in xrange(histo_nbins) :
      f.write('%d %f %d\n' % (i, (i * histo_max) / histo_nbins, histogram[i]))
    f.close()
    gpcfile.write('plot \'%s\' using 2:3 with boxes\n' % fname)
    gpcfile.write('pause -1 \'Hit return\'\n')
