# $Id$

# $Log$
# Revision 1.2  2005/06/15 10:56:12  jtk
# changed gene structure, introduced regexp use for start/end finding
#
# Revision 1.1  2005/06/14 18:35:42  jtk
# transdecode for decoding transsys programs from "DNA" sequences
#

import re

import transsys


def baseToInt(b) :
  i = 'acgt'.find(b)
  if i == -1 :
    raise StandardError, 'bad base "%s"' % b
  return i


class BindingSite :

  def __init__(self, protein, position, strength) :
    self.protein = protein
    self.position = position
    self.strength = strength


  def arrowLine(self) :
    return (' ' * self.position) + ('^' * len(self.protein.matrix)) + (' %s @ %d (%4.2f)' % (self.protein.name, self.position, self.strength))


class DNABinder :

  def __init__(self, name, dna) :
    self.name = name
    self.matrix = []
    self.threshold = 0.0
    i = 0
    while len(dna) - i >= 5 :
      c = map(baseToInt, dna[i:i + 4])
      if c[0] == c[1] and c[0] == c[2] and c[0] == c[3] :
        i = i + 4
        break
      self.matrix.append(c)
      self.threshold = self.threshold + baseToInt(dna[i + 4]) + 1.0
      i = i + 5
    self.encoding_sequence = dna[:i]


  def __str__(self) :
    
    sa = 'a  '
    sc = 'c  '
    sg = 'g  '
    st = 't  '
    for c in self.matrix :
      sa = sa + ('%4.1f' % c[0])
      sc = sc + ('%4.1f' % c[1])
      sg = sg + ('%4.1f' % c[2])
      st = st + ('%4.1f' % c[3])
    return 'DNABinder %s, threshold = %5.1f\n%s\n%s\n%s\n%s\nencoded by: %s' % (self.name, self.threshold, sa, sc, sg, st, self.encoding_sequence)


  def bindingEnergy(self, seq) :
    e = 0.0
    if len(seq) < len(self.matrix) :
      return e
    for i in xrange(len(self.matrix)) :
      e = e + self.matrix[i][baseToInt(seq[i])]
    return e


  def bindingSite(self, seq, position = 0) :
    if len(self.matrix) == 0 :
      return None
    e = self.bindingEnergy(seq[position:])
    if e < self.threshold :
      return None
    else :
      if self.threshold > 0.0 :
        strength = e / self.threshold
      else :
        strength = 1.0
      return BindingSite(self, position, strength)


def find_binding_sites(proteome, seq) :
  blist = []
  for i in xrange(len(seq) - 1) :
    for p in proteome :
      bs = p.bindingSite(seq, i)
      if bs is not None :
        blist.append(bs)
  return blist


class RawDNAGene :

  def __init__(self, position, gene_name, product_name, geneStart, geneEnd, activatorArea, repressorArea, structuralArea) :
    self.position = position
    self.gene_name = gene_name
    self.product_name = product_name
    self.activatorArea = activatorArea
    self.repressorArea = repressorArea
    self.structuralArea = structuralArea
    self.geneStart = geneStart
    self.geneEnd = geneEnd


  def promoterArea(self) :
    return self.activatorArea + self.repressorArea


  def rawSequence(self) :
    return self.geneStart + self.promoterArea() + self.structuralArea + self.geneEnd


class TranssysDNADecoder :

  def __init__(self) :
    self.transsysName = None
    self.factorNameTemplate = None
    self.geneNameTemplate = None
    self.geneStartRE = None
    self.geneEndRE = None
    self.repressorAreaLength = None
    self.activatorAreaLength = None
    self.decay = None
    self.diffusibility = None
    self.a_spec = None
    self.a_max = None
    self.constitutive = None


  def setGeneStartRE(self, r) :
    self.geneStartRE = re.compile(r)


  def setGeneEndRE(self, r) :
    self.geneEndRE = re.compile(r)


  def rawDNAGenes(self, genome) :
    # print 'gene_index_list(%s, %s)' % (genome, gene_refpoint_tag)
    raw_gene_list = []
    m = self.geneStartRE.search(genome)
    while m :
      position = m.start()
      geneStart = m.group()
      a_start = position + len(geneStart)
      a_end = a_start + self.activatorAreaLength
      r_start = a_end
      r_end = r_start + self.repressorAreaLength
      s_start = r_end
      m = self.geneEndRE.search(genome, s_start)
      if m :
        s_end = m.start()
        geneEnd = m.group()
      else :
        s_end = len(genome)
        geneEnd = '*'
      activatorArea = genome[a_start:a_end]
      repressorArea = genome[r_start:r_end]
      structuralArea = genome[s_start:s_end]
      gene_name = self.geneNameTemplate % position
      product_name = self.factorNameTemplate % position
      raw_gene = RawDNAGene(position, gene_name, product_name, geneStart, geneEnd, activatorArea, repressorArea, structuralArea)
      raw_gene_list.append(raw_gene)
      m = self.geneStartRE.search(genome, position + 1)
    return raw_gene_list


  def decode_transsys(self, genome) :
    raw_gene_list = self.rawDNAGenes(genome)
    proteome = []
    factor_list = []
    for rg in raw_gene_list :
      p = DNABinder(rg.product_name, rg.structuralArea)
      proteome.append(p)
      decay_expr = transsys.ExpressionNodeValue(self.decay)
      diffusibility_expr = transsys.ExpressionNodeValue(self.diffusibility)
      factor = transsys.Factor(rg.product_name, decay_expr, diffusibility_expr)
      factor.comments.append('decoded to DNABinder:')
      for l in str(p).split('\n') :
        factor.comments.append(l)
      factor_list.append(factor)
      # print '%d: struct. area = "%s"' % (i, structArea)
      # print p
    gene_list = []
    for rg in raw_gene_list :
      activators = find_binding_sites(proteome, rg.activatorArea)
      repressors = find_binding_sites(proteome, rg.repressorArea)
      promoter = [transsys.PromoterElementConstitutive(transsys.ExpressionNodeValue(self.constitutive))]
      for a in activators :
        a_spec = transsys.ExpressionNodeValue(self.a_spec)
        a_max = transsys.ExpressionNodeValue(self.a_max)
        promoter.append(transsys.PromoterElementActivate(a_spec, a_max, [a.protein.name]))
      for r in repressors :
        r_spec = transsys.ExpressionNodeValue(self.a_spec)
        r_max = transsys.ExpressionNodeValue(self.a_max)
        promoter.append(transsys.PromoterElementRepress(r_spec, r_max, [r.protein.name]))
      gene = transsys.Gene(rg.gene_name, rg.product_name, promoter)
      gene.comments.append('gene:            %s' % rg.rawSequence())
      gene.comments.append('promoter area:   %s' % (' ' * len(rg.geneStart)) + rg.promoterArea())
      gene.comments.append('structural area: %s' % (' ' * (len(rg.geneStart) + len(rg.promoterArea())) + rg.structuralArea))
      
      gene.comments.append('activator area:')
      gene.comments.append(rg.activatorArea)
      for a in activators :
        gene.comments.append(a.arrowLine())
      gene.comments.append('repressor area:')
      gene.comments.append(rg.repressorArea)
      for r in repressors :
        gene.comments.append(r.arrowLine())
      gene_list.append(gene)
    return transsys.TranssysProgram(self.transsysName, factor_list, gene_list, resolve = True)
