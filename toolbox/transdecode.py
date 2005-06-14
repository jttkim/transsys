# $Id$

# $Log$
# Revision 1.1  2005/06/14 18:35:42  jtk
# transdecode for decoding transsys programs from "DNA" sequences
#


import transsys


def baseToInt(b) :
  i = 'acgt'.find(b)
  if i == -1 :
    raise StandardError, 'bad base "%s"' % b
  return i


def gene_index_list(genome, gene_refpoint_tag) :
  # print 'gene_index_list(%s, %s)' % (genome, gene_refpoint_tag)
  l = []
  i = genome.find(gene_refpoint_tag)
  while i > -1 :
    l.append(i)
    i = genome.find(gene_refpoint_tag, i + 1)
  return l


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
      self.threshold = self.threshold + baseToInt(dna[i + 4])
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
    e = self.bindingEnergy(seq[position:])
    if e <= self.threshold :
      return None
    else :
      strength = e / self.threshold
      return BindingSite(self, position, strength)


def find_binding_sites(proteome, seq) :
  blist = []
  for i in xrange(len(seq) - 1) :
    for p in proteome :
      bs = p.bindingSite(seq, i)
      if bs is not None :
        blist.append(bs)
  return blist


class TranssysDNADecodeParameters :

  def __init__(self) :
    self.transsysName = None
    self.factorNameTemplate = None
    self.geneNameTemplate = None
    self.polymeraseBindingWord = None
    self.polymeraseTerminator = None
    self.repressorAreaLength = None
    self.activatorAreaLength = None
    self.decay = None
    self.diffusibility = None
    self.a_spec = None
    self.a_max = None
    self.constitutive = None


def decode_transsys(genome, dp) :
  gi_list = gene_index_list(genome, dp.polymeraseBindingWord)
  proteome = []
  factor_list = []
  for i in gi_list :
    factor_name = dp.factorNameTemplate % i
    structArea = genome[i + dp.repressorAreaLength:]
    # print genome
    # print i + dp.repressorAreaLength
    p = DNABinder(factor_name, structArea)
    proteome.append(p)
    decay_expr = transsys.ExpressionNodeValue(dp.decay)
    diffusibility_expr = transsys.ExpressionNodeValue(dp.diffusibility)
    factor = transsys.Factor(factor_name, decay_expr, diffusibility_expr)
    factor.comments.append('decoded to DNABinder:')
    for l in str(p).split('\n') :
      factor.comments.append(l)
    factor_list.append(factor)
    # print '%d: struct. area = "%s"' % (i, structArea)
    # print p
  gene_list = []
  for i in gi_list :
    gene_name = dp.geneNameTemplate % i
    product_name = dp.factorNameTemplate % i
    promoterArea = genome[max(0, i - dp.activatorAreaLength):i + len(dp.polymeraseBindingWord) + dp.repressorAreaLength]
    activatorArea = genome[max(0, i - dp.activatorAreaLength):i]
    repressorArea = genome[i + len(dp.polymeraseBindingWord):i + len(dp.polymeraseBindingWord) + dp.repressorAreaLength]
    # print 'gene %s: act. area = "%s", rep. area = "%s"' % (gene_name, activatorArea, repressorArea)
    activators = find_binding_sites(proteome, activatorArea)
    repressors = find_binding_sites(proteome, repressorArea)
    promoter = [transsys.PromoterElementConstitutive(transsys.ExpressionNodeValue(dp.constitutive))]
    for a in activators :
      a_spec = transsys.ExpressionNodeValue(dp.a_spec)
      a_max = transsys.ExpressionNodeValue(dp.a_max)
      promoter.append(transsys.PromoterElementActivate(a_spec, a_max, [a.protein.name]))
    for r in repressors :
      r_spec = transsys.ExpressionNodeValue(dp.a_spec)
      r_max = transsys.ExpressionNodeValue(dp.a_max)
      promoter.append(transsys.PromoterElementRepress(r_spec, r_max, [r.protein.name]))
    gene = transsys.Gene(gene_name, product_name, promoter)
    gene.comments.append('promoter area: %s' % promoterArea)
    gene.comments.append('activator area:')
    gene.comments.append(activatorArea)
    for a in activators :
      gene.comments.append(a.arrowLine())
    gene.comments.append('repressor area:')
    gene.comments.append(repressorArea)
    for r in repressors :
      gene.comments.append(r.arrowLine())
    gene_list.append(gene)
  return transsys.TranssysProgram(dp.transsysName, factor_list, gene_list, resolve = True)
