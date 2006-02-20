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
  """Represents a DNA binding factor using a scoring matrix and a threshold."""

  def __init__(self, name, dna, thresholdOffset) :
    self.name = name
    self.decode(dna, thresholdOffset)


  def decode(self, dna, thresholdOffset) :
    """Decode a string of dna bases to a matrix and a threshold.
The matrix is decoded by taking words of length 5 from the encoding
sequence. If the first four characters are all identical, decoding
stops. Otherwise, they are mapped to integer values in {0, 1, 2, 3}
and appended to the matrix as a further column. The 5th character is
mapped to an integer, which is and added to the threshold. The threshold
is further controlled by the threshold offset, which is also added to
the threshold for each position in the matrix, thus allowing some
control over the binder's specificity."""
    self.threshold = 0.0
    self.matrix = []
    i = 0
    while len(dna) - i >= 5 :
      c = map(baseToInt, dna[i:i + 4])
      if c[0] == c[1] and c[0] == c[2] and c[0] == c[3] :
        i = i + 4
        break
      self.matrix.append(c)
      self.threshold = self.threshold + baseToInt(dna[i + 4]) + thresholdOffset
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
  """Represents a raw gene, i.e. the sequences of the gene's elements, which
are the gene start (promoter), the activator binding area, the repressor binding
area, the structural area, and the end (transcriptional terminator). The gene is
'raw' in the sense that this class holds just the sequences, not any data
structures derived from them, such as DNABinder etc."""
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
  """Decodes genome sequences into transsys programs.
The overall gene structure is:
...|----activator area-----|----repressor area ----|-geneStart-|--structural area--|-geneEnd...
   |<-activatorAreaLength->|<-repressorAreaLength->|
The structural area length is variable. Gene start and end can also be variable,
as enabled by regular expressions.

This decoder uses the following constants:
    * geneStartRE, geneEndRE: Regular expressions to determine start and end of
      a gene
    * repressorAreaLength: length of area where binding represses
    * activationAreaLength: length of area where binding activates
    * decay, diffusibility: constants for factor construction
    * a_spec, a_max: constants for promoter element construction
      (used for both activate and repress elements)
    * constitutive: constant for constitutive promoter element construction
Notice that strength of activation / repression does not depend on binding
strength in this decoder."""

  savefile_magic = 'TranssysDNADecoderParameters-1.1'

  def __init__(self) :
    self.thresholdOffset = None
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


  def __str__(self) :
    s = 'thresholdOffset: %g\n' % self.thresholdOffset
    s = s + 'factorNameTemplate: %s\n' % self.factorNameTemplate
    s = s + 'geneNameTemplate: %s\n' % self.geneNameTemplate
    s = s + 'geneStartRE: %s\n' % self.geneStartRE.pattern
    s = s + 'geneEndRE: %s\n' % self.geneEndRE.pattern
    s = s + 'repressorAreaLength: %d\n' % self.repressorAreaLength
    s = s + 'activatorAreaLength: %d\n' % self.activatorAreaLength
    s = s + 'decay: %g\n' % self.decay
    s = s + 'diffusibility: %g\n' % self.diffusibility
    s = s + 'a_spec: %g\n' % self.a_spec
    s = s + 'a_max: %g\n' % self.a_max
    s = s + 'constitutive: %g' % self.constitutive
    return s


  def write(self, f) :
    f.write('%s\n' % self.savefile_magic)
    f.write(str(self))
    f.write('\n')


  def parse(self, f) :
    line = f.readline()
    if line.strip() != self.savefile_magic :
      raise StandardError, 'TranssysDNADecode::parse: bad magic "%s"' % line.strip()
    self.thresholdOffset = transsys.parse_float(f, 'thresholdOffset')
    self.factorNameTemplate = transsys.parse_string(f, 'factorNameTemplate').strip()
    self.geneNameTemplate = transsys.parse_string(f, 'geneNameTemplate').strip()
    self.setGeneStartRE(transsys.parse_string(f, 'geneStartRE').strip())
    self.setGeneEndRE(transsys.parse_string(f, 'geneEndRE').strip())
    self.repressorAreaLength = transsys.parse_int(f, 'repressorAreaLength')
    self.activatorAreaLength = transsys.parse_int(f, 'activatorAreaLength')
    self.decay = transsys.parse_float(f, 'decay')
    self.diffusibility = transsys.parse_float(f, 'diffusibility')
    self.a_spec = transsys.parse_float(f, 'a_spec')
    self.a_max = transsys.parse_float(f, 'a_max')
    self.constitutive = transsys.parse_float(f, 'constitutive')


  def setGeneStartRE(self, r) :
    self.geneStartRE = re.compile(r)


  def setGeneEndRE(self, r) :
    self.geneEndRE = re.compile(r)


  def rawDNAGenes(self, genome) :
    """separate a genome into raw genes."""
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


  def decode_transsys(self, transsys_name, genome) :
    """construct a transsys program by decoding a genome sequence. The
process is: (1) split the genome into sequence portions that represent
genes, based on the gene start and end regular expressions. (2) For each
Gene, translate the structural parts of a gene into a DNABinder; the set
of all DNABinders constitutes the proteome. (3) For each gene, determine
the multiset of DNABinders that bind in the activating and the repressing
regions, respectively. (4) Construct a transsys program of the genes and
the proteome found, with promoters constructed based on (3)."""
    raw_gene_list = self.rawDNAGenes(genome)
    proteome = []
    factor_list = []
    for rg in raw_gene_list :
      p = DNABinder(rg.product_name, rg.structuralArea, self.thresholdOffset)
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
    return transsys.TranssysProgram(transsys_name, factor_list, gene_list, resolve = True)
