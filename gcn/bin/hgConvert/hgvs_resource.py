"""
.. module:: converter
    :platform: Unix, Windows, MacOSX
    :synopsis: hg coordinate version convertor

.. moduleauthor:: changjin.hong@gmail.com

This load ref genome and refseq tx for HGVS nomenclature.

import pyhgvs
import pyhgvs.utils
from pygr.seqdb import SequenceFileDB

refseq_acc='NM_000492'
chrom='chr7'
offset=117307123
ref='AGAG'
alt='A'

genome_fn = 'hg19.fa'
refseq_fn = 'refGene.txt'

genome = SequenceFileDB(genome_fn)
fp = open(refseq_fn,'r')
refseq = pyhgvs.utils.read_transcripts(fp)
fp.close()

transcript = refseq.get(refseq_acc)
hgvs_name = pyhgvs.format_hgvs_name(chrom, offset, ref, alt, genome, transcript)

"""

import os
import pyhgvs
import pyhgvs.utils
from gcn.etc import fileconfig
from pygr.seqdb import SequenceFileDB
from gcn.lib.databases.refgene import get_ucsc_chrom

CHROMOSOMES, ICHROMOSOMES = get_ucsc_chrom()

class Hgvs2:

    def __init__(self):
        """ register hg19 reference genome sequence and NCBI RefSeq transcript coordinates"""
        self.genome_fn = fileconfig.FILECONFIG['REFGENOME_UCSC']
        self.genome = None
        if not os.path.exists(self.genome_fn):
          raise IOError('Reference Genome Sequence (UCSC format) for %s is not found'%self.genome_fn)

        self.refseq_fn = fileconfig.FILECONFIG['REFGENE']
        self.refseq = None
        if not os.path.exists(self.refseq_fn):
          raise IOError('NCBI RefSeq transcript for %s is not found'%self.refseq_fn)
        
    def get_transcript(self, name):
        return self.refseq.get(name)
  
    def load_resource(self):
        #load genome sequence
        print 'loading the genome sequence [%s] for HGVS...' % self.genome_fn
        self.genome = SequenceFileDB(self.genome_fn)
        print 'done.'
        
        #load refseq into dic
        print 'loading the refseq transcript [%s] for HGVS...' % self.refseq_fn
        fp = open(self.refseq_fn,'r')
        self.refseq = pyhgvs.utils.read_transcripts(fp)    
        fp.close()
        print 'done.'
    
    def close_resource(self):
      self.genome.close()
    
    def to_cDNA(self, chrom, offset, ref, alt, refseq_acc):
        """ convert to HGVS nomenclature """
        transcript = self.get_transcript(refseq_acc)
        
        if not chrom.startswith('chr'):
          chrom = 'chr%s'%chrom
          if chrom not in CHROMOSOMES:
            return None
          
        if not chrom in self.genome.keys():
          return None
          
        hgvs_name = pyhgvs.format_hgvs_name(chrom, offset, ref, alt, self.genome, transcript)
        if hgvs_name:
          itms = hgvs_name.split(':')
          if len(itms)>1:
            return itms[1]
          else:
            return hgvs_name
        else:
          return None

if __name__ in "__main__":
    hgvs2 = Hgvs2()
    
    chrom, offset, ref, alt, refseq_acc = ['chr1', 897325,'G','C','NM_198317']
    locus,nomenclature = hgvs2.to_cDNA(chrom, offset, ref, alt, refseq_acc)
    
    if locus:
      print locus
      print nomenclature

    chrom, offset, ref, alt, refseq_acc = ['chr1', 236899899,'TC','T','NM_001278344']
    locus,nomenclature = hgvs2.to_cDNA(chrom, offset, ref, alt, refseq_acc)
    
    if locus:
      print locus
      print nomenclature
   
