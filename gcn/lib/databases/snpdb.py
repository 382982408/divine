#
# COPYRIGHT (C) 2012-2013 TCS Ltd
#
"""
.. module:: snpdb
    :platform: Unix, Windows, MacOSX
    :synopsis: Class for accessing SNP information from an sqlite
               database

.. moduleauthor:: Rajgopal Srinivasan (rajgopal.srinivasan@gmail.com)

Class for accessing variant information from an sqlite database.  This is
primarily intended for reference databases of variants such as dbSNP,
1000 Genomes, etc.
"""
from gcn.etc import dbconfig
from gcn.lib.io import db
import os
from html5lib.treewalkers._base import UNKNOWN

BENIGN,UNKNOWN2,PATHOGENIC = range(-1,2)
clitag = {BENIGN:'benign',UNKNOWN2:'vus',PATHOGENIC:'pathogenic'}
        
def classify_clnsig(clinvar_sigs):
  
  if '2' in clinvar_sigs:
    if '5' in clinvar_sigs:
      return UNKNOWN2
    else:
      return BENIGN
  elif '5' in clinvar_sigs:
    return PATHOGENIC
  elif '3' in clinvar_sigs:
    if '4' in clinvar_sigs:
      return UNKNOWN2
    else:
      return BENIGN
  elif '4' in clinvar_sigs:
    return PATHOGENIC
  else:
    return UNKNOWN2

class SNPDB(db.DB):
    """Class to retrieve SNP information from a SNP database

    The class makes the following assumptions

    1. The name of the table is snps
    2. The table defines the following columns
        a) chrom
        b) pos
        c) id
        d) ref
        e) alt

    The class is initialized with the name of the database to connect to.
    The name should be defined in `gcn.etc.dbconfig.DBCONFIG`. Alternatively
    the full path qualified name of the database may be passed
    """

    def __init__(self, name):
        """Class initialization

        Argument

            name (string): Name of the database or filename for the database
        """
        self.name = name
        super(SNPDB, self).__init__()
        if name in dbconfig.DBCONFIG:
            self.load(name=name)
        elif os.path.exists(name):
            self.load(db=name)
        else:
            raise ValueError('No such database %s' % name)

    def snp_by_location(self, chromosome, position):
        """Retrieve details of a SNP associated with a given
        chromosomal location.

        Args:
            chromosome (str):   Chromosome name in the form chrN where
                                N=1-23,X,Y
            position (integer): Location on the chromosome, using 0-based
                                indexing

        Returns:
            SNP associated with the position as an instance of a `SNP`
            namedtuple
        """
        stmt = 'SELECT * FROM snps WHERE chrom=(?) AND pos=(?)'
        if chromosome[:3] == 'chr':
            chromosome = chromosome[3:]
        if hasattr(self,'SNP'):
          C = self.SNP
        else:
          C = self.SNP = self.make_namedtuple('snps', 'SNP')

        return [C._make(r) for r in self.execute(stmt, (chromosome, position))]

    def has_snp(self, chromosome, position, ref=None, alt=None):
        """Check if a SNP is present in the database given postion
        and location. Optionally the ref and alt allele may also be
        specified.

        Args:
            chromosome(str): Chromosome name in the form chrN where
                             N=1-23,X,Y
            position (integer): Location on the chromosome, using 0-based
                                indexing
            ref (str): Reference Allele (Optional)
            alt (str): Alternate Allele (Optional)
        Returns:
            True if the SNP is present and False otherwise
        """
        if chromosome[:3] == 'chr':
            chromosome = chromosome[3:]

        args = (chromosome, position)
        stmt = 'SELECT * FROM snps WHERE chrom=(?) AND pos=(?)'
        if ref is not None:
            stmt += ' AND ref=(?)'
            args = args + (ref,)
        if alt is not None:
            stmt += ' AND alt=(?) limit 1'
            args = args + (alt,)

        results = self.execute(stmt, args).fetchone()
        if results:
            return True
        else:
            return False

    def get_meta_info(self, attrib_id):
        """Returns the meta data for the given attribute id.

        Args:
            attrib_id(str): VCF Attribute Identifier
                     Eg. AF, GMAF, OM etc
        Returns:
            Named Tuple comprising of the fields -
            id, count, type, description
        """
        TYPEMAP = {'int': 'Integer',
                   'float': 'Float',
                   'text': 'String'}

        meta_info = self.make_namedtuple('info', 'metainfo')
        stmt = 'SELECT * FROM info WHERE id=(?)'
        c = self.conn.cursor()
        r = c.execute(stmt, ('info' + '_' + attrib_id, )).fetchall()
        if r:
            idfr, idfr_type, idfr_count, idfr_desc = r[0]
            ctype = TYPEMAP[idfr_type]
            if ctype == 'Integer':
                if idfr_count == '0':
                    ctype = 'Flag'
            return meta_info._make((attrib_id, ctype, idfr_count, idfr_desc))
        else:
            return

    def iterate_db(self, chromosome):
        """Retrieve the SNPs present in the given chromosome.
        Args:
            chromosome (str):   Chromosome name in the form chrN where
                                N=1-23,X,Y
        Returns:
            SNP associated with the position as an instance of a `SNP`
            namedtuple
        """
        stmt = "SELECT * FROM snps WHERE chrom='%s'"
        if chromosome[:3] == 'chr':
            chromosome = chromosome[3:]

        if hasattr(self,'SNP'):
          C = self.SNP
        else:
          C = self.SNP = self.make_namedtuple('snps', 'SNP')

        for r in self.execute(stmt % chromosome):
            yield C._make(r)

    def is_clinically_associated(self, rsid):
        """This method is valid only for querying Clinvardb.
        Check if a SNP is present in the database and annotated with
        disease causing annotation. For this the 'CLNDBN' tag is checked
        if its value is not '.' and instead hold a disease name.

        Args:
            rsid(str): dbSNP ID
        Returns:
            True if the SNP is disease associated and False otherwise
        """
        stmt = "SELECT * FROM snps WHERE id='%s' and info_CLNDBN != '.'"\
                                                                 % (rsid)
        results = self.execute(stmt).fetchone()
        if results:
            return True
        else:
            return False
      
    def select_all_hgmd(self):
      """ This method is valid only for querying Clinvardb.
      select variant key and clnsig (benign or pathogenic only)
      """
      
      if self.name!='HGMDDB':
        raise EnvironmentError('this method only supports for hgmddb')
      outs = []
   
      stmt = "SELECT chrom, pos, id, ref, alt, info_VC from snps"
      
      for chrm,pos,id,ref,alt,clisig in self.execute(stmt).fetchall():
        if clisig in ['FP','R']:
          sig='benign_%s'%self.name
        elif clisig in ['DM','DM?']:
          sig='pathogenic_%s'%self.name
        else:
          sig='vus_%s'%self.name
        outs.append([chrm,pos,id,ref,alt,sig])
      return outs
            
    def select_all_clinvar(self):
      """ This method is valid only for querying Clinvardb.
      select variant key and clnsig (benign or pathogenic only)
      """
      outs = []
      if self.name!='CLINVARDB':
        raise EnvironmentError('this method only supports for clinvardb')

      stmt = "SELECT chrom, pos, id, ref, alt, info_CLNSIG from snps"

      for chrm,pos,id,ref,alt,clisig in self.execute(stmt).fetchall():
        sigclass = classify_clnsig(clisig.split('|'))
        #print '%s|%s'%(sigclass,self.name) #debug
        sig='%s_%s'%(clitag[sigclass],self.name)
        outs.append([chrm,pos,id,ref,alt,sig])
      return outs
    
    def select_snps(self,min_maf,max_maf,sample_rate=1.0,snp_tag='.'):
      import random
      sampling = False
      
      num2 = {}
      if sample_rate < 1.0:
        sampling = True
        stmt = "select count(*) from snps where info_AF > %g and info_AF < %g" % (min_maf,max_maf)
        outs = self.execute(stmt).fetchone()
        R = int(outs[0])
        R1 = int(round(sample_rate*R))
        print 'sampling %d out of %d'%(R1,R)
        nums = range(R)
        random.shuffle(nums)
        nums = nums[:R1]
        num2 = {}
        for j in nums:
          num2[j]=True
      
      stmt = "select chrom, pos, id, ref, alt, info_AF from snps where info_AF > %g and info_AF < %g" % (min_maf,max_maf)
      r = 0
      outs = []
      print 'selecting ...'
      for chrom,pos,id,ref,alt,maf in self.execute(stmt).fetchall():
        if not sampling or (sampling and r in num2):
          outs.append([chrom,pos,id,ref,alt,snp_tag])
        r += 1
      print 'done.'
      return outs
    
    def get_hgmd_genes(self, clnsig2sels=['DM','DM?']):
        """ This method return (genes) for HGMD
        Args:
          type (str): types of clinical significance
        Returns:
          genes (tuple): gene1,gene2,...
        """
        if self.name!='HGMDDB':
         raise EnvironmentError('this method only supports for hgmddb')
        
        genes = []
        stmt = "SELECT info_Gene, info_VC FROM snps"
        results = self.execute(stmt).fetchall()
        if results:
          for gene,vc in results:
            if vc in clnsig2sels:
              if gene not in genes:
                genes.append(gene)
        return tuple(genes)

    def get_clinvar_genes(self, clnsig2sels=['pathogenic']):
        """ This method return (genes) for Clinvar
          Args:
          type (str): types of clinical significance
         Returns:
          genes (tuple): gene1,gene2,...
        """

        if self.name!='CLINVARDB':
            raise EnvironmentError('this method only supports for clinvardb')
        
        genes = []
        stmt = "SELECT info_CLNSIG, info_GENEINFO FROM snps"
        results = self.execute(stmt).fetchall()
        if results:
            for clisig, geneinfo in results:
                if clisig and geneinfo:
                    sigclass = classify_clnsig(clisig.split('|'))
                    if clitag[sigclass] in clnsig2sels:
                       gene = geneinfo.split(':')[0]
                       if gene not in genes:
                          genes.append(gene)
        return tuple(genes)
       
    def get_snp_coord(self, snpid):
        """This method return (chrom, pos, ref, alts) for the given
        snpid.

        Args:
            snpid(str): dbSNP ID
        Returns:
            snp_coord(tuple): (chrom, pos, ref, [alt1, alt2, ])
        """
        stmt = "SELECT chrom, pos, ref, alt FROM snps WHERE id='%s' limit 1"\
                                                                % (snpid)
        results = self.execute(stmt).fetchone()
        if results:
            alts = results[3].split(',')
            results = list(results)
            results[3] = alts
            return tuple(results)
        else:
            return None


# for now the following are implemented as functions. If in the future
# there is a need to define methods specific to these databases
# these can be converted to classes while still supporting the current
# API


def KGDB():
    """Return an instance of `SNPDB` configured to search 1000 Genomes
    Variations"""
    return SNPDB('KGDB')


def DBSNP():
    """Return an instance of `SNPDB` configured to search the NCBI dbSNP
    database"""
    return SNPDB('DBSNP')


def ClinvarDB():
    """Return an instance of `SNPDB` configured to search the NCBI
    clinvar variation database"""
    return SNPDB('CLINVARDB')

def HgmdDB():
    """Return an instance of `SNPDB` configured to search HGMD professional database"""
    return SNPDB('HGMDDB')

def ClinvitaeDB():
    """Return an instance of `SNPDB` configured to search CLINVITAE database"""
    return SNPDB('CLINVITAE')
   
def EspDB():
    """Return an instance of `SNPDB` configured to search the Exome
    Sequencing Project's(ESP) database"""
    return SNPDB('ESP')

def ExACDB():
    """Return an instance of 'SNPDB' configured to search the ExAC database"""
    return SNPDB('EXAC')

if __name__ == '__main__':
    clnvar = ClinvarDB()
    print clnvar.get_version()
    print clnvar.has_snp('16', 730490, 'C', 'T')
    print clnvar.is_clinically_associated('rs41307846')
    
    dbsnp = DBSNP()
    print dbsnp.get_version()
    for rec in dbsnp.iterate_db('1'):
        print rec
        break
    print dbsnp.snp_by_location('chr1', 2337967)
    print dbsnp.get_snp_coord('rs59861892')
    for rec in dbsnp.snp_by_location('1', 14907):
        print rec.dbSNPBuildID
    espdb = EspDB()
    print espdb.get_version()
    print espdb.has_snp('1', 861357, 'C', 'G')
    print espdb.snp_by_location('1', 861357)[0].MAF
    
    #TODO: test ExAC
    
    kgdb = KGDB()
    print kgdb.get_version()
    print kgdb.snp_by_location('14', 94847262)
