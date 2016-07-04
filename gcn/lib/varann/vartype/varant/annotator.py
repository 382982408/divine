"""
.. module:: annotator
    :platform: Unix, Windows, MacOSX
    :synopsis: SNP Annotation

.. moduleauthor:: Kunal Kundu (kunal.kundu@tcs.com); modified by changjin.hong@gmail.com

This modules is the base module for annotating a given SNP
(Chr, Position, Ref allele, Obs allele). It calls other supporting
modules (annotateRegion, annotateSequence, annotateSite, annotateMutation)
to generate the annotaion.
"""

from annotateRegion import RegionAnnotation
from annotateSequence import SequenceAnnotation
from gcn.data.complement import COMPLEMENT
from annotateMutation import MutationAnnotation
from annotateLCR import LCR
from gcn.lib.databases.refgene import Refgene
from gcn.lib.databases.refgene import get_ucsc_chrom, to_ucsc_chrom
from bisect import bisect
from gcn.lib.databases.omim import Omim
from gcn.lib.databases.clnphe import Clnsnp
from annotateGerp import Gerp
from gcn.lib.databases.nsfp import Effpred
from gcn.lib.databases.utr import UTRdb
from gcn.lib.databases.regulomedb import Regulomedb
from gcn.lib.databases.mirna import MiRna
from gcn.lib.databases.splicedb import Splicedb
from gcn.lib.databases.snpdb import ClinvarDB
from gcn.lib.databases.snpdb import ClinvitaeDB
from gcn.lib.databases.snpdb import HgmdDB
from gcn.lib.databases.snpdb import DBSNP
from gcn.lib.databases.snpdb import KGDB
from gcn.lib.databases.snpdb import EspDB
from gcn.lib.databases.snpdb import ExACDB
from annotateNimblegenCap import NimblegenCapture
from gcn.data.pseudoautosomal_genes import PSEUDO_AUTO_GENES
from gcn.etc.dbconfig import DBCONFIG
from gcn.lib.io.vcf import VCFParser
from gcn.lib.io.vcfutils import normalize_variant
from gcn.bin.hgConvert.hgvs_resource import Hgvs2
from collections import namedtuple
import os,sys
import time

CHROMOSOMES, ICHROMOSOMES = get_ucsc_chrom()

class SNPAnnotation():

    def __init__(self,capture_kit_name='SeqCapEZ_Exome',probe_flanking_bp=50,hgmd_on=False):
        self.refgene = Refgene()
        self.mimdb = Omim()
        self.omim = self.mimdb.load_omim()
        self.clnsnpdb = Clnsnp()
        self.clnsnp = self.clnsnpdb.load_clnsnp()
        self.gerp = Gerp()
        self.genelist = self.load_refgene()
        self.effpred = Effpred()
        self.utrdb = UTRdb()
        self.utrdb.load_utrdb()
        self.mirnadb = MiRna()
        self.splicedb = Splicedb()
        self.regulomedb = Regulomedb()
        self.clinvardb = ClinvarDB()
        self.hgmd_genes = []
        self.hgmddb = None
        
        if hgmd_on:
          if os.path.exists(DBCONFIG['HGMDDB']['name']):
            self.hgmddb = HgmdDB() # added by CJ
            self.hgmd_on = True
          else:
            print '[Warning] user attempts to use HGMD but it seems not configured properly (a license requires!)'
            self.hgmd_on = False
        else:
          self.hgmd_on = False
        
        self.clinvtdb = ClinvitaeDB()
        
        self.NO_HGMD,self.HGMD_GENE,self.HGMD_POS,self.HGMD_ALT = range(4)
        self.kgdb = KGDB()
        self.dbsnp = DBSNP()
        self.espdb = EspDB()
        self.exacdb = ExACDB() # added by CJ
        self.nimbcap = NimblegenCapture(capture_kit_name,probe_flanking_bp)
        self.lcr = LCR()
        
        #setting for HGVS nomenclature
        self.hgvs =  Hgvs2()

    def _getpos(self, pos, ref, alt):
        """Internal Method. Since in VCF file the coordinate for the indels
        are mentioned with reference to the previous neucleotide, so before
        annoation this correction needs to be made. e.g. -
        For Insertion case -
        In VCF : 1       15903   .       G       GCAC
        is reframed as :  pos1=15903, pos2=15903, ref='-', alt='CAC'
        For Deletion case -
        In VCF : 1       17961   .       TGAGA      T
        is reframed as : pos1=17962, pos2=17965, ref='GAGA',alt='-'
        Args:
            - pos:    Genomic Position
            - ref:    Reference Allele as mentioned in VCF
            - alt:    Alternate Allele as mentioned in VCF
        Returns:
            - spos:    Start genomic position for the variant
            - epos:    End genomic position for the variant
            - c_ref:    Corrected Reference Allele
            - c_alt:    Corrected Alternate Allele
        """
        c_ref = None
        c_alt = None
        spos = None
        epos = None
        pos = int(pos)
        # Deletion
        if alt == '<DEL>':
            c_ref = ref
            c_alt = '-'
            spos = pos
            epos = spos + len(c_ref) - 1
        # Deletion
        elif len(ref) > len(alt):
            c_ref = ref[len(alt):]
            c_alt = '-'
            spos = pos + len(alt)
            epos = spos + len(c_ref) - 1
        # Insertion
        elif len(ref) < len(alt):
            c_ref = '-'
            c_alt = alt[len(ref):]
            spos = pos + len(ref) - 1
            epos = spos + 1
        # SNP
        else:
            c_ref = ref
            c_alt = alt
            spos = pos
            epos = spos + len(c_alt) - 1
        return spos, epos, c_ref, c_alt

    def load_refgene(self):
        genelist = [[] for c in CHROMOSOMES]
        for record in self.refgene.iterate_db():
            if record.chrom in ICHROMOSOMES:
                cid = ICHROMOSOMES[record.chrom]
            else:
                cid = len(CHROMOSOMES)
                ICHROMOSOMES[record.chrom] = cid
                CHROMOSOMES.append(record.chrom)
                genelist.append([])

            try:
                genelist[cid].append((int(record.str_transcript),
                                      int(record.end_transcript), record))
            except:
                print record.chrom, cid, CHROMOSOMES
                raise
        for key in genelist:
            key.sort()
        return genelist
    
    def retrieve_hgmdgenes(self):
        self.hgmd_genes = self.hgmddb.get_hgmd_genes()
      
    def retrieve_refgene(self, chrom, spos, epos):
        refgene_list = set([])
        chrom = to_ucsc_chrom(chrom)
        pos = int(spos)
        if chrom in ICHROMOSOMES:
            chrom = ICHROMOSOMES[chrom]
            sposidx = bisect(self.genelist[chrom], (spos, spos + 1))
            if sposidx != 0:
                for j in range(sposidx, -1, -1):
                    i1, i2, record = self.genelist[chrom][j - 1]
                    if pos >= i1 and pos <= i2:
                        refgene_list.add(record)
            if spos != epos:
                eposidx = bisect(self.genelist[chrom], (epos, epos + 1))
                if eposidx != 0:
                    for j in range(eposidx, -1, -1):
                        i1, i2, record = self.genelist[chrom][j - 1]
                        if pos >= i1 and pos <= i2:
                            refgene_list.add(record)
        refgene_list = list(refgene_list)
        return refgene_list

    def extract_updwn_genes(self, chrom, spos, epos, updwnrange=5000):
        upgene = 'NONE(dist=NONE)'
        dwngene = 'NONE(dist=NONE)'
        if not chrom.startswith('chr'):
                chrom = 'chr' + chrom
        if chrom in ICHROMOSOMES:
            chrom = ICHROMOSOMES[chrom]
            idx = bisect(self.genelist[chrom], (spos, spos + 1))
            if idx != 0:
                ui1, ui2, urecord = self.genelist[chrom][idx - 1]
                if idx == len(self.genelist[chrom]):
                    udist = spos - ui2
                    if udist <= updwnrange:
                        upgene = urecord.gene + '(dist=' + str(spos - ui2)\
                        + ')'
                else:
                    di1, di2, drecord = self.genelist[chrom][idx]
                    if ui2 < spos and epos < di1:
                        udist = spos - ui2
                        if udist <= updwnrange:
                            upgene = urecord.gene + '(dist=' + str(udist) + ')'
                        ddist = di1 - epos + 1
                        if ddist <= updwnrange:
                            dwngene = drecord.gene + '(dist=' + \
                                                str(di1 - epos + 1) + ')'

            elif spos < self.genelist[chrom][0][0]:
                di1, di2, drecord = self.genelist[chrom][0]
                ddist = di1 - spos + 1
                if ddist <= updwnrange:
                    dwngene = drecord.gene + '(dist=' + \
                                                str(di1 - spos + 1) + ')'

        updwngenes = ':'.join([upgene, dwngene])

        return updwngenes

    def get_revcmplmnt(self, allele):
        '''Returns reverse complemented allele'''
        if allele != '-':
            comp_allele = ''
            for neu in allele:
                comp_allele += COMPLEMENT[neu]
            comp_allele = comp_allele.upper()
            comp_allele = comp_allele[::-1]
        else:
            comp_allele = '-'
        return comp_allele

    def annotate(self, chrom, pos, ref, alt):
        '''Returns a list of annotation for the given Variant position and
        dbSNP Id'''

        annotations = []
        intron_entry = {}
        
        #back up the original input for HGVS representation input
        pos0 = pos
        ref0 = ref
        alt0 = alt
        
        self.chrom = chrom
        self.spos, self.epos, self.ref, self.alt = self._getpos(pos, ref, alt)

        # dbSNP
        dbsnp_version = self.dbsnp.get_version()[0].version[1:]
        dbsnp_ant = ['DB%s' % dbsnp_version, False]
        dbsnp_build = ['dbSNPBuildID', []]
        dbsnpids = []
        flag = False
        if self.dbsnp.has_snp(chrom, pos, ref):
            for rec in self.dbsnp.snp_by_location(chrom, pos):
                if rec.ref == ref and alt in rec.alt.split(','):
                    dbsnp_ant[1] = True
                    dbsnp_build[1] = str(rec.dbSNPBuildID)
                    dbsnpids = rec.id.split(';')[-1]
                    flag = True
                    break
        if flag is False:
            varlist = normalize_variant(chrom, pos, ref, alt)
            for var in varlist:
                if var[-1] in ['COMPLEX_DEL', 'COMPLEX_INS']:
                    dbsnp_build[1].append('.')
                    dbsnpids.append('.')
                    continue
                flag = False
                if self.dbsnp.has_snp(chrom, var[1], var[2]):
                    for rec in self.dbsnp.snp_by_location(chrom, var[1]):
                        if rec.ref == var[2] and var[3] in rec.alt.split(','):
                            dbsnp_ant[1] = True
                            dbsnp_build[1].append(str(rec.dbSNPBuildID))
                            dbsnpids.append(rec.id.split(';')[-1])
                            flag = True
                            break
                if flag is False:
                    dbsnp_build[1].append('.')
                    dbsnpids.append('.')
            if len([e for e in dbsnp_build[1] if e == '.']) == \
                                            len(dbsnp_build[1]):
                dbsnp_build[1] = '.'
            else:
                dbsnp_build[1] = ':'.join(dbsnp_build[1])
            if len([e for e in dbsnpids if e == '.']) == len(dbsnpids):
                dbsnpids = '.'
            else:
                dbsnpids = ':'.join(dbsnpids)

        # 1000Genome
        kgdb_ant = ['KGDB', False]
        kgaf = ['KGAF', []]
        flag = False
        if self.kgdb.has_snp(chrom, pos, ref):
            for rec in self.kgdb.snp_by_location(chrom, pos):
                if rec.ref == ref and alt in rec.alt.split(','):
                    kgdb_ant[1] = True
                    kgaf[1] = str(rec.AF)
                    flag = True
                    break
        if flag is False:
            varlist = normalize_variant(chrom, pos, ref, alt)
            for var in varlist:
                if var[-1] in ['COMPLEX_DEL', 'COMPLEX_INS']:
                    kgaf[1].append('0.0')
                    continue
                flag = False
                if self.kgdb.has_snp(chrom, var[1], var[2]):
                    for rec in self.kgdb.snp_by_location(chrom, var[1]):
                        if rec.ref == var[2] and var[3] in rec.alt.split(','):
                            kgdb_ant[1] = True
                            kgaf[1].append(str(rec.AF))
                            flag = True
                            break
                if flag is False:
                    kgaf[1].append('0.0')
            if len([e for e in kgaf[1] if e == '0.0']) == len(kgaf[1]):
                kgaf[1] = '0.0'
            else:
                kgaf[1] = ':'.join(kgaf[1])

        # Exome Sequencing Project (ESP)
        esp_ant = ['ESPDB', False]
        espaf = ['ESPAF', []]
        flag = False
        if self.espdb.has_snp(chrom, pos, ref):
            for rec in self.espdb.snp_by_location(chrom, pos):
                if rec.ref == ref and alt in rec.alt.split(','):
                    esp_ant[1] = True
                    espaf[1] = str(float(rec.MAF.split(',')[-1]) / 100)
                    flag = True
                    break
        if flag is False:
            varlist = normalize_variant(chrom, pos, ref, alt)
            for var in varlist: #0:chrom, 1:pos, 2:ref, 3:alt, 4.variant_type
                if var[-1] in ['COMPLEX_DEL', 'COMPLEX_INS']:
                    espaf[1].append('0.0')
                    continue
                flag = False
                if self.espdb.has_snp(chrom, var[1], var[2]):
                    for rec in self.espdb.snp_by_location(chrom, var[1]):
                        if rec.ref == var[2] and var[3] in rec.alt.split(','):
                            esp_ant[1] = True
                            espaf[1].append(str(float(
                                            rec.MAF.split(',')[-1]) / 100))
                            flag = True
                            break
                if flag is False:
                    espaf[1].append('0.0')
            if len([e for e in espaf[1] if e == '0.0']) == len(espaf[1]):
                espaf[1] = '0.0'
            else:
                espaf[1] = ':'.join(espaf[1])
        
        #ExAc data
        exac_ant = ['EXACDB', False]
        exacaf = ['EXACAF', []]
        flag = False
        if self.exacdb.has_snp(chrom,pos,ref):
          for rec in self.exacdb.snp_by_location(chrom,pos):
            alts = rec.alt.split(',')
            if rec.ref == ref and alt in alts:
              exac_ant[1] = True
              idx = alts.index(alt)
              denom = float(rec.AN_Adj)
              nom = float(rec.AC_Adj.split(',')[idx])
              if nom>0. and denom>0.:
                exacaf[1] = str(nom/denom)
              else:
                exacaf[1] = rec.AF.split(',')[idx]
              flag = True
              break
        if flag is False:
            varlist = normalize_variant(chrom, pos, ref, alt) #update primary key for variant locations
            for var in varlist:
                if var[-1] in ['COMPLEX_DEL', 'COMPLEX_INS']:
                    exacaf[1].append('0.0')
                    continue
                flag = False
                if self.exacdb.has_snp(chrom, var[1], var[2]):
                    for rec in self.exacdb.snp_by_location(chrom, var[1]):
                        alts = rec.alt.split(',')
                        if rec.ref == var[2] and var[3] in alts:
                            exac_ant[1] = True
                            idx=alts.index(var[3])
                            denom = float(rec.AN_Adj)
                            nom = float(rec.AC_Adj.split(',')[idx])
                            if nom>0. and denom>0.:
                              exacaf[1].append(str(nom/denom))
                            else:
                              exacaf[1].append(rec.AF.split(',')[idx])
                            flag = True
                            break
                if flag is False:
                    exacaf[1].append('0.0')
            if len([e for e in exacaf[1] if e == '0.0']) == len(exacaf[1]):
                exacaf[1] = '0.0'
            else:
                exacaf[1] = ':'.join(exacaf[1])

        # ClinvarDB
        cols = ['CLN' + ele for ele in ['ACC', 'DBN', 'DSDB', 'DSDBID', 'HGVS',
                                      'ORIGIN', 'SIG', 'SRC', 'SRCID']]
        clnacc, clndbn, clndsdb, clndsdbid, clnhgvs, clnorigin,\
        clnsig, clnsrc, clnsrcid = [[ele, []] for ele in cols]
        clnant = [clnacc, clndbn, clndsdb, clndsdbid, clnhgvs, clnorigin,
                  clnsig, clnsrc, clnsrcid]
        flag = False
        if self.clinvardb.has_snp(chrom, pos, ref):
            for rec in self.clinvardb.snp_by_location(chrom, pos):
                if '-1' in rec.CLNALLE.split(','):
                    continue
                clnalt = [rec.alt.split(',')[int(e) - 1] for e in
                            rec.CLNALLE.split(',')]
                if rec.ref == ref and alt in clnalt:
                    for cln in zip(*[e.split(',') for e in
                                        [','.join(clnalt), rec.CLNACC,
                                        rec.CLNDBN.replace(' ', '_'),
                                        rec.CLNDSDB.replace(' ', '_'),
                                        rec.CLNDSDBID, rec.CLNHGVS,
                                        rec.CLNORIGIN, rec.CLNSIG,
                                        rec.CLNSRC.replace(' ', '_'),
                                        rec.CLNSRCID]]):
                        if alt == cln[0]:
                            for f, v in zip(clnant, cln[1:]):
                                f[1] = v
                            flag = True
                            break
        if flag is False:
            varlist = normalize_variant(chrom, pos, ref, alt)
            for var in varlist:
                if var[-1] in ['COMPLEX_DEL', 'COMPLEX_INS']:
                    for f in clnant:
                        f[1].append('.')
                    continue
                flag = False
                if self.clinvardb.has_snp(chrom, var[1], var[2]):
                    for rec in self.clinvardb.snp_by_location(chrom, var[1]):
                        if '-1' in rec.CLNALLE.split(','):
                            continue
                        clnalt = [rec.alt.split(',')[int(e) - 1] for e in
                            rec.CLNALLE.split(',')]
                        if rec.ref == var[2] and var[3] in clnalt:
                            for cln in zip(*[e.split(',') for e in
                                        [','.join(clnalt), rec.CLNACC,
                                        rec.CLNDBN.replace(' ', '_'),
                                        rec.CLNDSDB.replace(' ', '_'),
                                        rec.CLNDSDBID, rec.CLNHGVS,
                                        rec.CLNORIGIN, rec.CLNSIG,
                                        rec.CLNSRC.replace(' ', '_'),
                                        rec.CLNSRCID]]):
                                if var[3] == cln[0]:
                                    for f, v in zip(clnant, cln[1:]):
                                        f[1].append(v)
                                    flag = True
                                    break
                if flag is False:
                    for f in clnant:
                        f[1].append('.')
            for f in clnant:
                if len([e for e in f[1] if e == '.']) == len(f[1]):
                    f[1] = '.'
                else:
                    f[1] = '__'.join(f[1])

        # Check for Nimblegen Capture Array
        nimb_capture = self.nimbcap.get_info(self.chrom, self.spos, self.epos)

        # Gerp Conservation
        gerp_consv = self.gerp.get_conserve(self.chrom, self.spos, self.epos)

        # LCR Annotation
        lcr_ant = self.lcr.get_info(self.chrom, self.spos, self.epos)

        # RegulomeDB
        regulome_score = []
        for rsid in dbsnpids.split(':'):
            score_info = self.regulomedb.retrieve(rsid)
            if score_info:
                key = score_info[0]
                regulome_score.append(score_info[1])
        if len(regulome_score) > 0:
            regulome_score = [key, ':'.join(regulome_score)]
        else:
            regulome_score = None

        # NHGRI-GWAS
        gwas_phen = []
        for rsid in dbsnpids.split(':'):
            if rsid in self.clnsnp['NHGRI']:
                gwas_phen.append('__'.join(self.clnsnp['NHGRI'][rsid]))

        if gwas_phen:
            gwas_phen = ':'.join(gwas_phen)

        # ------------------------------------
        # HGMD        
        hgmd_ant = ['HGMDDB', False]
        hgmd = ['HGMD_DESC', [], []] #Header, string_to_printout, hgmd_acc
        NotAvail = '.'
        
        if self.hgmd_on:  
          hgmd_alt_match_cnt = 0
          if self.hgmddb.has_snp(chrom, pos, ref):
            for rec in self.hgmddb.snp_by_location(chrom, pos):
              #to store hgmd hits per chrom_pos
              #SNP(chrom='1', pos=1959699, id='rs41307846', ref='G', alt='A', Gene='GABRD', HGMD_ACC='CM041768', PHENO=None, PMID='15115768,16452673', VC='DFP')
              match_type = self.HGMD_POS
              hgmd_alt_match, hgmd_desc = self._hgmd_has_alt(ref,alt,rec)
              if hgmd_alt_match:
                match_type = self.HGMD_ALT
                hgmd_alt_match_cnt += 1
              else:
                hgmd_desc = self._concat_hgmd_desc(rec)
              hgmd_desc = '%d(%s)'%(match_type, hgmd_desc)
              hgmd[1].append(hgmd_desc)
              hgmd[2].append(rec.HGMD_ACC)
  
          if hgmd_alt_match_cnt == 0:
            varlist = normalize_variant(chrom, pos, ref, alt)
            
            for var in varlist:#0:chrom, 1:pos, 2:ref, 3:alt, 4.variant_type
                flag = False
                if self.hgmddb.has_snp(chrom, var[1], var[2]):
                    for rec in self.hgmddb.snp_by_location(chrom, var[1]):
                        if not rec.HGMD_ACC in hgmd[2]: #was it printed prevly?
                          match_type = self.HGMD_POS
                          hgmd_alt_match, hgmd_desc = self._hgmd_has_alt(ref,alt,rec)
                          if hgmd_alt_match:
                            match_type = self.HGMD_ALT
                            hgmd_alt_match_cnt += 1
                          else:
                            hgmd_desc = self._concat_hgmd_desc(rec)
                          hgmd_desc = '%d(%s)'%(match_type, hgmd_desc)
                          hgmd[1].append(hgmd_desc)
                          hgmd[2].append(rec.HGMD_ACC)
            
          hgmd_pos_match = False
          if len(hgmd[2])>0:
            hgmd_ant[1] = True
            hgmd_pos_match = True
            hgmd[1] = ','.join(hgmd[1])

        clinvt_pathogenic = False
        clinvt = ['CLINVT',None]
        
        if self.clinvtdb.has_snp(chrom, pos, ref):
            for rec in self.clinvtdb.snp_by_location(chrom, pos):
                clinvt[1] = self._clinvt_has_alt(ref,alt,rec)
                if clinvt[1] and 'pathogenic' in clinvt[1]:
                    clinvt_pathogenic = True
                    break

        if not clinvt_pathogenic:
            varlist = normalize_variant(chrom, pos, ref, alt)
            
            for var in varlist:
                if clinvt_pathogenic: break
                if self.clinvtdb.has_snp(chrom, var[1], var[2]):
                    for rec in self.clinvtdb.snp_by_location(chrom, var[1]):
                        clinvt[1] = self._clinvt_has_alt(ref,alt,rec)
                        if clinvt[1] and 'pathogenic' in clinvt[1]:
                            clinvt_pathogenic = True
                            break

        # CADD raw and phred score
        cadd_raw, cadd_phred, cadd_aa = (None, None, None)
        if self.ref != '-' and self.alt != '-':
            rl = []
            pl = []
            ral = []
            for idx, vpos in enumerate(range(self.spos, self.epos + 1)):
                raw, phred, refaa, altaa = self.effpred.get_cadd(self.chrom, vpos,
                                    self.ref[idx], self.alt[idx])
                if raw is not None:
                    rl.append(str(raw))
                    pl.append(str(phred))
                    ral.append('%s/%s'%(refaa,altaa))

            if rl:
                cadd_raw = ','.join(rl)
                cadd_phred = ','.join(pl)
                cadd_aa = ','.join(ral)

        refgene_entries = self.retrieve_refgene(self.chrom, self.spos,
                                                self.epos)

        # Pseudoautosomal region annotation variable
        par = ['PAR', False]

        if refgene_entries:
            hgmd_gene_printed = {} #to avoid duplicated output
            for entry in refgene_entries:

                warning = ''

                # Allele changes based on strands
                if entry.strand == '+':
                    alt = self.alt
                    ref = self.ref
                elif entry.strand == '-':
                    alt = self.get_revcmplmnt(self.alt)
                    ref = self.get_revcmplmnt(self.ref)

                # Transcript and Gene
                refseq_acc = entry.refseq_acc
                gene = entry.gene
                error_state = int(entry.error)

                # Checks if gene is a pseudoautosomal gene
                if not par[1]:
                    if gene in PSEUDO_AUTO_GENES:
                        par[1] = True

                # Gene -Level Disease Association
                if gene in self.omim:
                    omim_phen, omim_ids = self.omim[gene]
                else:
                    omim_phen, omim_ids = [None] * 2
                if gene in self.clnsnp['GAD']:
                    site_gad = self.clnsnp['GAD'][gene]
                    site_gad = '__'.join(site_gad)
                else:
                    site_gad = None

                # Interpro Domain
                domains = self.effpred.get_domains(self.chrom, self.spos,
                                                   self.epos)

                # Region
                region_annot = RegionAnnotation(self.spos, self.epos,
                                                self.ref, self.alt, entry)
                region_annot.annotate()
                region = region_annot.region
                splice_site = region_annot.splice_site
                exons = region_annot.exons
                if region == 'Intergenic':
                    if entry.strand == '+':
                        gene = self.extract_updwn_genes(self.chrom,
                                        self.spos - 1, self.spos - 1)
                    else:
                        gene = self.extract_updwn_genes(self.chrom,
                                            self.spos, self.spos)
                
                if self.hgmd_on and hgmd_pos_match is False:
                    if gene in self.hgmd_genes:
                      if gene not in hgmd_gene_printed:
                        hgmd_ant[1] = True
                        hgmd[1] = '%s(%s)'%(self.HGMD_GENE,gene)
                        hgmd_gene_printed[gene] = True
              
                #if region_annot.cdna_poslist:
                if 'Intergenic' not in region:
                    #cdna_pos = '_'.join([str(e) for e in region_annot.cdna_poslist])
                    if alt0 == '<DEL>':
                      cdna_pos = None
                    else:
                      cdna_pos = self.hgvs.to_cDNA(self.chrom, pos0, ref0, alt0, refseq_acc)
                else:
                    cdna_pos = None

                protein_len = None
                # For Intronic Variants retaining only one transcript
                if region[-8:] == 'Intronic':
                    if '__' in splice_site:
                        sddist = int(splice_site.split('__')[0].split('_')[-1])
                        sadist = int(splice_site.split('__')[1].split('_')[-1])
                    else:
                        sddist, sadist = [0, 0]
                    if refseq_acc.startswith('NM'):
                        rrank = 1
                    else:
                        rrank = 2
                    annt = [gene, refseq_acc, region, exons, cdna_pos,
                            splice_site] + \
                           [None] * 9 + [omim_phen, omim_ids, site_gad,
                           gerp_consv, lcr_ant, regulome_score, gwas_phen,
                           [dbsnpids, dbsnp_ant, dbsnp_build,
                           kgdb_ant, kgaf, esp_ant, espaf, exac_ant, exacaf, hgmd_ant, hgmd, clinvt,
                           clnacc, clndbn, clndsdb, clndsdbid, clnhgvs,
                           clnorigin, clnsig, clnsrc, clnsrcid],
                           nimb_capture, par, cadd_raw, cadd_phred, cadd_aa, domains]
                    if gene not in intron_entry:
                        intron_entry[gene] = [annt, rrank, sddist, sadist]
                    else:
                        if rrank < intron_entry[gene][1]:
                            intron_entry[gene] = [annt, rrank, sddist, sadist]
                        elif sddist < intron_entry[gene][2] or \
                                sadist < intron_entry[gene][3]:
                            intron_entry[gene] = [annt, rrank, sddist, sadist]

                # For the refgene entries whose CDS is not multiple of 3
                elif error_state == 1:
                    annt = [gene, refseq_acc, region, exons, cdna_pos,
                    splice_site] + [None] * 8 + ['CDS_NOT_MULTIPLE_OF_3'] + \
                    [omim_phen, omim_ids, site_gad, gerp_consv,
                    lcr_ant, regulome_score, gwas_phen, [dbsnpids, dbsnp_ant,
                    dbsnp_build, kgdb_ant, kgaf, esp_ant, espaf, exac_ant, exacaf, hgmd_ant, hgmd, clinvt,
                    clnacc, clndbn, clndsdb, clndsdbid,
                    clnhgvs, clnorigin, clnsig, clnsrc,
                    clnsrcid], nimb_capture, par, cadd_raw, cadd_phred, cadd_aa,
                     domains]
                    annotations.append(annt)

                else:
                    # For Exonic and UTR variants
                    # Sequence and Mutation
                    seq_annot = SequenceAnnotation(region, entry,
                                    self.spos, self.epos, ref, alt, warning)
                    status, warning = seq_annot.annotate()
                    if status is True:
                        if seq_annot.cds:
                            protein_len = str((len(seq_annot.cds) / 3) - 1)
                        mut_annot = MutationAnnotation(region, seq_annot)
                        mut_annot.annotate()
                        mut_type = mut_annot.mut_type
                        codon_usage = mut_annot.codon_usage
                        if type(seq_annot.wt_codon) == list:
                            aa_change, codon_change = ([], [])
                            for waa, maa, wcdn, mcdn, aapos in zip(
                                         seq_annot.wt_aa, seq_annot.mut_aa,
                                         seq_annot.wt_codon,
                                         seq_annot.mut_codon,
                                         seq_annot.wt_aa_position):
                                codon_change.append(wcdn + '/' + mcdn)
                                aa_change.append(waa + str(aapos) + maa)
                                if mut_type in ['SynStop', 'StopLoss']:
                                    if int(aapos) <= int(protein_len):
                                        warning = 'NOT_ACTUAL_STOP_CODON__' \
                                        + 'TRANSCRIPT_WITH_MULTIPLE_STOP_CODON'
                            aa_change = '__'.join(aa_change)
                            codon_change = '__'.join(codon_change)
                        else:
                            if seq_annot.wt_codon:
                                codon_change = seq_annot.wt_codon +\
                                            '/' + seq_annot.mut_codon
                                aa_change = seq_annot.wt_aa + \
                                            str(seq_annot.wt_aa_position) \
                                            + seq_annot.mut_aa
                            else:
                                codon_change, aa_change = (None, None)
                    else:
                        mut_type, codon_usage, codon_change, aa_change,\
                                        wt_aa, mut_aa, aa_pos = [None] * 7

                    # UTRdb
                    utrannt = set([])
                    if 'UTR5' in region or 'UTR3' in region:
                        for pos in range(self.spos, self.epos + 1):
                            utra = self.utrdb.retrieve(self.chrom,
                                            self.spos, refseq_acc)
                            for e in list(utra):
                                utrannt.add(e)

                    # miRNA Binding Site Annotation
                    mirnas = []
                    if region == 'UTR3':
                        mirnas = self.mirnadb.get_annot(self.chrom,
                                                self.spos, refseq_acc)

                    # UTR Functional Site : Combine annotation from utrdb
                    # and mirna database
                    utr_funsites = None
                    func_sites = list(utrannt) + mirnas
                    if func_sites:
                        utr_funsites = '__'.join(func_sites)

                    # Splice Enhancer/Silencer Site Annotation
                    if region in ['CodingExonic', 'NonCodingExonic']:
                        splice_reg_site = self.splicedb.get_annot(self.chrom,
                                                    self.spos, refseq_acc)
                        if splice_reg_site:
                            splice_site = splice_reg_site

                    # Effect Prediction
                    sift_pred, pphen_pred = (None, None)
                    if type(seq_annot.wt_aa) == list:
                        sp = []
                        pp = []
                        for wt_aa, mut_aa, aa_pos in zip(seq_annot.wt_aa,
                                                seq_annot.mut_aa,
                                                seq_annot.wt_aa_position):
                            if wt_aa != mut_aa:
                                if self.spos == self.epos:
                                    self.effpred(self.chrom, self.spos, wt_aa,
                                                 mut_aa.replace('*', 'X'),
                                                 aa_pos, self.ref, self.alt)
                                else:
                                    self.effpred(None, None, wt_aa,
                                                 mut_aa.replace('*', 'X'),
                                                 aa_pos, None, None, gene)

                                # Sift
                                if self.effpred.sift_pred:
                                    sp.append(self.effpred.sift_pred + '_' +
                                            str(self.effpred.sift_score))

                                # Polypen2
                                if self.effpred.pp2_pred:
                                    pp.append(self.effpred.pp2_pred + '_' +
                                                str(self.effpred.pp2_score))
                        sift_pred = '__'.join(sp)
                        pphen_pred = '__'.join(pp)

                    annt = [gene, refseq_acc, region, exons, cdna_pos,
                            splice_site,
                            utr_funsites, mut_type, codon_change,
                            aa_change, protein_len, codon_usage, sift_pred,
                            pphen_pred, warning, omim_phen,
                            omim_ids, site_gad, gerp_consv, lcr_ant,
                            regulome_score, gwas_phen,
                            [dbsnpids, dbsnp_ant, dbsnp_build, kgdb_ant,
                             kgaf, esp_ant, espaf, exac_ant, exacaf, 
                             hgmd_ant, hgmd, clinvt, clnacc, clndbn, clndsdb,
                             clndsdbid, clnhgvs, clnorigin, clnsig, clnsrc,
                             clnsrcid], nimb_capture, par, cadd_raw,
                            cadd_phred, cadd_aa, domains]

                    annotations.append(annt)
            if intron_entry:
                for ele in intron_entry.values():
                    annotations.append(ele[0])
        else:
            updwngenes = self.extract_updwn_genes(chrom, self.spos, self.epos)
            annt = [updwngenes, None, 'Intergenic'] + 15 * [None]\
                                + [gerp_consv, lcr_ant, regulome_score,
                                   gwas_phen, [dbsnpids, dbsnp_ant,
                                   dbsnp_build, kgdb_ant, kgaf, esp_ant,
                                   espaf, exac_ant, exacaf, hgmd_ant, hgmd, clinvt,
                                   clnacc, clndbn, clndsdb, clndsdbid,
                                   clnhgvs, clnorigin, clnsig, clnsrc,
                                   clnsrcid], nimb_capture, par, cadd_raw, cadd_phred, cadd_aa,
                                   None]
            annotations.append(annt)
 
        return annotations
    
    def _concat_hgmd_desc(self,rec):
        NotAvail = '.'
        tmp = NotAvail
        if rec.HGMD_ACC: tmp = rec.HGMD_ACC
        hgmd_desc = '%s'%tmp
        
        tmp = NotAvail
        if rec.VC: tmp = rec.VC
        hgmd_desc += '|%s'%tmp
        
        tmp = NotAvail
        if rec.PMID: tmp = '_'.join(rec.PMID.split(','))
        hgmd_desc += '|%s'%tmp
        
        tmp = NotAvail
        if rec.PHENO: tmp = rec.PHENO.replace(',','.')
        hgmd_desc += '|%s'%tmp
        
        return hgmd_desc
              
    def _hgmd_has_alt(self,ref,alt,rec):
      """
      input1: ref (nt base in reference sequence)
      input2: alt (nt base in query sample)
      input3: rec (a row in the VCF record)
      input4: hgmd (hgmd_desc ['HGMD_DESC', hgmd annotation, matching type]
      
      output1: flag (whether both ref and alt are matched exactly?)
      output2: hgmd_desc ()
      """
      flag = False
      hgmd_desc = None
      
      alts = rec.alt.split(',')
      if rec.ref == ref and alt in alts:
        hgmd_desc = self._concat_hgmd_desc(rec)
        flag = True
      return flag,hgmd_desc
    
    def _concat_clinvt_desc(self,rec):
       vc_in_db = rec.VC.lower()
       cvt_vc = 'vus'
       if 'conflicting' in vc_in_db:
          cvt_vc = 'vus'
       elif 'benign' in vc_in_db:
          cvt_vc = 'benign'
       elif 'pathogenic' in vc_in_db:
          cvt_vc = 'pathogenic'
       return "%s(%s|%s)"%(cvt_vc,rec.cDNA,rec.URL)
        
    def _clinvt_has_alt(self,ref,alt,rec):

      clinvt_desc = None
      alts = rec.alt.split(',')
      if rec.ref == ref and alt in alts:
        clinvt_desc = self._concat_clinvt_desc(rec)
        
      return clinvt_desc
        

    def get_header(self):
        """ Intializing vcf formatted header information for populating VARANT
        annotation data"""

        Varant_h = namedtuple('Varant_h', 'id count type desc')

        header = [Varant_h('VARANT_INTERGENIC', 1, 'String', "Variant in "
                           "intergenic region and so reporting upstream and "
                           "downstream genes which are 5000bp away from "
                           "variant position and their exact (TSS/TES) "
                           "distance from variant position. Format - "
                           "VARANT_INTERGENIC"
                           "=UpstreamGene(dist=XYZ):DownstreamGene(dist=XYZ)"),
                  Varant_h('VARANT_GENIC', '.', 'String', "Variant in genic "
                      "region and so reporting the effects on genes. "
                      "Format - VARANT_GENIC="
                      "Gene1(Transcript1_id|Region|Exon_Number|AltId|cDNAPos|"
                      "SpliceSite|UTRSignal|Mutation|Codon_Change|"
                      "Amino_Acid_Change|Ref_Protein_Length|Codon_Usage|"
                      "SIFT(pred_score)|Polyphen2(pred_score)|Warning: "
                      "OMIM_Disease:OMIM_Ids:GAD_Disease). If there are "
                      "more than one transcripts for the gene, the effect "
                      "on them are appended after the annotation of first "
                      "transcript by \':\'"),
                  Varant_h('CADD_raw', '.', 'String', 'CADD raw score for '
                           'funtional prediction of a SNP. Please '
                           'refer to Kircher et al.(2014) Nature Genetics '
                           '46(3):310-5 for details. The larger the score '
                           'the more likely the SNP has damaging effect. '
                           'Scores are reported in the order in which ALT '
                           'alleles are reported.'),
                  Varant_h('CADD_phred', '.', 'String', 'CADD phred-like score.'
                           ' This is phred-like rank score based on whole '
                           'genome CADD raw scores. Please refer to Kircher '
                           'et al. (2014) Nature Genetics 46(3):310-5 for '
                           'details. The larger the score the more likely '
                           'the SNP has damaging effect. Scores are reported'
                           ' in the order in which ALT alleles are reported'),
                  Varant_h('CADD_aa', '.', 'String', 'Amino acid change in CADD prediction.'
                           ' This is change of amino acid based on whole '
                           'genome CADD raw scores. Please refer to Kircher '
                           'et al. (2014) Nature Genetics 46(3):310-5 for '
                           'details. The larger the score the more likely '
                           'the SNP has damaging effect. Scores are reported'
                           ' in the order in which ALT alleles are reported'),                  
                  Varant_h('Interpro_domains', '.', 'String', 'domain or '
                           'conserved site on which the variant locates. '
                           'Domain annotations come from Interpro database.')
                  ]
        return header


def add_meta(vcf, snpa):
    '''Add Annotation related meta data to VCF'''
    dbobjects = [snpa.clinvardb, snpa.clnsnpdb, snpa.dbsnp, snpa.effpred,
                 snpa.espdb, snpa.exacdb, snpa.gerp, snpa.lcr, snpa.kgdb, snpa.mirnadb,
                 snpa.nimbcap, snpa.mimdb, snpa.refgene,
                 snpa.regulomedb, snpa.utrdb, snpa.splicedb]
    if snpa.hgmd_on:
      dbobjects.append(snpa.hgmddb)
      
    dbversion = None
    for obj in dbobjects:
        for ele in obj.get_version():
            if dbversion:
                dbversion += ',' + ele.dbname + '(' + ele.version + ')'
            else:
                dbversion = ele.dbname + '(' + ele.version + ')'
    vcf.add_meta('VARANTDB_VERSIONS', '"' + dbversion + '"')
    # vcf.add_meta('VARANT_VERSION', '"v1.0b (Build 2013-DEC-10)"')
    for header in NimblegenCapture(snpa.nimbcap.kit_symbol,snpa.nimbcap.ext_bp).get_header():
        vcf.add_meta_info(header.id, header.count, header.type, header.desc)
    for header in Gerp().get_header():
        vcf.add_meta_info(header.id, header.count, header.type, header.desc)
    for header in snpa.lcr.get_header():
        vcf.add_meta_info(header.id, header.count, header.type, header.desc)
    for header in Clnsnp().get_header():
        vcf.add_meta_info(header.id, header.count, header.type, header.desc)
    for header in Regulomedb().get_header():
        vcf.add_meta_info(header.id, header.count, header.type, header.desc)
    for header in snpa.get_header():
        vcf.add_meta_info(header.id, header.count, header.type, header.desc)
    dbsnp_version = snpa.dbsnp.get_version()[0].version[1:]
    vcf.add_meta_info('PAR', 0, 'Flag', "Variant is in Pseudoautosomal" +
                      " Region")
    vcf.add_meta_info('DB%s' % dbsnp_version, 0, 'Flag', "Variant" +
                      " is present in dbSNP%s" % dbsnp_version)

    vcf.add_meta_info('dbSNPBuildID', 'A', 'String', "First dbSNP " +
                      "Build for RS ordered by ALT allele")

    vcf.add_meta_info('KGDB', 0, 'Flag', "Variant" +
                      " is present in 1000Genome Database")

    vcf.add_meta_info('KGAF', 'A', 'String', "Minor Allele Frequency" +
                      " reported in 1000Genome Database" +
                      " ordered by ALT allele")

    vcf.add_meta_info('ESPDB', 0, 'Flag', "Variant" +
                            " is present in Exome Sequencing Project" +
                            "(ESP) Database")
    vcf.add_meta_info('ESPAF', 'A', 'String', "The minor-allele frequency "
                      "computed based on 6503 exomes samples in Exome "
                      "Sequencing Project Database ordered by ALT allele")
    vcf.add_meta_info('EXACDB', 0, 'Flag', "Variant" +
                            " is present in Exome Aggregation Consortium" +
                            "(EXAC) Database")
    vcf.add_meta_info('EXACAF', 'A', 'String', "The minor-allele frequency "
                      "computed based on 60,706 unrelated individuals in Exome "
                      "Aggregation Consortium Database ordered by ALT allele")

    #add meta information regardless of the license
    vcf.add_meta_info('HGMDDB', 0, 'Flag', "Variant is present in HGMD database")
    vcf.add_meta_info('HGMD_DESC', '.', 'String',"[0:NO_HGMD,1:GENE_MATCH,2:POS_MATCH,3:POS_ALT_MATCH](ACC_ID_or_GeneSymbol|Variant_Class|PubMed|Phenotype)")
    
    vcf.add_meta_info('CLINVT','.','String',"class(cDNA|source_url)")

    clnvir = ClinvarDB()
    for field in ['CLNACC', 'CLNDBN', 'CLNDSDB',
                  'CLNDSDBID', 'CLNHGVS', 'CLNORIGIN',
                  'CLNSIG', 'CLNSRC', 'CLNSRCID']:
        meta_field = clnvir.get_meta_info(field)
        vcf.add_meta_info(meta_field.id, meta_field.count,
                          meta_field.type, meta_field.description)

def add_annotation(vcf, rec, annotations):
    '''Add annotation to the info field of the vcf record'''
    if len(annotations) > 0:
        for i in range(0, len(annotations), 2):
            vcf.add_info(rec, annotations[i], annotations[i + 1])


def main(vcffile, outfile, capture_kit_name, probe_flanking_bp, hgmd_on, logger):
    vcf = VCFParser(vcffile)
    outstream = open(outfile, 'w')

    snpa = SNPAnnotation(capture_kit_name, probe_flanking_bp, hgmd_on)
    logger.info('Instantiated SNP annotation class..')

    add_meta(vcf, snpa)
    vcf.writeheader(outstream)
    snpa.hgvs.load_resource()
    
    # print refgene_entries
    if hgmd_on:
        snpa.retrieve_hgmdgenes()
    
    snpcnt = 0
    timecnt = 0
    t1 = time.time()
    for rec in vcf:
        #print 'chrm[%s],pos[%d],ref[%s]'%(rec['chrom'],rec['pos'],rec['ref']) #debug
        vcf.parseinfo(rec)

        annotations = []
        chrom, pos, ref, alts = (rec['chrom'], rec['pos'], rec['ref'],
                                     rec['alt'])
        #if pos == '916549':
        #  catch_me = 'x'
        snpcnt += 1
        annot_info = {'Intergenic': None, 'Genic': {}}
        snpdb_ant = {}
        snpids = []
        regulome_score = []
        gwas_phens = []
        cadd_raw_list = []
        cadd_phred_list = []
        cadd_aa_list = []
        domainlist = []
        par_ant = None
        alt_cnt = 0
        for alt in alts:
            alt_cnt += 1
            for annt in snpa.annotate(chrom, pos, ref, alt):
                gene = annt.pop(0)
                region = annt[1]
                domains = annt.pop(27)
                if domains:
                    domainlist += domains
                cadd_aa = annt.pop(26)
                cadd_phred = annt.pop(25)
                cadd_raw = annt.pop(24)
                par = annt.pop(23)
                if par[1]:
                    par_ant = par
                nimb_capture = annt.pop(22)
                snpdbinfo = annt.pop(21)
                dbsnpids = snpdbinfo.pop(0)
                gwas = annt.pop(20)
                regu_score = annt.pop(19)
                lcr_ant = annt.pop(18)
                conserve = annt.pop(17)
                gad = annt.pop(16)
                omim_ids = annt.pop(15)
                omim_phen = annt.pop(14)
                data = ["" if ele is None else str(ele) for ele in annt]
                data.insert(3, str(alt_cnt))
                if region == 'Intergenic':
                    annot_info['Intergenic'] = gene
                else:
                    if gene in annot_info['Genic']:
                        ele = annot_info['Genic'][gene]
                        ele.append('|'.join(data))
                    else:
                        annot_info['Genic'][gene] = []
                        annot_info['Genic'][gene] = ["" if ele is None else
                                                     str(ele) for ele in
                                                     [omim_phen, omim_ids,
                                                      gad, '|'.join(data)]]
            snpids.append(dbsnpids)
            if cadd_raw is not None:
                cadd_raw_list.append(cadd_raw)
            else:
                cadd_raw_list.append('.')
            if cadd_phred is not None:
                cadd_phred_list.append(cadd_phred)
            else:
                cadd_phred_list.append('.')
            if cadd_aa is not None:
                cadd_aa_list.append(cadd_aa)
            else:
                cadd_aa_list.append('.')
            if not regu_score:
                regulome_score.append('.')
            else:
                regulome_score.append(regu_score[1])
            if not gwas:
                gwas_phens.append('.')
            else:
                gwas_phens.append(gwas)
            for ele in snpdbinfo:
                if ele[0] not in snpdb_ant:
                    snpdb_ant[ele[0]] = [ele[1]]
                else:
                    snpdb_ant[ele[0]].append(ele[1])

        if annot_info['Intergenic']:
            annotations += ['VARANT_INTERGENIC', annot_info['Intergenic']]
        if annot_info['Genic']:
            eff = []
            for gene, trans in annot_info['Genic'].iteritems():
                eff.append(gene + '(' + ':'.join(trans[3:] + trans[:3]) + ')')
            annotations += ['VARANT_GENIC', eff]
        if conserve:
            annotations += conserve
        if lcr_ant:
            annotations += lcr_ant
        if len([e for e in regulome_score if e == '.']) != len(regulome_score):
            annotations += ['RegulomeScore', regulome_score]
        if len([e for e in gwas_phens if e == '.']) != len(gwas_phens):
            annotations += ['GWASPhenotype', gwas_phens]
        if nimb_capture:
            annotations += nimb_capture
        if par_ant:
            annotations += par_ant
        dbsnp_version = snpa.dbsnp.get_version()[0].version[1:]
        if True in snpdb_ant['DB%s' % dbsnp_version]:
            annotations += ['DB%s' % dbsnp_version, True]
            annotations += ['dbSNPBuildID', snpdb_ant['dbSNPBuildID']]
            rec['id'] = snpids
        if True in snpdb_ant['KGDB']:
            annotations += ['KGDB', True]
            annotations += ['KGAF', snpdb_ant['KGAF']]
        if True in snpdb_ant['ESPDB']:
            annotations += ['ESPDB', True]
            annotations += ['ESPAF', snpdb_ant['ESPAF']]
        if True in snpdb_ant['EXACDB']:
            annotations += ['EXACDB', True]
            annotations += ['EXACAF', snpdb_ant['EXACAF']]

        if hgmd_on and True in snpdb_ant['HGMDDB']:
            annotations += ['HGMDDB', True]
            annotations += ['HGMD_DESC', snpdb_ant['HGMD_DESC']]
        
        if snpdb_ant['CLINVT'][0]:
            annotations += ['CLINVT',snpdb_ant['CLINVT']]
            
        if len([e for e in snpdb_ant['CLNACC'] if e == '.']) != \
                                            len(snpdb_ant['CLNACC']):
            cols = ['CLN' + ele for ele in ['ACC', 'DBN', 'DSDB', 'DSDBID',
                                            'HGVS', 'ORIGIN', 'SIG', 'SRC',
                                            'SRCID']]
            for key in cols:
                annotations += [key, snpdb_ant[key]]
        if len([e for e in cadd_raw_list if e == '.']) != len(cadd_raw_list):
            annotations += ['CADD_raw', cadd_raw_list]
        if len([e for e in cadd_phred_list if e == '.']) != len(cadd_phred_list):
            annotations += ['CADD_phred', cadd_phred_list]
        if len([e for e in cadd_aa_list if e == '.']) != len(cadd_aa_list):
            annotations += ['CADD_aa', cadd_aa_list]
        if domainlist:
            domainlist = list(set(domainlist))
            annotations += ['Interpro_domains', domainlist]

        # print annotations
        add_annotation(vcf, rec, annotations)
        vcf.write(outstream, rec, False)

        if time.time() - t1 > 60:
            timecnt += time.time() - t1
            d = 'Annotated - ' + str(snpcnt) + ' variants in ' + \
                  '%.2f' % ((timecnt) / 60) + 'min,' + \
                  ' Currently processing chrom-%s' % chrom.strip('chr')
            logger.info(d)
            t1 = time.time()
    
    snpa.hgvs.close_resource()
    vcf.stream.close()
    
    if time.time() - t1 < 60:
        timecnt += time.time() - t1
        d = 'Annotated - ' + str(snpcnt) + ' variants in ' + \
                  '%.2f' % ((timecnt) / 60) + 'min,' + \
                  ' Currently processing chrom-%s' % chrom.strip('chr')
        logger.info(d)
    if outstream != sys.stdout:
        outstream.close()


if __name__ == '__main__':
    # variant = ['15', 25425643, 'CGG', 'TGG']
    # variant = ['17', 2541615, 'G', 'A']
    # variant = ['1', 878230, 'G', 'C']
    # variant = ['1', 246703860, 'ATCT', 'A']
    # variant = ['7', 117307123, 'AGAG', 'A']
    # variant  = ['1',145103947,'TTTTATTTA','TTTTA,TTTTATTTATTTA,T']
    variant = ['1', 980824, 'G', 'A,C']
    variant = ['2', 152417133, 'GGC', 'GTT']
    info = {'Genic': {}, 'Intergenic': None}
    p = time.time()
    snpa = SNPAnnotation()
    snpa.hgvs.load_resource()
    
    q = time.time()
    print 'Time Taken to instantiate SNPAnnotation=', q - p, 'sec'

    dbobjects = [snpa.clinvardb, snpa.clnsnpdb, snpa.dbsnp, snpa.effpred,
                 snpa.espdb, snpa.exacdb, snpa.gerp, snpa.kgdb, snpa.mirnadb,
                 snpa.nimbcap, snpa.mirnadb, snpa.mimdb, snpa.refgene,
                 snpa.regulomedb, snpa.utrdb, snpa.splicedb]
    dbversion = None
    for obj in dbobjects:
        for ele in obj.get_version():
            if dbversion:
                dbversion += ', ' + ele.dbname.upper() + '(' \
                             + ele.version + ')'
            else:
                dbversion = ele.dbname.upper() + '(' + \
                            ele.version + ')'
    print dbversion

    alt_cnt = 1
    par_ant = False
    for annt in snpa.annotate(*variant):
        # print annt
        gene = annt.pop(0)
        region = annt[1]
        domains = annt.pop(27)
        cadd_aa = annt.pop(26)
        cadd_phred = annt.pop(25)
        cadd_raw = annt.pop(24)
        par = annt.pop(23)
        if par[1]:
            par_ant = par
        nimb_capture = annt.pop(22)
        kwnsnp_dbinfo = annt.pop(21)
        gwas = annt.pop(20)
        regulome_score = annt.pop(19)
        lcr_ant = annt.pop(18)
        conserve = annt.pop(17)
        gad = annt.pop(16)
        omim_ids = annt.pop(15)
        omim_phen = annt.pop(14)
        data = ["" if ele is None else str(ele) for ele in annt]
        data.insert(3, str(alt_cnt))

        if region == 'Intergenic':
            info['Intergenic'] = gene
        else:
            if gene in info['Genic']:
                ele = info['Genic'][gene]
                ele.append('|'.join(data))
            else:
                info['Genic'][gene] = [], True
                info['Genic'][gene] = ["" if ele is None else str(ele) for
                                       ele in [omim_phen, omim_ids, gad,
                                               '|'.join(data)]]

    for key, val in info.iteritems():
        if key == 'Genic':
            temp = []
            for gene, trans in info['Genic'].iteritems():
                temp.append(gene + '(' + ':'.join(trans[3:] + trans[:3]) + ')')
            info['Genic'] = ','.join(temp)
            
    snpa.hgvs.close_resource()
    
    print info
    print cadd_raw
    print cadd_phred
    print cadd_aa
    print domains
    print kwnsnp_dbinfo
    print nimb_capture
    print gwas
    print regulome_score
    print lcr_ant
    print conserve
    print par_ant
