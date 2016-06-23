#!/usr/bin/env python
'''
COPYRIGHT (C) 2016 changjin.hong@gmail.com
author: changjin.hong@gmail.com
'''
import os, sys, argparse, re, math, datetime, time, logging
from ConfigParser import SafeConfigParser
import subprocess as sp
import cPickle as pickle
import damaging_model, pagerank
import gcn.lib.io.vcf as vcf
from gcn.lib.io import anyopen
from gcn.lib.utils import lib_utils
from gcn.lib.varann.vartype.varant import varant_parser as vp
from gcn.lib.varann.vartype.varant import annotateRegion, vcf_mask
from gcn.lib.io.vcfutils import normalize_variant
from gcn.lib.databases import geneontology
 
VERSION = '0.1.1'
author_email = 'changjin.hong@gmail.com'

class PerturbedGene:
	def __init__(self):
		self.hpos = []
		self.diseases = []
		self.genes = {}
		self.pdmg = 0.
		self.gdmg = 0.
		self.score = 0.

class Divine:
	'''
	collect program configuration, input parameters, and computational resource
	'''
	def __init__(self, uargs):
		#transferring user input arguments to class member variables

		self.exp_tag = uargs.exp_tag
		self.vknown = uargs.vknown
		
		if uargs.no_cadd: self.cadd = False
		else: self.cadd = True
		
		self.excl_non_coding = False
		self.sparser = SafeConfigParser()
		
		self.pheno_dmg = {}
		self.genetic_dmg = {}
		self.gene_dmg = {}

		self.hpo2disease_fn = None
		self.pheno_dmg_fn = None
		self.hpo_query = None
		self.vcf = None
		self.xls = None
		
		lib_utils.msgout('notice','initializing Divine ...','Divine')
		
		divine_root_dir = os.environ.get("DIVINE")
		if not divine_root_dir:
			raise EnvironmentError("set DIVINE variable properly!")
		
		config_fn = os.path.join(divine_root_dir,'gcn','bin','prioritize','configs','divine.conf')
		
		if not lib_utils.check_if_file_valid(config_fn):
			raise IOError("check if the configuration file[%s] is valid!" % config_fn)
		
		self.config_fn = config_fn
		self.entries = {'divine_root':divine_root_dir}
		self._set_args(uargs)
		
		# damage factor w.r.t the location of variant within the transcript
		self.dm = damaging_model.DmgCoeff(uargs.indel_mode,uargs.top_seed_rate,self.logger)
		
		if uargs.wes_mask:
			msg = 'VCF will be masked by RefGene coding region'
			lib_utils.msgout('notice',msg);self.logger.info(msg)

		self.wes_mask = uargs.wes_mask

		lib_utils.msgout('notice','done. initialization')
	
	def _assign_out_fn(self,fbase,fext='tsv'):
		if self.exp_tag:
			fn = os.path.join(self.out_dir,'%s_%s.%s'%(fbase,self.exp_tag,fext))
		else:
			fn = os.path.join(self.out_dir,'%s.%s'%(fbase,fext))
		return fn
	
	def _set_args(self, uargs):
		'''
		-objective: checking input parameters, reading config, and storing user command line
		-input: uargs (args from main())
		-output: class initialization 
		'''
		job_name = '_set_args'
		lib_utils.msgout('notice','storing input condition ...',job_name)
		
		if not uargs.hpo_query_fn and not uargs.cf_fn:
			raise RuntimeError('either VCF (-v) or query phenotype (-q) file should be provided!')

		# check sanity of the input files
		if uargs.hpo_query_fn:
			if lib_utils.check_if_file_valid(uargs.hpo_query_fn):
				self.hpo_query = uargs.hpo_query_fn
			else:
				raise IOError('check if [%s] is valid' % uargs.hpo_query_fn)
		
		if uargs.vcf:
			if lib_utils.check_if_file_valid(uargs.vcf):
				self.vcf = uargs.vcf
			else:
				raise IOError('check if [%s] is valid' % uargs.vcf)

		if uargs.capkit in ['SureSelect_V6', 'SeqCapEZ_Exome']:
			self.capkit = uargs.capkit
		else:
			raise RuntimeError("revise capture kit symbol[%s]" % uargs.capkit)
		
		# check input condition
		if uargs.out_dir is None:
			if self.vcf:
				uargs.out_dir = os.path.join(os.path.dirname(self.vcf), 'divine')
			else:
				uargs.out_dir = os.path.join(os.path.dirname(self.hpo_query), 'divine')
		
		#create the output directory user specifies
		if uargs.out_dir.endswith('/'):
			uargs.out_dir = uargs.out_dir[:-1]
		self.out_dir = uargs.out_dir
		lib_utils.ensure_dir(self.out_dir)
		
		#prepare output file name
		self.rank_fn = self._assign_out_fn('rank','tsv')

		self.log_dir = os.path.join(self.out_dir, 'logs')
		lib_utils.ensure_dir(self.log_dir)
		
		msg = 'prepared log directory[%s]  ...'%self.log_dir
		lib_utils.msgout('notice',msg,job_name)
				
		#prepare loggig handler
		ts = datetime.datetime.fromtimestamp(time.time()).strftime('%Y%m%d_%H%M%S')

		FORMAT = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
		logging.basicConfig(filename=os.path.join(self.log_dir, 'divine_%s.log' % ts),\
											 filemode="w", level=logging.DEBUG, format=FORMAT)
		
		# ------------------------
		self.logger = logging.getLogger('divine')
		# ------------------------
		self.logger.info(msg)
		
		#read configuration file containing 3rd parties s/w path and database locations
		self._read_config()

		#record user command line
		self.record_commandline()
		msg = 'Divine initialization completed [%s]'%job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)

	def _set_config(self, section, entry):
		'''
		objective: read item in the configuration file
		input:
			-section: part/section in the configuration file
			-entry: item to search in the section
		output:
			-entries: member var updated
		'''
		try:
			self.entries[entry] = self.sparser.get(section, entry)
			if section in ['program_paths', 'database', 'config']:
				if not self.entries[entry].startswith('/'):
					if self.entries['divine_root'] is None:
						raise ValueError('[divine_root] should be defined first!')
					self.entries[entry] = '%s/%s'%(self.entries['divine_root'],self.entries[entry])
		except:
			raise IOError('check if [%s] exists in %s[%s]' % (entry, self.config_fn, section))

	def _read_config(self):
		'''
		objective: read configuration file
		'''
		job_name = '_read_config'
		msg = 'reading configuration file [%s;%s] ...'%(job_name,self.config_fn)
		lib_utils.msgout('notice',msg);self.logger.info(msg)

		self.sparser.read(self.config_fn)

		self._set_config('program_paths', 'temp_dir')
		self._set_config('program_paths', 'varant')
		self._set_config('program_paths', 'hposim')
		self._set_config('program_paths', 'vcf2xls')
		
		self._set_config('config', 'vcf_filter_conf')
		
		self._set_config('database', 'ext_disease_to_gene')
		
		self._set_config('database', 'hpo_obo')
		self._set_config('database', 'beta_fit')
		self._set_config('database', 'string_link')
		
		'''
		to access to UCSC mysql database(hg19)
		select e2g.value, gtp.protein from ensGtp as gtp
		inner join ensemblToGeneName as e2g on e2g.name=gtp.transcript;
		'''
		self._set_config('database', 'esp_to_gene')

		# check if the file or directory all exists before long journey!
		for key, path2 in self.entries.iteritems():
			if not lib_utils.check_if_file_valid(path2):
				raise IOError('check [%s = %s] in the file [%s]' %\
										(key, path2, self.config_fn))

		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		return self.entries
	
	def record_commandline(self):
		'''
		objective: record the divine run condition into logger
		'''
		import socket
		job_name = 'record_commandline'
		msg='capturing user command line [%s] ...'%job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		try:
			host_name = socket.gethostname()
		except:
			host_name = 'N/A'
		self.logger.info('host:%s'%host_name)
		
		try:
			user = os.environ.get('USER')
		except:
			user = 'N/A'
		self.logger.info('user:%s'%user)
		
		try:
			pwd = os.environ.get('PWD')
		except:
			pwd = 'N/A'
		self.logger.info('pwd:%s'%pwd)

		self.logger.info('cmd:%s'%(' '.join(sys.argv)))
		self.logger.info("divine configuration file:%s" % self.config_fn)
		
		self.logger.info('exclude_non_coding:%s'%self.excl_non_coding)
		
		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
	def hpo_to_diseases(self):
		'''
		objective: match HPO IDs from a given patint phenotype to known disease database
		input: hpo_query, hpo database
		method: system call
		output: phenotype matching score w.r.t disease
		'''
		
		job_name = 'hpo_to_diseases'
		# prepare output file
		self.hpo2disease_fn = self._assign_out_fn(job_name,'tsv')
		
		msg = 'matching query phenotypes to diseases in semantic HPO ontology[%s;%s]'%(job_name,self.hpo2disease_fn)
		lib_utils.msgout('notice',msg);self.logger.info(msg)

		# run hpo similarity
		cmd = ["python", self.entries['hposim'], \
					"-q", self.hpo_query, \
					"-b", self.entries['hpo_obo'], \
					"-f", self.entries['ext_disease_to_gene'], \
					"-S", self.dm.hpo_max_queries, \
					"--normalize", \
					"-o", self.hpo2disease_fn]

		self.run_cmd(cmd, job_name)
		
		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)
	
	def disease_to_genes_sum(self,gene_norm=True):
		'''
		objective: from hit scores of query hpo to disease, associate disease to genes
		input:
		  -hpo2disease_fn: a file generated by hpo_to_disease()
		    '#query(file_name)\tomim\tgenes\tscore\n'
		  -gene_norm: want to normalize accumulated phenotype score per gene? [True]
		'''

		job_name = 'disease_to_genes'
		if self.hpo2disease_fn is None:
			self.hpo_to_diseases()

		msg = 'aggregating HPO hit scores of disease to each gene [%s;%s]...' % \
							(job_name,self.hpo2disease_fn)
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		fp = lib_utils.open2(self.hpo2disease_fn, 'r')

		#accumulating phenotye-matching score into genes associated with the disease
		pheno_genes = {}
		pheno_genes_cnt = {}
		for i in fp:  # for each disease
			if i.startswith('#'): continue
			i = i.rstrip()
			_,omim,geneStr,funsimMatAvg = i.rstrip().split('\t')
			genes = geneStr.split(',')
			funsimMatAvg = float(funsimMatAvg)
			for gene in genes:  # for each gene
				if funsimMatAvg > 0.:
					if gene not in pheno_genes:
						pheno_genes[gene] = 0.
						pheno_genes_cnt[gene] = 0.
					pheno_genes[gene] += funsimMatAvg
					pheno_genes_cnt[gene] += 1.
		fp.close()
		
		if gene_norm:
			msg = 'normalizing a bipartite graph between diseases and genes...'
			lib_utils.msgout('notice',msg,job_name);self.logger.info(msg)
			for gene in pheno_genes.keys():
				pheno_genes[gene] /= pheno_genes_cnt[gene]

		self.pheno_dmg = lib_utils.normalize_dic(pheno_genes, 'sum')
		
		#print phenotypic damage scores
		self.rank_pheno_gene()

		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)

		#clean up variables
		pheno_genes = None
		pheno_genes_cnt = None

	def disease_to_genes(self):
		'''
		objective: from hit scores of query hpo to disease, associate disease to genes
		input:
		  -hpo2disease_fn: a file generated by hpo_to_disease()
		    '#omim\tgenes\tscore\n'
		  -gene_norm: want to normalize accumulated phenotype score per gene? [True]
		'''

		job_name = 'disease_to_genes'
		if self.hpo2disease_fn is None:
			self.hpo_to_diseases()

		msg = 'aggregating HPO hit scores of disease to each gene [%s;%s]...' % \
							(job_name,self.hpo2disease_fn)
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		fp = lib_utils.open2(self.hpo2disease_fn, 'r')

		#accumulating phenotye-matching score into genes associated with the disease
		pheno_genes = {}
		
		for i in fp:  # for each disease
			if i.startswith('#'): continue
			i = i.rstrip()
			omim,geneStr,funsimMatAvg = i.rstrip().split('\t')
			genes = geneStr.split(',')
			funsimMatAvg = float(funsimMatAvg)
			for gene in genes:  # for each gene
				if funsimMatAvg > 0.:
					if gene not in pheno_genes:
						pheno_genes[gene] = 0.
					
					if funsimMatAvg>pheno_genes[gene]: #keep only maximum
						pheno_genes[gene] = funsimMatAvg

		fp.close()

		self.pheno_dmg = lib_utils.normalize_dic(pheno_genes, 'sum')
		
		#print phenotypic damage scores
		self.rank_pheno_gene()

		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)

		#clean up variables
		pheno_genes = None
		
	def vannotate(self,reuse=False):
		'''
		objective: run varant (GCN) annotator
		input: self.vcf
		output: annotated vcf
		'''
		job_name = 'vannotate'
		msg = 'annotating VCF file[%s;%s] ...'%(job_name,self.vcf)
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		# prepare output file
		varant_vcf = os.path.join(self.out_dir,'divine.vcf')
		
		# if necessary, masking the raw vcf file
		coding_vcf = None
		if self.wes_mask:
			if not lib_utils.check_if_file_valid(varant_vcf) or not reuse:
				cRef = annotateRegion.RefGeneUcscTB(work_dir=self.out_dir,logger=self.logger)
				coding_bed_fn = cRef.create_bed(ext_bp=20,reuse=False)
				
				msg = 'extracting variants in coding region from [%s] @ %s ...'%(self.vcf,job_name)
				lib_utils.msgout('notice',msg);self.logger.info(msg)
				
				coding_vcf = os.path.join(self.out_dir,'refgene_e20.vcf')
				self.vcf = vcf_mask.by_bed(self.vcf,coding_bed_fn,coding_vcf,logger=self.logger)
				
				msg = 'done. @ %s'%job_name
				lib_utils.msgout('notice',msg);self.logger.info(msg)

		if not lib_utils.check_if_file_valid(varant_vcf) or not reuse:
			self.logger.info('annotating [%s,%s] ...' % (job_name, self.vcf))
			
			cmd = ["python", self.entries['varant'], \
						"-i", self.vcf, \
						"-o", varant_vcf, \
						"-l", self.log_dir]
			if self.capkit:
				cmd.extend(["-c", self.capkit, "-e", "180"])

			self.run_cmd(cmd,job_name)
		self.vcf = varant_vcf
		
		if coding_vcf:
			os.unlink(coding_vcf)
		
		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)
	
	def vfilter(self):
		'''
		objective:apply a standard filter to VCF file and classify variants
		input: annotated vcf from varant (GCN) annotator
		output: filtered vcf
		'''
		job_name = 'vfilter'
		msg = 'filtering the annotated VCF [%s;%s] ...'%(job_name,self.vcf)
		lib_utils.msgout('notice',msg);self.logger.info(msg)

		filtered_vcf = self._assign_out_fn(job_name,'vcf')

		msg='applying a standard filter/class tagging [%s]' % self.vcf
		lib_utils.msgout('notice',msg,job_name);self.logger.info(msg)
		
		gcn_filter = os.path.join(self.entries['divine_root'], 'gcn', 'lib', 'utils', 'filter_cj.py')
		
		cmd = ["python", gcn_filter, \
					"-i", self.vcf, \
					"-o", filtered_vcf]
		
		filter_conf = self.entries['vcf_filter_conf']
		cmd.extend(["-f", filter_conf])
		
		self.logger.info('filter config [%s] is applied' % filter_conf)
		self.run_cmd(cmd, job_name)
		self.vcf = filtered_vcf
		
		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)
	
	def _store_variants(self,beta_fits):
		'''
		collect essential info on each variant
		'''
		job_name = '_store_variants'
		msg='collecting variant information and class label to determine genetic damage [%s;%s]...'%(job_name,self.vcf)
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		mutation_info = []
		v = vcf.VCFParser(self.vcf)
		for rec in v:
			v.parseinfo(rec)
			
			varlist = normalize_variant(rec.chrom, rec.pos, rec.ref, rec.alt)
			mut_type = varlist[0][-1]

			if ':' in rec.id[0]:
				mut_type = 'mnp'
			
			# collect conservation prediction score (CADD and GERP++) 
			cadd_aa = './.'
			px_cadd = None
			if rec.info.CADD_raw:
				# to get CADD_raw (average)
				px_cadd, cadd_aa = vcf.get_CADD_scores(mut_type, rec.info.CADD_aa, rec.info.CADD_raw, beta_fits)
			
			# to get GERP++ score
			px_gerp = None
			if rec.info.GerpConserve:
				px_gerp = vcf.get_GERP_scores(mut_type, cadd_aa, rec.info.GerpRSScore, beta_fits)
			
			# which score can be chosen
			px = 0.5
			if self.cadd and px_cadd is not None:
				px = px_cadd
			elif px_gerp is not None:
				px = px_gerp

			vpop = vp.parse(rec.info)
			genes = []
			
			# to get MAF in the order of ExAC, ESP, and 1K
			if rec.info.EXACDB:
				maf = float(rec.info.EXACAF[0])
			elif rec.info.ESPDB:
				maf = float(rec.info.ESPAF[0])
			elif rec.info.KGDB:
				maf = float(rec.info.KGAF[0])
			else:
				maf = 0.
			
			# to compute a significance of MAF
			maf_offset = 0.
			if maf > 0:
				maf_offset = (1. - self.dm.beta1 * math.exp(1000.*maf)) / self.dm.beta2
				if maf_offset < 0.:
					maf_offset = 0.
			
			# to get transcript length
			for altnum, val in vpop.items():
				# for each gene involved with the variant
				for gene, gd in val.items():
					protein_len = self.dm.avg_protein_len
					if gd:
						for t in gd['TRANSCRIPTS']:
							if t.protein_len:
								protein_len = float(t.protein_len)
								break

					# store a set of essential annotation to be used for genetic damage
					if gene not in genes:
						mutation_info.append([gene, rec.info.INDEL, rec.info.CLASS_TAG, protein_len, px, maf_offset])
						genes.append(gene)

		# done reading filterd VCF file
		v.stream.close()
		msg = 'done. [%s]'%job_name
		lib_utils.msgout('notice',msg); self.logger.info(msg)

		return mutation_info
	
	def _predict_dmg_score(self,mutation_info):
		# 
		job_name = '_predict_dmg_score'
		msg = 'combining impact by variant location in tx and conservation pred score [%s] ...'%job_name
		lib_utils.msgout('notice',msg); self.logger.info(msg)
		
		#[max for non-known pathogenic variants, min for known pathogenic variants]
		Mm_score = [-1., 1.]
		
		maf_offset_max = 0.
		printed_gene_cli = {}
		genetic_dmg2 = [{}, {}] #[unpublished, clinvar(HGMD)]
		
		# estimate genetic damage scores per gene
		for gene, indel, tag, protein_len, px, maf_offset in mutation_info:
			vreg_dmg = 0.
			# filter out frequent or known-benign
			if '1' in tag or '2' in tag: continue
			
			# to locate variant location
			if 'n' in tag:
				if self.excl_non_coding: continue
				else:
					vreg_dmg += self.dm.vncoding
			elif 'i' in tag:
				vreg_dmg += self.dm.vintronic
			elif 'e' in tag:
				vreg_dmg += self.dm.vexonic
			elif 'S' in tag:
				vreg_dmg += self.dm.vsplice
			elif 's' in tag:
				vreg_dmg += self.dm.vsplice_syn
			elif 'w' in tag:
				vreg_dmg += self.dm.warning

			j = 0
			if self.vknown:
				if 'c' in tag:  # previously known pathogenic (from ClinVar or HGMD)?
					j = 1
					if gene in printed_gene_cli: continue
					vreg_dmg += self.dm.cexonic   
					printed_gene_cli[gene]=True
					
				elif 'I' in tag: # was it in intronic/intergenic/non-coding?
					j = 1
					vreg_dmg += self.dm.cintronic
				elif 'g' in tag: # pathogenic gene?
					vreg_dmg += self.dm.cgene

			if not indel:
				if '3' in tag: # rare
					vreg_dmg += maf_offset
					if maf_offset>maf_offset_max:
						maf_offset_max = maf_offset
				elif '4' in tag: # de-novo
					vreg_dmg += maf_offset_max
					
			if '3' in tag:  # rare
				if 'D' in tag: # is it homozygous/compound het/known AD disease?
					dmg_allele_wt = self.dm.hom_rare
				else:
					dmg_allele_wt = self.dm.het_rare
			elif '4' in tag:  # de-novo
				if indel:
					if 'D' in tag:
						dmg_allele_wt = self.dm.hom_denovo_indel  
					else:
						dmg_allele_wt = self.dm.het_denovo_indel
				else:
					if 'D' in tag:
						dmg_allele_wt = self.dm.hom_denovo_snp
					else:
						dmg_allele_wt = self.dm.het_denovo_snp

			if j == 1:
				if 'D' in tag:
					dmg_allele_wt = self.dm.hom_denovo_snp
				else:
					dmg_allele_wt = self.dm.het_denovo_snp

			if gene not in genetic_dmg2[j]:
				genetic_dmg2[j][gene] = 0.
			
			#combine dmg predicted by variant within tx and conservation pred dmg
			genetic_dmg2[j][gene] += ((1. - self.dm.prwt) * vreg_dmg \
															+ self.dm.prwt * px) * dmg_allele_wt / protein_len
			
			if j == 0:
				if genetic_dmg2[j][gene] > Mm_score[j]:
					Mm_score[j] = genetic_dmg2[j][gene]
			else:  # prev known deleterious
				if genetic_dmg2[j][gene] < Mm_score[j]:
					Mm_score[j] = genetic_dmg2[j][gene]

		msg = 'done. [%s]'%job_name
		lib_utils.msgout('notice',msg); self.logger.info(msg)
		
		return genetic_dmg2, Mm_score
	
	def variants_to_genedmg(self):
		'''
		-objective:depending on the location of the variant, type, predicted effect function in similar species conservation, we compute the genetic damage score  
		-input: filtered vcf from vfilter()
		-output: dictionary {gene:genetic damaged score}
		'''
		
		job_name = 'variants_to_genedmg'
		msg='start to predict genetic damage score from variants in the provided VCF [%s]' % (job_name)
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		# 3:novel, 4: rare, 5:known_pathogenic variant, e:exonic, i:intronic/splice_effect, g: known pathogenic gene, D
		try:
			beta_fit_pyv = self.entries['beta_fit']
			msg='loading beta fit cdf[%s] for conservation score w.r.t. AA'%beta_fit_pyv
			lib_utils.msgout('notice',msg); self.logger.info(msg)
			fp = open(beta_fit_pyv, 'rb')
			beta_fits = pickle.load(fp)
			fp.close()
		except:
			beta_fits = [None, None, None]
		
		# to store essential info on each variants to a list
		mutation_info = self._store_variants(beta_fits)
		
		# combine variant location and conservation pred dmg
		genetic_dmg2, Mm_score = self._predict_dmg_score(mutation_info)
		
		if self.vknown:
			msg = 'adjust clipatho gene score s.t. the min score > non-clipatho score.'
			lib_utils.msgout('notice',msg); self.logger.info(msg)
			gamma2 = 1.1
			delta2 = Mm_score[0] * gamma2 - Mm_score[1]
			for gene in genetic_dmg2[1].keys():
				genetic_dmg2[1][gene] += delta2
			msg = 'done.'
			lib_utils.msgout('notice',msg); self.logger.info(msg)

		# normalize damage prediction score
		msg = 'normalizing genetic damage scores ...'
		lib_utils.msgout('notice',msg); self.logger.info(msg)
		
		genes_to_infer2 = {}
		for j in range(2):
			for gene, score in genetic_dmg2[j].iteritems():
				if gene not in genes_to_infer2:
					genes_to_infer2[gene] = 0.
				genes_to_infer2[gene] += score

		self.genetic_dmg = lib_utils.normalize_dic(genes_to_infer2, 'sum')
		
		msg = 'done. [%s]' % job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		#clean up
		genes_to_infer2 = None
		genetic_dmg2 = None
		printed_gene_cli = None
		mutation_info = None
		
		return self.genetic_dmg
			
	def run_cmd(self, cmd, job_name=None):

		cmd_str = lib_utils.joined(cmd,' ')
		lib_utils.msgout('notice',cmd_str) #debug
		self.logger.info('running [%s] ...' % cmd_str)
		
		if job_name:
			stdofp, stdefp = self.get_process_msg_handler(job_name)
		else:
			stdofp = sp.PIPE
			stdefp = sp.PIPE
		proc = sp.Popen(cmd_str, stdout=stdofp, stderr=stdefp, shell=True)
		retcode = proc.wait()
		
		if job_name:
			stdofp.close()
			stdefp.close()
			
		if retcode > 0:
			self.logger.error('[%s] failed' % cmd_str)
			raise RuntimeError('[%s] failed' % cmd_str)
		self.logger.info('done. [%s]' % job_name)

	def get_process_msg_handler(self, step_name):
		
		stdout_fn = os.path.join(self.log_dir, '%s.out'%step_name)
		stdofp = open(stdout_fn, 'w')
		stderr_fn = os.path.join(self.log_dir, '%s.err'%step_name)
		stdefp = open(stderr_fn, 'w')
		return stdofp, stdefp
	
	def get_GO_seeds(self,top_seed_rate):
		'''
		to collect genes associated a disease whose matching score to HPO is relatively high
		'''
		job_name = 'get_GO_seeds'
		msg = 'collecting genes associated with diseases [%s] showing high HPO matching'%self.hpo2disease_fn
		lib_utils.msgout('notice',msg); self.logger.info(msg)
		
		#count the total number of disease hit whose score > 0.
		fp = anyopen.openfile(self.hpo2disease_fn)
		num_omim = 0
		for i in fp:
			if i[0] == '#': continue
			omim,genes,score = i.rstrip().split('\t')
			score = float(score)
			if score>0.:
				num_omim += 1
		fp.close()
		
		t = 0
		T = round(num_omim*top_seed_rate)
		fp = anyopen.openfile(self.hpo2disease_fn)
		go_seeds = []
		for i in fp:
			if i[0] == '#': continue
			if t>T: break
			omim,genes,score = i.rstrip().split('\t')
			go_seeds.extend(genes.split(','))
			t += 1
		fp.close()
		go_seeds = list(set(go_seeds))
		
		msg = 'total [%d] genes are chosen for GO seeds in [%d] out of [%d] diseases\n'%(len(go_seeds),T,num_omim)
		msg += 'done. [%s]'%job_name
		lib_utils.msgout('notice',msg); self.logger.info(msg)
		
		return go_seeds
	
	def gene_ontology_enrichment(self):
		'''
		Objective:Gene-ontology enrichment (select private members of purturbed gene that highly matched with phenotypic-scored genes and assign predicted phenotypic score instead of assigning de-novo prior)
		Input:
			-pheno_dmg = {gene1:0.2,gene2:0.9,...} #e.g. phenotype score
			-genetic_dmg = {gene2:0.4,gene3:0.3,...} #e.g. genetic score
		'''
		job_name = 'gene_ontology_enrichment'
		msg = 'enriching perturbed genes with GO semantic similarity [%s] ...'%job_name
		lib_utils.msgout('notice',msg); self.logger.info(msg)
		
		# collect genes from both phenotype and genotype perturbation
		pgenes = list(self.pheno_dmg.keys())
		P = len(pgenes)
		
		msg = 'total phenotypic genes before enrichment:%d' % P
		lib_utils.msgout('notice',msg,job_name); self.logger.info(msg)

		ggenes = list(self.genetic_dmg.keys())
		msg = 'total perturbed genes:%d' % len(ggenes)
		lib_utils.msgout('notice',msg,job_name); self.logger.info(msg)
		
		# draw a venn diagram and get genes not reported by phenotypes among genetic dmg genes
		priv_ggenes = lib_utils.difference(pgenes, ggenes)
		msg ='the number of genes not associated with the given phenotypes:%d' % len(priv_ggenes)
		lib_utils.msgout('notice',msg,job_name); self.logger.info(msg)

		# to collect genes highly matched to do GO enrichment
		pgenes2 = self.get_GO_seeds(self.dm.top_seed_rate) #update self.go_seeds

		#query high-scored phenotype genes against private genetic-perturbed genes and bring high-matched ones
		msg='quering total [%d] seed phenotype genes into SQL ...' % len(pgenes2)
		lib_utils.msgout('notice',msg,job_name); self.logger.info(msg)
		go = geneontology.Geneontology()
		goSimScores = go.get_funsim(pgenes2,priv_ggenes,min_score=self.dm.gosim_min)
	
		# updating the original phenotype damage score
		# weighting enriched phenotype matching score to the gene not reported in the original phenotypes
		pheno_delta = []
		for pair, go_sc in goSimScores.iteritems():
			
			#search for a gene matched to seed pheno gene
			if pair[0] not in self.pheno_dmg:
				gene_enriched = pair[0]
				seed_sc = self.pheno_dmg[pair[1]]
			elif pair[1] not in self.pheno_dmg:
				gene_enriched = pair[1]
				seed_sc = self.pheno_dmg[pair[0]]
				
			#initialize score
			if gene_enriched not in self.pheno_dmg:
				self.pheno_dmg[gene_enriched] = 0.
			
			#keep only maximum
			indirect_sc = seed_sc*go_sc*self.dm.go_penalty
			if indirect_sc > self.pheno_dmg[gene_enriched]:
				self.pheno_dmg[gene_enriched]=indirect_sc
				if gene_enriched not in pheno_delta:
					pheno_delta.append(gene_enriched)
		
		P_delta = len(self.pheno_dmg.keys()) - P
		
		msg='Total %d perturbed genes are added by phenotype gene enrichment!\ndone. [%s]'%(P_delta,job_name)
		lib_utils.msgout('notice',msg); self.logger.info(msg)
		msg='genes enriched by GO:[%s]'%lib_utils.joined(pheno_delta,',')
		lib_utils.msgout('notice',msg); self.loger.info(msg)

	def combine_damage_scores(self):

		# Gene-ontology enrichment (select private members of purturbed gene that highly matched with phenotypic-scored genes and assign predicted phenotypic score instead of assigning de-novo prior)
		job_name = 'combine_damage_scores'
		
		msg='combining both phenotypes[%s] and geneotype[%s] damage scores ... [%s]' %\
			(self.hpo_query, self.vcf, job_name)
		lib_utils.msgout('notice',msg);self.logger.info(msg)
		
		if self.dm.top_seed_rate>0.:
			#to select perturbed genes whose GO is highly similar to phenotype genes 
			self.gene_ontology_enrichment()

		# to obtain min damage score for both pheno and genetic perturb
		pdmg_min = lib_utils.get_stat_dic(self.pheno_dmg, 'min')
		if pdmg_min == 0.:
			raise ValueError('pheno has 0 dmg score[self.pheno_dmg]')
			
		gdmg_min = lib_utils.get_stat_dic(self.genetic_dmg, 'min')
		if gdmg_min == 0.:
			raise ValueError('genetic has 0 dmg score[self.genetic_dmg]')
		
		msg = 'calculating damage scores in a Bayesian framework...'
		lib_utils.msgout('notice',msg,job_name);self.logger.info(msg)
		
		for gene in self.genetic_dmg.keys():
			self.gene_dmg[gene] = PerturbedGene()
			pdmg = pdmg_min * self.dm.prior
			if gene in self.pheno_dmg:
				pdmg = self.pheno_dmg[gene]
			pdmg *= self.dm.ptwt
			
			gdmg = self.genetic_dmg[gene]
			self.gene_dmg[gene].gdmg = gdmg

			gdmg *= self.dm.gtwt
			self.gene_dmg[gene].score = pdmg * gdmg / (pdmg * gdmg \
																					+ (1. - pdmg) * (1. - gdmg))

			#skip normalization since it will be done in heat_diffusion
			
		msg = 'done. [%s]'%job_name
		lib_utils.msgout('notice',msg);self.logger.info(msg)

		return self.gene_dmg

	def run_heatdiffusion(self):
		job_name = 'heat_diffusion'

		pagerank.heat_diffusion(self.entries['string_link'],\
													self.entries['esp_to_gene'],\
													self.gene_dmg,self.rank_fn,self.logger)

	def run_vcf2xls(self):
		job_name = 'run_vcf2xls'
		msg = 'converting vcf file to excel file [%s] ...'%job_name
		lib_utils.msgout('notice',msg); self.logger.info(msg)
		
		rank_fn_tmp = self.rank_fn+'.tmp'
		
		cmd = ["cut","-f1,2",self.rank_fn,"|","grep","-v","'#'",">",rank_fn_tmp]
		self.run_cmd(cmd, "extract_pred_rank")
		
		self.xls = self._assign_out_fn('divine','xls')
		
		cmd = ["python", self.entries['vcf2xls'], \
					"-i", self.vcf, \
					"-o", self.xls, \
					"-l", self.log_dir, \
					"-g", rank_fn_tmp]
		self.run_cmd(cmd, job_name)
		
		msg = 'done. [%s]'%job_name
		lib_utils.msgout('notice',msg); self.logger.info(msg)
		
		os.unlink(rank_fn_tmp)
	
	def rank_pheno_gene(self):
		job_name = 'rank_pheno_gene'
		
		msg = 'selecting genes matched by patient phenotypes ... [%s;%s]'%(job_name,self.hpo_query)
		lib_utils.msgout('notice',msg); self.logger.info(msg)

		self.pheno_dmg_fn = self._assign_out_fn(job_name,'tsv')
		tmp_fn = '%s.tmp' % self.pheno_dmg_fn
		fp2=open(tmp_fn,'w')
		fp2.write('#gene\tphenotypic_score\n')
		for gene,dmg_score in self.pheno_dmg.iteritems():
			fp2.write('%s\t%g\n'%(gene,dmg_score))
		fp2.close()
		
		lib_utils.sort_tsv_by_col2(tmp_fn,[2],['gr'],False,self.pheno_dmg_fn)
		msg = 'done. [%s]'%job_name
		os.unlink(tmp_fn)
		lib_utils.msgout('notice',msg); self.logger.info(msg)
	
	def ranking_vcf(self):
		'''
		this function is obsolete and replaced by vcf2xls_varant()
		'''
		import gcn.lib.io.vcf as vcf
		job_name = 'ranking_vcf'
		
		msg = 'annotating Divine prediction score into filtered VCF ... [%s;%s]'%(job_name,self.vcf)
		lib_utils.msgout('notice',msg); self.logger.info(msg)
		
		ranked_vcf = '%s.ranked'%self.vcf
		ostream = open(ranked_vcf, 'w')
		v = vcf.VCFParser(self.vcf)

		v.add_meta_info("DVN", "1", "Float",\
			"Gene damage score predicted by Divine:%s"%self.command)

		v.writeheader(ostream)

		for rec in v:
			v.parseinfo(rec)
			vpop = vp.parse(rec.info)
			max_dmg_sc = 0.
			for altnum, val in vpop.items():
				for gene, gd in val.items():
					if gene in self.gene_dmg:
						if self.gene_dmg[gene]>max_dmg_sc:
							max_dmg_score = self.gene_dmg[gene]
			rec.info.DVN = max_dmg_score
			v.write(ostream, rec)
			
		ostream.close()
		v.stream.close()
		
		os.rename(ranked_vcf,self.vcf)
		msg= 'done. [%s]'%job_name
		lib_utils.msgout('notice',msg); self.logger.info(msg)

def main():
	parser = argparse.ArgumentParser(description="Divine (v%s) [author:%s]"%(VERSION,author_email))
	parser.add_argument('-q', dest='hpo_query_fn', required=False, default=None, help='Input patient HPO file. A file contains HPO IDs (e.g., HP:0002307), one entry per line. Refer to http://compbio.charite.de/phenomizer or https://mseqdr.org/search_phenotype.php')
	parser.add_argument('-v', dest='vcf', required=False, default=None, help='input vcf file')
	parser.add_argument('-o', action='store', dest='out_dir', required=False, default=None, help='output directory without white space. If not exist, it will create the directory for you.')
	parser.add_argument('-d', action='store', dest='exp_tag', required=False, default=None, help='specify experiment tag without white space. The tag will be contained in the output file name.')
	parser.add_argument('-I', action='store', dest='indel_mode', required=False, type=int, default=1, help='the level of fidelity of indell call in VCF [1]:low (e.g., samtools), 2:high (GATK haplotype caller)')
	
	parser.add_argument('-r', action='store', dest='top_seed_rate', required=False, type=float, default=0.0015, help='the rate of choosing matched diseases from top for gene enrichment [0.0015]; set to 0. to disable')
	
	parser.add_argument('--wes_mask', action='store_const', dest='wes_mask', required=False, default=False, const=True, help='to make the annotation process faster; the annotation process only runs on RefSeq coding regions [False]')
	parser.add_argument('--no_cadd', action='store_const', dest='no_cadd', required=False, default=False, const=True, help='disable CADD [False]')
	parser.add_argument('--hgmd', action='store_const', dest='hgmd', required=False, default=False, const=True, help='enable HGMD (requires a license) [False]')
	parser.add_argument('-k', action='store', dest='vknown', required=False, default=1, type=int, help='apply variant-level pathogenic annotation (e.g., either ClinVar or HGMD) to prioritization strategy [1:Yes], 0:No')

	parser.add_argument('--reuse', action='store_const', dest='reuse', required=False, default=False, const=True, help='Reuse previous annotation file (divine.vcf) if it is available [False]')
	parser.add_argument('-c', dest='capkit', required=False, default='SureSelect_V6', help='capture kit symbol [SureSelect_V6]')
		
	args = parser.parse_args()
	
	lib_utils.msgout('banner','Divine (v%s) | contact to %s for any question/error report/feedback'%(VERSION,author_email))
	
	# get a Divine instance/configure program condition and inputs 
	dv = Divine(args)
	
	# analyze phenotype if avaiable
	if dv.hpo_query:
		msg = 'analyzing query phenotypes on [%s] ...'%dv.hpo_query
		lib_utils.msgout('notice',msg); dv.logger.info(msg)
		
		# to get disease score associated with the given phenotypes (syscall)
		dv.hpo_to_diseases()

		# to get gene scores associated with the disease
		_ = dv.disease_to_genes()
		
		msg = 'done. [phenotype analysis]'
		lib_utils.msgout('notice',msg); dv.logger.info(msg)
		 
	# analyze genotype if available
	if dv.vcf:
		msg = 'analyzing variants on [%s] ...'%dv.vcf
		lib_utils.msgout('notice',msg); dv.logger.info(msg)

		# to create an instance of varant
		dv.vannotate(args.reuse)
		
		# to apply basic filter/tagging variant class
		dv.vfilter()
		
		# to predict genetic damage score
		dv.variants_to_genedmg()
		
		# combine two damage scores
		_ = dv.combine_damage_scores()
		
		# pagerank
		dv.run_heatdiffusion()
		
		# to generate an excel file report with a ranking score
		dv.run_vcf2xls()
		msg = 'done [check %s;%s].'%(dv.vcf, dv.xls)
		lib_utils.msgout('notice',msg,'Divine'); dv.logger.info(msg)
	else:
		# print ranked phenotype score per gene
		dv.rank_pheno_gene()
		msg = 'done [check %s;%s].'%(dv.hpo2disease_fn,dv.pheno_dmg_fn)
		lib_utils.msgout('notice',msg,'Divine'); dv.logger.info(msg)

if __name__ == '__main__':
	main()
