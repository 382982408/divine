#!/usr/bin/env python
'''
COPYRIGHT (C) 2016 changjin.hong@gmail.com
author: changjin.hong@gmail.com
'''

class DmgCoeff:
	
	def __init__(self,indel_mode=1,top_seed_rate=0.0015,logger=None):
		msg='setting up damaging factors ...'
		print msg
		if logger: logger.info(msg)
		
		self.vncoding = 0.050
		self.vintronic = 0.055
		self.vexonic = 0.35
		self.vsplice = 0.45
		self.vsplice_syn = 0.060
		self.warning = 0.070
		
		# damage factor for known pathogenic variants 
		self.cgene = 0.04
		self.cintronic = 0.25
		self.cexonic = 0.50
		
		# logistic regression coeff for damage w.r.t MAF
		self.beta1 = 0.25
		self.beta2 = 1.5

		# hetero/homozygous/compound het allele impact for rare SNVs
		self.het_rare = 0.05
		self.hom_rare = 0.85
		
		self.het_denovo_snp = 0.75
		self.hom_denovo_snp = 1.0
		
		# hetero/homozyous/compound het allele impact w.r.t. the mutation type for denovo SNVs 
		
		if indel_mode==1:
			self.het_denovo_indel = 0.1
			self.hom_denovo_indel = 0.25
		elif indel_mode==2:
			self.het_denovo_indel = 0.50
			self.hom_denovo_indel = 1.0

		# default average transcript length when its length is unknown
		# https://www.quora.com/Whats -the-average-size-of-a-human-protein-in-kDa
		self.avg_protein_len = 480
		
		self.prior = 1.25
		self.top_seed_rate = top_seed_rate
		self.go_penalty = 0.4
		self.gosim_min = 0.95
		self.ptwt = 0.5
		self.gtwt = 0.5
		self.prwt = 0.3
		self.hpo_max_queries = 15
