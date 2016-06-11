#!/usr/bin/env python
'''
COPYRIGHT (C) 2016 changjin.hong@gmail.com
author: changjin.hong@gmail.com
'''

import os, re, math
import numpy as np
from scipy.sparse import coo_matrix
from sklearn.preprocessing import normalize
from gcn.lib.utils import lib_utils

def extract_ensembl_protein(protein):
	mObj=re.search(r'\d+\.(\w+)',protein)
	return mObj.group(1)

def get_sparse_elements(proteinLinkFile,min_edge_weight):
	'''
	to store ppi network
	input: dProtein2gene, dGenes(whether the gene is in ppi or not)- protein-gene relation; proteinLinkFile- ppi link
	output: update dProtein2gene, dGenes when add_dangled is enabled. Store ppi and lnkProteins
	'''
	#read string DB and assign an integer to each protein symbol
	fp = lib_utils.open2(proteinLinkFile,'r')
	
	nNodes = 0
	linked = [-1,-1]
	dProtein2num = {}
	lnkProteins = []
	ppi = [[],[],[]] #from protein, to protein, link weight
	
	lib_utils.msgout('notice', 'preparing a genetic network matrix. Please, be patient......', 'pagerank|heat_diffusion')
	#store col,row,weight from ppi file
	fp.next()
	for i in fp:
		#print '%s'%i #debug
		linked[0],linked[1],weight=i.rstrip().split(' ')
		weight = float(weight)
		if weight < min_edge_weight: continue
		for c in range(2):
			protein=extract_ensembl_protein(linked[c])
			
			#to register a protein node
			if not protein in dProtein2num:
				dProtein2num[protein]=nNodes
				lnkProteins.append(protein) #item index corresponds to a node number of the protein
				nNodes+=1

			ppi[c].append(dProtein2num[protein])
		ppi[2].append(weight)
	fp.close()

	dProtein2num = None
	
	return nNodes,ppi,lnkProteins

def parse_gene_prior(genePriorFile,mode):
	#assume that genePriorFile is likes
	#gene\trank\tp-value
	
	dGene ={}
	fp=lib_utils.open2(genePriorFile,'r')
	scoreSum = 0.
	maxScore = 0.
	fp.next()
	for i in fp:
		gene,rank,score=i.rstrip().split('\t')
		if mode == 1:
			dGene[gene] = float(rank)
		elif mode == 2:
			dGene[gene] = float(score)
			
		scoreSum+=dGene[gene]
		
		if dGene[gene]>maxScore:
			maxScore = dGene[gene]
	fp.close()
	
	#normalize
	if mode == 1:
		for gene, val in dGene.iteritems():
			dGene[gene]=(maxScore-val+1.)/scoreSum

	return dGene

def protein_to_gene(esp2geneFile):
	
	dProtein2gene = {}
	
	fp = lib_utils.open2(esp2geneFile,'r')
	fp.next()
	for i in fp:
		gene,protein = i[:-1].split('\t')
		if gene and protein:
			if protein not in dProtein2gene:
				dProtein2gene[protein]=gene
			
	fp.close()
	return dProtein2gene

def assign_node_prior(cPerturbedGene,dProtein2gene,nNodes,lnkProteins,verbose=False):
	
	Y = np.zeros(shape=(nNodes,1),dtype=float)
	perturbedGeneOnLnk = lib_utils.list_to_dic(list(cPerturbedGene.keys()),False)
	for n,protein in enumerate(lnkProteins):
		if protein in dProtein2gene:
			gene=dProtein2gene[protein]
			if gene in cPerturbedGene:
				Y[n] = cPerturbedGene[gene].score
				perturbedGeneOnLnk[gene] = True
		else:
			if verbose:
				print 'protein[%s] is not found from a gene-product mapping database!'%protein
	
	#count the number of perturbed genes not on the network
	
	nNodes0 = len([1 for isLnk in perturbedGeneOnLnk.values() if isLnk is False])
	Y0 = np.zeros(shape=(nNodes0,1),dtype=float)
	Y0_sum = 0.
	
	dangledGenes = []
	n = 0
	#print nNodes0
	for gene,isLnk in perturbedGeneOnLnk.iteritems():
		#print gene
		if not isLnk:
			dangledGenes.append(gene)
			Y0[n] = cPerturbedGene[gene].score
			n += 1
			
	prior_raw_sum = np.sum(Y) + np.sum(Y0)
	#normalize in a range of 0 to 1.
	Y = Y/prior_raw_sum
	Y0 = Y0/prior_raw_sum
	return Y,Y0,dangledGenes

def convert_node2gene(FinalNodeScores,PerturbedGenes,dProtein2gene,lnkProteins,rank_fn):
	
	nodeScores,dangledScores = FinalNodeScores
	cPerturbedGenes,dangledGenes = PerturbedGenes
	
	rank_fn2 = lib_utils.file_tag2(rank_fn, 'tmp', None)
	fp2=lib_utils.open2(rank_fn2,'w')
	fp2.write('#gene\tpredicted_score[-1/log10(x)]\tseed_score\n')	
	for n,protein in enumerate(lnkProteins):
		seed_score = 0.
		gene = protein
		genetic_dmg_score = 0.
		if protein in dProtein2gene:
			gene = dProtein2gene[protein]
			if gene in cPerturbedGenes:
				seed_score = cPerturbedGenes[gene].score
				genetic_dmg_score = cPerturbedGenes[gene].gdmg
				
		pred_score = 0.
		if nodeScores[n]>0:
			pred_score = -1./math.log10(nodeScores[n])
		if genetic_dmg_score>0.:
			fp2.write('%s\t%g\t%g\n'%(gene,pred_score,seed_score))
	
	#add dangled node score
	for n,gene in enumerate(dangledGenes):
		pred_score = 0.
		if dangledScores[n]>0:
			pred_score = -1./math.log10(dangledScores[n])
		if cPerturbedGenes[gene].gdmg>0.:
			fp2.write('%s\t%g\t%g\n'%(gene,pred_score,cPerturbedGenes[gene].score))

	fp2.close()
	
	#sort by score
	lib_utils.sort_tsv_by_col2(rank_fn2, [2], ['gr'], False, rank_fn)
	os.unlink(rank_fn2)

def cal_array_distance(npVec1, npVec2):
	return sum(abs(npVec1-npVec2))[0]

def heat_diffusion_core(gamma,M,alpha,A,Y,Y0,maxIter,logger=None):
	job_name = 'pagerank|heat_diffusion_core'
	
	N=len(Y)
	s = np.zeros(shape=(N,1))
	
	N0=len(Y0)
	s0 = np.zeros(shape=(N0,1))
	
	epsilon = 1e-4
	iter = 0

	msg = 'running heat diffusion on [%dx%d, gamma=%g, alpha=%g, max_iter=%d, M=%d]. Please, be patient......' % (N,N,gamma,alpha,maxIter,M)
	lib_utils.msgout('notice', msg, job_name)
	if logger: logger.info(msg)
		
	e = 1.
	while (e>=epsilon and iter<maxIter):
		#heat diffusion
		s_new  = (1.-gamma/M)*s  + (gamma/M)*(alpha*A.dot(s)+(1.-alpha)*Y)
		s0_new = (1.-gamma/M)*s0 + (gamma/M)*(1.-alpha)*Y0
		
		#normalize
		s_new = s_new/(np.sum(s_new)+np.sum(s0_new))
		s0_new = s0_new/(np.sum(s_new)+np.sum(s0_new))
		
		e = cal_array_distance(s_new,s)+cal_array_distance(s0_new,s0)
		
		s  = s_new
		s0 = s0_new
		iter+=1
	
	msg = 'done. [iteration:%d/%d,e:%g]'%(iter,maxIter,e)
	lib_utils.msgout('notice', msg, job_name)
	if logger: logger.info(msg)
	
	return s,s0


def heat_diffusion(proteinLinkFile,esp2geneFile,cPerturbedGene,rank_fn,logger,\
									min_edge_weight=350.,M=100,gamma=2.,alpha=0.9):
	
	job_name = 'pagerank.heat_diffusion'
	msg = 'transferring gene product to a matrix [%s;%s]'%(proteinLinkFile,esp2geneFile)
	lib_utils.msgout('notice',msg,job_name);logger.info(msg)
	
	nNodes,ppi,lnkProteins=get_sparse_elements(proteinLinkFile,min_edge_weight)
	
	msg = 'creating a genetic network matrix[%d x %d].'%(nNodes,nNodes)
	lib_utils.msgout('notice',msg,job_name);logger.info(msg)
	
	A = coo_matrix((ppi[2], (ppi[0],ppi[1])), dtype=np.float, shape=(nNodes, nNodes))
	ppi = None
	
	#convert to csr_matrix for faster/reliable matrix operation
	msg = 'reformatting the genetic network matrix.'
	lib_utils.msgout('notice',msg,job_name);logger.info(msg)
	A = A.tocsr()
	
	#normalize PPI matrix
	msg = 'normalizing the genetic network matrix.'
	lib_utils.msgout('notice',msg,job_name);logger.info(msg)
	A = normalize(A, norm='l1', axis=0)

	#to get a map between protein and gene
	dProtein2gene = protein_to_gene(esp2geneFile)
	
	#read phenotype_gene_ranks_file
	Y,Y0,dangledGenes = assign_node_prior(cPerturbedGene,dProtein2gene,nNodes,lnkProteins)
	
	#core heat diffusion in recursion
	msg = 'running heat diffusion on genetic networks labeled by perturbed genes.'
	lib_utils.msgout('notice',msg,job_name);logger.info(msg)
	s,s0 = heat_diffusion_core(gamma,M,alpha,A,Y,Y0,maxIter=150,logger=logger)
	
	A = None
	
	#annotate gene to node
	msg = 'reporting ranked genes.'
	lib_utils.msgout('notice',msg,job_name);logger.info(msg)
	convert_node2gene([s,s0],[cPerturbedGene,dangledGenes],dProtein2gene,lnkProteins,rank_fn)
	
	msg = 'done.'
	lib_utils.msgout('notice',msg,job_name);logger.info(msg)