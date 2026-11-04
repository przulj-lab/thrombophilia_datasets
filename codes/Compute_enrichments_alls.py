###CODES in python 2!
import sys
import math
import os.path
import operator
import numpy as np
import Gnuplot as gp
import networkx as nx
import scipy.stats as ss
import matplotlib as plt
from scipy.stats import hypergeom
import scipy.cluster.hierarchy as sch
from scipy.cluster.vq import vq, kmeans, kmeans2, whiten

###Functions
def load_clusters(fname):
	clusters = {}
	nodes = set()
	
	ifile=open(fname,'r')
	for line in ifile.readlines():
		lspt = line.strip().split()
		if len(lspt)>1:
			gene = lspt[0].lower()
			cid = int(lspt[1])
			if cid not in clusters:
				clusters[cid]=set()
			clusters[cid].add(gene)
			nodes.add(gene)
	ifile.close()
	carray = []
	for id in clusters.keys():
		carray.append(list(clusters[id]))
	return nodes, carray, 

# Benjamini-Hochberg p-value correction for multiple hypothesis testing
def p_adjust_bh(p):
    p = np.asfarray(p)
    by_descend = p.argsort()[::-1]
    by_orig = by_descend.argsort()
    steps = float(len(p)) / np.arange(len(p), 0, -1)
    q = np.minimum(1, np.minimum.accumulate(steps * p[by_descend]))
    return q[by_orig]
    #return p


# Read GO terms annotations
def load_GO(fname, genes):
	geneset=set(genes)
	GO_counter = {}
	gene2GO = {}
	#loading raw gene & go relationships
	fRead = open(fname, "r")
	for line in fRead.readlines():
		lspt = line.strip().split(';')
		if len(lspt)>1:
			geneid = lspt[0].lower()
			term = lspt[1]
			if geneid in geneset:
				if term not in GO_counter:
					GO_counter[term] = 0
				GO_counter[term] += 1
				if geneid not in gene2GO:
					gene2GO[geneid] = set()
				gene2GO[geneid].add(term)
	fRead.close()
	#filter out GO terms that annotates only one gene
	GO_counter2 = {}
	gene2GO2 = {}
	removed = set()
	for term in GO_counter:
		if GO_counter[term] > 1:
			GO_counter2[term] = GO_counter[term]
		else:
			removed.add(term)
	for gene in gene2GO:
		genego = gene2GO[gene].difference(removed)
		if len(genego)>0:
			gene2GO2[gene]=genego
	return [GO_counter2, gene2GO2]



def go_enrichment(clusters, gene2go, go_counts):
	M = len(gene2go) # total number of annotated genes in the dataset
	data = []
	i = -1
	enrichment = [[] for j in range(len(clusters))]
	cts = 0
	
	cls = []
	gos = []
	pvals = []
	
	NE_list = []
	
	for cluster in clusters:
		i+=1
		annotated_genes = []
		for gene in cluster:
			if gene in gene2go:
				annotated_genes.append(gene)
		N = len(annotated_genes) # number of annotated genes in the cluster
		#print N
		annotation_set = set()
		for gene in annotated_genes:
			for term in gene2go[gene]:
				annotation_set.add(term)
		#print len(annotation_set)
		for term in annotation_set:
			K = go_counts[term] # number of genes annotated with the given term in all the data
			X = 0	# number of gene in the clusters that are annotated with the go term
			for gene in annotated_genes:
				if term in gene2go[gene]:
					X += 1
			pval = hypergeom.sf(X-1, M, K, N) # faster and more accurate than 1-cdf
			cls.append(i)
			gos.append(term)
			pvals.append(pval)
			if pval <= 0.05:
				cts += 1
				#print "Cluster %i, term %s: X=%i, K=%i, pval = %s"%(i,term,X,K,str(pval))
				enrichment[i].append([term,pval])
			#print "%s %s"%(term,prb)
		
		#Mean normalized entropy:
		d = float(len(annotation_set))	# number of different annotations in cluster c
		nec = 0.
		if d>1.:
			Cl = float(len(cluster))	# number of gene in cluster c
			sum_pi = 0.
			#counting map
			counting_map = {}
			for gene in cluster:
				if gene in gene2go:
					for go in gene2go[gene]:
						if go not in counting_map:
							counting_map[go] = 0.
						counting_map[go] += 1.
			for go in counting_map.keys():
				pi = counting_map[go]/Cl	# percentage of genes in c annotated with the considered go term
				sum_pi += pi*math.log(pi)
			nec = (-1./(math.log(d)))*sum_pi
		NE_list.append( nec )
	
	#applying BH correction
	BHpvals = p_adjust_bh(pvals)
	BHcts = 0
	BH_enrichment = [[] for j in range(len(clusters))]
	enriched_genes = []
	for i in range(len(cls)):
		#print pvals[i], BHpvals[i]
		if BHpvals[i]<=0.05:
			BHcts += 1
			BH_enrichment[cls[i]].append([gos[i],BHpvals[i]])
	for i in range(len(clusters)):
			cluster_set = set()
			enriched_gos = BH_enrichment[i]
			for gene in clusters[i]:
				for go, pval in enriched_gos:
					if gene in gene2go:
						if go in gene2go[gene]:
							cluster_set.add(gene)
			enriched_genes.append(list(cluster_set))
	
	#print cts, BHcts
	MNE = sum(NE_list)/float(len(NE_list))
	#print "Mean Normalized Entropy = %s"%(str(MNE))
	#print BH_enrichment
	enr_cluster = 0
	total_cluster = 0
	for i in range(len(clusters)):
		if len(clusters[i])>0:
			total_cluster += 1.
			if len(BH_enrichment[i])>0:
				enr_cluster += 1
	perc_enr_cluster = 100.*enr_cluster/total_cluster
	perc_genes = sum([len(enriched_genes[i]) for i in range(len(clusters))])
	perc_genes = 100.*perc_genes/float(len(gene2go))
	return [BH_enrichment, MNE, perc_genes, perc_enr_cluster]


def output_enrichment(fname, C, BPbh, MFbh, CCbh, BPmne, MFmne, CCmne, go2name):
	#summary of cluster enrichments
	ofile = open(fname+"_Summary_clusters.csv", 'w')
	ofile.write("Cluster_ID\t#Genes\t#BP\t#MF\t#CC\n")
	lenC = len(C)
	cBP = [len(BPbh[i]) for i in range(lenC)]
	cMF = [len(MFbh[i]) for i in range(lenC)]
	cCC = [len(CCbh[i]) for i in range(lenC)]
	Perc_enr = [0., 0., 0., 0.]
	Real_nbclust = 0.
	for i in range(len(C)):
		if len(C[i]) == 0:
			ofile.write("%i\t 0 \t - \t - \t - \n"%(i+1))
		else:
			Real_nbclust += 1.
			if cBP[i] > 0:
				Perc_enr[0] += 1.
			if cMF[i] > 0:
				Perc_enr[1] += 1.
			if cCC[i] > 0:
				Perc_enr[2] += 1.
			if cBP[i] + cMF[i] + cCC[i] > 0:
				Perc_enr[3] += 1.
			ofile.write("%i\t%i\t%i\t%i\t%i\n"%(i+1, len(C[i]), cBP[i], cMF[i], cCC[i] ) )
	Perc_enr[0] = Perc_enr[0] / Real_nbclust
	Perc_enr[1] = Perc_enr[1] / Real_nbclust
	Perc_enr[2] = Perc_enr[2] / Real_nbclust
	Perc_enr[3] = Perc_enr[3] / Real_nbclust
	
	ofile.write("TOTAL\t \t%i\t%i\t%i\n"%( sum(cBP), sum(cMF), sum(cCC) ))
	ofile.write("\t \t \t \t%i\n"%( sum([sum(cBP), sum(cMF), sum(cCC)]) ))
	ofile.write("MNE\t \t%.3f\t%.3f\t%.3f\n"%( BPmne, MFmne, CCmne ))
	ofile.close()
	#details for GO-BP, GO-MF and GO-CC
	ofile = open(fname+"_Enrichments-GOBP.csv", 'w')
	for i in range(len(BPbh)):
		ofile.write("Cluster %i: Enriched_terms (BH corrected)\n"%(i+1))
		for pair in BPbh[i]:
			#GO version
			if pair[0] in go2name:
				ofile.write("%s, %.3e, %s\n"%(pair[0],pair[1], go2name[pair[0]]))
			else:
				ofile.write("%s, %.3e, ?\n"%(pair[0],pair[1]))
	ofile.close()
	ofile = open(fname+"_Enrichments-GOMF.csv", 'w')
	for i in range(len(MFbh)):
		ofile.write("Cluster %i: Enriched_terms (BH corrected)\n"%(i+1))
		for pair in MFbh[i]:
			#GO version
			if pair[0] in go2name:
				ofile.write("%s, %.3e, %s\n"%(pair[0],pair[1], go2name[pair[0]]))
			else:
				ofile.write("%s, %.3e, ?\n"%(pair[0],pair[1]))
	ofile.close()
	ofile = open(fname+"_Enrichments-GOCC.csv", 'w')
	for i in range(len(CCbh)):
		ofile.write("Cluster %i: Enriched_terms (BH corrected)\n"%(i+1))
		for pair in CCbh[i]:
			#GO version
			if pair[0] in go2name:
				ofile.write("%s, %.3e, %s\n"%(pair[0],pair[1], go2name[pair[0]]))
			else:
				ofile.write("%s, %.3e, ?\n"%(pair[0],pair[1]))
	ofile.close()
	return [sum(cBP), sum(cMF), sum(cCC), sum( [sum(cBP), sum(cMF), sum(cCC)] ), Perc_enr]


def output_gene_enrichment(fname, C, BPbh, MFbh, CCbh):
	#summary of cluster enrichments
	ofile = open(fname+"_Summary_gene.csv", 'w')
	ofile.write("Cluster_ID\t#Genes\t#BP\t#MF\t#CC\n")
	lenC = len(C)
	cBP = [len(BPbh[i]) for i in range(lenC)]
	cMF = [len(MFbh[i]) for i in range(lenC)]
	cCC = [len(CCbh[i]) for i in range(lenC)]
	for i in range(len(C)):
		if len(C[i]) == 0:
			ofile.write("%i\t 0 \t - \t - \t - \n"%(i+1))
		else:
			ofile.write("%i\t%i\t%i\t%i\t%i\n"%(i+1, len(C[i]), cBP[i], cMF[i], cCC[i] ) )
	ofile.write("TOTAL\t \t%i\t%i\t%i\n"%( sum(cBP), sum(cMF), sum(cCC) ))
	ofile.close()
	return [sum(cBP), sum(cMF), sum(cCC), sum( [sum(cBP), sum(cMF), sum(cCC)] )]


def Process_Enr(idir, goname):
	enr_clusters = []
	genes, cluster0 = load_clusters(idir)
	go_counts,gene2go = load_GO(goname, genes)
	bh, mne, eg, perc_cluster = go_enrichment(cluster0, gene2go, go_counts)
	
	return bh, mne, eg, perc_cluster


def Process_Enr_Rand(idir, goname):
	enr_clusters = []
	genes, cluster0 = load_clusters(idir)
	flat_list = [item for sublist in cluster0 for item in sublist]
	flat_list=[random.choice(flat_list) for _ in range(len(flat_list))]

	n=[]
	for i in cluster0:
		n.append(len(i))
	it = iter(flat_list)
	list_new=[[next(it) for _ in range(size)] for size in n]
	go_counts,gene2go = load_GO(goname, genes)
	bh, mne, eg, perc_cluster = go_enrichment(list_new, gene2go, go_counts)
	
	return bh, mne, eg, perc_cluster

###Enrichments of Genes in Reactome pathways and Genome Ontology terms
	
exps = ['RR', 'RP', 'BP', 'MF', 'CC']
enr_cp = []
enr_gp = []

#enrichments all
for exp in exps:
	goname = "./Annotations/Human_%s.lst"%(exp)
	bh, mne, eg, perc_cluster = Process_Enr("Gene_clusters_All.txt", goname)
	enr_cp.append(perc_cluster)
	enr_gp.append(eg)


enr_cc = []
enr_gc = []

#enrichments in all, fixed P
for exp in exps:
	goname = "./Annotations/Human_%s.lst"%(exp)
	bh, mne, eg, perc_cluster = Process_Enr("Gene_clusters_All_P.txt", goname)
	enr_cc.append(perc_cluster)
	enr_gc.append(eg)
	
    
	ofile = open("ENR_%s_All_P_prova.csv"%(exp), 'w')
	for c in range(len(bh)):
		for gt, p in bh[c]:
			ofile.write("%s\t%s\t%s\n"%(c, gt, p))
	ofile.close()

enr_cr = []
enr_gr = []

#Random enrichments
for exp in exps:
	goname = "./Annotations/Human_%s.lst"%(exp)
	bh, mne, eg, perc_cluster = Process_Enr_Rand("Gene_clusters_All_P.txt", goname)
	enr_cr.append(perc_cluster)
	enr_gr.append(eg)

#Plots

xp = [0,5,10,15,20]
xc = [1,6,11,16,21]
xr = [2,7,12,17,22]

matplotlib.rc('xtick', labelsize=22) 
matplotlib.rc('ytick', labelsize=22) 
plt.figure(figsize=(8,10), dpi=120)
plt.bar(xp, enr_cp, width=1, color="r", edgecolor="black")
plt.bar(xc, enr_cc, width=1, color="b", edgecolor="black")
plt.bar(xr, enr_cr, width=1, color="g", edgecolor="black")
plt.legend(["Default", "Fized P", "Random"], fontsize=20,bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(xc, ["RR", "RP", "BP", "MF", "CC"])
plt.ylabel("Enriched clusters (%)",fontsize=22)
plt.savefig("Enrichments_clusters_Alls_P_Random.png",bbox_inches='tight')


matplotlib.rc('xtick', labelsize=22) 
matplotlib.rc('ytick', labelsize=22) 
plt.figure(figsize=(8,10), dpi=120)
plt.bar(xp, enr_gp, width=1, color="r", edgecolor="black")
plt.bar(xc, enr_gc, width=1, color="b", edgecolor="black")
plt.bar(xr, enr_gr, width=1, color="g", edgecolor="black")
plt.legend(["Default", "Fized P", "Random"], fontsize=20,bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xticks(xc, ["RR", "RP", "BP", "MF", "CC"])
plt.ylabel("Enriched genes (%)",fontsize=22)
plt.savefig("Enrichments_genes_Alls_P_Random.png",bbox_inches='tight')
