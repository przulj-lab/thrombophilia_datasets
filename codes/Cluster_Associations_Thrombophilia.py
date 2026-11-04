# -*- coding: UTF-8 -*-

""" Non-negative matrix tri-factorization (numpy)"""
# Author: NMD

import networkx as nx
import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns


def load_diseasegene(fname):
	genes=set()
	ifile=open(fname, 'r')
	ifile.readline()
	for line in ifile.readlines():
		lspt = line.strip().split('\t')
		if len(lspt)>2:
			gene = lspt[0][:-1].lower()
			genes.add(gene)
	ifile.close()
	return genes



def Plot_Mat(matrix, fname, limit=False):
	if limit == True:
		ax = sns.heatmap(matrix, linewidth=0., vmax=30)
	else:
		ax = sns.heatmap(matrix, linewidth=0.)
	plt.savefig(fname)
	plt.clf()


def Load_Clusters(fname):
	Clusters = {}
	ifile = open(fname, 'r')
	for line in ifile.readlines():
		lspt = line.strip().split('\t')
		if len(lspt)> 1:
			g = lspt[0]
			c = int(lspt[1])
			if c not in Clusters:
				Clusters[c] = []
			Clusters[c].append(g)
	ifile.close()
	return Clusters


#Load network
def Load_Network(fname):
	print "Loading Molecular Network"
	net = nx.read_edgelist(fname)
	nodes = [n for n in net.nodes()]
	nodeset = set(nodes)
	nb_nodes = len(nodes)
	nodes2ind = {}
	for n in range(nb_nodes):
		nodes2ind[nodes[n]]=n
	print " - %i genes, %i interactions"%(net.number_of_nodes(), net.number_of_edges())
	return net, nodes, nodes2ind





#Computing laplacian matrix, L
def Make_Laplacian(net, net_nodes, net_n2i):
	print "Computing laplacian matrix"
	nb_nodes = len(net_nodes)
	A = np.zeros((nb_nodes,nb_nodes))
	D = np.zeros((nb_nodes,nb_nodes))
	for n1 in range(nb_nodes):
		D[n1][n1] = float(net.degree(net_nodes[n1]))
	
	for e in net.edges():
		n1=net_n2i[ e[0] ]
		n2=net_n2i[ e[1] ]
		A[n1][n2] = 1.
		A[n2][n1] = 1.
	
	L = D-A
	return L

#Computing laplacian matrix, L
def Make_Adj(net, net_nodes, net_n2i):
	print "Computing laplacian matrix"
	nb_nodes = len(net_nodes)
	A = np.zeros((nb_nodes,nb_nodes))
	
	for e in net.edges():
		if e[0] in net_n2i and e[1] in net_n2i:
			n1=net_n2i[ e[0] ]
			n2=net_n2i[ e[1] ]
			A[n1][n2] = 1.
			A[n2][n1] = 1.
	
	return A


def Load_Annotations(fname, net_nodes, net_n2i):
	print "Loading patient mutation profiles"
	patients = []
	M = []
	nb_nodes = len(net_nodes)
	with open(fname, 'r') as f:
		line = f.readline()
		genes = line.strip().split('\t')
		#print genes
		while True:
			line = f.readline()
			if line: 
				lspt = line.strip().split('\t')
				if len(lspt)>10:
					patients.append(lspt[0])
					row = [0. for i in range(nb_nodes)]
					for t in range(len(genes)):
						val = int(lspt[t+1])
						gene = genes[t]
						if (gene in net_n2i) and (val == 1):
							row[ net_n2i[gene] ] = 1.
					M.append( row )
					print " - Patient %s: %i mutated gene in the network"%(lspt[0], sum(row))
			else:
				break
	#print patients
	return np.array(M), patients


#Loading and preparing all the data
PPI_net = nx.read_edgelist("./Nets/Human_PPI_Biogrid.net")


#Considering the union of genes from the two networks
nodes_p = set([n for n in PPI_net.nodes()])
nodeset = nodes_p
nodes = list(nodeset)
nb_nodes = len(nodes)
nodes2ind = {}
for n in range(nb_nodes):
		nodes2ind[nodes[n]]=n


#Formating matrices for NMTF
R1, patients = Load_Annotations("./Mutations/Mutations.txt", nodes, nodes2ind) 

#thrombolphilia genes
TG = load_diseasegene("Thrombophilia_gene.txt")


#PPI clusters
Clusters = Load_Clusters('Gene_clusters_All.txt')
k = len(Clusters)
cluster_names = Clusters.keys()
cluster_names.sort()



#COEX clusters
Clusters = Load_Clusters('Gene_clusters_All_P.txt')
k = len(Clusters)
cluster_names = Clusters.keys()

ofile = open("Cluster_Association_Thombophilia_All_P.txt", 'w')

for c in range(len(cluster_names)):
	for gene in Clusters[cluster_names[c]]:
		idg = nodes2ind[gene]
		if gene in TG:
			ofile.write("%s\t%s\t%s\t%i\t%i\t%i\t%i\t%i\n"%(str(c), str(cluster_names[c]), gene, R1[0][idg], R1[1][idg], R1[2][idg], R1[3][idg], R1[4][idg]))

ofile.close()