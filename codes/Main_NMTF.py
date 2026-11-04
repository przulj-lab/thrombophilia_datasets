# -*- coding: UTF-8 -*-

""" Non-negative matrix tri-factorization (numpy)"""
# Author: NMD

import networkx as nx
import numpy as np
import MSNMTF2p
import math
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns



def Plot_Mat_P(matrix, fname, limit=False):
	ax = sns.heatmap(matrix, linewidth=0.)
	
	#plt.ylabel("Patients", fontsize=20)
	#plt.xlabel("Clusters", fontsize=20)
	plt.xticks([0.5,1.5], ["1", "2"])
	plt.yticks([0.5,1.5,2.5,3.5,4.5], ["B1", "B2", "D1", "S1", "S2"])
	plt.tick_params(axis='both', which='major', labelsize=18)
	plt.savefig(fname)
	plt.clf()



def Plot_Mat(matrix, fname, limit=False):
	if limit == True:
		ax = sns.heatmap(matrix, linewidth=0., vmax=0.01)
	else:
		ax = sns.heatmap(matrix, linewidth=0.)
	plt.savefig(fname)
	plt.clf()


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
	idx = net_n2i["fgb"]
	print M[0][idx], M[1][idx], M[2][idx], M[3][idx], M[4][idx]
	return np.array(M), patients


#Loading and preparing all the data
PPI_net = nx.read_edgelist("./Nets/Human_PPI_Biogrid.net")
COEX1_net= nx.read_edgelist("./Nets/Human_COEX_v7_1.net")
GI_net= nx.read_edgelist("./Nets/Human_GI_new.net")


#Considering the union of genes from the two networks
nodes_p = set([n for n in PPI_net.nodes()])
#nodes_c = set([n for n in COEX_net.nodes()])
#nodeset = nodes_p & nodes_c
nodeset = nodes_p
nodes = list(nodeset)
nb_nodes = len(nodes)
nodes2ind = {}
for n in range(nb_nodes):
		nodes2ind[nodes[n]]=n


#Formating matrices for NMTF
R1, patients = Load_Annotations("./Mutations/Mutations.txt", nodes, nodes2ind) 
R2 = Make_Adj(PPI_net, nodes, nodes2ind)
R3 = Make_Adj(COEX1_net, nodes, nodes2ind)
R4 = Make_Adj(GI_net, nodes, nodes2ind)


P = np.array( [  [0.5,0.],[0.,1.],[0.5,0.],[0.5,0.],[0.5,0.] ] )
Plot_Mat_P(P, 'P_All_P.png')

#print R

#Running NMTF for PPI + COEX
Solver = MSNMTF2p.MSNMTF(max_iter=100, verbose = 10)
k_p = 2
k_g = int( math.sqrt(nb_nodes/2.) )

print "Clustering parameters: %i, %i patient/gene clusters"%(k_p,k_g)
P, S, G, U, V, W = Solver.Solve_MUR(R1, R2, R3, R4, P, k_p, k_g)

print P

Clusters_P = Solver.Get_Clusters(P, patients)
print Clusters_P

#print G

Clusters_G = Solver.Get_Clusters(G, nodes)
ofile = open("Gene_clusters_All_P.txt", 'w')
for i in range(len(Clusters_G)):
	for j in range(len(Clusters_G[i])):
		gene = Clusters_G[i][j]
		idx = nodes2ind[gene]
		val = G[idx][i]
		ofile.write("%s\t%i\t%f\n"%(gene, i, val))
ofile.close()


Plot_Mat_P(P, 'P_All_P.png')
Plot_Mat(S, 'S_All_P.png')
Plot_Mat(G, 'G_All_P.png', True)


