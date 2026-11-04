# -*- coding: UTF-8 -*-

""" Non-negative matrix tri-factorization (numpy)"""
# Author: NMD

import numpy as np
import scipy.sparse.linalg as ssl

#  Input: R, a matrix relating n patients to m genes. Entry R_ij = 1 if patient i has a germline mutation on gene j, 0 otherwise
#         L, the laplacian matrix of a molecular network capturing interactions between the genes  
#
#                                                            
#           n×m                n×k_1                                  
#     ┌───────────────┐        ┌───┐                             
#     │               │        │   │     k_1×k_2       k_2×m                   
#     │               │        │   │     ┌───┐   ┌──────────────┐
#     │      R1       │    ≈   │ P │  x  │ S │ x │      G'      │
#     │               │        │   │     └───┘   └──────────────┘
#     │               │        │   │   
#     └───────────────┘        └───┘   
#
#           m×m                m×k_2                                  
#     ┌───────────────┐        ┌───┐                             
#     │               │        │   │     k_2×k_2       k_2×m                   
#     │               │        │   │     ┌───┐   ┌──────────────┐
#     │      R2       │    ≈   │ G │  x  │ U │ x │      G'      │
#     │               │        │   │     └───┘   └──────────────┘
#     │               │        │   │   
#     │               │        │   │   
#     └───────────────┘        └───┘   
#
# 
#           m×m                m×k_2                                  
#     ┌───────────────┐        ┌───┐                             
#     │               │        │   │     k_2×k_2       k_2×m                   
#     │               │        │   │     ┌───┐   ┌──────────────┐
#     │      R3       │    ≈   │ G │  x  │ V │ x │      G'      │
#     │               │        │   │     └───┘   └──────────────┘
#     │               │        │   │   
#     │               │        │   │   
#     └───────────────┘        └───┘   

#  With: S, G, U, V >= 0
#  With U given by the user
#  With the added contraints of bi-orthogonality: G'G = I, P'P = I
#
#
#





class MSNMTF:
	"""Compute Non-negative Matrix Tri-Factorization (NMTF)"""
	def __init__(self, max_iter=1000, verbose = 10):
		self.max_iter = max_iter
		self.verbose = verbose
		
	
	def Score(self, R1, R2, R3, R4, P, S, G, U, V, W, norm_R1, norm_R2, norm_R3, norm_R4):
		GT = np.transpose(G)
		
		Err1 = R1 - np.matmul(P, np.matmul(S, GT))
		norm_err1 = np.linalg.norm(Err1, ord='fro')
		
		Err2 = R2 - np.matmul(G, np.matmul(U, GT))
		norm_err2 = np.linalg.norm(Err2, ord='fro')
		
		Err3 = R3 - np.matmul(G, np.matmul(V, GT))
		norm_err3 = np.linalg.norm(Err3, ord='fro')
		
		Err4 = R4 - np.matmul(G, np.matmul(W, GT))
		norm_err4 = np.linalg.norm(Err4, ord='fro')
		
		rel_err1 = norm_err1 / norm_R1
		rel_err2 = norm_err2 / norm_R2
		rel_err3 = norm_err3 / norm_R3
		rel_err4 = norm_err4 / norm_R4
		
		return norm_err1+norm_err2+norm_err3+norm_err4, rel_err1, rel_err2, rel_err3, rel_err4
	
	
	def Get_Clusters(self, M, nodes):
		n,k = np.shape(M)
		Clusters = [[] for i in range(k)]
		for i in range(n):
			idx = np.argmax(M[i])
			Clusters[idx].append(nodes[i])
		return Clusters
		
	
	
	def Solve_MUR(self, R1, R2, R3, R4, P, k1, k2):
		print "Starting NMTF"
		n,m=np.shape(R1)
		norm_R1 = np.linalg.norm(R1, ord='fro')
		norm_R2 = np.linalg.norm(R2, ord='fro')
		norm_R3 = np.linalg.norm(R3, ord='fro')
		norm_R4 = np.linalg.norm(R4, ord='fro')
		
		print " * Using random initialization"
		S = np.random.rand(k1,k2)+1e-4
		G = np.random.rand(m,k2)+1e-4
		U = np.random.rand(k2,k2)+1e-4
		V = np.random.rand(k2,k2)+1e-4
		W = np.random.rand(k2,k2)+1e-4
		
		
		OBJ, REL1, REL2, REL3, REL4= self.Score(R1, R2, R3, R4, P, S, G, U, V, W, norm_R1, norm_R2, norm_R3, norm_R4)
		print " - Init:\t OBJ:%.4f\t REL1:%.4f\t REL2:%.4f\t REL3:%.4f\t REL4:%.4f"%(OBJ, REL1, REL2, REL3, REL4)
		
		
		#Begining M.U.R.
		for it in range(1,self.max_iter+1):
			
			R1T = np.transpose(R1)
			R2T = np.transpose(R2)
			R3T = np.transpose(R3)
			R4T = np.transpose(R4)
			PT = np.transpose(P)
			ST = np.transpose(S)
			GT = np.transpose(G)
			
			#update rule for G
			R1T_P_S = np.matmul(np.matmul(R1T, P), S)
			G_GT_R1T_P_S = np.matmul(G, np.matmul(GT, R1T_P_S))
			
			R2T_G_U = np.matmul(np.matmul(R2T, G), U)
			G_GT_R2T_G_U = np.matmul(G, np.matmul(GT, R2T_G_U))
			
			R3T_G_V = np.matmul(np.matmul(R3T, G), V)
			G_GT_R3T_G_V = np.matmul(G, np.matmul(GT, R3T_G_V))
			
			R4T_G_W = np.matmul(np.matmul(R4T, G), W)
			G_GT_R4T_G_W = np.matmul(G, np.matmul(GT, R4T_G_W))
			
			G_mult = np.sqrt( np.divide(R1T_P_S + R2T_G_U + R3T_G_V + R4T_G_W, G_GT_R1T_P_S + G_GT_R2T_G_U + G_GT_R3T_G_V + G_GT_R4T_G_W+ 1e-4) )
			
			
			#update rule for S
			PT_R1_G = np.matmul(PT, np.matmul(R1, G))
			PT_P = np.matmul(PT, P)
			GT_G = np.matmul(GT, G)
			PT_P_S_GT_G = np.matmul( PT_P, np.matmul(S, GT_G))
			
			S_mult = np.sqrt( np.divide(PT_R1_G,PT_P_S_GT_G + 1e-4))
			
			#update rule for U
			GT_R2_G = np.matmul(GT, np.matmul(R2, G))
			GT_G_U_GT_G = np.matmul( GT_G, np.matmul(U, GT_G))
			
			U_mult = np.sqrt( np.divide(GT_R2_G,GT_G_U_GT_G + 1e-4))
			
			#update rule for V
			GT_R3_G = np.matmul(GT, np.matmul(R3, G))
			GT_G_V_GT_G = np.matmul( GT_G, np.matmul(V, GT_G))
			
			V_mult = np.sqrt( np.divide(GT_R3_G,GT_G_V_GT_G + 1e-4))
			
			#update rule for W
			GT_R4_G = np.matmul(GT, np.matmul(R4, G))
			GT_G_W_GT_G = np.matmul( GT_G, np.matmul(W, GT_G))
			
			W_mult = np.sqrt( np.divide(GT_R4_G,GT_G_W_GT_G + 1e-4))
			
			# Applying M.U.R.
			G = np.multiply(G, G_mult) + 1e-4
			S = np.multiply(S, S_mult) + 1e-4
			U = np.multiply(U, U_mult) + 1e-4
			V = np.multiply(V, V_mult) + 1e-4
			W = np.multiply(W, W_mult) + 1e-4
			
			if (it%self.verbose == 0) or (it==1):
				OBJ, REL1, REL2, REL3, REL4 = self.Score(R1, R2, R3, R4, P, S, G, U, V, W, norm_R1, norm_R2, norm_R3, norm_R4)
				print " - It %i:\t OBJ:%.4f\t REL1:%.4f\t REL2:%.4f\t REL3:%.4f\t REL4:%.4f"%(it, OBJ, REL1, REL2, REL3, REL4)
		return P, S, G, U, V, W



# Test code
#R = np.random.rand(20,8000)
#L = np.zeros((8000,8000))
#
#Solver = NMTF(max_iter=100, verbose = 10)
#P, S, G = Solver.Solve(R, L, 4, 50, 0.1)
#
#print P
#print S
#print G

