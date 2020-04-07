#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Dec 15 15:11:32 2019

@author: fangsa
"""

import pandas as pd
import numpy as np
import networkx as nx
import os

from networkx.algorithms.dag import descendants

def Creatdir(path):
    if not os.path.exists(path):
        os.mkdir(path)

def inter(path):
    #TS
    #im = pd.read_csv(path, index_col = 0)
    im = pd.read_table(path, header=None)
    return im


def BuildEstGraph(k, dir = "MDSINE", count_type = "RA"):
    try:
        if(dir in ["ANM", "CDS", "IGCI"]): #"ANM", "BivariateFit", "CDS", "RECI", "eLSA", "TE"
            am = pd.read_csv(dir + "\\" + count_type + "_adjacent_matrix_" + str(k) + ".csv" , index_col = 0)
            p = pd.read_csv(dir + "\\" + count_type + "_p_" + str(k) + ".csv" , index_col = 0)
        elif(dir in ["SparCC"]):
            try:
                am = pd.read_table(dir + "\\" + str(k) + "_sparcc.txt", index_col = 0)
                p = pd.read_table(dir + "\\"  + str(k) + "_pvals.two_sided.txt", index_col = 0)
            except:
                try:
                    am = pd.read_table(dir + "\\RA_am_" + str(k) + ".txt", index_col = 0)
                    p = pd.read_table(dir + "\\p_"  + str(k) + ".txt", index_col = 0)
                except:
                    am = pd.read_table(dir + "\\" + str(k) + "_.txt", index_col = 0)
                    p = pd.read_table(dir + "\\"  + str(k) + "_pvals.two_sided.txt", index_col = 0)
        elif(dir in ["MDSINE"]):
            am = pd.read_csv(dir + "\\" + count_type + "_adjacent_matrix_" + str(k) + ".csv", index_col = 0)
            p = pd.read_csv(dir + "\\" + count_type + "_p_" + str(k) + ".csv", index_col = 0)
        elif(dir in ["beem"]):
            am = pd.read_csv(dir + "\\" +  count_type + "_beta_" + str(k) + ".csv", index_col = 0)
            p = "None"
            #p = pd.read_csv(dir + "//RA_p_" + str(k) + ".csv", index_col = 0)
        elif(dir in ["LiNGAM"]):
            am = pd.read_csv(dir + "\\" + count_type + "_adjacent_matrix_" + str(k) + ".csv", index_col = 0)
            p = "None"
        return am, p
    except:
        return "None"

def G2intervention(G):
    intervention_matrix = pd.DataFrame(np.zeros([len(G.nodes()), len(G.nodes())]), index = list(G.nodes()), columns = list(G.nodes()))
    for i in list(G.nodes()):
        for j in descendants(G, i):
            intervention_matrix.loc[j, i] = 1
    return intervention_matrix

def intervention_distenct(true_graph, estimated_graph):
    return sum((true_graph - estimated_graph) != 0)
    

for x in ["R", "SF"]:
    for y in ["100", "200", "300", "400", "500", "600", "700", "800"]:
        os.chdir(r"D:\Fangsa\causal_compare\Data-CS\CS_" + x + "_" + y)

        out_file = open("D:\\Fangsa\\causal_compare\\Data-CS\\evaluation\\intervention\\" + x + "_" + y + "_auc_am.txt", "w")

        out_file.write("\t".join(["Method", "Count_type", "index", "AUC", "AUPR"])+"\n")
        for i in ["ANM", "CDS", "IGCI", "SparCC", "beem", "LiNGAM"]:#, "MDSINE"
            s = ["AA", "RA"]
            #if(i == "beem"):
            #    s = ["RA"]
            for j in s:
                for k in range(100):

                        out_file.write("\t".join([i, j, str(k+1), str(auc), str(aupr)])+"\n")                    
        out_file.close()

#Adjacency Matrix
G_mat = np.array([[0, 1, 1, 1, 0, 1, 0, 0, 0, 0],
                  [0, 0, 0, 1, 0, 0, 1, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 1, 0, 1, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                  [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])
G_mat = pd.DataFrame(G_mat, index = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'],
                     columns = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j'])
G = nx.DiGraph(G_mat)

for i in ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j']:
    print(i)
    print(descendants(G, i))
    #ReturnsThe descendants ofsourceinG
    print(len(descendants(G, i)))

nx.draw_networkx(G)

def fetch_parent_nodes(G, nodes):
    #寻找一个节点的直接相连的子代
    neighbor = set(G.neighbors(nodes))
    descendant = descendants(G, nodes)
    return neighbor.intersection(descendant)

def fetch_parent_step(G, nodes, step=2):
    #在食物链上一定步长上的子代
    parents =  set(fetch_parent_nodes(G, nodes))#保留最终的结果
    parents_layer = set(fetch_parent_nodes(G, nodes)) #保留每一层的结果
    step -= 1
    while step:
        parent_temp = set() #保留下次轮寻子代的结果
        for i in parents_layer:
            parent_temp = parent_temp | fetch_parent_nodes(G, i)   
        parents_layer = parent_temp
        parents = parents | parent_temp
        step -= 1
    return parents
    
