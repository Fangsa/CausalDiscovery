#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 23 10:21:22 2019

@author: fangsa
"""

import os
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import spearmanr

def Creatdir(path):
    if not os.path.exists(path):
        os.mkdir(path)

def inter(path):
    #TS
    im = pd.read_csv(path, index_col = 0)
    #im = pd.read_table(path, header=None)
    return im

def BuildEstGraph(k, dir, count_type = "AA"):
    try:
        if(dir in ["ANM", "CDS", "IGCI"]): #"ANM", "BivariateFit", "CDS", "RECI", "eLSA", "TE"
            am = pd.read_csv(dir + "//" + count_type + "_adjacent_matrix_" + str(k) + ".csv" , index_col = 0)
            p = pd.read_csv(dir + "//" + count_type + "_p_" + str(k) + ".csv" , index_col = 0)
            estGraph = am * (p <= 0.05)
        elif(dir in ["SparCC"]):
            try:
                am = pd.read_table(dir + "//" + str(k) + "_.txt", index_col = 0)        
            except:
                am = pd.read_table(dir + "//" + str(k) + "_sparcc.txt", index_col = 0)
            p = pd.read_table(dir + "//"  + str(k) + "_pvals.two_sided.txt", index_col = 0)
            estGraph = am * (p <= 0.01)
        elif(dir in ["MDSINE"]):
            am = pd.read_csv(dir + "//" + count_type + "_adjacent_matrix_" + str(k) + ".csv", index_col = 0)
            p = pd.read_csv(dir + "//" + count_type + "_p_" + str(k) + ".csv", index_col = 0)
            estGraph = am*(p>=10)
        elif(dir in ["CCDr", "GES", "GIES", "LiNGAM"]):
            am = pd.read_csv(dir + "//" + count_type + "_adjacent_matrix_" + str(k) + ".csv", index_col = 0)
            estGraph = am
        elif(dir in ["GS", "iamb", "mmhc", "mmpc", "PC"]):
            am = pd.read_csv(dir + "//" + count_type + "_adjacent_matrix_" + str(k) + ".csv", index_col = 0)
            estGraph = am
        elif(dir in ["beem"]):
            am = pd.read_csv(dir + "//RA_adjacent_matrix_" + str(k) + ".csv", index_col = 0)
            p = pd.read_csv(dir + "//RA_p_" + str(k) + ".csv", index_col = 0)
            estGraph = am*(p>=1)
        elif(dir in ["SpiecEasi_GL", "SpiecEasi_MB"]):
            am = pd.read_csv(dir + "//" + count_type + "_" + str(k) + ".csv", index_col = 0)
            estGraph = am
        return estGraph
    except:
        return "None"


def topology(G):
    ##degree
    G_degree = G.degree()
    ##clustering
    G_clustering = nx.clustering(G)
    ##average_clustering
    G_average_clustering = nx.average_clustering(G)
    ##diameter
    #G_diameter = nx.diameter(G)
    ####average_shortest_path_length 
    #G_average_shortest_path_length = nx.average_shortest_path_length(G)
    ##degree_centrality (The degree centrality for a node v is the fraction of nodes it is connected to)
    #G_degree_centrality = nx.degree_centrality(G)
    ##closeness_centrality (average shortest path distance to u over all n-1 reachable nodes)
    G_closeness_centrality = nx.closeness_centrality(G)
    ##betweenness_centrality (Betweenness centrality of a node v is the sum of the fraction of all-pairs shortest paths that pass through v)
    G_betweenness_centrality = nx.betweenness_centrality(G)
    ##current_flow_closeness_centrality
    #G_current_flow_closeness_centrality = nx.current_flow_closeness_centrality(G)
    ##current_flow_betweenness_centrality
    #G_current_flow_betweenness_centrality = nx.current_flow_betweenness_centrality(G)
    ##eigenvector_centrality
    #G_eigenvector_centrality = nx.eigenvector_centrality(G)
    #G_pagerank = nx.pagerank(G)
    G_pagerank = nx.pagerank_numpy(G)
    return dict(G_degree), G_clustering, G_closeness_centrality, G_betweenness_centrality, G_pagerank

def mre(true_value, est_value):
    true_value = np.array(true_value)
    est_value = np.array(est_value)
    #return np.mean(abs((true_value - est_value)/true_value))
    return np.mean(abs(true_value - est_value))

    
def rank_error(true_value, est_value):
    true_value = pd.Series(true_value)
    est_value = pd.Series(est_value)
    true_rank = true_value.rank()
    est_rank = est_value.rank()
    return np.mean(abs(true_rank - est_rank))/len(true_rank)

def Spearman(true_value, est_value):
    true_value = pd.Series(true_value)
    est_value = pd.Series(est_value)
    return spearmanr(true_value, est_value)

for s in ["R", "SF"]:#
    for t in ["S10", "S20", "S30", "S40", "S50", "S60"]:#
       os.chdir(r"/Users/fangsa/compare_causal/causal/Data/TS_"+ s +"_"+ t)
       out_file = open("/Users/fangsa/compare_causal/causal/Fig-TS/evaluation-TS/mre_"+ s +"_"+ t + ".txt", "w")
       out_file.write("\t".join(["Method", "Count_type", "Index", "Degree", "Clustering", "Closeness_Centrality", "Betweenness_Centrality", "PageRank"]) + "\n")

       out_file1 = open("/Users/fangsa/compare_causal/causal/Fig-TS/evaluation-TS/dist_"+ s +"_"+ t + "_r.txt", "w")
       out_file1.write("\t".join(["Method", "Count_type", "Index", "Degree", "Clustering", "Closeness_Centrality", "Betweenness_Centrality", "PageRank"]) + "\n")

       out_file2 = open("/Users/fangsa/compare_causal/causal/Fig-TS/evaluation-TS/dist_"+ s +"_"+ t + "_p.txt", "w")
       out_file2.write("\t".join(["Method", "Count_type", "Index", "Degree", "Clustering", "Closeness_Centrality", "Betweenness_Centrality", "PageRank"]) + "\n")

       for k in range(1, 101):
           true_adj = inter("beta_" + str(k) + ".csv")
           #true_adj = inter("inter_" + str(k) + ".txt")
           np.fill_diagonal(true_adj.values, 0)
           #true_adj.index = ["sp" + str(i) for i in range(10)]
           #true_adj.columns = ["sp" + str(i) for i in range(10)]
           true_graph = nx.from_numpy_matrix(np.matrix(true_adj))
           true_degree, true_clustering, true_closeness_centrality, true_betweenness_centrality, true_pagerank = topology(true_graph)
           true_degree_list = [true_degree[l] if(l in true_degree.keys()) else 0 for l in range(10)]
           true_clustering_list = [true_clustering[l] if(l in true_clustering.keys()) else 0 for l in range(10)]
           true_closeness_centrality_list = [true_closeness_centrality[l] if(l in true_closeness_centrality.keys()) else 0 for l in range(10)]
           true_betweenness_centrality_list = [true_betweenness_centrality[l] if(l in true_betweenness_centrality.keys()) else 0 for l in range(10)]
           true_pagerank_list = [true_pagerank[l] if(l in true_pagerank.keys()) else 0 for l in range(10)]
           df = pd.DataFrame({"Index" : range(10),
                       "Degree" : true_degree_list,
                       "Clustering" : true_clustering_list,
                       "Closeness_Centrality" : true_closeness_centrality_list,
                       "Betweenness_Centrality" : true_betweenness_centrality_list,
                       "PageRank": true_pagerank_list})
           Creatdir("topology")
           df.to_csv("topology//top_" + str(k) + ".txt", sep = "\t")
    #for key in true_degree.keys():
        #out_file.write("\t".join([str(key), str(true_degree[key]), str(true_clustering[key]),\
                                  #str(true_closeness_centrality[key]), str(true_betweenness_centrality[key])])+"\n")
    #out_file.close()
           for i in ["SparCC", "SpiecEasi_GL", "SpiecEasi_MB", 
              "GS", "PC", "iamb", "mmpc", "CCDr", "GES", "GIES", "mmhc",
              "beem", "MDSINE", "ANM", "LiNGAM", "CDS", "IGCI"]:#"TE",
               Creatdir(i + "//topology")
               for j in ["AA", "RA"]:
                   if(i == "beem" and j == "AA"):
                       continue
                   if(i == "SparCC" and j == "AA"):
                       continue
                   est_adj = BuildEstGraph(k, i, j)
                   if(isinstance(est_adj, str)):
                       continue
                   np.fill_diagonal(est_adj.values, 0)
                   print("\t".join([i, j, str(k)]))
                   #est_adj.index = ["sp" + str(i) for i in range(10)]
                   #est_adj.columns = ["sp" + str(i) for i in range(10)]
                   est_graph = nx.from_numpy_matrix(np.matrix(est_adj))
                   est_degree, est_clustering, est_closeness_centrality, est_betweenness_centrality, est_pagerank = topology(est_graph)
                   est_degree_list = [est_degree[l] if(l in est_degree.keys()) else 0 for l in range(10)]
                   est_clustering_list = [est_clustering[l] if(l in est_clustering.keys()) else 0 for l in range(10)]
                   est_closeness_centrality_list = [est_closeness_centrality[l] if(l in est_closeness_centrality.keys()) else 0 for l in range(10)]
                   est_betweenness_centrality_list = [est_betweenness_centrality[l] if(l in est_betweenness_centrality.keys()) else 0 for l in range(10)]
                   est_pagerank_list = [est_pagerank[l] if(l in est_pagerank.keys()) else 0 for l in range(10)]
                   df = pd.DataFrame({"Index" : range(10),
                              "Degree" : est_degree_list,
                              "Clustering" : est_clustering_list,
                              "Closeness_Centrality" : est_closeness_centrality_list,
                              "Betweenness_Centrality" : est_betweenness_centrality_list,
                              "PageRank": est_pagerank_list})
                   df.to_csv(i + "//topology//topo_" + j + "_" + str(k) + ".txt", sep = "\t")
                   
                   mre_d = mre(true_degree_list, est_degree_list)
                   mre_c = mre(true_clustering_list, est_clustering_list)       
                   mre_cc = mre(true_closeness_centrality_list, est_closeness_centrality_list)
                   mre_bc = mre(true_betweenness_centrality_list, est_betweenness_centrality_list)
                   mre_p = mre(true_pagerank_list, est_pagerank_list)
                   out_file.write("\t".join([i, j, str(k), str(mre_d), str(mre_c), str(mre_cc), str(mre_bc), str(mre_p)]) + "\n")
       
                   dist_d_r = Spearman(true_degree_list, est_degree_list)[0]
                   dist_c_r = Spearman(true_clustering_list, est_clustering_list)[0]      
                   dist_cc_r = Spearman(true_closeness_centrality_list, est_closeness_centrality_list)[0]
                   dist_bc_r = Spearman(true_betweenness_centrality_list, est_betweenness_centrality_list)[0]
                   dist_p_r = Spearman(true_pagerank_list, est_pagerank_list)[0]
                   out_file1.write("\t".join([i, j, str(k), str(dist_d_r), str(dist_c_r), str(dist_cc_r), str(dist_bc_r), str(dist_p_r)]) + "\n")

                   dist_d_p = Spearman(true_degree_list, est_degree_list)[1]
                   dist_c_p = Spearman(true_clustering_list, est_clustering_list)[1]   
                   dist_cc_p = Spearman(true_closeness_centrality_list, est_closeness_centrality_list)[1]
                   dist_bc_p = Spearman(true_betweenness_centrality_list, est_betweenness_centrality_list)[1]
                   dist_p_p = Spearman(true_pagerank_list, est_pagerank_list)[1]
                   out_file2.write("\t".join([i, j, str(k), str(dist_d_p), str(dist_c_p), str(dist_cc_p), str(dist_bc_p), str(dist_p_p)]) + "\n")

       out_file.close()
       out_file1.close() 
       out_file2.close() 


            


            
            
