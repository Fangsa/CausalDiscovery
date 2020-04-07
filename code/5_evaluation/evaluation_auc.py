#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 23:24:30 2019

@author: fangsa
"""

import numpy as np
import pandas as pd
import os
from sklearn import metrics
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score

def BuildEstGraph(k, dir = "MDSINE", count_type = "RA"):
    try:
        if(dir in ["ANM", "CDS", "IGCI"]): #"ANM", "BivariateFit", "CDS", "RECI", "eLSA", "TE"
            am = pd.read_csv(dir + "//" + count_type + "_adjacent_matrix_" + str(k) + ".csv" , index_col = 0)
            p = pd.read_csv(dir + "//" + count_type + "_p_" + str(k) + ".csv" , index_col = 0)
        elif(dir in ["SparCC"]):
            try:
                am = pd.read_table(dir + "//" + str(k) + "_.txt", index_col = 0)        
            except:
                am = pd.read_table(dir + "//" + str(k) + "_sparcc.txt", index_col = 0)
            p = pd.read_table(dir + "//"  + str(k) + "_pvals.two_sided.txt", index_col = 0)
        elif(dir in ["MDSINE"]):
            am = pd.read_csv(dir + "//" + count_type + "_adjacent_matrix_" + str(k) + ".csv", index_col = 0)
            p = pd.read_csv(dir + "//" + count_type + "_p_" + str(k) + ".csv", index_col = 0)
        elif(dir in ["beem"]):
            am = pd.read_csv(dir + "//RA_adjacent_matrix_" + str(k) + ".csv", index_col = 0)
            p = pd.read_csv(dir + "//RA_p_" + str(k) + ".csv", index_col = 0)
        elif(dir in ["LiNGAM"]):
            am = pd.read_csv(dir + "//" + count_type + "_adjacent_matrix_" + str(k) + ".csv", index_col = 0)
            p = "None"
        return am, p
    except:
        return "None"

def CalAUC(k, pre_matrix):
    #pre_matrix = pd.read_csv("beem/RA_adjacent_matrix_"+k+".csv", index_col=0)
    true_matrix = pd.read_csv("beta_"+str(k)+".csv", index_col=0)
    pre_matrix.fillna(0, inplace=True)
    pre_list = [abs(pre_matrix.iloc[i, j]) for i in range(pre_matrix.shape[0]) for j in range(pre_matrix.shape[1]) if i!=j]
    true_list = [true_matrix.iloc[i, j] for i in range(true_matrix.shape[0]) for j in range(true_matrix.shape[1]) if i!=j]
    true_list = [1 if i!=0 else 0 for i in true_list]
    #print(pre_list)
    #print(true_list)
    return(roc_auc_score(true_list, pre_list), average_precision_score(true_list, pre_list))
   

for x in ["R", "SF"]:
    for y in ["S10", "S20", "S30", "S40", "S50", "S60"]:
        os.chdir(r"/Users/fangsa/compare_causal/causal/Data/TS_" + x + "_" + y)

        out_file = open("/Users/fangsa/compare_causal/causal/Fig-TS/evaluation-TS/" + x + "_" + y + "_auc_am.txt", "w")
        out_file1 = open("/Users/fangsa/compare_causal/causal/Fig-TS/evaluation-TS/" + x + "_" + y + "_auc_p.txt", "w")
        
        out_file.write("\t".join(["Method", "Count_type", "index", "AUC", "AUPR"])+"\n")
        out_file1.write("\t".join(["Method", "Count_type", "index", "AUC", "AUPR"])+"\n")
        for i in ["ANM", "CDS", "IGCI", "SparCC", "MDSINE", "beem", "LiNGAM"]:#
            s = ["AA", "RA"]
            if(i == "beem"):
                s = ["RA"]
            for j in s:
                for k in range(100):
                    res = BuildEstGraph(str(k+1), i, j)
                    if(res == "None"):
                        continue
                    else:
                        am, p = res
                        if(type(p) == str):                   
                            continue
                        else:
                            auc, aupr = CalAUC(str(k+1), p)
                            out_file1.write("\t".join([i, j, str(k+1), str(auc), str(aupr)])+"\n")  
                        auc, aupr = CalAUC(str(k+1), am)
                        out_file.write("\t".join([i, j, str(k+1), str(auc), str(aupr)])+"\n")                    
        out_file.close()
        out_file1.close()