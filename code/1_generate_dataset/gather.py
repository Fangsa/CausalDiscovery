# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 14:17:23 2019

@author: fangsa
"""

import os, shutil

for i in ["SF"]:
    for j in ['S10', "S20", 'S30', 'S40', 'S50', 'S60']:
        path = "//scratch//ose-taojc//Fang//TS_"+i+"_"+j 
        os.chdir(path+"//beem")
        if(not os.path.exists(path+"//filtered")):
            os.mkdir(path+"//filtered")
            os.mkdir(path+"//filtered//beem")
        files = os.listdir(os.getcwd())
        files = [k for k in files if os.path.isfile(k)]
        c = 1
        for f in files:
            if(f.startswith("RA_p_")):
                s = f.strip(".csv").split("_")[2]
                shutil.copy(path+"//AA_"+s+".txt", path+"//filtered//AA_"+str(c)+".txt")
                shutil.copy(path+"//Abundance_trajectroy_"+s+".pdf", path+"//filtered//Abundance_trajectroy_"+str(c)+".pdf")
                shutil.copy(path+"//alpha_"+s+".csv", path+"//filtered//alpha_"+str(c)+".csv")
                shutil.copy(path+"//beta_"+s+".csv", path+"//filtered//beta_"+str(c)+".csv")
                shutil.copy(path+"//biomass_"+s+".txt", path+"//filtered//biomass_"+str(c)+".txt")     
                shutil.copy(path+"//metadata_"+s+".txt", path+"//filtered//metadata_"+str(c)+".txt")             
                shutil.copy(path+"//RA_"+s+".txt", path+"//filtered//RA_"+str(c)+".txt")

                shutil.copy(path+"//beem//RA_adjacent_matrix_"+s+".csv", path+"//filtered//beem//RA_adjacent_matrix_"+str(c)+".csv")             
                shutil.copy(path+"//beem//RA_p_"+s+".csv", path+"//filtered//beem//RA_p_"+str(c)+".csv")
                
                c+=1
                
            
        
        
