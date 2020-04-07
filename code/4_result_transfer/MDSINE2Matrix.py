import pandas as pd
import numpy as np
import os

def GetPred(path):
    file = open(path)
    label = ["sp" + str(i) for i in range(10)]
    am = pd.DataFrame(np.zeros([10,10]), columns=label, index=label)
    p = pd.DataFrame(np.ones([10,10]), columns=label, index=label)
    for line in file:
        line = line.strip().split("\t")
        if(True):
            if((line[1] == line[2]) or (line[0] == "growth_rate") or ( line[0] == "parameter_type")):
                pass
            else:
                #source = int(line[1].strip("Strain"))-1; target = int(line[2].strip("Strain"))-1
                source = int(line[1].strip("sp"))-1; target = int(line[2].strip("sp"))-1
                value = float(line[3]); pvalue = float(line[4])
                am.ix[target, source] = value
                p.ix[target, source] = pvalue
    file.close()
    return am, p

if __name__=="__main__":
    for i in ["MDSINE"]:
        #os.chdir("C:\\Users\\Administrator\\Desktop\\causal_compare\\data_e\\" + i + "_TS_dense\\MDSINE")
        os.chdir(r"C:\Users\Administrator\Desktop\causal_compare\beem_e\beem_data1\MDSINE")
        for j in range(100):
            #path = "C:\\Users\\Administrator\\Desktop\\causal_compare\\beem\\beem_data3\\MDSINE\\" \
            #      "output_"+str(j + 1)+"_AA\\BVS.results.parameters.txt"
            #path = str(j + 1) + "\\output_RA\\BVS.results.parameters.txt"
            #path = "paramBEEM_"+str(j+1)+ ".csv"
            path = "output_"+str(j + 1)+"_RA\\BVS.results.parameters.txt"
            try:
                am, p = GetPred(path)
                am.to_csv("RA_adjacent_matrix_" + str(j+1) + ".csv", index=True, header=True)
                p.to_csv("RA_p_" + str(j +1) + ".csv", index=True, header=True)

                #am, p = GetPred(path)
                #am.to_csv("RA_adjacent_matrix_" + str(j+1) + ".csv", index=True, header=True)
                #p.to_csv("RA_p_" + str(j +1) + ".csv", index=True, header=True)
            except:
                print(j)
                continue
