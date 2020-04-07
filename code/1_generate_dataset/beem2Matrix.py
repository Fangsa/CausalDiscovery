import pandas as pd
import numpy as np
import os

def GetPred(path):
    file = open(path)
    label = ["sp" + str(i) for i in range(10)]
    am = pd.DataFrame(np.zeros([10,10]), columns=label, index=label)
    p = pd.DataFrame(np.ones([10,10]), columns=label, index=label)
    for line in file:
        if(True):
            line = line.strip().split(",")
            if((line[2] == line[3]) or (line[1] == '"growth_rate"') or ( line[1] == '"parameter_type"')):
                pass
            else:
                #print(line)
                source = int(eval(line[2]).strip("sp"))-1; target = int(eval(line[3]).strip("sp"))-1
                value = float(line[4]); pvalue = float(line[5])
                am.ix[target, source] = value
                p.ix[target, source] = pvalue
    file.close()
    return am, p

if __name__=="__main__":
    for i in ["R"]:#, "hubbell", "soi"
        for k in ["S10", "S20"]:#"S5",
            os.chdir( "D:\\Fangsa\\TS_" + i + "_" + k+ "\\output")
            for j in range(400):
                path = "paramBEEM_" + str(j + 1) + ".csv"
                try:
                    am, p = GetPred(path)
                    path = "//scratch//ose-taojc//Fang//TS_" + i + "_" + k+ "//beem"
                    if(not os.path.exists(path)):
                        os.mkdir(path)
                    am.to_csv(path+"//RA_adjacent_matrix_" + str(j + 1) + ".csv", index=True, header=True)
                    p.to_csv(path+"//RA_p_" + str(j + 1) + ".csv", index=True, header=True)
                except:
                    print(j)
                    continue
