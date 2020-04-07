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
                #source = int(line[1].strip("Strain"))-1; target = int(line[2].strip("Strain"))-1
                #source = int(line[1].strip("sp"))-1; target = int(line[2].strip("sp"))-1
                #value = float(line[3]); pvalue = float(line[4])
                print(line)
                source = int(eval(line[2]).strip("sp"))-1; target = int(eval(line[3]).strip("sp"))-1
                value = float(line[4]); pvalue = float(line[5])
                am.ix[target, source] = value
                p.ix[target, source] = pvalue
    file.close()
    return am, p

if __name__=="__main__":
    for i in ["R", "SF"]:#, "hubbell", "soi"
        #os.chdir("C:\\Users\\Administrator\\Desktop\\causal_compare\\data_e\\" + i + "_TS_dense\\MDSINE")
        for k in ["S10", "S20", "S40"]:#"S5",
            #    continue
            os.chdir(r"C:\Users\Administrator\Desktop\causal_compare\beem_e\TS_" + i + "_" + k+"\\output")
            for j in range(800):
            #if(i == "SF" and k == "S20"):
                # path = str(j + 1) + "\\output_AA\\BVS.results.parameters.txt"
                # path = str(j + 1) + "\\output_RA\\BVS.results.parameters.txt"
                path = "paramBEEM_" + str(j + 1) + ".csv"
                try:
                    am, p = GetPred(path)
                    am.to_csv("beem\\RA_adjacent_matrix_" + str(j + 1) + ".csv", index=True, header=True)
                    p.to_csv("beem\\RA_p_" + str(j + 1) + ".csv", index=True, header=True)

                    # am, p = GetPred(path)
                    # am.to_csv("RA_adjacent_matrix_" + str(j+1) + ".csv", index=True, header=True)
                    # p.to_csv("RA_p_" + str(j +1) + ".csv", index=True, header=True)
                except:
                    print(j)
                    continue
