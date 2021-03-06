import os, sys

#sys.path.append(r"C:\Users\Administrator\Downloads\dowhy-master-0514")
sys.path.append(r"C:\Users\Administrator\Downloads\dowhy-master")

import numpy as np
import pandas as pd
import networkx as nx
import random as rd
from sklearn.utils import shuffle
import multiprocessing as mp
import time
from sklearn.metrics import roc_auc_score


import dowhy
from dowhy.do_why import CausalModel

def BuildAM(par, t, data_type):
    k = par["k"]
    path = par["path"]
    os.chdir(path)
    if(t == "cg"):
        G = nx.complete_graph(10)
        am = nx.to_pandas_adjacency(G)
    elif(t == "half"):
        am = pd.read_table("C:\\Users\\Administrator\\Desktop\\causal_compare\\data_e\\inter_" + str(k) + ".txt", header=None)
        am = CreatPesudoAM(am)
    elif(t == "mmhc"):
        am = pd.read_csv("C:\\Users\\Administrator\\Desktop\\causal_compare\\beem_e\\" + par["i"] + "_" + par["j"] + "\\mmhc\\"+ data_type+ "_adjacent_matrix_" + str(k) + ".csv", index_col=0)
    elif(t == "sparcc"):
        p = pd.read_table("C:\\Users\\Administrator\\Desktop\\causal_compare\\data_e\\" + par["i"] + "_" + par["j"] + "\\SparCC\\p_" + str(k) + ".txt", index_col=0)
        am = pd.DataFrame(np.zeros([p.shape[0], p.shape[0]]))
        am.columns = ["sp" + str(i) for i in range(p.shape[0])]
        am.index = ["sp" + str(i) for i in range(p.shape[0])]
        for s in range(p.shape[0]):
            for t in range(p.shape[1]):
                if(p.ix[s, t] <= 0.05):
                    am.ix[s, t] = 1
    elif(t == "true"):
        am = pd.read_table("C:\\Users\\Administrator\\Desktop\\causal_compare\\data_e\\inter_" + str(k) + ".txt", header=None)
    label = ["sp" + str(i) for i in range(len(am))]
    am.index = label; am.columns = label
    #am = pd.DataFrame(am, index=label, columns=label)
    return am

def CreatPesudoAM(TrueAM, proba=0.5):
    PesudoAM = np.matrix(TrueAM) + np.diag([0.5] * len(TrueAM))
    TrueAM = np.matrix(TrueAM) + np.diag([0.5] * len(TrueAM))
    index_positive = np.argwhere(TrueAM != 0)
    index_negative = np.argwhere((TrueAM + np.diag([1] * len(TrueAM))) == 0)
    for x, y in rd.sample(list(index_positive), int(len(index_positive) * proba)):
        loc = rd.choice(index_negative)
        PesudoAM[x, y], PesudoAM[loc[0], loc[1]] = 0, PesudoAM[x, y]
        index_negative = np.delete(index_negative, loc, axis=0)
    return pd.DataFrame(PesudoAM)


def getAllEdges(am, vertices):
    edges_array = []
    for i in range(len(am)):
        for j in range(len(am)):
            if am.ix[i, j] != 0:
                edges_array.append((vertices[j], vertices[i],))  # edges_list : [(source, target),()]
    return edges_array


def CreatGML(am, vertices, treatment, outcome):
    edges_list = getAllEdges(am, vertices)
    gml_graph = 'graph[directed 1 '
    gml_graph = gml_graph + " ".join(['node[ id "{0}" label "{0}"]'.format(n) for n in vertices])
    gml_graph = gml_graph + " ".join(['edge[ source "{0}" target "{1}"]'.format(s, t) for s, t in edges_list])
    gml_graph = gml_graph + ']'
    return gml_graph


def CreatData(dat, am, treatment="sp0", outcome="sp1"):
    treatment = treatment
    outcome = outcome
    vertices = list(am.index)
    gml_graph = CreatGML(am, vertices, treatment, outcome)

    ret_dict = {
        "df": dat,
        "treatment_name": treatment,
        "outcome_name": outcome,
        # "common_causes_names": common_causes,
        # "instrument_names": instruments,
        "gml_graph": gml_graph,
    }
    return ret_dict

def CalDoWhy(dat):
    model = CausalModel(
        data=dat["df"],
        treatment=dat["treatment_name"],
        outcome=dat["outcome_name"],
        graph=dat["gml_graph"]
    )

    # Identification
    identified_estimand = model.identify_effect()

    # Estimation
    causal_estimate = model.estimate_effect(identified_estimand,
                                            method_name="backdoor.linear_regression")

    return causal_estimate

def CalPvalue(dat, value):
    nullmodel = []
    c = 0
    while(True):
        try:
            dat["df"] = shuffle(dat["df"])
            causal_estimate = CalDoWhy(dat)
            nullmodel.append(causal_estimate.value)
            print(causal_estimate.value)
        except ValueError:
            c -= 1
        finally:
            c += 1
            if c == 1000:
                break
    p_value = sum([1 if nullmodel[s] > abs(value) else 0 for s in range(len(nullmodel))]) / len(nullmodel)
    print(p_value)
    return p_value

def main(count, am):
    label = ["sp" + str(i) for i in range(len(am))]
    value = pd.DataFrame(np.zeros([am.shape[0], am.shape[0]]), index=label, columns=label)
    p = pd.DataFrame(np.ones([am.shape[0], am.shape[0]]), index=label, columns=label)
    for i in range(am.shape[0]):
        for j in range(am.shape[0]):
            if (i != j):
                print(i); print(j)
                start_time = time.time()
                treatment = label[i]
                outcome = label[j]
                dat = CreatData(count, am, treatment, outcome)
                try:
                    causal_estimate = CalDoWhy(dat)
                    value.ix[j, i] = causal_estimate.value
                    print(causal_estimate.value)
                except ValueError:
                    if(am.ix[j, i]):
                        value.ix[j, i]=1; p.ix[j, i]=0
                    else:
                        value.ix[j, i] = 0
                    continue
                print(causal_estimate)
                pvalue = CalPvalue(dat, causal_estimate.value)
                print(pvalue)
                p.ix[j, i] = pvalue
                print("***************************"+str(time.time()-start_time)+"***********")
    return value, p


def StandardOutput(am, t, method, k, path, data_type="RA", matrix_type="am"):
    if(not os.path.exists(path +"\\DoWhy")):
        os.mkdir(path +"\\DoWhy")
        os.chdir(path +"\\DoWhy")
    print(os.getcwd())
    if (not os.path.exists(method + "_" + t + str(1))):
        os.makedirs(method + "_" + t + str(1))
    am = pd.DataFrame(am)
    am.columns = ["sp" + str(i) for i in range(am.shape[0])]
    am.index = ["sp" + str(i) for i in range(am.shape[0])]
    if (matrix_type == "am"):
        am.to_csv(method + "_" + t + str(1) + "\\" + data_type + "_adjacent_matrix_" + str(k) + ".csv")
    elif (matrix_type == "p"):
        am.to_csv(method + "_" + t + str(1) +"\\" + data_type + "_p_" + str(k) + ".csv")

def matrix2array(m):
    a = []
    for i in range(m.shape[0]):
        for j in range(m.shape[1]):
            if(i != j):
                a.append(m.ix[i, j])
    return a

def threshold(prior_am, est_am):
    label = matrix2array(prior_am)
    prob = matrix2array(est_am)
    return roc_auc_score(label, prob)

def func(par):
    k = par["k"]
    path = par["path"]
    os.chdir(path)

    aa = pd.read_table("AA_" + str(k) + ".txt", index_col=0).T
    ra = pd.read_table("RA_" + str(k) + ".txt", index_col=0).T

    for t in ["mmhc"]:#half, "mmhc"
        for data_type in ["RA"]:#, "AA"
            start_time = time.time()
            am_prior = BuildAM(par, t, data_type)
            print(am_prior)
            if (data_type == "AA"):
                am, p = main(aa, am_prior)
            elif (data_type == "RA"):
                am, p = main(ra, am_prior)
            StandardOutput(am, t, method="DoWhy", k=k, data_type=data_type, matrix_type="am", path=path)
            StandardOutput(p, t, method="DoWhy", k=k, data_type=data_type, matrix_type="p", path=path)
            print("running time is " + str(time.time() - start_time))

if __name__ == "__main__":
    pool = mp.Pool(processes=15)
    for k in range(1, 101):
        for i in ["gLV_CS_R"]:#
            #if(i == "gLV"):
             #   s = ["TS", "TS_dense", "CS"]#
            #else:
             #   s = ["TS", "CS"]#
            s = ["100"]
            for j in s:
                path = "C:\\Users\\Administrator\\Desktop\\causal_compare\\beem_e\\" + i + "_" + j
                par = {
                    "i": i,
                    "j": j,
                    "k": k,
                    "path": path
                }
                func(par)
                #pool.apply_async(func, (par,))
    #pool.close()
    #pool.join()