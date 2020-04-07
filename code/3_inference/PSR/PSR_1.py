import os, sys

sys.path.append(r"C:\Users\Administrator\Downloads\dowhy-master")

import numpy as np
import pandas as pd
import networkx as nx
import random as rd
from sklearn.utils import shuffle
from sklearn.linear_model import LinearRegression
import multiprocessing as mp
import time
from sklearn.metrics import roc_auc_score

import dowhy
from dowhy.do_why import CausalModel


def getAllEdges(am, vertices):
    edges_array = []
    for i in range(len(am)):
        for j in range(len(am)):
            if am.ix[i, j] != 0:
                edges_array.append((vertices[j], vertices[i],))
    return edges_array

def CreatGML(am, vertices):
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
    gml_graph = CreatGML(am, vertices)

    ret_dict = {
        "df": dat,
        "treatment_name": treatment,
        "outcome_name": outcome,
        "gml_graph": gml_graph,
    }
    return ret_dict


def CalPSR(dat):
    model = CausalModel(
        data=dat["df"],
        treatment=dat["treatment_name"],
        outcome=dat["outcome_name"],
        graph=dat["gml_graph"]
    )

    treatment_name = model._treatment
    outcome_name = model._outcome
    common_causes_name = model._graph.get_common_causes(treatment_name, outcome_name)

    data = dat["df"]
    treatment = data[treatment_name]
    outcome = data[outcome_name]
    if("U" in common_causes_name):
        common_causes_name.remove("U")
    common_causes = data[common_causes_name]
    reg_ps = LinearRegression().fit(common_causes, treatment)

    ps = reg_ps.predict(common_causes)

    X = pd.DataFrame({"Treatment": treatment, "PS": ps})

    psr = LinearRegression().fit(X, outcome)

    return psr.coef_[0]


def CalPvalue(dat, value):
    nullmodel = []
    for i in range(1000):
        dat["df"] = shuffle(dat["df"])
        nullmodel.append(CalPSR(dat))
    p_value = sum([1 if nullmodel[s] > abs(value) else 0 for s in range(len(nullmodel))]) / len(nullmodel)
    print(p_value)
    return p_value

def main(count, am):
    label = ["sp" + str(i) for i in range(len(am))]
    value = pd.DataFrame(np.zeros([am.shape[0], am.shape[0]]))
    p = pd.DataFrame(np.ones([am.shape[0], am.shape[0]]))
    for i in range(am.shape[0]):
        for j in range(am.shape[0]):
            if (i != j):
                start_time = time.time()
                treatment = label[i]
                outcome = label[j]
                dat = CreatData(count, am, treatment, outcome)
                try:
                    causal_estimate = CalPSR(dat)
                    value.ix[j, i] = causal_estimate
                except ValueError:
                    value.ix[j, i] = 0
                    continue
                print(causal_estimate)
                pvalue = CalPvalue(dat, causal_estimate)
                print(pvalue)
                p.ix[j, i] = pvalue
                print("***************************"+str(time.time()-start_time)+"***********")
    return value, p


def StandardOutput(am, t, s, method, k, path, data_type="RA", matrix_type="am"):
    os.chdir(path +"\\PSR")
    print(os.getcwd())
    if (not os.path.exists(method + "_" + t + str(s))):
        os.makedirs(method + "_" + t + str(s))
    am = pd.DataFrame(am)
    am.columns = ["sp" + str(i) for i in range(am.shape[0])]
    am.index = ["sp" + str(i) for i in range(am.shape[0])]
    if (matrix_type == "am"):
        am.to_csv(method + "_" + t + str(s) + "\\" + data_type + "_adjacent_matrix_" + str(k) + ".csv")
    elif (matrix_type == "p"):
        am.to_csv(method + "_" + t + str(s) +"\\" + data_type + "_p_" + str(k) + ".csv")

def BuildAM(par, t, s, data_type):
    k = par["k"]
    path = par["path"]
    os.chdir(path)
    p = pd.read_csv(
        "C:\\Users\\Administrator\\Desktop\\causal_compare\\data\\" + par["i"] + "_" + par["j"] + "\\PSR\\PSR_" + t + str(s)+"\\" + data_type + "_p_" + str(k) + ".csv", index_col=0)
    am = pd.DataFrame(np.zeros([p.shape[0], p.shape[0]]))
    am.columns = ["sp" + str(i) for i in range(p.shape[0])]
    am.index = ["sp" + str(i) for i in range(p.shape[0])]
    for i in range(p.shape[0]):
        for j in range(p.shape[1]):
            if(p.ix[i, j] <= 0.05):
                am.ix[i, j] = 1
    return am

def matrix2array(m):
    a = []
    for i c range(m.shape[0]):
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

    for t in ["mmhc", "cg", "half", "sparcc"]:
        for data_type in ["AA", "RA"]:
            if(t == "sparcc" and data_type == "AA"):
                break
            step = 2
            while(step <= 50):
                start_time = time.time()
                am_prior = BuildAM(par, t, step-1, data_type)
                if(data_type == "AA"):
                    am, p = main(aa, am_prior)
                elif(data_type == "RA"):
                    am, p = main(ra, am_prior)
                StandardOutput(am, t, step, method="PSR", k=k, data_type=data_type, matrix_type="am", path=path)
                StandardOutput(p, t, step, method="PSR", k=k, data_type=data_type, matrix_type="p", path=path)
                print("running time is "+str(time.time()-start_time))
                step += 1
                try:
                    th = threshold(am_prior, p)
                    if(th > 0.9):
                        print(th)
                        break
                except ValueError:
                    break

if __name__ == "__main__":
    pool = mp.Pool(processes=30)
    for k in range(1, 11):
        for i in ["gLV", "hubbell", "soi"]:
            for j in ["CS", "TS"]:
                path = "C:\\Users\\Administrator\\Desktop\\causal_compare\\data\\" + i + "_" + j
                par = {
                    "i": i,
                    "j": j,
                    "k": k,
                    "path": path
                }
                #func(par)
                pool.apply_async(func, (par,))
    pool.close()
    pool.join()