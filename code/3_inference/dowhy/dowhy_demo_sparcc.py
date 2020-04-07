import os, sys

sys.path.append(r"C:\Users\Administrator\Downloads\dowhy-master")

import numpy as np
import pandas as pd
import networkx as nx
import random as rd
from sklearn.utils import shuffle
import multiprocessing as mp
import time

import dowhy
from dowhy.do_why import CausalModel

'''
class Graph_Matrix:
    """
    Adjacency Matrix
    """
    def __init__(self, vertices=[], matrix=[]):
        """
        :param vertices:a dict with vertex id and index of matrix , such as {vertex:index}
        :param matrix: a matrix
        """
        self.matrix = matrix
        self.vertices = nodes
        self.edges_array = []  # (tail, head, weight)
        self.G = nx.Graph()

        self.edge_list = self.getAllEdges()
        self.G = self.add_edges_from_list(self.edge_list)

    def add_edges_from_list(self, edges_list):  # edges_list : [(source, target, weight),()]
        self.G.add_edges_from(edges_list)
        return self.G

    def getAllEdges(self):
        for i in range(len(self.matrix)):
            for j in range(len(self.matrix)):
                if 0 < self.matrix[i][j] < float('inf'):
                    self.edges_array.append((self.vertices[j],self.vertices[i],))
        return self.edges_array

def create_directed_matrix(am, p=None):
    nodes = ['sp' + str(i) for i in range(len(am))]
    if(p==None):
        matrix = np.array(am)
    else:
        matrix = np.array(am)*np.array(p<=0.05)
    my_graph = Graph_Matrix(nodes, matrix)
    return my_graph.G
'''


def CreatPesudoAM(TrueAM, proba=0.5):
    PesudoAM = np.matrix(TrueAM) + np.diag([0.5] * len(TrueAM))
    TrueAM = np.matrix(TrueAM) + np.diag([0.5] * len(TrueAM))
    index_positive = np.argwhere(TrueAM != 0)
    index_negative = np.argwhere((TrueAM + np.diag([1] * len(TrueAM))) == 0)
    for x, y in rd.sample(list(index_positive), int(len(index_positive) * proba)):
        loc = rd.choice(index_negative)
        PesudoAM[x, y], PesudoAM[loc[0], loc[1]] = 0, PesudoAM[x, y]
        # PesudoAM[loc[0], loc[1]] = 1
        index_negative = np.delete(index_negative, loc, axis=0)
    return PesudoAM


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
    '''
    common_causes = []; instruments = []
    for i in list(am.index):
        if(i == treatment or i == outcome):
            pass
        else:
            if(am.ix[treatment,i] != 0 and am.ix[outcome, i] != 0):
                common_causes.append(i)

    current_node = [treatment, outcome]
    current_node.extend(common_causes)
    for i in list(am.index):
        if((i not in current_node) and am.ix[treatment, i]!=0):
            current_node.remove(treatment)
            for j in current_node:
                if(am.ix[i, j] == 0 and am.ix[i,j] == 0):
                    instruments.append(i)
                else:
                    pass

    gml_graph = ('graph[directed 1'
                 'node[ id "{0}" label "{0}"]'
                 'node[ id "{1}" label "{1}"]'
                 'node[ id "{2}" label "{2}"]'
                 'edge[source "{0}" target "{1}"]'
                 'edge[source "{2}" target "{0}"]'
                 'edge[source "{2}" target "{1}"]'
                 ).format(treatment, outcome, "Unobserved Confounders")

    gml_graph = gml_graph + " ".join(['node[ id "{0}" label "{0}"] edge[ source "{0}" target "{1}"]'.format(v, treatment) for v in common_causes])
    gml_graph = gml_graph + " ".join(['edge[ source "{0}" target "{1}"]'.format(v, outcome) for v in common_causes])
    gml_graph = gml_graph + " ".join(['node[ id "{0}" label "{0}"] edge[ source "{0}" target "{1}"]'.format(v, treatment) for v in instruments])
    gml_graph = gml_graph + ']'
    '''
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
                                            #test_significance=False)
    return causal_estimate
    ##Refute
    # Adding a random common cause variable
    # res_random=model.refute_estimate(identified_estimand, causal_estimate, method_name="random_common_cause")

    # Replacing treatment with a random (placebo) variable
    # res_placebo=model.refute_estimate(identified_estimand, causal_estimate,
    #    method_name="placebo_treatment_refuter", placebo_type="permute")

    # Removing a random subset of the data
    # res_subset=model.refute_estimate(identified_estimand, causal_estimate,
    #    method_name="data_subset_refuter", subset_fraction=0.9)

#from tqdm import tqdm
def CalPvalue(dat, value):
    nullmodel = []
    c = 0
    while(True):
        try:
            dat["df"] = shuffle(dat["df"])
            causal_estimate = CalDoWhy(dat)
            nullmodel.append(causal_estimate.value)
        except ValueError:
            c -= 1
        finally:
            c += 1
            if c == 100:
                break
    '''            
    for i in range():
        dat["df"] = shuffle(dat["df"])
        causal_estimate = CalDoWhy(dat)
        nullmodel.append(causal_estimate.value)
        '''
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
                    causal_estimate = CalDoWhy(dat)
                    value.ix[j, i] = causal_estimate.value
                except ValueError:
                    value.ix[j, i] = 0
                    continue
                print(causal_estimate)
                pvalue = CalPvalue(dat, causal_estimate.value)
                print(pvalue)
                p.ix[j, i] = pvalue
                print("***************************"+str(time.time()-start_time)+"***********")
    return value, p


def StandardOutput(am, method, k, path, data_type="RA", matrix_type="am"):
    os.chdir(path)
    print(os.getcwd())
    if (not os.path.exists(method+"_sparcc1")):
        os.makedirs(method+"_sparcc1")
    am = pd.DataFrame(am)
    am.columns = ["sp" + str(i) for i in range(am.shape[0])]
    am.index = ["sp" + str(i) for i in range(am.shape[0])]
    if (matrix_type == "am"):
        am.to_csv(method + "_sparcc1\\" + data_type + "_adjacent_matrix_" + str(k) + ".csv")
    elif (matrix_type == "p"):
        am.to_csv(method + "_sparcc1\\" + data_type + "_p_" + str(k) + ".csv")


def func(par):
    start_time = time.time()
    k = par["k"]
    #am = par["am"]
    path = par["path"]
    os.chdir(path)

    p_aa = pd.read_csv(
        "C:\\Users\\Administrator\\Desktop\\causal_compare\\data\\" + par["i"] + "_" + par["j"] + "\\DoWhy_sparcc\\AA_p_" + str(
            k) + ".csv", index_col=0)
    am_aa = pd.DataFrame(np.zeros([p_aa.shape[0], p_aa.shape[0]]))
    am_aa.columns = ["sp" + str(i) for i in range(p_aa.shape[0])]
    am_aa.index = ["sp" + str(i) for i in range(p_aa.shape[0])]
    for s in range(p_aa.shape[0]):
        for t in range(p_aa.shape[1]):
            if(p_aa.ix[s, t] <= 0.05):
                am_aa.ix[s, t] = 1

    p_ra = pd.read_csv(
        "C:\\Users\\Administrator\\Desktop\\causal_compare\\data\\" + par["i"] + "_" + par["j"] + "\\DoWhy_sparcc\\RA_p_" + str(
            k) + ".csv", index_col=0)
    am_ra = pd.DataFrame(np.zeros([p_ra.shape[0], p_ra.shape[0]]))
    am_ra.columns = ["sp" + str(i) for i in range(p_ra.shape[0])]
    am_ra.index = ["sp" + str(i) for i in range(p_ra.shape[0])]
    for s in range(p_ra.shape[0]):
        for t in range(p_ra.shape[1]):
            if(p_ra.ix[s, t] <= 0.05):
                am_ra.ix[s, t] = 1

    aa = pd.read_table("AA_" + str(k) + ".txt", index_col=0).T
    ra = pd.read_table("RA_" + str(k) + ".txt", index_col=0).T

    # DoWhy
    am_aa, p_aa = main(aa, am_aa)
    print(am_aa)
    print(p_aa)
    am_ra, p_ra = main(ra, am_ra)
    StandardOutput(am_aa, method="DoWhy", k=k, data_type="AA", matrix_type="am", path=path)
    StandardOutput(am_ra, method="DoWhy", k=k, data_type="RA", matrix_type="am", path=path)
    StandardOutput(p_aa, method="DoWhy", k=k, data_type="AA", matrix_type="p", path=path)
    StandardOutput(p_ra, method="DoWhy", k=k, data_type="RA", matrix_type="p", path=path)
    print("running time is "+str(time.time()-start_time))

if __name__ == "__main__":
    pool = mp.Pool(processes=30)
    for k in range(1, 11):
        for i in ["gLV", "hubbell", "soi"]:#
            for j in ["CS", "TS"]:#
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