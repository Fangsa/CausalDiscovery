"""
import pandas as pd
import numpy as np
import statsmodels.api as sm
from linearmodels.iv import IV2SLS
from linearmodels.datasets import wage
data = wage.load()
dependent = np.log(data.wage)
exog = sm.add_constant(data.exper)
endog = data.educ
instruments = data.sibs

mod = IV2SLS(dependent, exog, endog, instruments)
res = mod.fit(cov_type='unadjusted')
res
"""
import pandas as pd
import numpy as np
import random as rd
from sklearn.utils import shuffle
import multiprocessing as mp
import os


def CalPvalue(iv_est, treatment, outcome, instrument):
    nullmodel = []
    for i in range(1000):
        treatment = shuffle(treatment)
        outcome = shuffle(outcome)
        instrument = shuffle(instrument)

        num_yz = np.cov(outcome, instrument)[1, 0]
        deno_xz = np.cov(treatment, instrument)[1, 0]
        iv_pesudo = num_yz / deno_xz

        nullmodel.append(iv_pesudo)
    p_value = sum([1 if nullmodel[s] > abs(iv_est) else 0 for s in range(len(nullmodel))]) / len(nullmodel)
    print(p_value)
    return p_value


def main(dat):
    label = ["sp" + str(i) for i in range(dat.shape[0])]
    value = pd.DataFrame(np.zeros([dat.shape[0], dat.shape[0]]))
    p = pd.DataFrame(np.ones([dat.shape[0], dat.shape[0]]))
    for i in range(dat.shape[0]):
        for j in range(dat.shape[0]):
            if (i != j):
                treatment = dat.ix[label[i],]
                outcome = dat.ix[label[j],]
                instrument = rd.uniform(-1, 1) * treatment

                num_yz = np.cov(outcome, instrument)[1, 0]
                deno_xz = np.cov(treatment, instrument)[1, 0]
                iv_est = num_yz / deno_xz

                value.ix[j, i] = iv_est

                pvalue = CalPvalue(iv_est, treatment, outcome, instrument)
                print(pvalue)
                p.ix[j, i] = pvalue
    return value, p


def StandardOutput(am, method="IV", k=0, path=os.getcwd(), data_type="RA", matrix_type="am"):
    os.chdir(path)
    print(os.getcwd())
    if (not os.path.exists(method)):
        os.makedirs(method)
    am = pd.DataFrame(am)
    am.columns = ["sp" + str(i) for i in range(am.shape[0])]
    am.index = ["sp" + str(i) for i in range(am.shape[0])]
    if (matrix_type == "am"):
        am.to_csv(method + "\\" + data_type + "_adjacent_matrix_" + str(k) + ".csv")
    elif (matrix_type == "p"):
        am.to_csv(method + "\\" + data_type + "_p_" + str(k) + ".csv")


def func(par):
    k = par["k"]
    path = par["path"]
    os.chdir(path)
    aa = pd.read_table("AA_" + str(k) + ".txt", index_col=0)
    ra = pd.read_table("RA_" + str(k) + ".txt", index_col=0)

    # IV
    am_aa, p_aa = main(aa)
    am_ra, p_ra = main(ra)
    StandardOutput(am_aa, method="IV", k=k, data_type="AA", matrix_type="am", path=path)
    StandardOutput(am_ra, method="IV", k=k, data_type="RA", matrix_type="am", path=path)
    StandardOutput(p_aa, method="IV", k=k, data_type="AA", matrix_type="p", path=path)
    StandardOutput(p_ra, method="IV", k=k, data_type="RA", matrix_type="p", path=path)

if __name__ == "__main__":
    pool = mp.Pool(processes=30)
    for k in range(1, 11):
        for i in ["gLV", "hubbell", "soi"]:
            for j in ["CS", "TS"]:
                print("\t".join([i,j,str(k)]))
                path = "C:\\Users\\Administrator\\Desktop\\causal_compare\\data\\" + i + "_" + j
                par = {
                    "k": k,
                    "path": path
                }
                #func(par)
                pool.apply_async(func, (par,))
    pool.close()
    pool.join()