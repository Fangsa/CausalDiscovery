def generate(path, i, j):
    for ID in range(1, 101):
        for dt in ["AA", "RA"]:
            cfg = ''.join(open(r'C:\Users\Administrator\Desktop\causal_compare\beem_e\parameters.cfg', 'r').readlines())
            cfg = cfg.replace('output_dir = output_cdiff', 'output_dir = ' + path + 'output_' + str(ID) + "_" + dt)
            cfg = cfg.replace("metadata_file = data_cdiff/metadata.txt",
                              "metadata_file = " + path + "metadata_1.txt")
            cfg = cfg.replace("counts_file = data_cdiff/counts.txt",
                              "counts_file = " + path + dt + "_" + str(ID) + ".txt")
            cfg = cfg.replace("biomass_file = data_cdiff/biomass.txt",
                              "biomass_file = " + path + "biomass_" + str( ID) + ".txt")
            print(ID)
            if (dt == "AA"):
                cfg = cfg.replace("normalize_counts = 0", "normalize_counts = 1")
            f = open("D:\Fangsa\causal_compare\Data\TS_" + i + "_" + j + '\MDSINE\parameters_' + str(ID) + "_" + dt + ".cfg", "w")
            f.write(cfg)
            f.close()

for i in ["R"]:#, "SF"
    for j in ["S10", "S20"]:#, "S30", "S40", "S50", "S60"
        path = "/scratch/ose-taojc/Fang/Data/TS_" + i + "_" + j + "/MDSINE/"
        generate(path, i, j)