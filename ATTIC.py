import argparse
import pycaret
import os
import re
import numpy as np
import pandas as pd
import iFeatureOmegaCLI
import read_fasta
from pycaret.classification import *


datapath = os.path.abspath('.')

class Fly(object):
    def __init__(self, filepath):
        self.filepath = filepath

    #def check_sequence_length(self,file):


    def generate_encoding(self):
        # generate features
        RNA = iFeatureOmegaCLI.iRNA(self.filepath)
        RNA.import_parameters(datapath + '/iFeatureOmegaCLI/parameters/RNA_parameters_setting.json')
        RNA.get_descriptor("DPCP")
        RNA.to_csv("./data/D/DPCP.csv", "index=False", header=True)
        RNA.get_descriptor("CKSNAP type 1")
        RNA.to_csv("./data/D/CKSNAP.csv", "index=False", header=True)
        RNA.get_descriptor("PseDNC")
        RNA.to_csv("./data/D/PseDNC.csv", "index=False", header=True)
        RNA.get_descriptor("ASDC")
        RNA.to_csv("./data/D/ASDC.csv", "index=False", header=True)
        RNA.get_descriptor("Kmer type 1")
        RNA.to_csv("./data/D/Kmer.csv", "index=False", header=True)
        # combine features
        DPCP = pd.read_csv("./data/D/DPCP.csv", header=0).iloc[:, 1:]
        CKSNAP = pd.read_csv("./data/D/CKSNAP.csv", header=0).iloc[:, 1:]
        PseDNC = pd.read_csv("./data/D/PseDNC.csv", header=0).iloc[:, 1:]
        ASDC = pd.read_csv("./data/D/ASDC.csv", header=0).iloc[:, 1:]
        Kmer = pd.read_csv("./data/D/Kmer.csv", header=0).iloc[:, 1:]
        D_train = pd.concat([DPCP, CKSNAP, PseDNC, Kmer, ASDC], axis=1)
        D_train.to_csv("./data/D/D_train.csv", header=True,index=True)
        return D_train

    def getOptimalfeatures(self):
        traindata = pd.read_csv("./data/D/D_train.csv").iloc[:, 1:]
        IFS_features = pd.read_csv("./Model/optimal features_93.txt", header=None).iloc[:, 0].to_list()
        col_name = traindata.columns
        out_list = []
        for i in col_name:
            i = i.replace("(RNA)", "")
            i = i.replace(" ", "")
            out_list.append(i)
        traindata.columns = out_list
        D_train = traindata.iloc[:, IFS_features]
        return D_train

    def predict(self,trainfile):
        D_final = load_model("./Model/D_final")
        label = pd.DataFrame(D_final.predict(trainfile))
        prob = [i[1] for i in D_final.predict_proba(trainfile)]
        prob = pd.DataFrame(prob)
        result = pd.concat([label, prob], axis=1)
        result.columns = ["label", "values"]
        #result = result[result["values"] > 0.5]
        return result

class Mouse(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.min_length = 41

    def generate_encoding(self):
        # generate features
        RNA = iFeatureOmegaCLI.iRNA(self.filepath)
        RNA.import_parameters(datapath+'/iFeatureOmegaCLI/parameters/RNA_parameters_setting.json')
        RNA.get_descriptor("NCP")
        RNA.to_csv("./data/M/NCP.csv", "index=False", header=True)
        RNA.get_descriptor("ENAC")
        RNA.to_csv("./data/M/ENAC.csv", "index=False", header=True)
        RNA.get_descriptor("binary")
        RNA.to_csv("./data/M/binary.csv", "index=False", header=True)
        RNA.get_descriptor("DBE")
        RNA.to_csv("./data/M/DBE.csv", "index=False", header=True)
        RNA.get_descriptor("DPCP type2")
        RNA.to_csv("./data/M/DPCP2.csv", "index=False", header=True)
        # combine features
        NCP = pd.read_csv("./data/M/NCP.csv", header=0).iloc[:, 1:]
        ENAC = pd.read_csv("./data/M/ENAC.csv", header=0).iloc[:, 1:]
        binary = pd.read_csv("./data/M/binary.csv", header=0).iloc[:, 1:]
        DBE = pd.read_csv("./data/M/DBE.csv", header=0).iloc[:, 1:]
        DPCP2 = pd.read_csv("./data/M/DPCP2.csv", header=0).iloc[:, 1:]
        DPCP2 = DPCP2.iloc[:, ~DPCP2.columns.str.contains('Tilt') & ~DPCP2.columns.str.contains('Twist')]
        M_train = pd.concat([NCP, ENAC, binary, DBE, DPCP2], axis=1)
        M_train.to_csv("./data/M/M_train.csv", header=True, index=True)
        return M_train

    def getOptimalfeatures(self):
        traindata = pd.read_csv("./data/M/M_train.csv").iloc[:, 1:]
        IFS_features = pd.read_csv("./Model/optimal features_163.txt", header=None).iloc[:, 0].to_list()
        M_train = traindata.iloc[:, IFS_features]
        return M_train

    def predict(self,trainfile):
        M_final = load_model("./Model/M_final")
        label = pd.DataFrame(M_final.predict(trainfile))
        prob = [i[1] for i in M_final.predict_proba(trainfile)]
        prob = pd.DataFrame(prob)
        result = pd.concat([label, prob], axis=1)
        result.columns = ["label", "values"]
        return result

class Human(object):
    def __init__(self, filepath):
        self.filepath = filepath
        self.min_length = 51

    def generate_encoding(self):
        # generate features
        RNA = iFeatureOmegaCLI.iRNA(self.filepath)
        RNA.import_parameters(datapath + '/iFeatureOmegaCLI/parameters/RNA_parameters_setting.json')
        RNA.get_descriptor("ENAC")
        RNA.to_csv("./data/H/ENAC.csv", "index=False", header=True)
        RNA.get_descriptor("binary")
        RNA.to_csv("./data/H/binary.csv", "index=False", header=True)
        RNA.get_descriptor("PS2")
        RNA.to_csv("./data/H/PS2.csv", "index=False", header=True)
        RNA.get_descriptor("NCP")
        RNA.to_csv("./data/H/NCP.csv", "index=False", header=True)
        # combine features
        ENAC = pd.read_csv("./data/H/ENAC.csv", header=0).iloc[:, 1:]
        binary = pd.read_csv("./data/H/binary.csv", header=0).iloc[:, 1:]
        PS2 = pd.read_csv("./data/H/PS2.csv", header=0).iloc[:, 1:]
        NCP = pd.read_csv("./data/H/NCP.csv", header=0).iloc[:, 1:]
        H_train = pd.concat([ENAC, binary, PS2, NCP], axis=1)
        H_train.to_csv("./data/H/H_train.csv", header=True, index=True)
        return H_train

    def getOptimalfeatures(self):
        traindata = pd.read_csv("./data/H/H_train.csv").iloc[:, 1:]
        IFS_features = pd.read_csv("./Model/optimal features_215.txt", header=None).iloc[:, 0].to_list()
        H_train = traindata.iloc[:, IFS_features]
        return H_train

    def predict(self,trainfile):
        H_final = load_model("./Model/H_final")
        label = pd.DataFrame(H_final.predict(trainfile))
        prob = [i[1] for i in H_final.predict_proba(trainfile)]
        prob = pd.DataFrame(prob)
        result = pd.concat([label, prob], axis=1)
        result.columns = ["label", "values"]
        return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(usage="it's usage tip.",
                                     description="running the ATTIC")
    parser.add_argument("-i", required=True, help="The input file with fasta format")
    parser.add_argument("-s", required=True,
                        choices=['H.sapiens','M.musculus','D.melanogaster'],
                        help="select species types")
    parser.add_argument('-o', required=True, help="output file")
    parser.add_argument('-t', required=True, help="threshold",type=float)
    args = parser.parse_args()
    sequence = read_fasta.Sequence(args.i)
    myfasta = sequence.read_fasta()
    samplename = sequence.output_fasta(myfasta,args.s)

    # Use different model to predict according to different species
    if args.s == "D.melanogaster":
        fly = Fly(datapath +'/data/D/D.txt')
        trainfile = fly.generate_encoding()
        optimalfile = fly.getOptimalfeatures()
        result = fly.predict(optimalfile)
        result.index = samplename
        result.to_csv("./data/D/Predict results.csv",index=True,header=True)
        result_threshold = result[result["values"] > args.t]
        print(result_threshold)
        result_threshold.to_csv("./data/D/Predict results-threshold.csv", index=True, header=True)
        print("Finish prediction.")

    elif args.s == "M.musculus":
        mouse = Mouse(datapath +'/data/M/M.txt')
        trainfile = mouse.generate_encoding()
        optimalfile = mouse.getOptimalfeatures()
        result = mouse.predict(optimalfile)
        result.index = samplename
        result.to_csv("./data/M/Predict results.csv", index=True, header=True)
        result_threshold = result[result["values"] > args.t]
        print(result_threshold)
        result_threshold.to_csv("./data/M/Predict results-threshold.csv", index=True, header=True)
        print("Finish prediction.")

    else:
        human = Human(datapath + '/data/H/H.txt')
        trainfile = human.generate_encoding()
        optimalfile = human.getOptimalfeatures()
        result = human.predict(optimalfile)
        result.index = samplename
        result.to_csv("./data/H/Predict results.csv", index=True, header=True)
        result_threshold = result[result["values"] > args.t]
        print(result_threshold)
        result_threshold.to_csv("./data/H/Predict results-threshold.csv", index=True, header=True)
        print("Finish prediction.")






