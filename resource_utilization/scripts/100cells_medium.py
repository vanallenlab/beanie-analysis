import pandas as pd
import numpy as np

import matplotlib.pyplot as plt
import seaborn as sns

import scanpy as sc

import os
import sys

from memory_profiler import memory_usage
from time import process_time,time

path = "/home/unix/sjohri/valab_sjohri/projects/beanie_private/"
sys.path.append("/home/unix/sjohri/valab_sjohri/projects/2_beanie_revision1/beanie/src")

import beanie_08112022.beanie as bn

def func(count_path, meta_path):
    bobj = bn.Beanie(counts_path = path+"/data/profiling_datasets/"+count_path, 
                     metad_path = path+"/data/profiling_datasets/"+meta_path,
                     sig_path = path+"/data/signatures/test_medium.gmt",
                     normalised = False, min_cells=50, bins=False, bin_size=20,
                     output_dir = path+"/profiling/results/beanie_out_medium/")
    bobj.SignatureScoring(scoring_method="beanie", no_random_sigs=1000)
    bobj.DifferentialExpression(test_name="mwu-test", subsamples=25)
    
paths_count = sorted([x for x in os.listdir(path+"/data/profiling_datasets/") if x[-20:]=="_100cellsperpat.h5ad"])
paths_metad = sorted([x for x in os.listdir(path+"/data/profiling_datasets/") if x[-19:]=="_100cellsperpat.csv"])

for i in range(len(paths_metad)):
    f = open(path+"/profiling/results/test_medium/"+paths_metad[i].split(".")[0]+'.txt', 'w')
    t1_start_wall = time()
    t1_start_cpu = process_time()
    mem_usage = memory_usage((func, [paths_count[i], paths_metad[i]]))
    t1_stop_cpu = process_time()
    t1_stop_wall = time()
    f.write("Wall time (in seconds):" + str(t1_stop_wall-t1_start_wall))
    f.write("\n")
    f.write("CPU time (in seconds):" + str(t1_stop_cpu-t1_start_cpu))
    f.write("\n")
    f.write(str(mem_usage))
    f.close()
    