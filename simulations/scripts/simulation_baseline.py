import pandas as pd
import numpy as np
import scanpy as sc

import os
import pickle as pkl

import beanie.beanie as bn

SIM_DIR_PATH = "/home/unix/sjohri/valab/projects/beanie-analysis/simulations/"
SIG_FILE_PATH = "/home/unix/sjohri/valab/projects/beanie-analysis/data/gene_signatures/hallmark.gmt"

def calculate_baseline(count_path, meta_path, sim_name):
    
    # BEANIE
    bobj = bn.Beanie(counts_path = f"{SIM_DIR_PATH}/counts/{count_path}", 
                     metad_path = f"{SIM_DIR_PATH}/metad/{meta_path}",
                     sig_path = SIG_FILE_PATH,
                     normalised = False, min_cells=50, bins=False,
                     output_dir = f"{SIM_DIR_PATH}/out/baseline_{sim_name}/")
    
    bobj.SignatureScoring(scoring_method="beanie", no_random_sigs=1000)
    bobj.signature_scores.to_csv(f"{SIM_DIR_PATH}/baseline/sig_scores/{meta_path}")
    pkl.dump(bobj._null_dist_scores, open(f"{SIM_DIR_PATH}/baseline/bg_scores/{meta_path[:-4]}.pkl","wb"))
    
    bobj.DifferentialExpression(test_name="mwu-test", subsamples=100)
    bobj.de_summary.to_csv(f"{SIM_DIR_PATH}/baseline/beanie_de/{meta_path}")
    

if __name__ == "__main__":
    
    sim_names = sorted([x[:-4] for x in os.listdir(f"{SIM_DIR_PATH}/metad/")])

    for sim in sim_names:
        calculate_baseline(f"{sim}.h5ad", f"{sim}.csv", sim)