import scanpy as sc
import pandas as pd
import numpy as np

import os
import gc
import random

from tqdm.auto import tqdm
import pickle as pkl

import beanie.beanie as bn

SIM_DIR_PATH = "/home/unix/sjohri/valab/projects/beanie-analysis/simulations/"
SIG_FILE_PATH = "/home/unix/sjohri/valab/projects/beanie-analysis/data/gene_signatures/hallmark.gmt"

def perturb_helper(pats,percent=0.5):
    
    # select which patients to perturb
    return random.sample(pats,max(1,int(len(pats)*percent)))


def create_perturbation(count_path, meta_path, sim_name):
    
    bobj = bn.Beanie(counts_path = f"{SIM_DIR_PATH}/counts/{count_path}", 
                     metad_path = f"{SIM_DIR_PATH}/metad/{meta_path}",
                     sig_path = SIG_FILE_PATH,
                     normalised = False, min_cells=50, bins=False,
                     output_dir = f"{SIM_DIR_PATH}/out/perturb_{sim_name}/")
    
    bobj._null_dist_scores = pkl.load(open(f"{SIM_DIR_PATH}/baseline/bg_scores/{meta_path[:-4]}.pkl", "rb"))
    unperturbed_sig_scores = pd.read_csv(f"{SIM_DIR_PATH}/baseline/sig_scores/{meta_path}", index_col=0)
    
    p_pct = [0.1, 0.2, 0.3, 0.4, 0.5, 1]
    p_mag = [1, 2, 3]
    
    
    for p1 in p_pct:
        for p2 in p_mag:            
            # perturb
            gr1_pats = bobj.t1_ids
            gr2_pats = bobj.t2_ids

            perturbed_sig_scores = unperturbed_sig_scores.copy()
            perturbed_pats = perturb_helper(gr1_pats, p1)

            for pat in gr1_pats:
                if pat in perturbed_pats:
                    flag=1
                else:
                    flag=0
                cells = bobj.d1_all[pat]
                stdevs1 = unperturbed_sig_scores.loc[cells,].std(axis=0)

                for sig_name in stdevs1.index:
                    perturbed_sig_scores.loc[cells,sig_name] = perturbed_sig_scores.loc[cells,sig_name] + flag*p2*stdevs1[sig_name]
            
            perturbed_sig_scores.to_csv(f"{SIM_DIR_PATH}/perturbations/{int(100*p1)}perc/{int(p2)}std/sig_scores/{meta_path}")
            
            bobj.signature_scores = perturbed_sig_scores
            
            # calculate DE
            bobj.DifferentialExpression(test_name="mwu-test", subsamples=100)
            bobj.de_summary.to_csv(f"{SIM_DIR_PATH}/perturbations/{int(100*p1)}perc/{int(p2)}std/beanie_de/{meta_path}")
            
            # reset object
            bobj.signature_scores = None
            bobj._differential_expression_run = False
    
    # cleanup
    del bobj
    gc.collect()
    
if __name__ == "__main__":
    
    random.seed(20)
    sim_names = sorted([x[:-4] for x in os.listdir(f"{SIM_DIR_PATH}/metad/")])

    for sim in sim_names:
        create_perturbation(f"{sim}.h5ad", f"{sim}.csv", sim)
        