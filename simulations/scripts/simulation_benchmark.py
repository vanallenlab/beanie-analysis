import pandas as pd
import numpy as np
import scipy as sp

import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.formula.api import glm
from statsmodels.stats.multitest import multipletests

import multiprocessing
from multiprocessing import Pool

SIM_DIR_PATH = "/home/unix/sjohri/valab/projects/beanie-analysis/simulations/"

def calc_mwu(sigs_df, t1, t2, corr_method = "fdr_bh"):
    mwu_res = pd.Series(index=sigs_df.columns, dtype=float)
    for x in sigs_df.columns:
        try:
            test = sp.stats.mannwhitneyu(sigs_df.loc[t1,x], sigs_df.loc[t2,x])
            mwu_res[x] = test[1]
        except ValueError:
            mwu_res[x] = 1

    mwu_res_bh = pd.Series(multipletests(mwu_res, method=corr_method)[1], index=sigs_df.columns)

    df1 = pd.DataFrame(mwu_res_bh, columns=["mwu_pval"])
    return df1


def calc_glm(sigs_df, t1, t2):
    temp = sigs_df.copy()
    temp.columns = temp.columns.str.replace(".","_")
    glm_res = pd.Series(index=temp.columns, dtype=float)
    temp["group_id"] = ["group1" if x in t1 else "group2" for x in temp.index]
    for x in temp.columns[:-1]:
        model = smf.glm(formula=f'group_id ~ 1 + {x}', data=temp, 
                        family=sm.families.Binomial())
        try:
            res = model.fit()
            glm_res[x] = (res.pvalues[1:]).values
        except:
            glm_res[x] = [1]
            
    df1 = pd.DataFrame(glm_res,columns=["glm_pval"])
    df1.index = sigs_df.columns
    return df1


def calc_pseudobulk_pval(sigs_df, t1, t2, corr_method = "fdr_bh"):
    mwu_res = pd.Series(index=sigs_df.columns, dtype=float)
    for x in sigs_df.columns:
        try:
            test = sp.stats.mannwhitneyu(sigs_df.loc[t1,x],sigs_df.loc[t2,x])
            mwu_res[x] = test[1]
        except ValueError:
            mwu_res[x] = 1

    mwu_res_bh = pd.Series(multipletests(mwu_res, method=corr_method)[1], index=sigs_df.columns)

    df1 = pd.DataFrame(mwu_res_bh,columns=["pseudobulk_pval"])
    return df1


def calc_mwu_loocv_pval(sigs_df, t1_map, t2_map, corr_method = "fdr_bh"):    
    all_res = []
    for pat in t1_map.keys():
        cells_to_drop = t1_map[pat]
        temp = sigs_df.drop(cells_to_drop, axis=0)
        t1_all = [x2 for x1 in sorted(t1_map.values()) for x2 in x1]
        t1_temp = sorted(set(t1_all) - set(cells_to_drop))
        t2 = [x2 for x1 in sorted(t2_map.values()) for x2 in x1]
        mwu_res = pd.Series(index=temp.columns, dtype=float)
        for x in temp.columns:
            try:
                test = sp.stats.mannwhitneyu(temp.loc[t1_temp,x],temp.loc[t2,x])
                mwu_res[x] = test[1]
            except ValueError:
                mwu_res[x] = 1

        mwu_res_bh = pd.Series(multipletests(mwu_res, method=corr_method)[1], index=temp.columns)

        df1 = pd.DataFrame(mwu_res_bh,columns=["mwu_pval"])
        all_res.append(df1)
        
    for pat in t2_map.keys():
        cells_to_drop = t2_map[pat]
        temp = sigs_df.drop(cells_to_drop, axis=0)
        t2_all = [x2 for x1 in sorted(t2_map.values()) for x2 in x1]
        t2_temp = sorted(set(t2_all) - set(cells_to_drop))
        t1 = [x2 for x1 in sorted(t1_map.values()) for x2 in x1]
        for x in temp.columns:
            try:
                test = sp.stats.mannwhitneyu(temp.loc[t1,x],temp.loc[t2_temp,x])
                mwu_res[x] = test[1]
            except ValueError:
                mwu_res[x] = 1

        mwu_res_bh = pd.Series(multipletests(mwu_res, method=corr_method)[1], index=temp.columns)

        df1 = pd.DataFrame(mwu_res_bh,columns=["mwu_pval"])
        all_res.append(df1)
    
    pval_df = pd.concat(all_res, axis=1)    
    df_final = pd.DataFrame([pval_df.median(axis=1),((pval_df<=0.05).sum(axis=1)!=len(all_res))], 
                            index=["mwu_loocv_pval","mwu_nonrobust"])
    
    return df_final.T


def calc_glm_loocv_pval(sigs_df, t1_map, t2_map, corr_method = "fdr_bh"):
    all_res = []
    for pat in t1_map.keys():
        temp = sigs_df.copy()
        temp.columns = temp.columns.str.replace(".","_")
        glm_res = pd.Series(index=temp.columns, dtype=float)
        
        cells_to_drop = t1_map[pat]
        temp = sigs_df.drop(cells_to_drop, axis=0)
        t1_all = [x2 for x1 in sorted(t1_dmap.values()) for x2 in x1]
        t1_temp = sorted(set(t1_all) - set(cells_to_drop))
        t2 = [x2 for x1 in sorted(t2_map.values()) for x2 in x1]
        temp["group_id"] = ["group1" if x in t1_temp else "group2" for x in temp.index]        
        for x in temp.columns[:-1]:
            model = smf.glm(formula=f'group_id ~ 1 + {x}', data=temp, 
                            family=sm.families.Binomial())
            try:
                res = model.fit()
                glm_res[x] = (res.pvalues[1:]).values
            except:
                glm_res[x] = [1]
            
        df1 = pd.DataFrame(glm_res,columns=["glm_pval"])
        df1.index = sigs_df.columns
        all_res.append(df1)
        
    for pat in t2_map.keys():
        temp = sigs_df.copy()
        temp.columns = temp.columns.str.replace(".","_")
        glm_res = pd.Series(index=temp.columns, dtype=float)
        
        cells_to_drop = t2_map[pat]
        temp = sigs_df.drop(cells_to_drop, axis=0)
        t2_all = [x2 for x1 in sorted(t2_dmap.values()) for x2 in x1]
        t2_temp = sorted(set(t2_all) - set(cells_to_drop))
        t1 = [x2 for x1 in sorted(t1_map.values()) for x2 in x1]
        temp["group_id"] = ["group2" if x in t2_temp else "group1" for x in temp.index]        
        for x in temp.columns[:-1]:
            model = smf.glm(formula=f'group_id ~ 1 + {x}', data=temp, 
                            family=sm.families.Binomial())
            try:
                res = model.fit()
                glm_res[x] = (res.pvalues[1:]).values
            except:
                glm_res[x] = [1]
            
        df1 = pd.DataFrame(glm_res,columns=["glm_pval"])
        df1.index = sigs_df.columns
        all_res.append(df1)
    
    pval_df = pd.concat(all_res, axis=1)    
    df_final = pd.DataFrame([pval_df.median(axis=1),((pval_df<0.05).sum(axis=1)!=len(all_res))], 
                            index=["glm_loocv_pval","glm_nonrobust"])
    
    return df_final.T


def calc_pseudobulk_loocv_pval(sigs_df, t1, t2, corr_method = "fdr_bh"):
    pats = sigs_df.index.tolist()
    
    all_res = []
    for pat in t1:
        temp = sigs_df.drop(pat, axis=0)
        t1_temp = [x for x in t1 if x != pat]
        mwu_res = pd.Series(index=temp.columns, dtype=float)
        for x in temp.columns:
            try:
                test = sp.stats.mannwhitneyu(temp.loc[t1_temp,x],temp.loc[t2,x])
                mwu_res[x] = test[1]
            except ValueError:
                mwu_res[x] = 1

        mwu_res_bh = pd.Series(multipletests(mwu_res, method=corr_method)[1], index=temp.columns)

        df1 = pd.DataFrame(mwu_res_bh,columns=["mwu_pval"])
        all_res.append(df1)
        
    for pat in t2:
        temp = sigs_df.drop(pat, axis=0)
        t2_temp = [x for x in t2 if x != pat]
        mwu_res = pd.Series(index=temp.columns, dtype=float)
        for x in temp.columns:
            try:
                test = sp.stats.mannwhitneyu(temp.loc[t1,x],temp.loc[t2_temp,x])
                mwu_res[x] = test[1]
            except ValueError:
                mwu_res[x] = 1

        mwu_res_bh = pd.Series(multipletests(mwu_res, method=corr_method)[1], index=temp.columns)

        df1 = pd.DataFrame(mwu_res_bh,columns=["mwu_pval"])
        all_res.append(df1)
    
    pval_df = pd.concat(all_res, axis=1)
    
    df_final = pd.DataFrame([pval_df.median(axis=1),((pval_df<0.05).sum(axis=1)!=len(all_res))], index=["pseudobulk_loocv_pval","pseudobulk_nonrobust"])
    
    return df_final.T

def process_simulation(sim):
    print(f"Processing sim {sim}.")
    
    mtd = pd.read_csv(f"{SIM_DIR_PATH}/metad/{sim}.csv", index_col=0)    
    groups = sorted(set(mtd.group_id))
    
    t1_ids = sorted(set(mtd.loc[mtd.group_id=="group1", "sample_id"]))
    t1_dmap = {}
    for t1 in t1_ids:
        t1_dmap[t1] = sorted(mtd[mtd.sample_id==t1].index)

    t2_ids = sorted(set(mtd.loc[mtd.group_id=="group2", "sample_id"]))
    t2_dmap = {}
    for t2 in t2_ids:
        t2_dmap[t2] = sorted(mtd[mtd.sample_id==t2].index)
        
      
    samples = t1_ids + t2_ids
    for s in samples:
        dmap[s] = sorted(mtd[mtd.sample_id==s].index)
    
    # process baseline

    s_score = pd.read_csv(f"{SIM_DIR_PATH}/baseline/sig_scores/{sim}.csv", index_col=0)         
    pseudobulk_df = pd.DataFrame(index = samples, columns=s_score.columns)
    for sid in samples:
        pseudobulk_df.loc[sid,:] = s_score.loc[dmap[sid],].mean(axis=0).values

    mw = calc_mwu_pval(s_score, t1_ids, t2_ids)
    gl = calc_glm_pval(s_score, t1_ids, t2_ids)
    ps = calc_pseudobulk_pval(pseudobulk_df, t1_ids, t2_ids)

    mw_loocv = calc_mwu_loocv_pval(s_score, t1_dmap, t2_dmap)
    gl_loocv = calc_glm_loocv_pval(s_score, t1_dmap, t2_dmap)
    ps_loocv = calc_pseudobulk_loocv_pval(pseudobulk, t1_ids, t2_ids)

    pd.concat([mw, gl, ps, mw_loocv, gl_loocv, ps_loocv], axis=1).to_csv(f"{SIM_DIR_PATH}/baseline/benchmark_de/{sim}.csv")

    # process perturbation
    perturbs = {}
    res_perturb = {}
        
    mag = [1, 2, 3]
    perc = [10, 20, 30, 40, 50, 100]
    
    for p in perc:
        if p not in perturbs.keys():
            perturbs[p] = {k:[] for k in mag}
            res_benchmark[p] = {k: [] for k in mag}
            
        for m in mag:            
            
            # perturbations
            s_score = pd.read_csv(f"{SIM_DIR_PATH}/perturbations/{int(p)}perc/{int(m)}std/sig_scores/{sim}.csv", index_col=0)         
            pseudobulk_df = pd.DataFrame(index = samples, columns=s_score.columns)
            for sid in samples:
                pseudobulk_df.loc[sid,:] = s_score.loc[dmap[sid],].mean(axis=0).values
            
            mw = calc_mwu_pval(s_score, t1_ids, t2_ids)
            gl = calc_glm_pval(s_score, t1_ids, t2_ids)
            ps = calc_pseudobulk_pval(pseudobulk_df, t1_ids, t2_ids)
            
            mw_loocv = calc_mwu_loocv_pval(s_score, t1_dmap, t2_dmap)
            gl_loocv = calc_glm_loocv_pval(s_score, t1_dmap, t2_dmap)
            ps_loocv = calc_pseudobulk_loocv_pval(pseudobulk, t1_ids, t2_ids)

            res_benchmark[p][m].append(pd.concat([mw, gl, ps, mw_loocv, gl_loocv, ps_loocv], axis=1))
            
    pkl.dump(res_benchmark, open(f"{SIM_DIR_PATH}/perturbations/benchmark_de/{sim}.pkl","wb"))
    print(f"Completed sim {sim}.")

if __name__ == "__main__":    
    
    sims = [f"sim_{i}" for i in range(1000)]

    num_cpus = multiprocessing.cpu_count()
    print(f"Parallelising over {num_cpus} cpus.")
    
    with Pool(processes=num_cpus) as p:
        p.map(process_simulation, [s for s in sims])
        