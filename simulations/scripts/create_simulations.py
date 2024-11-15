import random
import scanpy as sc
import pandas as pd

from tqdm.auto import tqdm

SIM_DIR_PATH = "/home/unix/sjohri/valab/projects/beanie-analysis/simulations/"

def create_simulations(adata, pats_list, num_pats=2, num_sims=1000, seed=20):
    random.seed(seed)
    pats = random.sample(pats_list,num_pats)

    for i in tqdm(range(num_sims)):
        cells_list1 = []
        cells_list2 = []
        for s in pats:
            cells = sorted(adata[adata.obs.patient_id==s].obs.index)
            if len(cells)>100:
                temp = random.sample(cells,100)
                cells_list1.extend(temp)
                cells_temp = sorted(set(cells).difference(temp))
                if len(cells_temp)>100:
                    cells_list2.extend(random.sample(cells_temp,100))
                elif len(cells_temp)==100:
                    cells_list2.extend(cells_temp)

        adata_subset1 = adata1[cells_list1]
        adata_subset1.obs["sample_id"] = [x + "_gr1" for x in adata_subset1.obs.patient_id]
        adata_subset1.obs["group_id"] = "group1"

        adata_subset2 = adata1[cells_list2]
        adata_subset2.obs["sample_id"] = [x + "_gr2" for x in adata_subset2.obs.patient_id]
        adata_subset2.obs["group_id"] = "group2"

        adata_concat = adata_subset1.concatenate(adata_subset2)
        
        # the adata.raw.X will be used for scoring by BEANIE, which contains log-normalized genes expression.
        adata_concat.write(f"{SIM_DIR_PATH}/counts/sim_{i}.h5ad")
        adata_concat.obs[["sample_id","group_id"]].to_csv(f"{SIM_DIR_PATH}/metad/sim_{i}.csv")
        

if __name__ == "__main__":
    
    adata1 = sc.read_10x_h5(path+"/data/brca/bassez_2021/1863-counts_cells_cohort1.h5")
    adata1.obs = pd.read_csv(path+"/data/brca/bassez_2021/raw_data/1872-BIOKEY_metaData_cohort1_web.csv", index_col=0)
    
    adata_gr1 = adata1[(adata1.obs.BC_type=="TNBC") & (adata1.obs.timepoint=="Pre") & (adata1.obs.cellType=="Cancer_cell")]
    
    sc.pp.normalize_total(adata_gr1, target_sum=1e4)
    sc.pp.log1p(adata_gr1)
    adata_gr1.raw = adata_gr1

    sc.pp.highly_variable_genes(adata_gr1, n_top_genes=2500)
    sc.pl.highly_variable_genes(adata_gr1)
    print("Highly variable genes: %d"%sum(adata_gr1.var.highly_variable))
    var_genes_all = adata_gr1.var.highly_variable
    adata_gr1 = adata_gr1[:,var_genes_all]
    
    sc.pp.scale(adata_gr1, max_value=10)
    sc.pp.pca(adata_gr1, random_state=initialization, svd_solver='arpack', n_comps=100)
    sc.pl.pca_variance_ratio(adata_gr1, n_pcs= 100, log=True, show = True)

    sc.external.pp.harmony_integrate(adata_gr1, key=["patient_id"],
                                 random_state=initialization, max_iter_harmony=50)
    adata_gr1.obsm['X_pca'] = adata_gr1.obsm['X_pca_harmony']
    
    npc = 50
    sc.pp.neighbors(adata_gr1, random_state=initialization, n_neighbors=20, n_pcs=npc)
    sc.tl.leiden(adata_gr1, random_state=initialization, resolution=0.1)
    sc.tl.umap(adata_gr1, random_state=initialization, min_dist=0.1)
    
    cell_counts = adata_gr1.obs[adata_gr1.obs.leiden=="0"].patient_id.value_counts()
    patients = sorted(cell_counts[cell_counts>=200].index)
    
    ### Create simulated dataset
    select_cells(adata_gr1[adata_gr1.obs.leiden=="0"], patients, 10)
    