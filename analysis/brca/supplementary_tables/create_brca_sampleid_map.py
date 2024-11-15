import pandas as pd
import scanpy as sc
import pickle as pkl

ad1 = sc.read_h5ad("./beanie_inputs/analysis_1/adata_full.h5ad")
ad2 = sc.read_h5ad("./beanie_inputs/analysis_2/adata_full.h5ad")

df1 = pd.DataFrame(ad1.obs.sample_id.value_counts())
df1.index = df1.index.astype(str)
df2 = pd.DataFrame(ad2.obs.sample_id.value_counts())
df2.index = df2.index.astype(str)
df_brca = pd.concat([df1,df2])
df_brca = df_brca.sort_values(by="sample_id", ascending=False)
df_brca = df_brca.reset_index()
df_brca.columns=["sample_id", "num_tumor_cells"]

# simplify sample ids
dmap = {x:f"BRCA_{idx}" for idx,x in enumerate(df_brca.sample_id)}

pkl.dump(dmap, open("brca_sampleid_map.pkl","wb"))