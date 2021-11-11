import matplotlib
import scanpy as sc
import pandas as pd
import numpy as np

CWD = 'F:/scanpy/ALL2'


filename = "{CWD}/all_scanpy.loom".format(CWD=CWD)

adata = sc.read("ALL_scanpy.h5ad", cache = True)
metadata = pd.read_csv("meta.csv", index_col = 0)
adata = sc.read("{CWD}/ALL_scanpy.h5ad".format(CWD=CWD))

adata = sc.read("{CWD}/seurat_obj_immu.mtx".format(CWD=CWD), cache = True).T

# load gene names
adata.var_names = pd.read_csv("{CWD}/immu_genenames.csv".format(CWD=CWD)).iloc[:, 1]

# load cell names
adata.obs_names = pd.read_csv("{CWD}/immu_cellnames.csv".format(CWD=CWD)).iloc[:, 1]

# add cell labels
#metadata = pd.read_csv("meta.csv", index_col = 0)
#metadata = pd.read_csv("{CWD}/seurat_obj_immu.csv".format(CWD=CWD), index_col = 0)
metadata = pd.read_csv("{CWD}/meta.csv".format(CWD=CWD), index_col = 0)

adata.obs = metadata

#umap=  pd.read_csv("{CWD}/immu_umap.csv".format(CWD=CWD), index_col = 0)
#hm =  pd.read_csv("{CWD}/immu_harmony.csv".format(CWD=CWD), index_col = 0)
#adata.obsm['X_umap'] = umap.values
#adata.obsm['X_pca'] = hm.values

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata = adata[:, adata.var.highly_variable]
#sc.pp.regress_out(adata, ['percent.mt'])
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata)
#batch correction
#sc.external.pp.bbknn(adata, batch_key='patient_id')
#harmony
import harmonypy as hm
ho = hm.run_harmony(adata.obsm['X_pca'], adata.obs, 'patient_id', theta=2)
res = pd.DataFrame(ho.Z_corr)
res = np.array(res.T)
adata.obsm['X_pca'] = res

#neighbors
sc.pp.neighbors(adata, n_neighbors=20, n_pcs=50)

#umap
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution = 1)
sc.pl.umap(adata, color=['leiden'], legend_loc = 'on data')
sc.pl.umap(adata, color=['CD8A', 'CD4', 'NKG7', 'GNLY'])

marker_genes = ['CD8A', 'CD4', 'NKG7', 'GNLY', 'CCR7','FOXP3','CXCL13','CCL5']
sc.pl.umap(adata, color=['CD8A', 'CD4', 'NKG7', 'GNLY', 'CCR7','FOXP3','CXCL13','CCL5'])

#diffusion map
sc.pp.neighbors(adata, n_neighbors = 30, n_pcs = 25, knn = True, random_state = 10, method = "gauss")


# compute diffusion map
sc.tl.diffmap(adata, n_comps = 20)
sc.pl.diffmap(adata, color = ['sample_type','leiden'])
#save UMAP
pd.DataFrame(adata.obsm['X_umap']).to_csv("{CWD}/all_scanpy_umap.csv".format(CWD=CWD))
#save meta_data
adata.obs.to_csv("{CWD}/resolution2.csv".format(CWD=CWD))

#dotplot marker
marker_genes = ['CD3D', 'CD4','CD8A','GNLY', 'NKG7', 'CD79A', 'CCR7',  'LYZ', 'KLRB1', 'CLDN4', 'DCN',
                'FCGR3A',  'CST3']
ax = sc.pl.dotplot(adata, marker_genes, groupby='leiden')
ax = sc.pl.dotplot(adata, marker_genes, groupby='leiden_label', dendrogram=True)

#leiden to label
new_label = ['CD4-C1', 'CD8-C2', 'CD4-C2', 'CD4-C3', 'CD4-C2', 'CD8-C3',
        'CD4-C5', 'CD4-C6', 'CD8-C4', 'CD4-C5', 'CD4-C4', 'CD4-C5',
        'CD8-C2', 'CD4-C2', 'CD4-C2', 'NK-C1', 'CD8-C1', 'CD8-C4', 'NK-C3',
        'CD4-C4', 'CD8-C5', 'NK-C2']
my_marker = ['CD14', 'LYZ','CD68', 'CCL3',  'S100A12', 'SPP1', 'CPA3',
            'CST6', 'MKI67',  'CHI3L1', 'LAMP3', 'PTGDS','CLEC9A','CD1C', 'CD7', 'TRBC2']
ax = sc.pl.dotplot(adata, my_marker, groupby='leiden', dendrogram=True)
label = [new_label[int(i)] for i in adata.obs.leiden]
adata.obs['leiden_label'] = label

#save file
save_file = 'F:/scanpy/ALL'
#save metafile
adata.obs.to_csv("{SAVE}/immue_harmony2_leiden2_2.csv".format(SAVE=save_file))

adata.obs.to_csv("{CWD}/meta.csv".format(CWD=CWD))
#save umap
pd.DataFrame(adata.obsm['X_umap']).to_csv("{CWD}/immue_scanpy_umap.csv".format(CWD=CWD))
#save Project
adata.write("{SAVE}/ALL_scanpy.h5ad".format(SAVE=SAVE))

adata.write("{CWD}/ALL_scanpy.h5ad".format(CWD=CWD))

adata.obs.to_csv("{CWD}/ALL_scanpy.csv".format(CWD=CWD))

#subset
