
sc.pp.normalize_total(adata_subset, target_sum=1e4)
sc.pp.log1p(adata_subset)
adata_subset.raw = adata_subset

adata.raw = adata
sc.pp.highly_variable_genes(adata_subset, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata_subset = adata_subset[:, adata_subset.var.highly_variable]
sc.pp.regress_out(adata_subset, ['nCount_RNA', 'percent.mt'])
sc.pp.scale(adata_subset, max_value=10)
sc.tl.pca(adata_subset)
#batch correction
sc.external.pp.bbknn(adata_subset, batch_key='patient_id')
#harmony
adata_subset.obs
import harmonypy as hm
ho = hm.run_harmony(adata_subset.obsm['X_pca'], adata_subset.obs, 'patient_id', theta=1)
res = pd.DataFrame(ho.Z_corr)
res = np.array(res.T)
adata_subset.obsm['X_pca'] = res

#neighbors
sc.pp.neighbors(adata_subset, n_neighbors=20, n_pcs=40)

#umap
sc.tl.umap(adata_subset)
sc.pl.umap(adata_subset, color=['sample_type', 'patient_id'])
sc.pl.umap(adata_subset, color=['main_label', 'CD3D'])
sc.tl.leiden(adata_subset, resolution = 1.5)
sc.pl.umap(adata_subset, color=['leiden'], legend_loc = 'on data')
sc.pl.umap(adata_subset, color=['CD8A', 'CD4', 'NKG7', 'CD3D'])
ax = sc.pl.dotplot(adata_subset, marker_genes, groupby='leiden')

sc.tl.rank_genes_groups(adata_subset), 'leiden', method='wilcoxon')

filename = "{CWD}/immu_scanpy.loom".format(CWD=CWD)
