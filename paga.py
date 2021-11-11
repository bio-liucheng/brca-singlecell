
adata_subset = sc.read("{CWD}/seurat_obj_immu.mtx".format(CWD=CWD), cache = False).T

# load gene names
adata_subset.var_names = pd.read_csv("{CWD}/immu_genenames.csv".format(CWD=CWD)).iloc[:, 1]

# load cell names
adata_subset.obs_names = pd.read_csv("{CWD}/immu_cellnames.csv".format(CWD=CWD)).iloc[:, 1]

# add cell labels
metadata = pd.read_csv("{CWD}/seurat_obj_immu.csv".format(CWD=CWD), index_col = 0)

adata_subset.obs = metadata


umap=  pd.read_csv("{CWD}/immu_umap.csv".format(CWD=CWD), index_col = 0)
hm =  pd.read_csv("{CWD}/immu_harmony.csv".format(CWD=CWD), index_col = 0)
adata_subset.obsm['X_umap'] = umap.values
adata_subset.obsm['X_pca'] = hm.values

sc.pp.neighbors(adata_subset, n_neighbors = 30, n_pcs = 20, knn = True, random_state = 10, method = "gauss")

# compute diffusion map
sc.tl.diffmap(adata_subset, n_comps = 20)
sc.pl.diffmap(adata_subset, color = ['sample_type','cell_label'])
#---------------------------------------------------
# set root
root_cell_type = 'CD8-C4'
adata_subset.uns['iroot'] = np.flatnonzero(adata_subset.obs['cell_label'] == root_cell_type)[0]

# compute dpt
print("computing sc.tl.dpt")
sc.tl.dpt(adata_subset, n_dcs = 20)
pdt = adata_subset.obs["dpt_pseudotime"].to_csv("{CWD}/pseudotime.csv".format(CWD=CWD))

0
sc.pl.dpt_groups_pseudotime(adata_subset)
sc.pl.dpt_timeseries(adata_subset)

# pdt is at scObj.obs["dpt_pseudotime"]
print("displaying pdt table stored in scObj")
print(adata_subset.obs["dpt_pseudotime"])
pdt = scObj.obs["dpt_pseudotime"].to_csv("{CWD}/material/pseudotime.csv".format(CWD=CWD))

# save the pseudotime
#PAGA
sc.pp.neighbors(adata_subset, n_neighbors=5, n_pcs=10)
sc.tl.draw_graph(adata_subset)
sc.pl.draw_graph(adata_subset, color='cell_label', legend_loc='on data')

#---------------------------------------------------
adata_subset.uns['velocity_graph'] = adata_vlc_select.uns['velocity_graph']
sc.tl.paga(adata_subset,groups='cell_label', use_rna_velocity=True)
sc.tl.paga(adata_subset,groups='cell_label')
sc.pl.paga(adata_subset, threshold=0.1)

sc.pl.paga(adata_subset, threshold=0.5,
    arrowsize = 10, edge_width_scale = 0.5,
    transitions = "transitions_confidence",
#    layout = 'rt',
    dashed_edges = "connectivities",
    color=['leiden'],
    frameon = False,
    labels = None)
