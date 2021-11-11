
adata.obs = imu_obs
adata_subset2 = adata[adata.obs['new_label'].isin(['CD4-C1', 'CD4-C2', 'CD4-C3','CD4-C4', 'CD4-C5', 'CD4-C6'])]
adata_subset2 = adata[adata.obs['new_label'].isin(['NK cells', 'CD8 T', 'CD4 T','Immue_unkown'])]
adata_subset2 = adata[adata.obs['new_label'].isin(['CAFs'])]
adata_subset2 = adata[adata.obs['new_label'].isin(['B', 'Plasma B'])]
adata_subset2 = adata[adata.obs['new_label'].isin(['Malignant cells'])]
adata_subset2 = adata[adata.obs['new_label'].isin(['Myeloid cells'])]
adata_subset2 = adata_subset2[~adata_subset2.obs['cell_label'].isin(['Tumor-CDT', 'Tumor-unkown','Tumor-B', 'Tumor-C11'])]
adata_subset2 = adata[adata.obs['orig.ident'].isin(['CA1', 'CA2', 'CA3', 'CA4', 'CA5',
        'CA6', 'CA7', 'CA8', 'LN1', 'LN3', 'LN4', 'LN5', 'LN6','LN7', 'LN8' ])]
adata_subset2 = adata[adata.obs['new_label'].isin(['TECs'])]
adata_subset2 = adata_subset2[~adata_subset2.obs['cell_label'].isin(['CAF-CDT', 'CAF-unkown','CAF-B'])]

adata_subset2 = adata[adata.obs['cell_label'].str.contains('CD4-|CD8-|^NK-|KI67-CD8')]
adata_subset2 = adata[adata.obs['cell_label'].str.contains('^Macro-|^DC-|Mast')]
adata_subset2 = adata_subset2[~adata_subset2.obs['cell_label'].isin(['Macro-CDT', 'Macro-B','CAF-B'])]
adata_subset2 = adata_subset[adata_subset.obs['cell_label'].str.contains('CD8-')]
#adata = adata[adata.obs_names.isin(cell_name.x)]
adata_subset2 = adata_subset2[adata_subset2.obs['patient_id'].isin(['patient_4'])]


adata_subset2 =adata_subset[adata_subset.obs['cell_label'].str.contains('Macro-')]
adata_subset2 =adata_subset[adata_subset.obs['cell_label'].str.contains('Macro-')]
adata_subset2 =adata_subset2[~adata_subset2.obs['cell_label'].isin(['Macro-C6'])]
meta = adata_subset2.obs
adata_subset2 = adata[adata_subset2.obs_names]
adata_subset2.obs = meta


sc.pp.normalize_total(adata_subset2, target_sum=1e4)
sc.pp.log1p(adata_subset2)
#adata_subset2.raw = adata_subset2

#adata.raw = adata
sc.pp.highly_variable_genes(adata_subset2, min_mean=0.0125, max_mean=3, min_disp=0.5)

adata_subset2 = adata_subset2[:, adata_subset2.var.highly_variable]
#sc.pp.regress_out(adata_subset2, ['nCount_RNA', 'percent.mt'])
sc.pp.scale(adata_subset2, max_value=10)
sc.tl.pca(adata_subset2)
#batch correction
#sc.external.pp.bbknn(adata_subset2, batch_key='patient_id')
#harmony
adata_subset2.obs
import harmonypy as hm
ho = hm.run_harmony(adata_subset2.obsm['X_pca'], adata_subset2.obs, 'patient_id', theta=2)
res = pd.DataFrame(ho.Z_corr)
res = np.array(res.T)
adata_subset2.obsm['X_pca'] = res

#neighbors
sc.pp.neighbors(adata_subset2, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata_subset2)
#sc.pl.umap(adata_subset2, color = ['FCER1G','TYROBP', 'cell_label'], legend_loc = 'on data')
sc.tl.leiden(adata_subset2, resolution = 1)
sc.pl.umap(adata_subset2, color = ['sample_type','cell_label', 'leiden'], legend_loc = 'on data')
#diffusion mapvalues
sc.pp.neighbors(adata_subset2, n_neighbors = 10, n_pcs = 10, knn = True, random_state = 10, method = "gauss")
sc.tl.diffmap(adata_subset2, n_comps = 20)
sc.pl.diffmap(adata_subset2, color = ['VCAN','THBS1'])
#set cd8 root CD8-C1
root_cell_type = 'CD8-C1'
adata_subset2.uns['iroot'] = np.flatnonzero(adata_subset2.obs['cell_label'] == root_cell_type)[0]
#comput dpt
sc.tl.dpt(adata_subset2, n_dcs = 20)
sc.pl.diffmap(adata_subset2, color = ['sample_type','cell_label', 'dpt_pseudotime'])


zeileis_colors = np.array(sc.pl.palettes.zeileis_28)
new_colors = np.array(adata_subset2.uns['louvain_anno_colors'])
adata_subset2.uns['louvain_anno_colors'] = new_colors


gene_names = [ 'CD248', 'PDCD1', 'CXCL13', 'GZMK', 'CCL5', 'NKG7','IL7R','ANXA1','ZNF683',
              'SELL', 'CCL4', 'GZMA']

paths = [('CD8T_path', ['CD8-C1', 'CD8-C2', 'CD8-C3', 'CD8-C4', 'CD8-C5'])]


adata_subset2.obs['distance'] = adata_subset2.obs['dpt_pseudotime']

adata_subset2.obs['clusters'] = adata_subset2.obs['cell_label']

adata_subset2.uns['clusters_colors'] = adata_subset2.uns['cell_label_colors']

!mkdir write


_, axs = pl.subplots(ncols=2, figsize=(10, 5), gridspec_kw={'wspace': 0.2, 'left': 0.12})
pl.subplots_adjust(left=0.05, right=0.98, top=0.82, bottom=0.2)
for ipath, (descr, path) in enumerate(paths):
    _, data = sc.pl.paga_path(
        adata_subset2, path, gene_names,
        show_node_names=False,
        ax=axs[ipath],
        ytick_fontsize=12,
        left_margin=0.15,
        n_avg=50,
        annotations=['distance'],
        show_yticks=True if ipath==0 else False,
        show_colorbar=False,
        color_map='Greys',
        groups_key='cell_label',
        color_maps_annotations={'distance': 'viridis'},
        title='{} path'.format(descr),
        return_data=True,
        show=False)
    data.to_csv('{SAVE}/paga_path.csv'.format(SAVE=SAVE2))

pl.savefig('{SAVE}/paga_path_paul15.pdf'.format(SAVE=SAVE2))
pl.show()

SAVE2 = 'F:/scanpy/ALL2/myeloid/macro'

adata_subset2 = sc.read("{SAVE}/myeloid.h5ad".format(SAVE=SAVE))
adata_subset2.write("{SAVE}/macro.h5ad".format(SAVE=SAVE2))
adata_subset2.obs.to_csv("{SAVE}/meta.csv".format(SAVE=SAVE2))
pd.DataFrame(adata_subset2.obsm['X_diffmap']).to_csv("{SAVE}/scanpy_diffmap.csv".format(SAVE=SAVE2))
adata_subset2.obs["dpt_pseudotime"].to_csv("{CWD}/material/pseudotime.csv".format(CWD=CWD))
filename = "{CWD}/immu_scanpy.loom".format(CWD=CWD)

meta = pd.read_csv("{SAVE}/meta.csv".format(SAVE=SAVE), index_col = 0)
adata_subset2.obs = meta
