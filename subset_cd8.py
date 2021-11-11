
adata.obs = imu_obs
adata_subset = adata[adata.obs['new_label'].isin(['CD4-C1', 'CD4-C2', 'CD4-C3','CD4-C4', 'CD4-C5', 'CD4-C6'])]
adata_subset = adata[adata.obs['new_label'].isin(['NK cells', 'CD8 T', 'CD4 T'])]
adata_subset = adata[adata.obs['new_label'].isin(['CAFs'])]
adata_subset = adata[adata.obs['new_label'].isin(['B', 'Plasma B'])]
adata_subset = adata[adata.obs['cell_label'].isin(['Epithelial Cells'])]
adata_subset = adata[adata.obs['new_label'].isin(['Myeloid cells'])]
adata_subset = adata_subset[~adata_subset.obs['cell_label'].isin(['Tumor-CDT', 'Tumor-unkown','Tumor-B', 'Tumor-C11'])]
adata_subset = adata[adata.obs['orig.ident'].isin(['CA1', 'CA2', 'CA3', 'CA4', 'CA5',
        'CA6', 'CA7', 'CA8', 'LN1', 'LN3', 'LN4', 'LN5', 'LN6','LN7', 'LN8' ])]
adata_subset = adata[adata.obs['new_label'].isin(['TECs'])]
adata_subset = adata_subset[~adata_subset.obs['leiden'].isin(['CAF-CDT', 'CAF-unkown','CAF-B'])]
adata_subset = adata_subset[~adata_subset.obs['leiden'].isin(['10', '13'])]
adata_subset = adata[adata.obs['cell_label'].str.contains('CD4-|CD8-')]
adata_subset = adata[adata.obs['cell_label'].str.contains('Macro') & ~adata.obs['cell_label'].isin(['Macro-C3-VCAN', 'Macro-C7-SCGB2A2'])]
adata_subset = adata_subset[~adata_subset.obs['cell_label'].isin(['Macro-CDT', 'Macro-B','Macro-C9'])]
adata_subset = adata_subset[~adata_subset.obs['cell_label'].isin(['Macro-C7'])]
adata_subset = adata[adata.obs['cell_label'].str.contains('Macro-')]
adata_subset = adata[adata.obs['cell_label'].isin(['Epithelial Cells'])]
adata_subset = adata[adata.obs['leiden'].isin(['0','9'])]
adata_subset = adata[adata.obs['cell_label'].str.contains('^Macro-') & ~adata.obs['cell_label'].str.contains('SCGB2A2')]
adata_subset = adata[adata.obs['cell_label'].str.contains('CAF-') & ~adata.obs['cell_label'].str.contains('-un')]
adata_subset = adata[adata.obs['cell_label'].str.contains('^CD4-|^CD8-|^CD-|^NK-') & ~adata.obs['cell_label'].str.contains('CD8-C6-IFI44L')]
adata_subset = adata[adata.obs['cell_label'].str.contains('CD4-C6-CXCL13')]
adata_subset = adata[adata.obs['cell_label'].str.contains('-un')]
adata_subset = adata[adata.obs['cell_label'].str.contains('CAF-C1|CAF-C2|CAF-C3')]

adata_subset = adata_subset[adata_subset.obs['cell_label'].str.contains('CAF-')]

adata_subset = adata[adata.obs['cell_label'].str.contains('Macro') & ~adata.obs['cell_label'].str.contains('VCAN')]


adata_subset = adata[adata.obs['sample_type'].str.contains('Tumor')]
#adata = adata[adata.obs_names.isin(cell_name.x)]
adata_subset = adata_subset[adata_subset.obs['sample_type'].isin(['Tumor'])]
adata_subset = adata_subset[adata_subset.obs['patient_id'].isin(['patient_7'])]
adata_subset = adata[adata_subset.obs_names]
adata_subset = adata[meta.index]

sc.pp.normalize_total(adata_subset, target_sum=1e4)
sc.pp.log1p(adata_subset)
adata_subset.raw = adata_subset

#adata.raw = adata
sc.pp.highly_variable_genes(adata_subset, min_mean=0.0125, max_mean=3, min_disp=0.5,n_top_genes =3000)

adata_subset = adata_subset[:, adata_subset.var.highly_variable]
#sc.pp.regress_out(adata_subset, ['nCount_RNA', 'percent.mt'])
sc.pp.scale(adata_subset, max_value=10)
sc.tl.pca(adata_subset)
#batch correction

sc.external.pp.bbknn(adata_subset, batch_key='patient_id', n_pcs = 20)
#harmony
#adata_subset.obs
import harmonypy as hm
ho = hm.run_harmony(adata_subset.obsm['X_pca'], adata_subset.obs, 'patient_id', theta=2)
res = pd.DataFrame(ho.Z_corr)
res = np.array(res.T)
adata_subset.obsm['X_pca'] = res

#mye_anntation = ['LYZ','CD68','C1QA','IGFBPC', 'CXCL9','PTGDS', 'APBEC3A',  'SCGB2A2']
#dc_annotation = ['LAMP3','CD1C', 'CLEC9A', 'CLEC4C']
#t_annotation = ['CD4', 'CD8A', 'GNLY', 'GZMK', 'FOXP3', 'PDCD1', 'CCR7', 'CD69']
#neighbors
sc.pp.neighbors(adata_subset, n_neighbors=20, n_pcs=30)
sc.tl.umap(adata_subset)
#sc.pl.umap(adata_subset, color = ['FCER1G','TYROBP', 'cell_label'], legend_loc = 'on data')
sc.tl.leiden(adata_subset, resolution = 0.6)
sc.pl.umap(adata_subset, color = ['cell_label'])

sc.pl.umap(adata_subset, color = ['chain_pairing'],
    groups=["Extra beta", "Two full chains"], palette = ['#fc8d62','#fc8d62', 'grey'], frameon = False,
    size=[10 if x == "Extra beta" else 3 for x in adata_subset.obs["chain_pairing"]])

'''
color = plt.cm.get_cmap('hsv', 300)
color_list = [color(i) for i in range(300)]


#生成随机颜色 按照clonotype 的数量
import seaborn as sns
#获得tab20 color list
color = sns.color_palette("tab20")
#获得colonal abundance 最多的 colonal id
top_clonotypes = adata.obs.clonotype.value_counts()[:8].index.values.tolist() # A better way might be needed especailly to take normalization into account
top_vgenes = adata.obs.TRB_1_v_gene.value_counts()[:8].index.values.tolist()

j = 0
color_list = []
for i, x in enumerate(adata_subset.obs["clonotype"].unique()):
    if x in top_clonotypes:
        color_list.append(color[j])
        j = j+1
    else:
        color_list.append("grey")

top_clonotypes = adata.obs.clonotype.value_counts()[:8].index.values.tolist() # A better way might be needed especailly to take normalization into account
top_vgenes = adata.obs.TRB_1_v_gene.value_counts()[:8].index.values.tolist()
#上色
sc.pl.umap(adata_subset, color = ['clonotype'],
    groups=top_clonotypes, palette = color_list, frameon = False,
    size=[20 if x in top_clonotypes else 3 for x in adata_subset.obs["clonotype"]])

'''
sc.pl.umap(adata_subset, color = ['chain_pairing'],
    groups=top_clonotypes, frameon = False,
    size=[20 if x in top_clonotypes else 3 for x in adata_subset.obs["clonotype"]])


sc.pl.umap(adata_subset, color = ['cell_label'])
sc.pl.umap(adata_subset, color = ['sample_type'])
sc.pl.umap(adata_subset, color = ['CYP2A6','CNTNAP2', 'CLEC3A', 'CRYAB'], legend_loc = 'on data')
sc.pl.umap(adata_subset, color = ['CD8A', 'CD4', 'GNLY', 'TYROBP'], frameon= False, save = "umap_CDT_marker.png", show = False, title = ["","","",""], size = 5, ncols = 2)



#diffusion mapvalues
sc.pp.neighbors(adata_subset, n_neighbors = 10, n_pcs = 15, knn = True, random_state = 10, method = "gauss")

sc.tl.diffmap(adata_subset, n_comps = 20)
sc.pl.diffmap(adata_subset, color = ['sample_type','cell_label', 'leiden'])

#set cd8 root CD8-C1
root_cell_type = 'CD8-C1-CD8B'
adata_subset.uns['iroot'] = np.flatnonzero(adata_subset.obs['new_label'] == root_cell_type)[0]
#comput dpt
sc.tl.dpt(adata_subset, n_dcs = 20)
sc.pl.diffmap(adata_subset, color = 'leiden', frameon= False, save = "3d_diffmap_cd4_cxcl13_cell_label.png", show = False, title = "", size = 20 )

sc.pl.diffmap(adata_subset, color = ['CCR7','TIGIT','GZMK','sample_type'], frameon= False, save = "diffmap_cd8T_sample_type.png", show = False, title = "", size = 5,ncols = 2)


marker_genes = {'T cell markers': ['CD8A','CD4'],
                'Naive markers': ['SELL', 'CCR7', 'LEF1'],
                'Cytotoxic markers': ['GNLY','NKG7','IFNG','GZMA','GZMB','GZMH','GZMK'],
                'Effector memory': ['S100A4','ANXA1'],
                'Tissue resident': ['CD69','ITGAE', 'CXCR6'],
                'Treg':['IL2RA', 'FOXP3'],
                'Exhausted markers':['PDCD1', 'CTLA4', 'LAG3','TIGIT'],
                'NK cell markers':['TYROBP','FGFBP2', 'XCL1', 'XCL2'],
                'Other':['HOPX','FOS','TRDV2',]
                }


marker_genes = ['CD8A', 'CD4','CD3D', 'CCR7', 'SELL','TCF7', 'HOPX','CCL5',
                'GNLY', 'NKG7','GZMA', 'GZMK', 'HSPA1A', 'CD69', 'CXCR6',
                'CXCL13', 'PDCD1','TIGIT', 'CTLA4','CD200','TOX2', 'MT1G','FOXP3', 'IL2RA',
                'FGFBP2', 'XCL1', 'TRDV2']

marker_genes = ['VSIR', 'PLA2G2A', 'OGN', 'FN1', 'DPP4']

sc.pl.stacked_violin(adata_subset, var_names = marker_genes, groupby = "cell_label",
                    cmap = "Reds", colorbar_title = "", standard_scale = 'var', save = 'violin_T.pdf')

SAVE = 'F:/scanpy/ALL2/T_cell/CD4-C1'

SAVE = 'F:/zhang zemin/GSE146771_scanpy'
SAVE = 'G:/scRNA-seq/scCancer_analysis/soup_seurat/t_cell/cd8'
SAVE = 'G:/scRNA-seq/scCancer_analysis/soup_seurat/t_cell/cd4/cd4-c6'
SAVE = 'G:/scRNA-seq/scCancer_analysis/soup_seurat/myeloid/macrophage'
SAVE = 'G:/scRNA-seq/scCancer_analysis/soup_seurat/tec_caf'
SAVE = 'G:/scRNA-seq/scCancer_analysis/soup_seurat/epi_cell/maliganat'

adata_subset = sc.read("{SAVE}/cd5-c1.h5ad".format(SAVE=SAVE))

adata_subset = sc.read("{SAVE}/scanpy.h5ad".format(SAVE=SAVE))


adata_subset.write("{SAVE}/scanpy.h5ad".format(SAVE=SAVE))
adata_subset.obs.to_csv("{SAVE}/meta.csv".format(SAVE=SAVE))
adata_subset.obs.to_csv("{SAVE}/meta_tcr.csv".format(SAVE=SAVE))
pd.DataFrame(adata_subset.obsm['X_pca']).to_csv("{SAVE}/scanpy_pca.csv".format(SAVE=SAVE))
pd.DataFrame(adata_subset.obsm['X_umap']).to_csv("{SAVE}/scanpy_umap.csv".format(SAVE=SAVE))


filename = "{CWD}/immu_scanpy.loom".format(CWD=CWD)

meta = pd.read_csv("{SAVE}/meta.csv".format(SAVE=SAVE), index_col = 0)
adata_subset.obs = meta
