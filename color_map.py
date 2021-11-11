
CONFIG = 'G:/scRNA-seq/scanpy/ALL2/configure'

color = pd.read_csv("{CWD}/cell_label_color_map3.csv".format(CWD=CONFIG))

# replace cell type
color_s = color['cell_label'].isin(adata_subset.obs['cell_label'].unique())
color_s = color[color_s]
c = pd.CategoricalDtype(categories=color_s['cell_label'], ordered = True)
s_cat = adata_subset.obs['cell_label'].astype(c)
adata_subset.obs['cell_label'] = s_cat
adata_subset.uns['cell_label_colors'] = color_s['cell_color'].unique()
sc.pl.umap(adata_subset, color = ['cell_label'], frameon = False, title = "")

# replace sample_type
color_s = color['cell_label'].isin(adata_subset.obs['sample_type'].unique())
color_s = color[color_s]
c = pd.CategoricalDtype(categories=color_s['cell_label'], ordered = True)
s_cat = adata_subset.obs['sample_type'].astype(c)
adata_subset.obs['sample_type'] = s_cat
adata_subset.uns['sample_type_colors'] = color_s['cell_color'].unique()
sc.pl.umap(adata_subset, color = ['sample_type'], frameon = False, title = "")

# replace patient_id
color_s = color['cell_label'].isin(adata_subset.obs['patient_id'].unique())
color_s = color[color_s]
c = pd.CategoricalDtype(categories=color_s['cell_label'], ordered = True)
s_cat = adata_subset.obs['patient_id'].astype(c)
adata_subset.obs['patient_id'] = s_cat
adata_subset.uns['patient_id_colors'] = color_s['cell_color'].unique()
sc.pl.umap(adata_subset, color = ['patient_id'], frameon = False, size = 3, title = "")

# replace her2_statu
color_s = color['cell_label'].isin(adata_subset.obs['Her2_statu'].unique())
color_s = color[color_s]
c = pd.CategoricalDtype(categories=color_s['cell_label'], ordered = True)
s_cat = adata_subset.obs['Her2_statu'].astype(c)
adata_subset.obs['Her2_statu'] = s_cat
adata_subset.uns['Her2_statu_colors'] = color_s['cell_color'].unique()
sc.pl.umap(adata_subset, color = ['Her2_statu'], frameon = False, size = 3, title = "")


# replace main_label
color_s = color['cell_label'].isin(adata_subset.obs['main_label'].unique())
color_s = color[color_s]
c = pd.CategoricalDtype(categories=color_s['cell_label'], ordered = True)
s_cat = adata_subset.obs['main_label'].astype(c)
adata_subset.obs['main_label'] = s_cat
adata_subset.uns['main_label_colors'] = color_s['cell_color'].unique()
sc.pl.umap(adata_subset, color = ['main_label'], frameon = False, size = 3, title = "")
#设置matlab字体
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif']=['Arial']

#更改pdf为可修改模式
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['figure.dpi'] = 350
matplotlib.rcParams.update({'font.size': 15})
matplotlib.rcParams['figure.figsize'] = 6, 6
