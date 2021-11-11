
CWD = 'F:/zhang zemin/GSE146771_scanpy'
CWD2 = 'F:/scCancer_analysis/soup_seurat'

# load gene names
adata_ref = sc.read("{CWD}/ALL_scanpy_BBKNN.h5ad".format(CWD=CWD2), cache = True)
adata = sc.read("{CWD}/seurat_obj_immu.mtx".format(CWD=CWD), cache = True).T

adata.var_names = pd.read_csv("{CWD}/immu_genenames.csv".format(CWD=CWD)).iloc[:, 1]

# load cell names
adata.obs_names = pd.read_csv("{CWD}/immu_cellnames.csv".format(CWD=CWD)).iloc[:, 1]


metadata = pd.read_csv("{CWD}/meta.csv".format(CWD=CWD2), index_col = 0)
adata_ref.obs = metadata

var_names = adata_ref.var_names.intersection(adata.var_names)
adata_ref = adata_ref[:, var_names]
adata = adata[:, var_names]

adata

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.tl.ingest(adata, adata_ref, obs='cell_label')
