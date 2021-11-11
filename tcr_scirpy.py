import scirpy as ir
import pandas as pd
import numpy as np
import scanpy as sc
import os
from matplotlib import pyplot as plt
TCR = 'G:/scRNA-seq/TCR/vdj_files'
sample_name = os.listdir(TCR)
paths = [TCR+ '/'+i + '/all_contig_annotations.json' for i in sample_name]

adatas = [ir.io.read_10x_vdj(filename) for filename in paths]

for i in range(len(adatas)):
    obs_names = adatas[i].obs_names
    obs_names = obs_names.str.replace('-1', '')
    obs_names = [ sample_name[i]+ '_' + j for j in obs_names]
    adatas[i].obs_names = obs_names

adatas_list = [adatas[i].obs for i in range(len(adatas))]

adata = pd.concat(adatas_list)
adata_all = sc.AnnData(obs = adata)
#adata = adatas[0].concatenate(adatas[1:])
ir.pp.merge_with_tcr(adata_subset, adata_all)


ir.tl.chain_pairing(adata_subset)
ir.tl.clonal_expansion(adata_subset)
ir.pp.tcr_neighbors(adata_subset, receptor_arms="all", dual_tcr="primary_only")
ir.tl.define_clonotypes(adata_subset)

ir.pl.group_abundance(adata_subset, groupby="chain_pairing", target_col="sample_type",)

#we visualize the Multichain cells on the UMAP plot and exclude them from downstream analysis
sc.pl.umap(adata_subset, color="multi_chain")

#Compute CDR3 neighborhood graph


ir.pp.tcr_neighbors(
    adata_subset,
    metric="alignment",
    sequence="aa",
    cutoff=15,
    receptor_arms="all",
    dual_tcr="all",
)
ir.tl.define_clonotype_clusters(
    adata_subset, partitions="connected", sequence="aa", metric="alignment"
)

ir.pp.tcr_neighbors(adata_subset, receptor_arms="all", dual_tcr="primary_only", cutoff=0)
ir.tl.define_clonotypes(adata_subset)

ir.tl.clonotype_network(adata_subset, min_size=10)

ir.pl.clonotype_network(adata_subset,
                        color=["clonotype", "patient_id"],
                        edges=False, size=50, ncols=2,
                        legend_fontoutline=2,
                        legend_loc=["on data", "right margin"])

ir.pl.clonotype_network(adata_subset,
                        color=["clonotype", "clonotype_nt"],
                        edges=False, size=50, ncols=2,
                        legend_fontoutline=2,
                        legend_loc=["on data", "none"])

nt_in_aa = adata_subset.obs.groupby(["clonotype_nt", "clonotype"]).size().reset_index().groupby(["clonotype"]).size().reset_index()
convergent_clonotypes = nt_in_aa.loc[nt_in_aa[0] > 1, "clonotype"]

nt_in_orig = adata_subset.obs.groupby(["clonotype_orig", "clonotype_nt"]).size().reset_index().groupby("clonotype_orig").size().reset_index()
nt_in_orig[nt_in_orig[0] > 1]


sc.pl.umap(adata_subset, color="is_convergent", groups=["convergent"],
    size=[10 if x == "convergent" else 3 for x in adata_subset.obs["is_convergent"]])


ir.pl.vdj_usage(adata_subset, full_combination=False, top_n=30)

ir.pl.spectratype(adata_subset, target_col="cell_label", viztype="bar", fig_kws={"dpi": 120})


df, dst, lk = ir.tl.repertoire_overlap(adata_subset, 'orig.ident', inplace=False)

ir.pl.repertoire_overlap(adata_subset, 'orig.ident')

ir.pl.repertoire_overlap(adata_subset, 'orig.ident', heatmap_cats=['patient_id', 'sample_type'])

ir.pl.repertoire_overlap(adata_subset, 'orig.ident', dendro_only=True, heatmap_cats=['sample_type'])

ir.pl.repertoire_overlap(adata_subset, 'orig.ident', pair_to_plot=['CA8', 'LN8'])

ir.tl.clonotype_imbalance(adata_subset, replicate_col='orig.ident', groupby='sample_type', case_label='Tumor')
ir.pl.clonotype_imbalance(adata_subset, replicate_col='orig.ident', groupby='sample_type', case_label='Tumor', plot_type='strip')


SAVE = 'F:/scanpy/TCR'
adata_subset.write("{SAVE}/immu_tcr_cd48.h5ad".format(SAVE=SAVE))
adata_subset = sc.read("{SAVE}/immu_tcr_cd48.h5ad".format(SAVE=SAVE))
adata_subset.obs.to_csv("{CWD}/ALL_scanpy.csv".format(CWD=SAVE))
