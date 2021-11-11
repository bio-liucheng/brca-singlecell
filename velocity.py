import scvelo as scv
VEL = 'G:/scRNA-seq/single_cell/velocity/loom_file'
adata_vlc = scv.read_loom("{VEL}/combine_all.loom".format(VEL=VEL))

#change obs_names
obs_names = adata_vlc.obs_names

obs_names = obs_names.str.replace('x', '')

obs_names = obs_names.str.replace('_trans5_outs:', '_')

adata_vlc.obs_names = obs_names


#subset
cell_select = adata_subset.obs_names
adata_vlc_select = adata_vlc[cell_select]
adata_vlc_select.obs = adata_subset.obs

#scevlocity

scv.pp.filter_genes(adata_vlc_select, min_shared_counts=10)
scv.pp.normalize_per_cell(adata_vlc_select)
scv.pp.filter_genes_dispersion(adata_vlc_select, n_top_genes=3000)
#scv.pp.filter_and_normalize(adata_vlc_select, min_shared_counts=30, n_top_genes=3000)
scv.pp.log1p(adata_vlc_select)

scv.pp.moments(adata_vlc_select, n_pcs=30, n_neighbors=20)
scv.tl.velocity(adata_vlc_select)
scv.tl.velocity_graph(adata_vlc_select)


#
adata_vlc_select.obsm['X_umap'] = adata_subset.obsm['X_umap']
adata_vlc_select.obsm['X_diffmap'] = adata_subset.obsm['X_diffmap']
#plot velocity
scv.pl.velocity_embedding_stream(adata_vlc_select, basis='umap', color = ['cell_label'])
scv.pl.velocity_embedding_grid(adata_vlc_select, basis='umap', alpha = 0.7, size =20,
        color = ['cell_label'], arrow_length=1.2, title = "",
        legend_fontsize = 8,
        legend_loc = 'right margin',frameon=False, dpi=150)
scv.pl.velocity(adata_vlc_select, ['PLA2G2A'])

adata_vlc_select.uns['cell_label_colors'] = adata_subset.uns['cell_label_colors']
scv.pl.velocity_embedding(adata_vlc_select, color = ['cell_label','patient_id'], basis='umap', arrow_length=1.2, arrow_size=1.2, dpi=150)


adata_vlc_select.write("{SAVE}/tumor_velocity.h5ad".format(SAVE=SAVE))
#next plot
scv.tl.recover_dynamics(adata_vlc_select)
scv.tl.velocity(adata_vlc_select, mode='dynamical')
scv.tl.velocity_graph(adata_vlc_select)

scv.tl.recover_latent_time(adata_vlc_select)

scv.pl.scatter(adata_vlc_select, color='latent_time', fontsize=24, size=100,
               color_map='gnuplot', perc=[2, 98], colorbar=True, rescale_color=[0,1])
top_genes = adata_vlc_select.var_names[adata_vlc_select,.var.fit_likelihood.argsort()[::-1]][:300]

top_genes = adata_vlc_select.var_names[adata_vlc_select.var.fit_likelihood.argsort()[::-1]][:300]
scv.pl.heatmap(adata_vlc_select, var_names=top_genes, tkey='latent_time', n_convolve=100, col_color='cell_label')

Copy to clipboard
scv.pl.scatter(adata, basis=top_genes[:10], legend_loc='none',
               size=80, frameon=False, ncols=5, fontsize=20)





adata_vlc_select.uns['neighbors']['distances'] = adata_subset.obsp['distances']
adata_vlc_select.uns['neighbors']['connectivities'] = adata_subset.obsp['connectivities']
adata_vlc_select.uns['diffmap_evals'] = adata_subset.uns['diffmap_evals']

scv.tl.paga(adata_vlc_select, groups='cell_label')
df = scv.get_df(adata_vlc_select, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')
scv.pl.paga(adata_vlc_select, basis='umap', size=50, alpha=0,
            min_edge_width=2, node_size_scale=1.5, threshold_root_end_prior = 0.1, minimum_spanning_tree = False)
