import pyvdj
import os

TCR = 'F:/TCR/vdj_files'
sample_name = os.listdir(TCR)
paths = [TCR+ '/'+i + '/filtered_contig_annotations.csv' for i in sample_name]
samples = dict(zip(paths, sample_name))

adata_subset = pyvdj.load_vdj(samples, adata_subset, obs_col='vdj_obs', cellranger=3)

adata = pyvdj.load_vdj(samples, adata, obs_col='vdj_obs', cellranger=3)

adata_subset = pyvdj.add_obs(adata_subset, obs=['clonotype'])
adata = pyvdj.add_obs(adata, obs=['clonotype', 'is_clone'])

adata_subset = pyvdj.add_obs(adata_subset, obs=['is_clone'])
adata_subset = pyvdj.add_obs(adata_subset, obs=['clone_count'])


sc.pl.umap(adata_subset, color=['vdj_has_vdjdata', 'vdj_is_clone'])

sc.pl.umap(adata_subset, color=['vdj_all_productive', 'vdj_any_productive'])

sc.pl.umap(adata_subset, color=['vdj_chain_TRA', 'vdj_chain_TRB'])

sc.pl.umap(adata_subset, color='vdj_clone_count')


#diversity index
#Diversity index of all clonotypes, for each metadata category, for each 10x sample (channel)
obs_col = adata_subset.uns['pyvdj']['obs_col']

vdjdf = adata_subset.uns['pyvdj']['df']
vdjdf.columns

shanno_list = []
vdjdf = vdjdf.loc[vdjdf['barcode_meta'].isin(adata_subset.obs[obs_col])]  # optional step: subset for cells in anndata
meta_cat = 'cell_label'
for m in adata_subset.obs[meta_cat].unique():
    m_cells = adata_subset.obs.loc[(adata_subset.obs[meta_cat] == m)][obs_col]
    vdjdf_m = vdjdf.loc[vdjdf['barcode_meta'].isin(m_cells)]
    for s in vdjdf_m['sample'].unique():
        sd = dict(vdjdf_m.loc[(vdjdf_m['sample'] == s)].groupby('clonotype_meta').size())
        shannon_index = pyvdj.shannon(sd)
        simpson_index = pyvdj.simpson(sd)
        shanno_list.append([s, m, shannon_index, simpson_index])
        print(s)
        print('The Shannon index for sample %s in category %s is %s.' % (s, m, shannon_index))
        print('The Simpson index for sample %s in category %s is %s.' % (s, m, simpson_index))
        print()
        print()

adat = pd.DataFrame(shanno_list)
adat.rename(columns = {0:'orig.ident',1:'cell_label',2:'shannon_index', 3:'simpson_index'}, inplace = True)


patient = adata_subset.obs.patient_id.unique()

sample_donor = dict(zip(adata_subset.uns['pyvdj']['df']['sample'].unique(), list(patient) + list(patient)))
adata_subset = pyvdj.find_clones(adata_subset, sample_donor)

shanno_list = []
for s in adata_subset.obs['patient_id'].unique():
    sd = dict(adata_subset.obs.loc[adata_subset.obs['patient_id'] == s]['vdj_donor_clonotype'].value_counts())
    shannon_index = pyvdj.shannon(sd)
    simpson_index = pyvdj.simpson(sd)
    shanno_list.append([s, shannon_index, simpson_index])
    print(s)
    print('The Shannon index for %s is %s.' % (s, shannon_index))
    print('The Simpson index for %s is %s.' % (s, simpson_index))
    print()
    print()


for s in vdjdf['sample'].unique():
    sd = dict(vdjdf.loc[(vdjdf['sample'] == s)].groupby('clonotype_meta').size())
    shannon_index = pyvdj.shannon(sd)
    simpson_index = pyvdj.simpson(sd)
    print(s)
    print('The Shannon index for %s is %s.' % (s, shannon_index))
    print('The Simpson index for %s is %s.' % (s, simpson_index))
    print()
    print()


#finding public clonotype

meta = 'cell_label'
adata_subset = pyvdj.stats(adata_subset, meta)
cdr3 = adata_subset.uns['pyvdj']['stats'][meta]['cdr3']

set(cdr3['C']) & set(cdr3['P1']) & set(cdr3['P2'])

set(cdr3['P1']) & set(cdr3['P2']) - set(cdr3['C'])

set(cdr3['C']) & set(cdr3['P1'])

order = ['CA1','LN1','CA2','LN2','CA3','LN3','CA4','LN4',
        'CA5','LN5','CA6','LN6','CA7','LN7','CA8','LN8']

cdr3_codes_rev = adata_subset.uns['pyvdj']['cdr3']['cdr3_codes_rev']

vdjdf = vdjdf.loc[vdjdf['barcode_meta'].isin(adata_subset.obs[obs_col])]
cdr3_coded = {}
for key, value in cdr3.items():
    print(key)
    cdr3_coded[key] = [cdr3_codes_rev[cdr3] for cdr3 in value]


cdr3_coded.keys()

order = ['CD8-C1', 'CD8-C2','CD8-C3', 'CD8-C4', 'CD8-C5', 'CD8-C6']

order = ['CD4-C1', 'CD4-C2','CD4-C3', 'CD4-C4', 'CD4-C5']
cdr3_coded_ord = {}
for key in order:
    cdr3_coded_ord[key] = cdr3_coded[key]


contents = upsetplot.from_contents(cdr3_coded_ord)
contents

upsetplot.plot(contents, show_counts = '%d', sort_by = 'cardinality')
plt.show()

upsetplot.UpSet(contents, subset_size='count')



#get cell
meta = 'orig.ident'
vdjdf = adata_subset.uns['pyvdj']['df']
cdr3_simple_dict = {}
for m in adata_subset.obs[meta].unique():
    print(m)
    cells = adata_subset.obs.loc[adata_subset.obs[meta] == m][obs_col].tolist()
    cdr3_m = vdjdf.loc[vdjdf['barcode_meta'].isin(cells)]['cdr3']
    cdr3_simple_dict[m] = set(cdr3_m)


n = 2  # clonotypes with at least 2 clones are considered expanded
adata_subset.obs['vdj_expanded_n_2'] = adata_subset.obs['vdj_clone_count']
adata_subset.obs['vdj_expanded_n_2'][adata_subset.obs['vdj_expanded_n_2'] < n] = 0
adata_subset.obs['vdj_expanded_n_2'][adata_subset.obs['vdj_expanded_n_2'] >= n] = 1


obs_vdj_df = adata_subset.obs[adata_subset.obs['vdj_has_vdjdata'] == True]
exp_counts = obs_vdj_df.groupby(['cell_label']).vdj_expanded_n_2.value_counts()
exp_counts.to_csv('expanded_proportions.csv', header=True)
