numiter=200 # Number of NMF replicates. Set this to a larger value ~200 for real data. We set this to a relatively low value here for illustration at a faster speed
numworkers=4 # Number of parallel factorization jobs to run. Set this to a value reflective of the number of cores on your computer.
numhvgenes=2000 ## Number of over-dispersed genes to use for running the factorizations
seed = 14

output_directory = 'example_PBMC\\cNMF'
countfn = 'example_PBMC\\counts.h5ad'
run_name = 'pbmc_cNMF'

prepare_cmd = 'python cnmf.py prepare --output-dir %s --name %s -c %s -k %s --n-iter %d --total-workers 1 --seed %d --numgenes %d --beta-loss frobenius' % (output_directory, run_name, countfn, K, numiter, seed, numhvgenes)
os.system(prepare_cmd)
factorize_cmd = 'python cnmf.py factorize --output-dir %s --name %s --worker-index 0' % (output_directory, run_name)
os.system(factorize_cmd)
cmd = 'python cnmf.py combine --output-dir %s --name %s' % (output_directory, run_name)
os.system(cmd)
worker_index = ' '.join([str(x) for x in range(numworkers)])

kselect_plot_cmd = 'python cnmf.py k_selection_plot --output-dir %s --name %s' % (output_directory, run_name)
os.system(kselect_plot_cmd)
Image(filename = "example_PBMC/cNMF2/pbmc_cNMF/pbmc_cNMF.k_selection.png", width=1000, height=1000)
