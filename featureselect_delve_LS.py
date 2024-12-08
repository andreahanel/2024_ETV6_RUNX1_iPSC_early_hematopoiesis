import delve_benchmark
import os
import scanpy as sc
import pandas as pd
import anndata as ad

#-----------------------------------
adata_directory = "/research/groups/sysgen/PROJECTS/sysgen_team/andrea_work/2023_Fischer_iPSC/revision_final/data/" 
adata = sc.read_h5ad(os.path.join(adata_directory, 'adata2_norm_PCAharmony_leiden.h5ad'))

adata.X=adata.layers["log1p_norm"]

#RPE will be the name of dir
feature_directory = os.path.join("/research/groups/sysgen/PROJECTS/sysgen_team/andrea_work/delve_benchmark/outs/", 'RPE', 'predicted_features')
delve_benchmark.pp.make_directory(feature_directory)
#-----------------------------------

n_selected = 4000
trial = 0
fs_method_list = [ [delve_benchmark.tl.laplacian_score_fs, {'k': 10}],
                    [delve_benchmark.tl.neighborhood_variance_fs, {}],
                    [delve_benchmark.tl.mcfs_fs, {'k': 10, 'n_selected_features': 30, 'n_clusters':4}],
                    [delve_benchmark.tl.variance_score, {}],
                    [delve_benchmark.tl.hvg, {'n_top_genes': n_selected, 'log': True}],
                    [delve_benchmark.tl.scmer_fs, {'k': 10, 'n_pcs': 50}],
                    [delve_benchmark.tl.hotspot_fs, {'k': 10, 'n_pcs': 50, 'model':'danb'}],
                   # [delve_benchmark.tl.run_delve_fs, {'num_subsamples': 1000, 'n_clusters': 5, 'k': 10, 'return_modules':True, 'n_random_state': 10, 'random_state': 0}]
]


for fs_method, fs_method_params in fs_method_list:
    fs_name = fs_method.__name__
    fs = delve_benchmark.tl.fs(adata = adata, 
                               fs_method = fs_method, 
                               fs_method_params = fs_method_params)

    if fs_name == 'run_delve_fs':
        delta_mean, modules, ranked_features = fs.select()
        ranked_features = ranked_features[:n_selected]
        ranked_features.to_csv(os.path.join(feature_directory, f'delve_fs_trial{trial}.csv'))
        delta_mean.to_csv(os.path.join(feature_directory, f'delta_mean_trial{trial}.csv'))
        modules.to_csv(os.path.join(feature_directory, f'modules_trial{trial}.csv'))
    elif fs_name == 'all_features':
        ranked_features = fs.select()
        ranked_features = pd.DataFrame(ranked_features)
        ranked_features.to_csv(os.path.join(feature_directory, f'{fs_name}.csv')) #all features
    else:
        ranked_features = fs.select()
        ranked_features = pd.DataFrame(ranked_features[:n_selected])
        ranked_features.to_csv(os.path.join(feature_directory, f'{fs_name}.csv')) #all features
    print(f'{fs_name}: selected {len(ranked_features)} features')