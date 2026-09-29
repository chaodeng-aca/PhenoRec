
import argparse

import pandas as pd
import scanpy as sc

from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection
import numpy as np
import anndata as ad
from sklearn.neighbors import NearestNeighbors


def benchmark_one(adata_path: str, embedding_obsm_keys, batch_key: str, label_key: str, n_jobs: int = 15):
    adata = sc.read(adata_path)

    bm = Benchmarker(
        adata,
        batch_key=batch_key,
        label_key=label_key,
        bio_conservation_metrics=BioConservation(
            isolated_labels=True,
            nmi_ari_cluster_labels_leiden=False,
            nmi_ari_cluster_labels_kmeans=True,
            silhouette_label=False,
            clisi_knn=False,
        ),
        batch_correction_metrics=BatchCorrection(
            silhouette_batch=False,
            ilisi_knn=True,
            kbet_per_label=True,
            graph_connectivity=False,
            pcr_comparison=False,
        ),
        pre_integrated_embedding_obsm_key="pca",
        embedding_obsm_keys=embedding_obsm_keys,
        n_jobs=n_jobs,
    )

    bm.benchmark()

    
    df = bm.get_results(min_max_scale=False).transpose()

    
    df.insert(0, "adata", adata_path)
    df.insert(1, "metrics", df.index.astype(str))

    return df.reset_index(drop=True)

  
def average(df_all, embed_cols):

  
    order_df = df_all[["metrics", "Metric Type"]].drop_duplicates()

    df_mean = (
        df_all.groupby(["metrics", "Metric Type"], sort=False, as_index=False)[embed_cols]
        .mean()
    )

   
    df_mean = df_mean.merge(
        order_df.assign(_order=range(len(order_df))),
        on=["metrics", "Metric Type"],
        how="left"
    ).sort_values("_order").drop(columns="_order")

    df_mean.insert(0, "adata", "mean")
    
    return df_mean



def fit_knn(mat_train, mat_holdout, n_neighbors, algorithm = 'kd_tree'):
    
    # fit knn using mat_train
    # return nn indices and distances in train set for holdout set
    knn = NearestNeighbors(n_neighbors = n_neighbors, algorithm = algorithm).fit(mat_train)
    distances, indices = knn.kneighbors(mat_holdout)
    indices = indices[:,1:]
    distances = distances[:,1:]

    return indices, distances


def calc_knn_prop(knn_indices, labels_train, label_categories):

    # knn_indices: shape = (n_holdout_samples, (knn-1)), np.array
    # labels_train: shape = (n_train_samples, ), pd.object
    # label_categories: shape = (n_label_categories, ), np.array
    n = knn_indices.shape[0]
    n_category = label_categories.shape[0]
    nn_prop = np.zeros(shape = (n, n_category))

    for i in range(n):
        knn_labels = labels_train[knn_indices[i,]]
        for k in range(n_category):
            nn_prop[i, k] = sum(knn_labels == label_categories[k]) 

    nn_prop = nn_prop / knn_indices.shape[1]
    return nn_prop


def calc_oobNN(adata_orig, batch_key, condition_key, dim_key, n_neighbors=15):
    ''' Compute out-of-batch k-nearest-neighbor composition 

    Parameters
    ----------
    adata_orig : anndata object
        Expression data stored in adata_orig.X, based on which to compute out of batch nearest neighbors.
    batch_key : str
        Variable name indicating batch. Should be a column name of adata.obs.
    condition_key : str
        Variable name indicating condition. Should be a column name of adata.obs. 
        We compute out-of-batch proportion of each condition level within each cell's neighborhood.
    n_neighbors : int, optional
        Number of k-nearest neighbors.
    
    Returns
    ----------
    res: anndata object
        One new attribute added. 
        res.obsm['knn_prop'] : pd.DataFrame, out-of-batch k-nearest-neighbor composition, cell-by-condition
    '''

    np.random.seed(123)
    
    list_holdout = []
    for holdout_idx in np.unique(adata_orig.obs[batch_key]):

        adata_train = adata_orig[~adata_orig.obs[batch_key].isin([holdout_idx])]
        adata_holdout = adata_orig[adata_orig.obs[batch_key].isin([holdout_idx])]
        num_cells = adata_train.obs[condition_key].value_counts().min()
    
        a_list = []
        for x in np.unique(adata_train.obs[condition_key]):
            a1 = adata_train[adata_train.obs[condition_key].isin([x])]
            random_indices = np.random.choice(a1.shape[0], size=num_cells, replace = False)
            a1 = a1[random_indices,:]
            a_list.append(a1)
 
        adata_train = ad.concat(a_list)
        adata = ad.concat([adata_train, adata_holdout])
    
        mat = adata.obsm[dim_key]
        mat_train = mat[~adata.obs[batch_key].isin([holdout_idx]),]
        mat_holdout = mat[adata.obs[batch_key].isin([holdout_idx]),]
    
        # fit knn
        indices, distances = fit_knn(mat_train=mat_train, mat_holdout=mat_holdout, n_neighbors=n_neighbors, algorithm = 'kd_tree')

        # compute proprotion
        labels_train = adata_train.obs[condition_key].astype('object')
        label_categories = np.unique(labels_train)
        result = calc_knn_prop(indices, labels_train, label_categories)
        knn_df = pd.DataFrame(data=result, 
                              index =  adata_holdout.obs_names,
                              columns = label_categories)
        adata_holdout.obsm['knn_prop'] = knn_df
        list_holdout.append(adata_holdout)

    res = ad.concat(list_holdout)
    
    df = res.obsm['knn_prop']
    df['condition'] = res.obs[condition_key]
    df = df.reset_index()
    df = pd.melt(df, id_vars=['index', 'condition'], var_name='neighbor', value_name='proportion')
    df = df.rename(columns={'index': 'cell'})

    return df


def calc_oobNN_for_multidims(adata_orig, batch_key, condition_key, dim_keys, n_neighbors=30):
    '''
    Compute out-of-batch k-nearest-neighbor composition for multiple dimensionality reductions.
    Parameters:
    ----------
    adata_orig : anndata object
        Expression data stored in adata_orig.X, based on which to compute out of batch nearest neighbors.
    batch_key : str
        Variable name indicating batch. Should be a column name of adata.obs.
    condition_key : str
        Variable name indicating condition. Should be a column name of adata.obs. 
        We compute out-of-batch proportion of each condition level within each cell's neighborhood.
    dim_keys : list of str
        List of variable names indicating the keys for the different dimensionality reductions in adata.obsm.
    n_neighbors : int, optional
        Number of k-nearest neighbors.

    Returns
    ----------
    result_df : pd.DataFrame
        Dataframe with concatenated out-of-batch k-nearest-neighbor composition for each dimension reduction.
    '''

    result_list = []

    # Loop through all provided dimensionality reductions
    for dim_key in dim_keys:
        print(f"Processing dimensionality space: {dim_key}")
        result_df = calc_oobNN(adata_orig, batch_key, condition_key, dim_key, n_neighbors)

        # Rename the 'proportion' column to include the dim_key for distinction
        result_df = result_df.rename(columns={'proportion': f'proportion_{dim_key}'})
        result_list.append(result_df)
        
    cell_order = result_list[0]['cell']
    for df in result_list:
        if not df['cell'].equals(cell_order):
            raise ValueError("The 'cell' column order is not consistent across dataframes!")
            
    # Merge all dataframes on 'cell', 'condition', and 'neighbor' columns
    final_result_df = result_list[0]
    for df in result_list[1:]:
        final_result_df = pd.merge(final_result_df, df, on=['cell', 'condition', 'neighbor'], how='outer')

        
    return final_result_df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--adata",
        required=True,
        nargs="+",
        help="One or more .h5ad files",
    )
    ap.add_argument(
        "--embeddings",
        required=True,
        nargs="+",
        help="Embedding keys in adata.obsm (e.g. pca harmony denoise)",
    )
    ap.add_argument(
        "--batch_key",
        required=True,
        help="Batch key in adata.obs",
    )
    ap.add_argument(
        "--label_key",
        required=True,
        help="Label key in adata.obs",
    )
    ap.add_argument(
        "--out_csv",
        required=True,
        help="Output CSV path",
    )
    ap.add_argument(
        "--n_jobs",
        type=int,
        default=15,
    )
    ap.add_argument("--adata_whole", required=False, default=None, help="Optional whole .h5ad file for oobnn only")
    ap.add_argument("--pheno_key", required=False, default=None, help="Pheno key in adata.obs")
    
    ap.add_argument("--calc_oobNN", action="store_true", help="Flag to compute out-of-batch nearest neighbors for multiple dimensions")

    args = ap.parse_args()

    all_dfs = []
    for adata_path in args.adata:
        df = benchmark_one(
            adata_path=adata_path,
            embedding_obsm_keys=args.embeddings,
            batch_key=args.batch_key,
            label_key=args.label_key,
            n_jobs=args.n_jobs,
        )
        all_dfs.append(df)

    df_all = pd.concat(all_dfs, ignore_index=True)
    embed_cols = args.embeddings
    df_all[embed_cols] = df_all[embed_cols].apply(pd.to_numeric, errors="coerce")  
    
    df_mean = average(df_all, embed_cols)
    
    df_out = pd.concat([df_all, df_mean], ignore_index=True)
    
    out_csv = args.out_csv
    df_out.to_csv(out_csv, index=False)
    
    print(f"Saved results to: {out_csv}")
    
    if args.calc_oobNN:
        print("Calculating out-of-batch nearest neighbors for multiple dimensions...")
        adata_whole = sc.read(args.adata_whole)
        result_oobNN = calc_oobNN_for_multidims(
            adata_orig=adata_whole,
            batch_key=args.batch_key,
            condition_key=args.pheno_key,  
            dim_keys=args.embeddings,
            n_neighbors=30,
        )
        out_oobNN_csv = args.out_csv.replace(".csv", "_oobNN.csv")
        result_oobNN.to_csv(out_oobNN_csv, index=False)
        print(f"Saved oobNN results to: {out_oobNN_csv}")


if __name__ == "__main__":
    main()
