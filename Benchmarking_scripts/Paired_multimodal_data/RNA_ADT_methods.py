import scvi
import scanpy as sc
import numpy as np
import pandas as pd
import muon
import multigrate as mtg
import anndata
import argparse

import h5py
from scipy import sparse
import sys
import subprocess


def Multigrate(rna, adt, batch_id, out_csv):
    
    rna.layers['counts'] = rna.X.copy()
    sc.pp.normalize_total(rna, target_sum=1e4)
    sc.pp.log1p(rna)
    sc.pp.highly_variable_genes(rna, n_top_genes=3000, flavor='seurat', batch_key=batch_id)
    rna_hvg = rna[:, rna.var.highly_variable].copy()
    rna_hvg
    
    
    # adt preprocess
    adt.layers['counts'] = adt.X.copy()
    muon.prot.pp.clr(adt)
    adt.layers["clr"] = adt.X.copy()
    
    
    adata = mtg.data.organize_multimodal_anndatas(
        adatas = [[rna_hvg], [adt]],           # a list of anndata objects per modality, RNA-seq always goes first
        layers = [['counts'], ['clr']],    # if need to use data from .layers, if None use .X
    )
    
    # Setup anndata
    mtg.model.MultiVAE.setup_anndata(
        adata,
        rna_indices_end=3000, # how many features in the rna-seq modality
        categorical_covariate_keys=[batch_id]
    )
    
    # Initialize Model
    model = mtg.model.MultiVAE(
        adata,
        losses=['nb', 'mse']
    )
    
    # Train model
    model.train(max_epochs=100)
    model.get_model_output()
    
    latent = pd.DataFrame(data=adata.obsm['X_multigrate'],
                          index=adata.obs_names)
    
    latent.to_csv(out_csv, index=True)

    
def totalVI(rna, adt, batch_id, out_csv):
    
    adata = rna
    adata.layers["counts"] = adata.X.copy()
    adata.obsm["protein_expression"]=adt.X.todense()
    adata.uns['protein_names']=adt.var_names
    
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.raw = adata
    
    sc.pp.highly_variable_genes(adata,n_top_genes=3000,flavor="seurat_v3",batch_key=batch_id,layer="counts")
    scvi.model.TOTALVI.setup_anndata(adata,protein_expression_obsm_key="protein_expression",layer="counts",batch_key=batch_id)
    
    # Prepare and run model
    model = scvi.model.TOTALVI(adata, latent_distribution="normal")
    model.train()
    
    latent = model.get_latent_representation() 
    latent = pd.DataFrame(latent, index=adata.obs_names)
    latent.to_csv(out_csv, index=True)
    
def scMDC(rna, adt, batch_id, out_csv):
    
    N_CLUSTERS = 10 
    batch_category = pd.Categorical(rna.obs[batch_id])
    batch = batch_category.codes.astype(np.int32)
    n_batch = len(batch_category.categories)
   
    h5_file = "scMDC_output/scMDC_input.h5"
    
    with h5py.File(h5_file, "w") as h5f:
        rna_x = rna.X.toarray() if sparse.issparse(rna.X) else np.asarray(rna.X)
        h5f.create_dataset(
            "X1", data=rna_x, dtype="f4", compression="gzip"
        )
        del rna_x
        
        adt_x = adt.X.toarray() if sparse.issparse(adt.X) else np.asarray(adt.X)
        h5f.create_dataset(
            "X2", data=adt_x, dtype="f4", compression="gzip"
        )
        del adt_x
        
        h5f.create_dataset(
            "Batch", data=batch, compression="gzip"
        )
    
    # 运行多批次版 scMDC
    cmd = [
        sys.executable,
        "-u",
        "/dssg/home/acct-clswt/clswt-chaodeng/phenotype_project/Multiome/benchmark_method/scMDC-master/src/run_scMDC_batch.py",
        "--device", "cpu",
        "--data_file", str(h5_file),
        "--save_dir", 'scMDC_output/',
        "--n_clusters", str(N_CLUSTERS),
        "--nbatch", str(n_batch),
        "--no_labels",
        "--embedding_file",
        "--filter1",
        "--filter2",
        "--f2", "10000",
        "-el", "256", "128", "64",
        "-dl1", "64", "128", "256",
        "-dl2", "64", "128", "256",
        "--phi1", "0.005",
        "--phi2", "0.005",
        "--sigma2", "2.5",
        "--tau", "0.1",
        "--pretrain_epochs", "200",
    ]
    subprocess.run(cmd, cwd='./', check=True)
    
    kept_idx = np.loadtxt("scMDC_output/kept_cell_indices.csv",dtype=int,delimiter=",",)
    latent = pd.read_csv("scMDC_output/1_embedding.csv", header=None)
    assert latent.shape[0] == len(kept_idx)
    
    kept_obs_names = rna.obs_names[kept_idx]
    latent.index = kept_obs_names
    
    latent.to_csv(out_csv, index=True)
    

def python_benchmark(rna_path, adt_path, batch_id, out_csv,
                     run_Multigrate=True, run_totalVI=True, run_scMDC=True):
    
        
    rna = sc.read_h5ad(rna_path)
    adt = sc.read_h5ad(adt_path)
    
    scvi.settings.seed = 0
    
    if run_Multigrate:
        Multigrate(rna=rna.copy(), adt=adt.copy(), batch_id=batch_id, out_csv=f"{out_csv}_Multigrate.csv")
    
    if run_totalVI:
        totalVI(rna=rna.copy(), adt=adt.copy(), batch_id=batch_id, out_csv=f"{out_csv}_totalVI.csv")
        
    if run_scMDC:    
        scMDC(rna=rna.copy(), adt=adt.copy(), batch_id=batch_id, out_csv=f"{out_csv}_scMDC.csv")
    

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--rna_path', type=str, required=True, help="AnnData")
    parser.add_argument('--adt_path', type=str, required=True, help="AnnData")
    parser.add_argument('--batch_id', type=str, required=True, help="batch")
    parser.add_argument('--out_csv', type=str, required=True, help="Path")
    parser.add_argument('--run_Multigrate', action='store_true', help="Multigrate")
    parser.add_argument('--run_totalVI', action='store_true', help="MultiVI")
    parser.add_argument('--run_scMDC', action='store_true', help="scMDC")
    args = parser.parse_args()
    
    python_benchmark(args.rna_path, args.adt_path, args.batch_id, args.out_csv, 
                     args.run_Multigrate, args.run_totalVI, args.run_scMDC)

if __name__ == "__main__":
    main()
    
    
    
