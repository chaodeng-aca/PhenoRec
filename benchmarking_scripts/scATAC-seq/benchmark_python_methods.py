import matplotlib.pyplot as plt
import scanpy as sc
import scvi
import torch

import anndata as ad
import pandas as pd
import numpy as np

import argparse

def PeakVI(adata, batch_id, out_csv):
    scvi.settings.seed = 0
    print("Last run with scvi-tools version:", scvi.__version__)
    
    adata.X = adata.X.toarray()
    scvi.model.PEAKVI.setup_anndata(adata, batch_key=batch_id)
    model = scvi.model.PEAKVI(adata, n_latent=50)
    model.train()
    latent = model.get_latent_representation()
    
    res_scvi = pd.DataFrame(latent)
    res_scvi.index = adata.obs.index
    res_scvi.to_csv(out_csv, index=True)
    
    

def python_benchmark(adata_path, batch_id, out_csv,
                     run_PeakVI=True):
    
    adata = ad.read_h5ad(adata_path)
    
    if run_PeakVI:
        PeakVI(adata.copy(), batch_id, out_csv = f"{out_csv}_PeakVI.csv")

    
def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--adata_path', type=str, required=True, help="AnnData")
    parser.add_argument('--batch_id', type=str, required=True, help="batch")
    parser.add_argument('--out_csv', type=str, required=True, help="Path")
    parser.add_argument('--run_PeakVI', action='store_true', help="run_scVI")
    args = parser.parse_args()
    
    python_benchmark(args.adata_path, args.batch_id, args.out_csv, 
                     args.run_PeakVI)

if __name__ == "__main__":
    main()
