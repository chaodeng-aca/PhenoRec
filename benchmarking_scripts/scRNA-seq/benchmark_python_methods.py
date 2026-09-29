import scanpy as sc
import scvi
import pandas as pd

import scipy.sparse

import cellanova as cnova
import anndata as ad
import numpy as np

import argparse

import bbknn



def scVI_benchmark(adata,batch_id, out_csv):
    scvi.settings.seed = 0
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key=batch_id)
    model = scvi.model.SCVI(adata, n_layers=2, n_latent=50, gene_likelihood="nb")
    model.train()
    res_scvi = pd.DataFrame(model.get_latent_representation())
    res_scvi.index = adata.obs.index
    res_scvi.to_csv(out_csv, index=True)
    
    
    
def cellANOVA_benchmark(adata,batch_id,condition_id,control_id,out_csv):
    
    np.random.seed(0)
    
    integrate_key = batch_id
    condition_key = condition_id
    control_name = control_id
    
    adata_prep = cnova.model.preprocess_data(adata, integrate_key=integrate_key)
    control_batches = list(set(adata_prep[adata_prep.obs[condition_key]==control_name,].obs[integrate_key]))
    control_dict = {'g1': control_batches,}
    
    adata_prep= cnova.model.calc_ME(adata_prep, integrate_key=integrate_key)
    adata_prep = cnova.model.calc_BE(adata_prep, integrate_key=integrate_key, control_dict=control_dict)
    adata_prep = cnova.model.calc_TE(adata_prep, integrate_key=integrate_key)
    
    integrated = ad.AnnData(adata_prep.layers['denoised'], dtype=np.float32)
    integrated.obs = adata_prep.obs.copy()
    integrated.var_names = adata_prep.var_names
    sc.pp.pca(integrated)
    embed = pd.DataFrame(integrated.obsm['X_pca'])
    embed.index = integrated.obs.index
    embed.to_csv(out_csv, index=True)

def BBKNN_benchmark(adata, batch_id, out_csv):
    sc.pp.pca(adata)
    bbknn.bbknn(adata, batch_key=batch_id)
    sc.tl.umap(adata)
    embed = pd.DataFrame(adata.obsm['X_umap'])
    embed.index = adata.obs.index
    embed.to_csv(out_csv, index=True)
    
def python_benchmark(adata_path,batch_id,condition_id,control_id,cell_type_id,out_csv,
                     run_scVI=True, run_cellANOVA=True, run_BBKNN=True):
    
    
    adata = sc.read(adata_path)
    #adata.X = adata.X.toarray()
    adata.layers["counts"] = adata.X.copy()
    
    if run_scVI:
        scVI_benchmark(adata.copy(), batch_id, out_csv = f"{out_csv}_scVI.csv")
    
    if run_cellANOVA:
        cellANOVA_benchmark(adata.copy(), batch_id, condition_id, control_id, out_csv = f"{out_csv}_CellANOVA.csv")
        
    if run_BBKNN:
        BBKNN_benchmark(adata.copy(), batch_id, out_csv = f"{out_csv}_BBKNN.csv")
    
def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--adata_path', type=str, required=True, help="AnnData")
    parser.add_argument('--batch_id', type=str, required=True, help="batch")
    parser.add_argument('--cell_type_id', type=str, required=True, help="cell_type")
    parser.add_argument('--condition_id', type=str, required=True, help="phenotype")
    parser.add_argument('--control_id', type=str, required=True, help="control_group")
    parser.add_argument('--out_csv', type=str, required=True, help="Path")
    parser.add_argument('--run_scVI', action='store_true', help="scVI")
    parser.add_argument('--run_cellANOVA', action='store_true', help="CellANOVA")
    parser.add_argument('--run_BBKNN', action='store_true', help="BBKNN")
    args = parser.parse_args()
    
    python_benchmark(args.adata_path, args.batch_id, args.condition_id, args.control_id, args.cell_type_id, args.out_csv, 
                     args.run_scVI, args.run_cellANOVA, args.run_BBKNN)

if __name__ == "__main__":
    main()
