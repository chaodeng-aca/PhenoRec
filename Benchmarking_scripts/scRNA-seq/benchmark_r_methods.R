library(Seurat)
library(rliger)

seurat_benchmark = function(seurat_data,batch_id){
  
  seurat_data[["RNA"]] <- split(seurat_data[["RNA"]], f = batch_id)
  
  seurat_data <- NormalizeData(seurat_data)
  seurat_data <- FindVariableFeatures(seurat_data)
  seurat_data <- ScaleData(seurat_data)
  seurat_data <- RunPCA(seurat_data)
  
  seurat_data <- IntegrateLayers(object = seurat_data, method = CCAIntegration, orig.reduction = "pca", 
                                 new.reduction = "integrated.cca",verbose = FALSE)
  seurat_data[["RNA"]] <- JoinLayers(seurat_data[["RNA"]])
  seurat_embed = seurat_data@reductions$integrated.cca@cell.embeddings
  
  return(seurat_embed)
}


liger_benchmark = function(seurat_data, batch_id){
  
  count_matrix = seurat_data@assays$RNA$counts
  ligerObj = as.liger(count_matrix, datasetVar = batch_id)
  
  ligerObj = ligerObj %>%
    normalize() %>%
    selectGenes() %>%
    scaleNotCenter()
  
  ligerObj = runIntegration(ligerObj, k = 50)
  ligerObj = alignFactors(ligerObj, method = "centroidAlign")
  liger_embed = ligerObj@H.norm
  
  rownames(liger_embed) <- sub(paste0("^(", paste(unique(batch_id), collapse = "|"), ")_"), "", rownames(liger_embed))
  
  return(liger_embed)
  
}


R_benchmark = function(seurat_data_path, batch_id, run_Seurat = TRUE, run_Liger = TRUE){
  
  seurat_data = readRDS(seurat_data_path)
  batch_id = seurat_data@meta.data[[batch_id]]
  
  result_list <- list()
  
  if (run_Seurat) {
    seurat_result = seurat_benchmark(seurat_data = seurat_data, batch_id = batch_id)
    result_list$Seurat = seurat_result
  }
  if (run_Liger) {
    liger_result = liger_benchmark(seurat_data = seurat_data, batch_id = batch_id)
    result_list$Liger = liger_result
  }
  return(result_list)
}


run__benchmark <- function(seurat_data_path, adata_path, batch_id, cell_type_id, condition_id, control_id, out_csv, 
                                 run_Seurat = TRUE, run_Liger = TRUE, run_scVI = TRUE, run_cellANOVA = TRUE, run_BBKNN = TRUE) {
  
  results = R_benchmark(seurat_data_path = seurat_data_path, batch_id = batch_id, run_Seurat = run_Seurat, run_Liger = run_Liger)

  python_path <- "/dssg/home/acct-clswt/clswt-chaodeng/miniconda3/envs/python310/bin/python"
  script_path <- "~/phenotype_project/scRNA/benchmark_methods/benchmark_python_methods.py"  
  
  args <- c(script_path,
            "--adata_path", adata_path,
            "--batch_id", batch_id,
            "--cell_type_id", cell_type_id,
            "--condition_id", condition_id,
            "--control_id", control_id,
            "--out_csv", out_csv)
  
  if (run_scVI) {
    args <- c(args, "--run_scVI")
  }
  if (run_cellANOVA) {
    args <- c(args, "--run_cellANOVA")
  }
  if (run_BBKNN) {
    args <- c(args, "--run_BBKNN")
  }
  system2(python_path, args = args, stdout = TRUE, stderr = TRUE)
  
  
  if (run_scVI) {
    results$scVI = as.matrix(read.csv(paste0(out_csv,'_scVI.csv'), row.names = 1))
  }
  if (run_cellANOVA) {
    results$CellANOVA = as.matrix(read.csv(paste0(out_csv,'_CellANOVA.csv'), row.names = 1))
  }
  if (run_BBKNN) {
    results$BBKNN = as.matrix(read.csv(paste0(out_csv,'_BBKNN.csv'), row.names = 1))
  }
  
  
  return(results)
  
}






