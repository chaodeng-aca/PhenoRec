library(Signac)
library(Seurat)

signac_benchmark = function(seurat_data,batch_id){
  
  seurat_data = RunTFIDF(seurat_data)
  seurat_data = RunSVD(seurat_data,features = rownames(seurat_data))
  
  
  obj = seurat_data
  obj.list <- SplitObject(obj, split.by = batch_id)
  
  obj.list <- lapply(obj.list, function(x) {
    
    x <- RunTFIDF(x)
    x <- RunSVD(
      object = x,
      n = 50, features = rownames(x),
      reduction.name = "lsi",
      reduction.key = "LSI_"
    )
    
    return(x)
  })
  
  
  integration.anchors <- FindIntegrationAnchors(
    object.list = obj.list,
    anchor.features = rownames(obj),
    reduction = "rlsi",
    dims = 2:50
  )
  
  obj.integrated <- IntegrateEmbeddings(
    anchorset = integration.anchors,
    reductions = obj[["lsi"]],
    new.reduction.name = "integrated_lsi",
    dims.to.integrate = 1:50,
    k.weight = 50
  )
  
  signac_embed = obj.integrated@reductions$integrated_lsi@cell.embeddings[,2:50]
  
  return(signac_embed)
  
}




R_benchmark = function(seurat_data_path, batch_id, run_Signac = TRUE){
  
  seurat_data = readRDS(seurat_data_path)
  
  result_list <- list()
  
  if (run_Signac) {
    seurat_result = signac_benchmark(seurat_data = seurat_data, batch_id = batch_id)
    result_list$Signac = seurat_result
  }

  return(result_list)
}


run__benchmark <- function(seurat_data_path, adata_path, batch_id, out_csv, 
                                 run_Signac = TRUE, run_PeakVI = TRUE) {
  
  results = R_benchmark(seurat_data_path = seurat_data_path, batch_id = batch_id, run_Signac = run_Signac)

  python_path <- "/dssg/home/acct-clswt/clswt-chaodeng/miniconda3/envs/python310/bin/python"
  script_path <- "~/phenotype_project/scATAC/benchmark_methods/python_benchmark.py"  
  
  args <- c(script_path,
            "--adata_path", adata_path,
            "--batch_id", batch_id,
            "--out_csv", out_csv)
  
  if (run_PeakVI) {
    args <- c(args, "--run_PeakVI")
  }


  system2(python_path, args = args, stdout = TRUE, stderr = TRUE)
  
  
  if (run_PeakVI) {
    results$PeakVI = as.matrix(read.csv(paste0(out_csv,'_PeakVI.csv'), row.names = 1))
  }
  
  return(results)
  
}






