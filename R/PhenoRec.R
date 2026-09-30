#' PhenoRec: phenotype-preserving integration of single-cell data
#'
#' PhenoRec provides functions for recovering phenotype-associated variation
#' that may be removed during conventional single-cell data integration.
#' The current implementation supports Gaussian measurements (for example,
#' normalized RNA expression) and Bernoulli measurements (for example, binarized
#' chromatin accessibility).
#'

#' @importFrom pracma pinv
#' @import Matrix
#' @import foreach
#' @import doMC
#' @import glmnet
#' @importFrom stats model.matrix coef
#' @import doParallel
#'
#' @keywords internal
"_PACKAGE"


#' Estimate phenotype effects
#'
#' This internal helper projects the estimated cell-state-dependent batch
#' coefficients onto the phenotype design space. It estimates the phenotype-associated coefficient matrix,
#' and remaining batch-contrast component.
#'
#' @param U A encoding batch design matrix with one row per cell.
#' @param Y A phenotype design matrix with one row per cell.
#' @param estimate_B_uz A numeric matrix of estimated coefficients for the
#'   interaction between batch and latent cell state.
#' @param d A positive integer giving the number of latent cell-state dimensions.
#'
#' @return A named list with two elements:
#' \itemize{
#'   \item \code{estimate_Gamma_yz}: the estimated phenotype effects matrix.
#'   \item \code{bulk_estimate_B_uz}: the residual batch-contrast coefficient matrix after removing the phenotype-associated component.
#' }
#'
#' @keywords internal
gamma_estimation = function(U, Y, estimate_B_uz, d){

  G = cbind(0,U)
  G[rowSums(U)==0,1] = 1
  G = t(G)
  G = G/rowSums(G)

  sample_y = G %*% Y
  sample_y = round(sample_y)

  # IMPORTANT:
  # The row order of sample_y must be identical to the sample order represented
  # by the columns of U. Constructing G directly from U guarantees this order.

  w = sample_y[-1, , drop = FALSE]

  if (all(Y==U%*%w)) {

    #print('sample_y is right')

    m_1 = nrow(w)
    I_d <- diag(d)

    y_intercept = kronecker(matrix(data = 1,nrow = m_1,ncol = 1), I_d) # (1_m ⊗ I_d),
    y_bulk = kronecker(w, I_d) # (GY ⊗ I_d),
    y_bulk <- cbind(y_intercept, y_bulk)

    estimate_Gamma_yz = pinv(t(y_bulk) %*% y_bulk) %*% t(y_bulk) %*% estimate_B_uz
    estimate_Gamma_yz = estimate_Gamma_yz[-(1:d), ,drop = FALSE]

    bulk_estimate_B_uz = estimate_B_uz - kronecker(w, I_d)%*%estimate_Gamma_yz

    return(list(estimate_Gamma_yz = estimate_Gamma_yz, bulk_estimate_B_uz = bulk_estimate_B_uz))

  }else{
    stop("w is wrong")
  }

}


#' Run PhenoRec for single-cell RNA measurements
#'
#' Fits the Gaussian version of PhenoRec using a previously estimated latent
#' cell-state representation. The function estimates cell-state-dependent
#' batch effects, identifies phenotype-associated effects, removes the
#' batch-contrast component, and reconstructs phenotype-preserving integrated data.
#'
#' @param input A named list containing:
#' \itemize{
#'   \item \code{rna_data}: a Seurat object. The current implementation expects
#'   scaled RNA measurements in \code{object@assays$RNA$scale.data} and a
#'   latent cell-state embedding in
#'   \code{object@reductions$estimate_z@cell.embeddings}.
#'   \item \code{encoding_data}: a list containing \code{batch_encoding} and
#'   \code{phenotype_encoding}, typically returned by
#'   \code{encoding_function()}.
#' }
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{estimate_eta}: fitted Gaussian values containing cell-state,
#'   phenotype, and batch-contrast components.
#'   \item \code{integrated_data}: phenotype-preserving integrated values after
#'   removal of the batch-contrast component.
#'   \item \code{estimate_Gamma_yz}: estimated phenotype effects.
#'   \item \code{estimate_phi}: estimated latent cell-state loadings.
#'   \item \code{estimate_B_uz}: estimated residual batch-contrast coefficients.
#'   \item \code{estimate_alpha}: intercept values repeated across cells.
#' }
#'
#' @details
#' The rows of the latent embedding, expression matrix, batch encoding, and
#' phenotype encoding must refer to the same cells in the same order.
#'
#' @export
PhenoRec_RNA = function(input){

  seurat_data = input$rna_data
  rna_data = seurat_data@assays$RNA$scale.data
  rna_data= t(rna_data)

  estimate_z = seurat_data@reductions$estimate_z@cell.embeddings
  X = rna_data; U = input$encoding_data$batch_encoding; Y = input$encoding_data$phenotype_encoding

  openblasctl::openblas_set_num_threads(60)

  UZ = do.call(cbind, lapply(1:ncol(U), function(b) U[, b] * estimate_z))
  YZ = do.call(cbind, lapply(1:ncol(Y), function(b) Y[, b] * estimate_z))

  #OLS
  Q = cbind(1, UZ, estimate_z)
  theta_matrix = pinv(t(Q) %*% Q) %*% t(Q) %*% X

  estimate_alpha = matrix(theta_matrix[1,],ncol = ncol(theta_matrix))
  estimate_B_uz = theta_matrix[2:(1+ncol(UZ)),]
  estimate_phi = theta_matrix[(2+ncol(UZ)):nrow(theta_matrix),]

  bulk_result = gamma_estimation(U = U, Y = Y, estimate_B_uz = estimate_B_uz, d = ncol(estimate_z))

  estimate_B_uz = bulk_result$bulk_estimate_B_uz
  estimate_Gamma_yz = bulk_result$estimate_Gamma_yz

  alpha_matrix = matrix(rep(estimate_alpha, nrow(X)), nrow = nrow(X), byrow = TRUE)

  estimate_eta = alpha_matrix + UZ%*%estimate_B_uz + YZ%*%estimate_Gamma_yz + estimate_z%*%estimate_phi
  rownames(estimate_eta) = rownames(X)

  denoise_eta = alpha_matrix + YZ%*%estimate_Gamma_yz + estimate_z%*%estimate_phi
  rownames(denoise_eta) = rownames(X)

  message("PhenoRec completed.")

  result = list(estimate_eta = estimate_eta, integrated_data = denoise_eta,
                estimate_Gamma_yz = estimate_Gamma_yz,
                estimate_phi = estimate_phi, estimate_B_uz = estimate_B_uz, estimate_alpha = alpha_matrix)

  return(result)
}


#' Construct batch and phenotype encoding matrices
#'
#' Creates reference-coded design matrices for sample/batch identity and
#' phenotype labels from the metadata of a Seurat object. The batch containing
#' the largest number of cells is used as the reference batch. The phenotype
#' associated with that reference batch is used as the reference phenotype.
#'
#' @param seurat_data A Seurat object whose \code{meta.data} contains the batch
#'   and phenotype variables.
#' @param batch_id A single character string giving the column name in
#'   \code{seurat_data@meta.data} that identifies sample or batch membership.
#' @param phenotype_id A single character string giving the column name in
#'   \code{seurat_data@meta.data} that identifies the phenotype of interest.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{batch_encoding}: a reference-coded batch design matrix.
#'   \item \code{phenotype_encoding}: a reference-coded phenotype design matrix.
#' }
#'
#' @details
#' PhenoRec assumes phenotype is defined at the sample level. Therefore, all
#' cells from the same sample/batch should have the same phenotype label.
#' The current implementation also assumes that the reference batch maps to a
#' single phenotype level.
#'
#' @export
encoding_function = function(seurat_data, batch_id, phenotype_id){
  cell_batch = seurat_data@meta.data[[batch_id]]
  cell_batch = factor(cell_batch)
  batch_count <- table(cell_batch)
  base_batch <- names(which.max(batch_count))
  batch_count;base_batch
  cell_batch <- factor(cell_batch, levels = c(base_batch, setdiff(levels(cell_batch), base_batch)))
  dummy_batch = model.matrix(~cell_batch)
  dummy_batch = dummy_batch[,-1,drop = FALSE]

  if(!is.null(phenotype_id )){
    cell_pheno = seurat_data@meta.data[[phenotype_id]]
    cell_pheno = factor(cell_pheno)
    base_pheno = cell_pheno[rowSums(dummy_batch)==0]
    length(base_pheno);unique(base_pheno)
    base_pheno <- as.character(unique(base_pheno))
    cell_pheno = factor(cell_pheno, levels = c(base_pheno, setdiff(levels(cell_pheno), base_pheno)))
    dummy_pheno = model.matrix(~cell_pheno)
    dummy_pheno = dummy_pheno[, -1, drop = FALSE]
  }
  return(list(batch_encoding = dummy_batch, phenotype_encoding = dummy_pheno))

}

#' Run PhenoRec for single-cell ADT measurements
#'
#' Fits the Gaussian version of PhenoRec using a previously estimated latent
#' cell-state representation. The function estimates cell-state-dependent
#' batch effects, identifies phenotype-associated effects, removes the
#' batch-contrast component, and reconstructs phenotype-preserving integrated data.
#'
#' @param input A named list containing:
#' \itemize{
#'   \item \code{adt_data}: a Seurat object. The current implementation expects
#'   scaled ADT measurements in \code{object@assays$ADT$scale.data} and a
#'   latent cell-state embedding in
#'   \code{object@reductions$estimate_z@cell.embeddings}.
#'   \item \code{encoding_data}: a list containing \code{batch_encoding} and
#'   \code{phenotype_encoding}, typically returned by
#'   \code{encoding_function()}.
#' }
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{estimate_eta}: fitted Gaussian values containing cell-state,
#'   phenotype, and batch-contrast components.
#'   \item \code{integrated_data}: phenotype-preserving integrated values after
#'   removal of the batch-contrast component.
#'   \item \code{estimate_Gamma_yz}: estimated phenotype effects.
#'   \item \code{estimate_phi}: estimated latent cell-state loadings.
#'   \item \code{estimate_B_uz}: estimated residual batch-contrast coefficients.
#'   \item \code{estimate_alpha}: intercept values repeated across cells.
#' }
#'
#' @details
#' The rows of the latent embedding, ADT matrix, batch encoding, and
#' phenotype encoding must refer to the same cells in the same order.
#'
#' @export
PhenoRec_ADT = function(input){

  seurat_data = input$adt_data
  adt_data = seurat_data@assays$ADT$scale.data
  adt_data= t(adt_data)

  estimate_z = seurat_data@reductions$estimate_z@cell.embeddings
  X = adt_data; U = input$encoding_data$batch_encoding; Y = input$encoding_data$phenotype_encoding

  openblasctl::openblas_set_num_threads(60)

  UZ = do.call(cbind, lapply(1:ncol(U), function(b) U[, b] * estimate_z))
  YZ = do.call(cbind, lapply(1:ncol(Y), function(b) Y[, b] * estimate_z))

  #OLS
  Q = cbind(1, UZ, estimate_z)
  theta_matrix = pinv(t(Q) %*% Q) %*% t(Q) %*% X

  estimate_alpha = matrix(theta_matrix[1,],ncol = ncol(theta_matrix))
  estimate_B_uz = theta_matrix[2:(1+ncol(UZ)),]
  estimate_phi = theta_matrix[(2+ncol(UZ)):nrow(theta_matrix),]

  bulk_result = gamma_estimation(U = U, Y = Y, estimate_B_uz = estimate_B_uz, d = ncol(estimate_z))

  estimate_B_uz = bulk_result$bulk_estimate_B_uz
  estimate_Gamma_yz = bulk_result$estimate_Gamma_yz

  alpha_matrix = matrix(rep(estimate_alpha, nrow(X)), nrow = nrow(X), byrow = TRUE)

  estimate_eta = alpha_matrix + UZ%*%estimate_B_uz + YZ%*%estimate_Gamma_yz + estimate_z%*%estimate_phi
  rownames(estimate_eta) = rownames(X)

  denoise_eta = alpha_matrix + YZ%*%estimate_Gamma_yz + estimate_z%*%estimate_phi
  rownames(denoise_eta) = rownames(X)

  message("PhenoRec completed.")

  result = list(estimate_eta = estimate_eta, integrated_data = denoise_eta,
                estimate_Gamma_yz = estimate_Gamma_yz,
                estimate_phi = estimate_phi, estimate_B_uz = estimate_B_uz, estimate_alpha = alpha_matrix)

  return(result)
}


#' Subsample negative cells for a binary feature
#'
#' Internal helper for the Bernoulli model. All positive cells are retained,
#' while a subset of negative cells is sampled according to a user-specified
#' negative-to-positive ratio.
#'
#' @param xj A binary numeric vector containing 0/1 measurements for one
#'   feature across cells.
#' @param neg_ratio A positive numeric value specifying the maximum number of
#'   negative cells sampled per positive cell.
#' @param seed An integer random seed used for reproducible negative sampling.
#'
#' @return A named list containing \code{use_cells}, the indices of retained
#'   cells, and \code{s0}, the fraction of all negative cells retained.
#'
#' @keywords internal
sample_case_control_cells <- function(xj, neg_ratio, seed) {

  set.seed(seed)

  pos <- which(xj == 1)
  neg <- which(xj == 0)

  n_neg_sample <- min(length(neg), neg_ratio * length(pos))

  neg_sub <- sample(neg, n_neg_sample)
  use_cells <- sort(c(pos, neg_sub))

  sampling_info <- list(
    use_cells = use_cells,
    s0 = n_neg_sample / length(neg)
  )

  return(sampling_info)
}

#' Correct a logistic-regression intercept after negative subsampling
#'
#' Internal helper that adjusts the fitted intercept to account for the
#' sampling fraction of negative observations used in the Bernoulli model.
#'
#' @param beta A numeric coefficient vector whose first element is the fitted
#'   intercept.
#' @param sampling_info_s0 A numeric value in \eqn{(0, 1]} giving the fraction
#'   of negative observations retained during subsampling.
#'
#' @return The coefficient vector \code{beta} with its intercept corrected for
#'   negative subsampling.
#'
#' @keywords internal
correct_intercept <- function(beta, sampling_info_s0) {

  s1 <- 1
  beta[1] <- beta[1] - log(s1 / sampling_info_s0)

  return(beta)
}

#' Fit logistic regression
#'
#' Internal helper that fits a binomial ridge-regression model using
#' \code{glmnet} in a forked process and terminates the process when the
#' requested timeout is exceeded.
#'
#' @param Q_use A numeric design matrix for the selected cells.
#' @param y_use A binary numeric response vector aligned with the rows of
#'   \code{Q_use}.
#' @param lambda0 A positive numeric value giving the ridge penalty parameter
#'   passed to \code{glmnet::glmnet()}.
#' @param timeout_sec A positive numeric value giving the maximum fitting time
#'   in seconds. Default is 300 seconds.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{ok}: logical indicator of whether fitting completed.
#'   \item \code{reason}: \code{"ok"} or \code{"timeout"}.
#'   \item \code{beta}: fitted coefficient vector, or \code{NULL} after timeout.
#'   \item \code{dev_ratio}: fitted deviance ratio, or \code{NA} after timeout.
#' }
#'
#' @details
#' This implementation relies on \code{parallel::mcparallel()} and is therefore
#' intended for Unix-like operating systems that support forked processes.
#'
#' @keywords internal

fit_glmnet_timeout <- function(Q_use, y_use, lambda0, timeout_sec = 300) {

  job <- parallel::mcparallel({

    fit <- glmnet::glmnet(
      x = Q_use,
      y = y_use,
      family = "binomial",
      alpha = 0,
      lambda = lambda0,
      intercept = TRUE,
      standardize = FALSE
    )

    list(
      beta = as.numeric(coef(fit)),
      dev_ratio = as.numeric(fit$dev.ratio)
    )

  }, silent = TRUE)

  res <- parallel::mccollect(job, wait = FALSE, timeout = timeout_sec)

  if (is.null(res)) {
    try(tools::pskill(job$pid, signal = 15), silent = TRUE)
    Sys.sleep(0.2)
    try(tools::pskill(job$pid, signal = 9), silent = TRUE)
    return(list(ok = FALSE, reason = "timeout", beta = NULL, dev_ratio = NA_real_))
  }

  res <- res[[1]]

  return(list(
    ok = TRUE,
    reason = "ok",
    beta = res$beta,
    dev_ratio = res$dev_ratio
  ))
}


#' Run PhenoRec for single-cell ATAC measurements
#'
#' Fits the Bernoulli version of PhenoRec to binarized single-cell chromatin
#' accessibility measurements. For each feature, the function performs
#' negative-cell subsampling, fits ridge logistic regression, identifies
#' phenotype-associated effects, removes the batch-contrast component, and
#' reconstructs a phenotype-preserving integrated data.
#'
#' @param input A named list containing:
#' \itemize{
#'   \item \code{atac_data}: a Seurat object. The current implementation expects
#'   raw ATAC counts in \code{object@assays$ATAC$counts}, total fragment/count
#'   information in \code{object$nCount_ATAC}, and a latent cell-state embedding
#'   in \code{object@reductions$estimate_z@cell.embeddings}.
#'   \item \code{encoding_data}: a list containing \code{batch_encoding} and
#'   \code{phenotype_encoding}, typically returned by
#'   \code{encoding_function()}.
#' }
#' @param ncores A positive integer specifying the number of worker processes
#'   used by the \code{doMC}/\code{foreach} backend.
#' @param neg_ratio A positive numeric value specifying the maximum number of
#'   negative cells sampled per positive cell for each feature. Default is 5.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{integrated_data}: phenotype-preserving integrated data after removing the batch-contrast component.
#'   \item \code{estimate_Gamma_yz}: estimated phenotype effects.
#'   \item \code{estimate_phi}: estimated latent cell-state loadings.
#'   \item \code{estimate_B_uz}: estimated residual batch-contrast coefficients.
#'   \item \code{estimate_alpha}: estimated feature-specific intercepts.
#'   \item \code{estimate_depth}: coefficients for standardized log ATAC depth.
#'   \item \code{dev_ratio}: feature-wise deviance ratios from \code{glmnet}.
#' }
#'
#' @details
#' ATAC counts are binarized before model fitting. The function currently processes 100 features per chunk, and applies
#' a 300-second timeout to each feature fit.
#'
#' Parallel fitting uses \code{doMC} and \code{parallel::mcparallel()}, so the
#' current implementation is intended for Unix-like systems. The function also
#' changes BLAS-related thread environment variables to one thread during
#' fitting.
#'
#' @export
PhenoRec_ATAC = function(input,ncores,neg_ratio=5){

  seurat_data = input$atac_data
  nFrag <- seurat_data$nCount_ATAC
  log_nFrag <- scale(log(nFrag))[, 1]

  atac_data = seurat_data@assays$ATAC$counts
  atac_data = drop0(atac_data)
  atac_data@x = rep(1, length(atac_data@x))
  X = t(atac_data)

  estimate_z = seurat_data@reductions$estimate_z@cell.embeddings

  U = input$encoding_data$batch_encoding; Y = input$encoding_data$phenotype_encoding

  UZ = do.call(cbind, lapply(1:ncol(U), function(b) U[, b] * estimate_z))
  YZ = do.call(cbind, lapply(1:ncol(Y), function(b) Y[, b] * estimate_z))

  Q = cbind(UZ, estimate_z, log_nFrag = log_nFrag)  ##glmnet includes an intercept internally
  print(dim(Q))

  lambda0 <- 1e-5  #Ridge penalty used to stabilize logistic regression under separation
  n_coef <- ncol(Q) + 1L

  #Number of peaks processed per chunk
  chunk_size <- 100
  peak_ids <- seq_len(ncol(X))
  peak_chunks <- split(peak_ids, ceiling(seq_along(peak_ids) / chunk_size))

  n_chunks <- length(peak_chunks)

  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1"
  )

  if (requireNamespace("openblasctl", quietly = TRUE)) {
    openblasctl::openblas_set_num_threads(1)
  }

  #Fork-based parallel
  registerDoMC(cores = ncores)

  print('parallel...........')

  theta_matrix <- foreach(
    i = seq_along(peak_chunks),
    .combine = cbind,
    .packages = c("glmnet", "Matrix", "parallel"),
    .inorder = TRUE,
    .multicombine = TRUE,
    .maxcombine = 8,
    .options.multicore = list(preschedule = FALSE)
  ) %dopar% {

    message(sprintf("Start chunk %d/%d", i, n_chunks))

    idx <- peak_chunks[[i]]
    out <- matrix(NA_real_, nrow = n_coef+1, ncol = length(idx))

    for (k in seq_along(idx)) {

      j <- idx[k]

      xj <- as.numeric(X[, j])

      sampling_info <- sample_case_control_cells(xj = xj, neg_ratio = neg_ratio, seed = j)
      use_cells <- sampling_info$use_cells

      fit <- fit_glmnet_timeout(
        Q_use = Q[use_cells, , drop = FALSE],
        y_use = xj[use_cells],
        lambda0 = lambda0,
        timeout_sec = 300
      )

      if (!fit$ok) {

        out[, k] <- NA_real_

        cat(sprintf(
          "Chunk %d/%d | peak j = %d | skipped: %s, n_used = %d\n",
          i, n_chunks, j, fit$reason, length(use_cells)
        ))

        rm(xj, fit)
        next
      }

      beta <- fit$beta

      beta <- correct_intercept(beta = beta,sampling_info_s0 = sampling_info$s0)

      out[, k] <- c(beta, fit$dev_ratio)

      if (i %% ncores == 0 && k %% 1 == 0) {cat(sprintf("Chunk %d/%d | k = %d, n_used = %d | dev.ratio = %.4f\n", i, n_chunks, k, length(use_cells), fit$dev_ratio))}

      rm(xj, fit, beta)
    }

    gc(FALSE)

    out
  }

  colnames(theta_matrix) <- colnames(X)

  ## Remove features whose fitted coefficients contain NA/NaN/Inf/-Inf
  valid_cols <- colSums(!is.finite(theta_matrix)) == 0
  theta_matrix <- theta_matrix[, valid_cols, drop = FALSE]
  cat("Kept peaks:", sum(valid_cols), "\n")

  openblasctl::openblas_set_num_threads(ncores)

  estimate_alpha = matrix(theta_matrix[1,],ncol = ncol(theta_matrix))
  estimate_B_uz = theta_matrix[2:(1+ncol(UZ)),]
  estimate_phi = theta_matrix[(2+ncol(UZ)):(1+ncol(UZ)+ncol(estimate_z)),]

  estimate_depth = theta_matrix[(2+ncol(UZ)+ncol(estimate_z)),]

  dev_ratio = theta_matrix[(3+ncol(UZ)+ncol(estimate_z)),]

  bulk_result = gamma_estimation(U = U, Y = Y, estimate_B_uz = estimate_B_uz, d = ncol(estimate_z))

  estimate_B_uz = bulk_result$bulk_estimate_B_uz
  estimate_Gamma_yz = bulk_result$estimate_Gamma_yz

  alpha_matrix = matrix(rep(estimate_alpha, nrow(X)), nrow = nrow(X), byrow = TRUE)

  estimate_eta = alpha_matrix + UZ%*%estimate_B_uz + YZ%*%estimate_Gamma_yz + estimate_z%*%estimate_phi
  rownames(estimate_eta) = rownames(X)

  denoise_eta = alpha_matrix + YZ%*%estimate_Gamma_yz + estimate_z%*%estimate_phi
  rownames(denoise_eta) = rownames(X)
  colnames(denoise_eta) = colnames(theta_matrix)
  #denoise_prob <- plogis(denoise_eta)

  message("PhenoRec completed.")

  result = list(integrated_data = denoise_eta,
                estimate_Gamma_yz = estimate_Gamma_yz,
                estimate_phi = estimate_phi, estimate_B_uz = estimate_B_uz, estimate_alpha = alpha_matrix[1,],
                estimate_depth = estimate_depth, dev_ratio = dev_ratio)

  return(result)
}





