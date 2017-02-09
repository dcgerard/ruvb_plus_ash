library(snow)
library(seqgendiff)

#' new_params and current_params are basically interchangeable.
#' They get combined into one big list. I separate them just so that I can
#' parallelize the implementation.
#'
#' mat contains a matrix of gene expression values. I think this is better than reading
#' in the matrix each iteration --- it should speed things up a tiny bit.
one_rep <- function(new_params, current_params, submat) {
  source("../code/adjustment_methods.R")
  args_val <- append(current_params, new_params)
  set.seed(new_params$current_seed)
  stopifnot(ncol(submat) == args_val$Ngene)
  dout <- seqgendiff::poisthin(mat = submat, nsamp = args_val$Nsamp,
                               ngene = args_val$Ngene, gselect = "custom",
                               gvec = rep(TRUE, args_val$Ngene),
                               skip_gene = 0,
                               signal_fun = stats::rnorm,
                               signal_params = list(mean = args_val$log2foldmean,
                                                    sd = args_val$log2foldsd),
                               prop_null = args_val$nullpi)

  which_null    <- abs(dout$beta) < 10^-6
  control_genes <- which_null
  nnull         <- sum(control_genes)
  control_genes[control_genes][sample(1:nnull, size = nnull - args_val$ncontrol)] <- FALSE

  beta_true <- dout$beta
  X <- dout$X
  colnames(X) <- c("Intercept", "Treatment")
  Y <- log2(as.matrix(dout$Y + 1))


  num_sv <- max(sva::num.sv(t(Y), mod = X, method = "be"), 1)

  method_list            <- list()
  method_list$ols        <- ols(Y = Y, X = X)

  ## control gene methods --------------------------------------------------
  method_list$ruv2         <- ruv2(Y = Y, X = X, num_sv = num_sv,
                                   control_genes = control_genes)
  method_list$ruv3_nomult  <- ruv3(Y = Y, X = X, num_sv = num_sv,
                                   control_genes = control_genes,
                                   multiplier = FALSE)
  method_list$ruv3_mult    <- ruv3(Y = Y, X = X, num_sv = num_sv,
                                   control_genes = control_genes,
                                   multiplier = TRUE)
  method_list$ruv4         <- ruv4(Y = Y, X = X, num_sv = num_sv,
                                   control_genes = control_genes)
  method_list$ruv4_rsvar   <- ruv4_rsvar_ebayes(Y = Y, X = X, num_sv = num_sv,
                                                control_genes = control_genes)
  method_list$catenc_nocal <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                      control_genes = control_genes,
                                      calibrate = FALSE)
  method_list$catenc_cal   <- cate_nc(Y = Y, X = X, num_sv = num_sv,
                                      control_genes = control_genes,
                                      calibrate = TRUE)
  method_list$ruv4v_norm   <- vruv4(Y = Y, X = X,
                                    num_sv = num_sv,
                                    control_genes = control_genes,
                                    likelihood = "normal")
  method_list$ruv4v_t      <- vruv4(Y = Y, X = X,
                                    num_sv = num_sv,
                                    control_genes = control_genes,
                                    likelihood = "t")

  ## non control gene methods -----------------------------------------------
  method_list$sva          <- sva(Y = Y, X = X, num_sv = num_sv)
  method_list$caterr_nocal <- cate_rr(Y = Y, X = X, num_sv = num_sv, calibrate = FALSE)
  method_list$caterr_cal   <- cate_rr(Y = Y, X = X, num_sv = num_sv, calibrate = TRUE)

  ## LEAPP SUCKS! Slow and always with the BUGS!
  ## method_list$leapp <- leapp(Y = Y, X = X, num_sv = num_sv)

  pout <- sapply(method_list, FUN = function(x) x$pvalue)

  ## Note: RUVB's p-values are actually the lfsr's.
  ## Need to calculate them from betahat and sebetahat
  method_list$ruvb         <- ruvb_bfa_gs_linked(Y = Y, X = X,
                                                 num_sv = num_sv,
                                                 control_genes = control_genes)
  ruvb_pvalues <- c(2 * stats::pnorm(-abs(method_list$ruvb$betahat / method_list$ruvb$sebetahat)))

  pout <- cbind(pout, ruvb_pvalues, which_null, control_genes)
  colnames(pout)[(ncol(pout) - 2):ncol(pout)] <- c("ruvb", "which_null", "control_genes")

  # ## Now get the AUC
  # get_auc <- function(pvalues, which_null){
  #   return(pROC::auc(predictor = pvalues, response =  which_null))
  # }
  # auc_vec <- apply(pout, 2, get_auc, which_null = which_null)
  #
  # ## Now get FDP
  # get_fdp <- function(pvalues, which_null, fdr_levels = c(0.01, 0.05, 0.1, 0.2)) {
  #   qvalues <- stats::p.adjust(pvalues, method = "BH")
  #   fdp_vec <- rep(NA, length = length(fdr_levels))
  #   for(findex in 1:length(fdr_levels)) {
  #     sig_genes <- qvalues <= fdr_levels[findex]
  #     if (sum(sig_genes, na.rm = TRUE) == 0) {
  #       fdp_vec[findex] <- 0
  #     } else {
  #       fdp_vec[findex] <- mean(which_null[sig_genes], na.rm = TRUE)
  #     }
  #   }
  #   return(fdp_vec)
  # }
  # qout <- apply(pout, 2, get_fdp, which_null = which_null)



  return(pout)
}

itermax <- 200
seed_start <- 2222

## these change
nullpi_seq   <- c(0.5, 0.9, 1)
Nsamp_seq    <- c(6, 10, 20, 40)
ncontrol_seq <- c(10, 100)

par_vals <- expand.grid(list((1 + seed_start):(itermax + seed_start),
                             nullpi_seq, Nsamp_seq, ncontrol_seq))
colnames(par_vals) <- c("current_seed", "nullpi", "Nsamp", "ncontrols")

par_list <- list()
for (list_index in 1:nrow(par_vals)) {
    par_list[[list_index]] <- list()
    for (inner_list_index in 1:ncol(par_vals)) {
        par_list[[list_index]][[inner_list_index]] <- par_vals[list_index, inner_list_index]
        names(par_list[[list_index]])[inner_list_index] <- colnames(par_vals)[inner_list_index]
    }
}

## these do not change
args_val              <- list()
args_val$log2foldsd   <- 0.8
args_val$Ngene        <- 1000
args_val$log2foldmean <- 0
args_val$skip_gene    <- 0

## read in gene expression matrix
mat <- as.matrix(read.csv("../output/gtex_tissue_gene_reads_v6p/muscle.csv")[, -(1:2)])
## select top genes
rowmed <- apply(mat, 1, median)
rowmean <- rowMeans(mat)
order_vec <- order(rowmed, rowmean, decreasing = TRUE)
submat <- t(mat[order_vec[1:args_val$Ngene], ])
rm(mat) # so that doesn't take up memory

## one_rep(par_list[[3]], args_val, submat)

## ## If on your own computer, use this
library(parallel)
cl <- makeCluster(detectCores() - 2)
sout <- t(snow::parLapply(cl = cl, par_list, fun = one_rep, current_params = args_val,
                          submat = submat))
stopCluster(cl)

sout <- t(rbind(sout, t(as.matrix(par_vals))))
colnames(sout)[1] <- "pvalues"

save(sout, file = "../output/sims_out/pvalue_matrices.Rd")

