## simulate under alpha = 0 and run when alpha = 0
## use_t_adjust = TRUE here

one_rep <- function(new_params, current_params) {
  args_val <- append(current_params, new_params)
  set.seed(new_params$current_seed)

  ## Choose all of the genes because already got top expressed
  stopifnot(args_val$Ngene == ncol(args_val$mat))
  d_out <- seqgendiff::poisthin(mat = args_val$mat,
                                nsamp = args_val$Nsamp,
                                ngene = args_val$Ngene,
                                skip_gene = args_val$skip_gene,
                                signal_params = list(mean = 0, sd = args_val$log2foldsd),
                                gvec = rep(TRUE, length(args_val$Ngene)),
                                gselect = "custom",
                                prop_null = args_val$nullpi,
                                alpha = 0)

  which_null <- abs(d_out$beta) < 10 ^ -6
  nnull         <- sum(which_null)
  control_genes <- which_null
  control_genes[control_genes][sample(1:nnull, size = nnull - args_val$ncontrol)] <- FALSE

  beta_true <- d_out$beta

  X <- d_out$X
  colnames(X) <- c("Intercept", "Treatment")
  Y <- log2(d_out$Y + 1)

  num_sv <- max(sva::num.sv(t(Y), mod = X, method = "be"), 1)


  ash_list <- list()
  ash_list$mouthwash_t_adjust_unif <- vicar::mouthwash(Y = Y, X = X, k = num_sv,
                                                       likelihood = "t",
                                                       mixing_dist = "sym_uniform",
                                                       cov_of_interest = 2,
                                                       include_intercept = FALSE,
                                                       use_t_adjust = TRUE,
                                                       scale_var = FALSE)
  ash_list$mouthwash_t_adjust_norm <- vicar::mouthwash(Y = Y, X = X, k = num_sv,
                                                       likelihood = "t",
                                                       mixing_dist = "normal",
                                                       cov_of_interest = 2,
                                                       include_intercept = FALSE,
                                                       use_t_adjust = TRUE,
                                                       scale_var = FALSE)
  ash_list$mouthwash_t_adjust_unif_sv <- vicar::mouthwash(Y = Y, X = X, k = num_sv,
                                                       likelihood = "t",
                                                       mixing_dist = "sym_uniform",
                                                       cov_of_interest = 2,
                                                       include_intercept = FALSE,
                                                       use_t_adjust = TRUE,
                                                       scale_var = TRUE)
  ash_list$mouthwash_t_adjust_norm_sv <- vicar::mouthwash(Y = Y, X = X, k = num_sv,
                                                       likelihood = "t",
                                                       mixing_dist = "normal",
                                                       cov_of_interest = 2,
                                                       include_intercept = FALSE,
                                                       use_t_adjust = TRUE,
                                                       scale_var = TRUE)


  get_mse <- function(args, beta_true) {
    mean((args$result$PosteriorMean - beta_true) ^ 2, na.rm = TRUE)
  }

  get_auc <- function(args, which_null) {
    pROC::roc(predictor = args$result$lfdr, response = which_null)$auc
  }

  ## pi0hat ----------------------------------------------------------
  ash_pi0 <- sapply(ash_list, FUN = function(args) args$pi0)
  names(ash_pi0) <- paste0("pi0_ash_", names(ash_pi0))
  pi0_vec <- ash_pi0

  ## auc ------------------------------------------------------------
  ash_auc <- sapply(ash_list, FUN = get_auc, which_null = which_null)
  names(ash_auc) <- paste0("auc_ash_", names(ash_auc))
  auc_vec <- ash_auc

  ## mse ------------------------------------------------------------
  ash_mse <- sapply(ash_list, FUN = get_mse, beta_true = beta_true)
  names(ash_mse) <- paste0("mse_ash_", names(ash_mse))
  mse_vec <- ash_mse

  return_vec <- c(mse_vec, auc_vec, pi0_vec)

  return(return_vec)
}

itermax <- 100 ## itermax should be 500
seed_start <- 2222

## these change
nullpi_seq   <- c(0.5)
Nsamp_seq    <- c(6, 40)
ncontrol_seq <- c(10)

par_vals <- expand.grid(list((1 + seed_start):(itermax + seed_start),
                             nullpi_seq, Nsamp_seq, ncontrol_seq))
colnames(par_vals) <- c("current_seed", "nullpi", "Nsamp", "ncontrols")
par_vals$poisthin <- TRUE
par_vals$poisthin[abs(par_vals$nullpi - 1) < 10 ^ -10] <- FALSE

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

## Create muscle_mat with most expressed genes
mat <- t(as.matrix(read.csv("../../reproduce_ruv3/Output/gtex_tissue_gene_reads_v6p/muscle.csv",
                            header = TRUE)[, -c(1,2)]))
args_val$mat <- mat[, order(apply(mat, 2, median), decreasing = TRUE)[1:args_val$Ngene]]
rm(mat)

start_time <- proc.time()
oout <- one_rep(par_list[[10]], args_val)
oout
end_time <- proc.time() - start_time
end_time

## If on your own computer, use this
library(snow)
library(parallel)
cl <- makeCluster(detectCores() - 2)
sout <- t(snow::parSapply(cl = cl, par_list, FUN = one_rep, current_params = args_val))
stopCluster(cl)

saveRDS(cbind(par_vals, sout), "../output/alpha_1_sims_out/sims_out_alpha0_adjust_by_t_mouth.RDS")






