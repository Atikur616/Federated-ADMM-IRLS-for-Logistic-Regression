## Packages 
library(pROC)
library(dplyr)

# Replications
total_rep <- 100

rep_start <- 1
rep_end   <- total_rep

cat("\nRunning replications", rep_start, "to", rep_end, "\n")


# parameters

# Number of features
p <- 50

# Number of samples per client
n <- 100

## Number of center
m <- 2

# Sequence of lambda values for tuning
LAMBDA_GRID <- exp(seq(log(0.05), log(8), length.out = 20))

# threshold
threshold <- 1e-3
# Absolute tolerance
eps_abs <- 1e-5
# Relative tolerance
eps_rel <- 1e-4

# ADMM penalty parameter rho
rho_values <- c(1, 1.5, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0)

### True beta
true_beta <- c(2, 3, 1, 4, -2, rep(0, p - 5))

# Helper functions
soft_threshold <- function(x, lambda) sign(x) * pmax(0, abs(x) - lambda)

simulate_data <- function(n, p, beta) {
  X <- matrix(rnorm(n * p), nrow = n)
  prob <- plogis(X %*% beta)
  y <- rbinom(n, 1, prob)
  list(X = X, y = y)
}

# Prob-MSE
val_prob_mse <- function(X, y, beta) {
  pr <- as.vector(plogis(X %*% beta))
  pr <- pmin(pmax(pr, 1e-6), 1 - 1e-6)
  mean((y - pr)^2)
}

# safe_auc
safe_auc <- function(y, score) {
  if (length(unique(y)) < 2) return(NA_real_)
  tryCatch(
    as.numeric(pROC::auc(pROC::roc(y, as.vector(score), quiet = TRUE))),
    error = function(e) NA_real_
  )
}

# validation AUC
val_auc_mean_clients <- function(clients_val, beta) {
  aucs <- sapply(seq_along(clients_val), function(k) {
    pr <- as.vector(plogis(clients_val[[k]]$X %*% beta))
    safe_auc(clients_val[[k]]$y, pr)
  })
  mean(aucs, na.rm = TRUE)
}

# ADMM-IRLS local update
admm_local_irls <- function(X, y, z, u, rho, irls_max_iter = 100, tol = 1e-6) {
  p_local <- ncol(X)
  beta <- z
  
  for (t in 1:irls_max_iter) {
    eta   <- as.vector(X %*% beta)
    p_hat <- pmin(pmax(plogis(eta), 1e-6), 1 - 1e-6)
    w     <- pmax(p_hat * (1 - p_hat), 1e-3)
    z_adj <- eta + (y - p_hat) / w
    
    WX   <- X * w
    XtWX <- crossprod(X, WX)
    XtWz <- crossprod(X, w * z_adj)
    
    beta_new <- solve(
      XtWX + rho * diag(p_local),
      XtWz + rho * (z - u)
    )
    
    if (max(abs(beta_new - beta)) < tol) break
    beta <- beta_new
  }
  
  as.vector(beta)
}

# ADMM residuals for convergence
admm_check <- function(beta_list, z_new, z_prev, u_list, rho, eps_abs, eps_rel, p, m) {
  r_norm <- sqrt(sum(sapply(beta_list, function(b) sum((b - z_new)^2))))
  s_norm <- rho * sqrt(m) * sqrt(sum((z_new - z_prev)^2))
  
  beta_norm <- sqrt(sum(sapply(beta_list, function(b) sum(b^2))))
  z_norm    <- sqrt(sum(z_new^2))
  u_norm    <- sqrt(sum(sapply(u_list, function(u) sum(u^2))))
  
  eps_primal <- sqrt(p * m) * eps_abs + eps_rel * max(beta_norm, sqrt(m) * z_norm)
  eps_dual   <- sqrt(p * m) * eps_abs + eps_rel * (rho * u_norm)
  
  list(
    r_norm = r_norm,
    s_norm = s_norm,
    eps_primal = eps_primal,
    eps_dual = eps_dual
  )
}

# Run ADMM
run_admm <- function(clients_train, lambda, rho, p, m, eps_abs, eps_rel,
                     z_init = NULL, u_init = NULL,
                     max_admm_iter = 3000,
                     irls_max_iter = 100, irls_tol = 1e-6) {
  
  z <- if (is.null(z_init)) rep(0, p) else z_init
  u_list <- if (is.null(u_init)) lapply(1:m, function(j) rep(0, p)) else u_init
  
  converged <- FALSE
  
  for (iter in 1:max_admm_iter) {
    z_prev <- z
    
    beta_list <- lapply(1:m, function(k)
      admm_local_irls(
        clients_train[[k]]$X,
        clients_train[[k]]$y,
        z,
        u_list[[k]],
        rho,
        irls_max_iter = irls_max_iter,
        tol = irls_tol
      )
    )
    
    avg_bu <- Reduce(`+`, Map(`+`, beta_list, u_list)) / m
    z_new  <- soft_threshold(avg_bu, lambda / (m * rho))
    u_list <- lapply(1:m, function(k) u_list[[k]] + beta_list[[k]] - z_new)
    
    chk <- admm_check(beta_list, z_new, z_prev, u_list, rho, eps_abs, eps_rel, p, m)
    z <- z_new
    
    if (chk$r_norm < chk$eps_primal && chk$s_norm < chk$eps_dual) {
      converged <- TRUE
      break
    }
  }
  
  list(
    z = z,
    u_list = u_list,
    iter = iter,
    converged = converged
  )
}


# Main simulation: Federated ADMM-IRLS 
out_list <- vector("list", length = length(rho_values) * (rep_end - rep_start + 1))
ptr <- 1

for (rho in rho_values) {
  for (rep in rep_start:rep_end) {
    
    set.seed(1000 + 100 * match(rho, rho_values) + rep)
    
    # train, validation and test dataset
    clients_train <- lapply(1:m, function(j) simulate_data(n, p, true_beta))
    clients_val   <- lapply(1:m, function(j) simulate_data(n, p, true_beta))
    clients_test  <- lapply(1:m, function(j) simulate_data(n, p, true_beta))
    
    # ADMM-IRLS: tune lambda by VALIDATION AUC
    best_lambda_admm <- LAMBDA_GRID[1]
    best_val_admm <- -Inf
    
    z_ws <- rep(0, p)
    u_ws <- lapply(1:m, function(j) rep(0, p))
    best_z <- z_ws
    best_u <- u_ws
    
    for (lam in LAMBDA_GRID) {
      fit_lam <- run_admm(clients_train,lam,rho,p,m,eps_abs,eps_rel,z_init = z_ws,u_init = u_ws,max_admm_iter = 3000)
      
      z_ws <- fit_lam$z
      u_ws <- fit_lam$u_list
      
      # VALIDATION AUC
      v <- val_auc_mean_clients(clients_val, z_ws)
      
      if (!is.na(v) && v > best_val_admm) {
        best_val_admm <- v
        best_lambda_admm <- lam
        best_z <- z_ws
        best_u <- u_ws
      }
    }
    
    # final ADMM at best lambda
    fit_admm <- run_admm(clients_train,best_lambda_admm,rho,p,m,eps_abs,eps_rel,z_init = best_z,u_init = best_u,max_admm_iter = 3000)
    
    beta_admm <- fit_admm$z
    admm_converged <- fit_admm$converged
    admm_iter <- fit_admm$iter
    
    # Evaluate ADMM-IRLS on test sets
    beta_est <- beta_admm
    best_lam <- best_lambda_admm
    
    est_active  <- abs(beta_est) > threshold
    true_active <- abs(true_beta) > threshold
    
    TP <- sum(est_active & true_active)
    FP <- sum(est_active & !true_active)
    FN <- sum(!est_active & true_active)
    TN <- sum(!est_active & !true_active)
    
    Sensitivity <- if ((TP + FN) > 0) TP / (TP + FN) else 0
    Specificity <- if ((TN + FP) > 0) TN / (TN + FP) else 0
    
    MSE <- mean((beta_est - true_beta)^2)
    MAE <- mean(abs(beta_est - true_beta))
    
    test_prob_mse <- mean(sapply(1:m, function(k)
      val_prob_mse(clients_test[[k]]$X, clients_test[[k]]$y, beta_est)))
    
    client_aucs <- sapply(1:m, function(k) {
      score <- as.vector(plogis(clients_test[[k]]$X %*% beta_est))
      safe_auc(clients_test[[k]]$y, score)
    })
    
    AUC <- mean(client_aucs, na.rm = TRUE)
    
    out_list[[ptr]] <- data.frame(Rho = rho,Replication = rep,Model = "ADMM-IRLS",Best_Lambda = best_lam,TP = TP,FP = FP,TN = TN,
      FN = FN,Sensitivity = Sensitivity,Specificity = Specificity,MSE = MSE,MAE = MAE,Test_ProbMSE = test_prob_mse,AUC = AUC,
      ADMM_Converged = admm_converged,ADMM_Iter = admm_iter)
    
    ptr <- ptr + 1
    
    cat("Finished rho =", rho, "| replication =", rep, "\n")
  }
}

results_df <- dplyr::bind_rows(out_list)

# Summary results
summary_results <- results_df %>%
  group_by(Rho) %>%
  summarise(Mean_TP = mean(TP),SD_TP = sd(TP),Mean_FP = mean(FP),SD_FP = sd(FP),Mean_TN = mean(TN),SD_TN = sd(TN),Mean_FN = mean(FN),
    SD_FN = sd(FN),Mean_SE = mean(Sensitivity),SD_SE = sd(Sensitivity),Mean_SP = mean(Specificity),SD_SP = sd(Specificity),
    Mean_MSE = mean(MSE),SD_MSE = sd(MSE),Mean_MAE = mean(MAE), SD_MAE = sd(MAE), Mean_Test_ProbMSE = mean(Test_ProbMSE),
    SD_Test_ProbMSE = sd(Test_ProbMSE),Mean_AUC = mean(AUC, na.rm = TRUE),SD_AUC = sd(AUC, na.rm = TRUE),
    Mean_Best_Lambda = mean(Best_Lambda),SD_Best_Lambda = sd(Best_Lambda),Convergence_Rate = mean(ADMM_Converged),
    Mean_ADMM_Iter = mean(ADMM_Iter),SD_ADMM_Iter = sd(ADMM_Iter),.groups = "drop")

print(summary_results)

write.csv(summary_results, "results/admm_irls_windows_summary_by_rho.csv", row.names = FALSE)

cat("\nSaved summary by rho to:\n")
cat("results/admm_irls_windows_summary_by_rho.csv\n")


