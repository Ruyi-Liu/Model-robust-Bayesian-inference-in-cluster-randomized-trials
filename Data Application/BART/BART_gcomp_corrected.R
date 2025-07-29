print("DA,BART,gform,antonelli,7.22")
library(dplyr)
library(haven)
library(MASS)
library(BART3)
library(mxBART)
library(data.table)
set.seed(999)


data <- read_sas("/Users/ruyi/Desktop/Bayesian CRTs Project/Data_analysis_sep21/ppact_public_bpi_062623.sas7bdat")

# View the data
#View(data)

#filtered_data <- data[data$TIMEPOINT == 12, ]

# Define the columns you want to keep
columns_to_keep <- c("PEGS","SID", "CLUST", "INTERVENTION",
                     "AGE", "FEMALE", "disable", "Current_Smoke", "BMI",
                     "Alcohol_Abuse","Drug_Abuse","comorbid","Depression",
                     "pain_count","BL_avg_daily","BL_avg_above90")

# Filter the data where TIMEPOINT is 12 and select only the specified columns
filtered_data <- data[data$TIMEPOINT == 12, columns_to_keep]

#View(filtered_data)

clean_data <- na.omit(filtered_data)

clean_data_new <- clean_data %>%
  group_by(CLUST) %>%
  mutate(N_freq = n()) %>%
  ungroup()

sorted_data <- clean_data_new %>%
  arrange(CLUST)



# Set new names for all columns
names(sorted_data) <- c("Y", "individual_id","cluster_id","A_i",
                        "X_ij1", "X_ij2", "X_ij3", "X_ij4", "X_ij5", "X_ij6",
                        "X_ij7", "X_ij8", "X_ij9", "X_ij10", "X_ij11", "X_ij12",
                        "N_freq")



data <- sorted_data
data$cluster_id <- match(data$cluster_id, unique(data$cluster_id))



# mc_num <- 10000
# burnin <- 5000
# thin_rate <- 1
# num_new_dataset = 100

mc_num <- 10000
burnin <- 5000
thin_rate <- 1
num_new_dataset = 100

BART_calculate_deltas_gform <- function(data_resample) {
  
  Y_1_pred <- data_resample$Y_1_pred_individual
  Y_0_pred <- data_resample$Y_0_pred_individual
  df_new <- data.frame(cluster_id=data_resample$new_cluster_id,Y_1_pred=Y_1_pred,Y_0_pred=Y_0_pred)
  
  means <- aggregate(cbind(Y_1_pred, Y_0_pred) ~ cluster_id, data = df_new, FUN = mean)
  # Print the resulting means
  mu_c1 <- means$Y_1_pred
  mu_c0 <- means$Y_0_pred
  Delta_C <- mean(mu_c1-mu_c0)
  
  sums <- aggregate(cbind(Y_1_pred, Y_0_pred) ~ cluster_id, data = df_new, FUN = sum)
  mu_I1 <- sum(sums$Y_1_pred)/(length(data_resample$Y))
  mu_I0 <- sum(sums$Y_0_pred)/(length(data_resample$Y))
  Delta_I <- mu_I1-mu_I0
  
  
  return(list(Delta_C_list = Delta_C, Delta_I_list = Delta_I))
}


resample_clusters <- function(df) {
  dt <- as.data.table(df)
  
  # Pre-split by cluster_id
  cluster_list <- split(dt, by = "cluster_id", keep.by = TRUE)
  
  # Sample cluster names with replacement
  cluster_ids <- names(cluster_list)
  sampled_ids <- sample(cluster_ids, length(cluster_ids), replace = TRUE)
  
  # Efficiently relabel and rbind
  resampled <- rbindlist(lapply(seq_along(sampled_ids), function(new_id) {
    cluster_copy <- copy(cluster_list[[sampled_ids[new_id]]])
    cluster_copy[, new_cluster_id := new_id]
    return(cluster_copy)
  }))
  
  return(resampled)
}


data <- data[order(data$cluster_id), ]

mxbart_func <- function(df,
                        outcome,      # outcome Y
                        treatment,    # treatment (character, 0/1)
                        covariates,   # covariates (character vector)
                        id,           # cluster ID column names (character)
                        ntree     = 50L,    # BART trees
                        nskip     = 5000L,  # MCMC burn‐in
                        ndpost    = 10000L,  # MCMC posterior sample
                        keepevery = 1L  
) {
  # 1. check columns
  stopifnot(all(c(outcome, treatment, id, covariates) %in% names(df)))
  
  # 2. sorting
  df        <- df[order(df[[id]]), ]
  y.train   <- df[[outcome]]
  x.train   <- as.matrix(df[, c(treatment, covariates)])
  id.train  <- as.integer(df[[id]])
  n         <- nrow(df)
  
  # 3. counterfactual matrix for S-learner
  df0       <- df; df0[[treatment]] <- 0
  df1       <- df; df1[[treatment]] <- 1
  x0        <- as.matrix(df0[, c(treatment, covariates)])
  x1        <- as.matrix(df1[, c(treatment, covariates)])
  x.test    <- rbind(x0, x1)   # A=0 then A=1
  
  
  fit <- mxbart(
    x.train   = x.train,
    y.train   = y.train,
    id.train  = id.train,
    x.test    = x.test,
    ntree     = ntree,
    nskip     = nskip,
    ndpost    = ndpost,
    keepevery = keepevery
  )
  
  # 5. posterior prediction
  posterior <- fit$fhat.test
  y0_post   <- posterior[,       1:n]      # A=0
  y1_post   <- posterior[, (n+1):(2*n)]   # A=1
  
  
  list(
    y0_post = t(y0_post),
    y1_post = t(y1_post)
  )
  
  
  
}


result <- mxbart_func(
  df         = data,
  outcome    = "Y",
  treatment  = "A_i",
  covariates = c("X_ij1", "X_ij2", "X_ij3", "X_ij4",
                 "X_ij5", "X_ij6", "X_ij7", "X_ij8",
                 "X_ij9", "X_ij10", "X_ij11", "X_ij12",
                 "N_freq"),
  id         = "cluster_id",
  ntree      = 50L,
  nskip      = burnin,
  ndpost     = mc_num - burnin,
  keepevery  = thin_rate
)


n_post_draws <- mc_num - burnin

# Preallocate result storage
Delta_C_list <- numeric(n_post_draws)
Delta_I_list <- numeric(n_post_draws)
Variance_matrix_delta_C <- matrix(NA, nrow = num_new_dataset, ncol = n_post_draws)
Variance_matrix_delta_I <- matrix(NA, nrow = num_new_dataset, ncol = n_post_draws)

# Step 1: Resample datasets once
resampled_datasets <- lapply(1:num_new_dataset, function(index_resampling) {
  resample_clusters(data)
})


# Step 2: Loop over posterior draws
for (draw_index in 1:n_post_draws) {
  
  cat("Now processing posterior draw:", draw_index, "of", n_post_draws, "\n")
  
  
  # Posterior predictions for this draw
  y1_pred <- result$y1_post[, draw_index]
  y0_pred <- result$y0_post[, draw_index]
  
  # Add predictions to original dataset for full-data delta
  data_w_Y_pred <- data
  data_w_Y_pred$new_cluster_id <- data$cluster_id  # or any appropriate cluster ID
  data_w_Y_pred$Y_1_pred_individual <- y1_pred
  data_w_Y_pred$Y_0_pred_individual <- y0_pred
  
  # Calculate full-data deltas
  deltas <- BART_calculate_deltas_gform(data_w_Y_pred)
  Delta_C_list[draw_index] <- deltas$Delta_C_list
  Delta_I_list[draw_index] <- deltas$Delta_I_list
  
  # Build a fast lookup table of posterior predictions
  pred_df <- data.frame(cluster_id = data$cluster_id,
                        individual_id = data$individual_id,
                        Y_1_pred_individual = y1_pred,
                        Y_0_pred_individual = y0_pred)
  
  id_key_pred <- paste(pred_df$cluster_id, pred_df$individual_id)
  
  # Loop over resampled datasets and inject predictions
  for (num_new_dataset_index in 1:num_new_dataset) {
    merged_df <- resampled_datasets[[num_new_dataset_index]]
    
    # Use match() to inject predictions
    id_key_merged <- paste(merged_df$cluster_id, merged_df$individual_id)
    idx <- match(id_key_merged, id_key_pred)
    
    merged_df$Y_1_pred_individual <- pred_df$Y_1_pred_individual[idx]
    merged_df$Y_0_pred_individual <- pred_df$Y_0_pred_individual[idx]
    
    # Optional: sort if your delta function depends on cluster order
    merged_df <- merged_df[order(merged_df$new_cluster_id), ]
    
    # Calculate deltas for the resampled dataset
    ret_delta_resampled <- BART_calculate_deltas_gform(merged_df)
    Variance_matrix_delta_C[num_new_dataset_index, draw_index] <- ret_delta_resampled$Delta_C_list
    Variance_matrix_delta_I[num_new_dataset_index, draw_index] <- ret_delta_resampled$Delta_I_list
  }
}



row_means_E_deltaC <- apply(Variance_matrix_delta_C, 1, mean)
row_means_E_deltaI <- apply(Variance_matrix_delta_I, 1, mean)

Var_deltaC <- var(row_means_E_deltaC)+var(Delta_C_list)
Var_deltaI <- var(row_means_E_deltaI)+var(Delta_I_list)






library(coda)
# Convert post-burnin samples to mcmc objects
cATE_mcmc <- as.mcmc(Delta_C_list)
iATE_mcmc <- as.mcmc(Delta_I_list)

# Plot traceplots
traceplot(cATE_mcmc, main = "Traceplot for c-ATE")
traceplot(iATE_mcmc, main = "Traceplot for i-ATE")

# Geweke diagnostic
geweke_c <- geweke.diag(cATE_mcmc)
geweke_i <- geweke.diag(iATE_mcmc)

# Print Z-scores
print("Geweke Z-scores for c-ATE:")
print(geweke_c$z)

print("Geweke Z-scores for i-ATE:")
print(geweke_i$z)


if (any(abs(geweke_c$z) > 1.96)) {
  cat("Warning: Some Geweke Z-scores for c-ATE are outside the 95% interval [-1.96, 1.96].\n")
} else {
  cat("All Geweke Z-scores for c-ATE are within the 95% interval [-1.96, 1.96].\n")
}

if (any(abs(geweke_i$z) > 1.96)) {
  cat("Warning: Some Geweke Z-scores for i-ATE are outside the 95% interval [-1.96, 1.96].\n")
} else {
  cat("All Geweke Z-scores for i-ATE are within the 95% interval [-1.96, 1.96].\n")
}


############# Performance Metric ###########################
upper_CI_cATE <- mean(Delta_C_list) + 1.96*sqrt(Var_deltaC)
lower_CI_cATE <- mean(Delta_C_list) - 1.96*sqrt(Var_deltaC)
print("--------------c-ATE---------------")
cat("Mean c-ATE: ", formatC(mean(Delta_C_list), format = "f", digits = 3), "\n")
cat("The credible interval for c-ATE is (", 
    formatC(lower_CI_cATE, format = "f", digits = 3), ", ", 
    formatC(upper_CI_cATE, format = "f", digits = 3), ")\n", sep = "")


##############################################################
upper_CI_iATE <- mean(Delta_I_list) + 1.96*sqrt(Var_deltaI)
lower_CI_iATE <- mean(Delta_I_list) - 1.96*sqrt(Var_deltaI)
print("--------------I-ATE---------------")
cat("Mean i-ATE: ", formatC(mean(Delta_I_list), format = "f", digits = 3), "\n")
cat("The credible interval for i-ATE is (", 
    formatC(lower_CI_iATE, format = "f", digits = 3), ", ", 
    formatC(upper_CI_iATE, format = "f", digits = 3), ")\n", sep = "")

