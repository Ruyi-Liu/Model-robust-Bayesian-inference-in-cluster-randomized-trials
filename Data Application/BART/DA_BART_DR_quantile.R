print("DA,BART,DR,quantile,9.27")
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



mc_num <- 10000
burnin <- 5000
thin_rate <- 1

BART_calculate_deltas_DR <- function(data_resample) {
  dt <- data.table(
    cluster_id = data_resample$new_cluster_id,
    Y_1_pred_individual = data_resample$Y_1_pred_individual,
    Y_0_pred_individual = data_resample$Y_0_pred_individual,
    Y_true_obs = data_resample$Y,
    trt_ori = data_resample$A_i,
    cluster_size = data_resample$N_freq
  )
  
  # Aggregate once per cluster
  means <- dt[, .(
    mean_Y_1_pred = mean(Y_1_pred_individual),
    mean_Y_0_pred = mean(Y_0_pred_individual),
    mean_Y_obs = mean(Y_true_obs),
    trt_ori = first(trt_ori),
    cluster_size_N = first(cluster_size)
  ), by = cluster_id]
  
  # Precompute indicators
  cluster_trt_indicator_1 <- as.numeric(means$trt_ori == 1)
  cluster_trt_indicator_0 <- 1 - cluster_trt_indicator_1
  
  # Cluster-level estimator
  mu_c1 <- mean(means$mean_Y_1_pred + 2 * cluster_trt_indicator_1 * (means$mean_Y_obs - means$mean_Y_1_pred))
  mu_c0 <- mean(means$mean_Y_0_pred + 2 * cluster_trt_indicator_0 * (means$mean_Y_obs - means$mean_Y_0_pred))
  Delta_C <- mu_c1 - mu_c0
  
  # Individual-level estimator
  mu_I1 <- sum((means$mean_Y_1_pred + 2 * cluster_trt_indicator_1 * (means$mean_Y_obs - means$mean_Y_1_pred)) * means$cluster_size_N) /
    sum(means$cluster_size_N)
  mu_I0 <- sum((means$mean_Y_0_pred + 2 * cluster_trt_indicator_0 * (means$mean_Y_obs - means$mean_Y_0_pred)) * means$cluster_size_N) /
    sum(means$cluster_size_N)
  Delta_I <- mu_I1 - mu_I0
  
  return(list(Delta_C_list = Delta_C, Delta_I_list = Delta_I))
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

Delta_C_list <- numeric(n_post_draws)
Delta_I_list <- numeric(n_post_draws)

for (draw_index in 1:n_post_draws) {
  data_w_Y_pred <- data
  data_w_Y_pred$new_cluster_id <- data$cluster_id
  data_w_Y_pred$Y_1_pred_individual <- result$y1_post[, draw_index]
  data_w_Y_pred$Y_0_pred_individual <- result$y0_post[, draw_index]
  
  deltas <- BART_calculate_deltas_DR(data_w_Y_pred)
  Delta_C_list[draw_index] <- deltas$Delta_C_list
  Delta_I_list[draw_index] <- deltas$Delta_I_list
}


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
upper_CI <- quantile(Delta_C_list, 0.975)
lower_CI <- quantile(Delta_C_list, 0.025)
print("--------------c-ATE---------------")
cat("Mean c-ATE: ", formatC(mean(Delta_C_list), format = "f", digits = 3), "\n")
cat("The credible interval for c-ATE is (", 
    formatC(lower_CI, format = "f", digits = 3), ", ", 
    formatC(upper_CI, format = "f", digits = 3), ")\n", sep = "")


##############################################################
upper_CI <- quantile(Delta_I_list, 0.975)
lower_CI <- quantile(Delta_I_list, 0.025)
print("--------------I-ATE---------------")
cat("Mean i-ATE: ", formatC(mean(Delta_I_list), format = "f", digits = 3), "\n")
cat("The credible interval for i-ATE is (", 
    formatC(lower_CI, format = "f", digits = 3), ", ", 
    formatC(upper_CI, format = "f", digits = 3), ")\n", sep = "")

