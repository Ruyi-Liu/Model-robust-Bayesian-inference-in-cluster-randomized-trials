print("da,DR,quantile,7.28")
library(dplyr)
library(haven)
library(MASS)
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

#View(sorted_data)

data <- sorted_data
data$cluster_id <- match(data$cluster_id, unique(data$cluster_id))


mc_num <- 2000
burnin <- 1000


calculate_deltas_DR <- function(data_resample, BETA_post, draw_index) {
  N <- nrow(data_resample)
  
  X_1<-cbind(rep(1,N), rep(1,N))
  
  X_0<-cbind(rep(1,N), rep(0,N))
  
  Y_1_pred_individual <- X_1 %*% BETA_post[draw_index, ]
  Y_0_pred_individual <- X_0 %*% BETA_post[draw_index, ]
  
  df_new <- data.frame(cluster_id = data_resample$cluster_id,
                       Y_1_pred_individual = Y_1_pred_individual,
                       Y_0_pred_individual = Y_0_pred_individual,
                       Y_true_obs = data_resample$Y,
                       trt_ori = data_resample$A_i,
                       cluster_size = data_resample$N_freq)
  
  means <- df_new %>%
    group_by(cluster_id) %>%
    summarise(mean_Y_1_pred = mean(Y_1_pred_individual),
              mean_Y_0_pred = mean(Y_0_pred_individual),
              mean_Y_obs = mean(Y_true_obs),
              trt_ori = first(trt_ori),
              cluster_size_N = first(cluster_size))
  
  cluster_trt_indicator_1 <- as.numeric(means$trt_ori == 1)
  cluster_trt_indicator_0 <- 1 - cluster_trt_indicator_1
  
  mu_c1 <- mean((means$mean_Y_1_pred) + 2 * cluster_trt_indicator_1 * (means$mean_Y_obs - means$mean_Y_1_pred))
  mu_c0 <- mean((means$mean_Y_0_pred) + 2 * cluster_trt_indicator_0 * (means$mean_Y_obs - means$mean_Y_0_pred))
  Delta_C <- mu_c1 - mu_c0
  
  mu_I1 <- sum(((means$mean_Y_1_pred) + 2 * cluster_trt_indicator_1 * (means$mean_Y_obs - means$mean_Y_1_pred)) *
                 (means$cluster_size_N)) / sum(means$cluster_size_N)
  mu_I0 <- sum(((means$mean_Y_0_pred) + 2 * cluster_trt_indicator_0 * (means$mean_Y_obs - means$mean_Y_0_pred)) *
                 (means$cluster_size_N)) / sum(means$cluster_size_N)
  Delta_I <- mu_I1 - mu_I0
  
  return(list(Delta_C_list = Delta_C, Delta_I_list = Delta_I))
}

Y<-data$Y
N<-length(Y)
# Design matrix
X<-cbind(rep(1,N), data$A_i)


p<-dim(X)[2]
id<-data$cluster_id
m<-length(unique(data$cluster_id))

# Design matrix of random effects
Z<-matrix(0,N,m)
for(i in 1:m){
  q<-which(data$cluster_id==i)
  Z[q,i]<-1
}   

# Initialization
sigmab2<-3
sigma2<- 10
B<-rep(0,m)

# priors
mu0<-rep(0,p)
L0<-diag(rep(100,p))
nu0<-2*0.001 # Inverse-Wishart degenerates at Inverse-Gamma
D<-2*matrix(0.001)
a<-b<-0.001



# MCMC

BETA.post <- matrix(NA, nrow = mc_num, ncol = p)
SIGMAb2 <- numeric(mc_num)
SIGMA2 <- numeric(mc_num)
B.post <- matrix(NA, nrow = mc_num, ncol = m)

for(s in 1:mc_num){
  # update beta-fixed effects
  Lm<-solve(solve(L0)+t(X)%*%X/sigma2)
  mum<-Lm%*%(solve(L0)%*%mu0+t(X)%*%(Y-Z%*%B)/sigma2)
  beta <- mvrnorm(n = 1, mu = mum, Sigma = Lm)
  
  ##
  
  # update b_i-random intercept
  for(j in 1:m){
    q<-which(id==j)
    Lamdam<-solve(solve(sigmab2)+t(Z[q,j])%*%Z[q,j]/sigma2)
    taum<-Lamdam%*%(t(Z[q,j])%*%(Y[q]-X[q,]%*%beta)/sigma2)
    B[j]<-c(mvrnorm(n = 1, mu = taum, Sigma = Lamdam))
  }
  
  # update random intercept variance (covariance matrix)
  D0<-t(B)%*%B
  sigmab2<-1/rgamma(1,shape=(nu0/2+m/2),rate=((D+D0)/2))
  
  # update sigma^2
  ss=t(Y-X%*%beta-Z%*%B)%*%(Y-X%*%beta-Z%*%B)
  sigma2<-1/rgamma(1,a+N/2,b+ss/2)
  ##
  
  # store some output
  BETA.post[s, ] <- beta
  SIGMAb2[s] <- sigmab2
  SIGMA2[s] <- sigma2
  B.post[s, ] <- B

}

Delta_C_list <- c()
Delta_I_list <- c()

for (draw_index in (burnin+1):mc_num){
  Delta_D_obs <- calculate_deltas_DR(data, BETA.post, draw_index) 
  Delta_C <- Delta_D_obs$Delta_C_list
  Delta_I <- Delta_D_obs$Delta_I_list
  Delta_C_list <- c(Delta_C_list,Delta_C)
  Delta_I_list <- c(Delta_I_list,Delta_I)
}



library(coda)
# Convert post-burnin samples to mcmc objects
iATE_mcmc <- as.mcmc(Delta_I_list)
traceplot(iATE_mcmc, main = "Traceplot for i-ATE")

# Geweke diagnostic
geweke_i <- geweke.diag(iATE_mcmc)

# Print Z-scores


print("Geweke Z-scores for i-ATE:")
print(geweke_i$z)


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

