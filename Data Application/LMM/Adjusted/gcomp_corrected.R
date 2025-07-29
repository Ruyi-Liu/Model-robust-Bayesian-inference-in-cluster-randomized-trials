print("da,gform,antonelli,9.21")
library(dplyr)
library(haven)
library(MASS)
set.seed(999)

#data <- read_sas("/Users/ruyi/Desktop/ppact_public_bpi_062623.sas7bdat")


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


#mc_num <- 2000
#burnin <- 1000
#num_new_dataset = 100
mc_num <- 2000
burnin <- 1000
num_new_dataset = 100


Variance_matrix_delta_C <- matrix(NA, nrow = num_new_dataset, ncol = mc_num-burnin)
Variance_matrix_delta_I <- matrix(NA, nrow = num_new_dataset, ncol = mc_num-burnin)


calculate_deltas_gform <- function(data_resample, BETA_post, draw_index) {
  N <- nrow(data_resample)
  
  X_1<-cbind(rep(1,N), rep(1,N), data_resample$X_ij1, data_resample$X_ij2,data_resample$X_ij3, data_resample$X_ij4,
           data_resample$X_ij5, data_resample$X_ij6,data_resample$X_ij7, data_resample$X_ij8,
           data_resample$X_ij9, data_resample$X_ij10,data_resample$X_ij11, data_resample$X_ij12,data_resample$N_freq, 
           rep(1,N)*data_resample$X_ij1,rep(1,N)*data_resample$X_ij2,
           rep(1,N)*data_resample$X_ij3,rep(1,N)*data_resample$X_ij4,
           rep(1,N)*data_resample$X_ij5,rep(1,N)*data_resample$X_ij6,
           rep(1,N)*data_resample$X_ij7,rep(1,N)*data_resample$X_ij8,
           rep(1,N)*data_resample$X_ij9,rep(1,N)*data_resample$X_ij10,
           rep(1,N)*data_resample$X_ij11,rep(1,N)*data_resample$X_ij12,
           rep(1,N)*data_resample$N_freq)
  
  X_0<-cbind(rep(1,N), rep(0,N), data_resample$X_ij1, data_resample$X_ij2,data_resample$X_ij3, data_resample$X_ij4,
           data_resample$X_ij5, data_resample$X_ij6,data_resample$X_ij7, data_resample$X_ij8,
           data_resample$X_ij9, data_resample$X_ij10,data_resample$X_ij11, data_resample$X_ij12,data_resample$N_freq, 
           rep(0,N)*data_resample$X_ij1,rep(0,N)*data_resample$X_ij2,
           rep(0,N)*data_resample$X_ij3,rep(0,N)*data_resample$X_ij4,
           rep(0,N)*data_resample$X_ij5,rep(0,N)*data_resample$X_ij6,
           rep(0,N)*data_resample$X_ij7,rep(0,N)*data_resample$X_ij8,
           rep(0,N)*data_resample$X_ij9,rep(0,N)*data_resample$X_ij10,
           rep(0,N)*data_resample$X_ij11,rep(0,N)*data_resample$X_ij12,
           rep(0,N)*data_resample$N_freq)
  
  Y_1_pred <- X_1%*%BETA.post[draw_index,]
  Y_0_pred <- X_0%*%BETA.post[draw_index,]
  df_new <- data.frame(cluster_id=data_resample$cluster_id,Y_1_pred=Y_1_pred,Y_0_pred=Y_0_pred)
  
  means <- aggregate(cbind(Y_1_pred, Y_0_pred) ~ cluster_id, data = df_new, FUN = mean)
  # Print the resulting means
  head(means)
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
  # Get a unique list of cluster IDs
  cluster_ids <- unique(df$cluster_id)
  
  # Sample cluster IDs with replacement
  sampled_cluster_ids <- sample(cluster_ids, size = length(cluster_ids), replace = TRUE)
  #print(sampled_cluster_ids)
  # Resample full clusters and create a list of data frames
  resampled_list <- lapply(seq_along(sampled_cluster_ids), function(new_id) {
    original_id <- sampled_cluster_ids[new_id]
    cluster_data <- df[df$cluster_id == original_id, ]
    cluster_data$cluster_id <- new_id  # Assign the new cluster ID directly
    return(cluster_data)
  })
  
  # Combine the list of data frames into one data frame
  resampled_df <- do.call(rbind, resampled_list)
  
  # Reset row names if necessary
  rownames(resampled_df) <- NULL
  
  return(resampled_df)
}

Y<-data$Y
N<-length(Y)
# Design matrix
X<-cbind(rep(1,N), data$A_i, data$X_ij1, data$X_ij2,data$X_ij3, data$X_ij4,
         data$X_ij5, data$X_ij6,data$X_ij7, data$X_ij8,
         data$X_ij9, data$X_ij10,data$X_ij11, data$X_ij12,data$N_freq, 
         data$A_i*data$X_ij1,data$A_i*data$X_ij2,
         data$A_i*data$X_ij3,data$A_i*data$X_ij4,
         data$A_i*data$X_ij5,data$A_i*data$X_ij6,
         data$A_i*data$X_ij7,data$A_i*data$X_ij8,
         data$A_i*data$X_ij9,data$A_i*data$X_ij10,
         data$A_i*data$X_ij11,data$A_i*data$X_ij12,
         data$A_i*data$N_freq)

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

resampled_datasets <- list()  # Initialize an empty list to store resampled datasets
for (index_resampling in 1:num_new_dataset) {
  resampled_datasets[[index_resampling]] <- resample_clusters(data)  # Resample and store
}

for (draw_index in (burnin+1):mc_num){
  Delta_D_obs <- calculate_deltas_gform(data, BETA.post, draw_index) 
  Delta_C <- Delta_D_obs$Delta_C_list
  Delta_I <- Delta_D_obs$Delta_I_list
  Delta_C_list <- c(Delta_C_list,Delta_C)
  Delta_I_list <- c(Delta_I_list,Delta_I)
  
  
  #stop("Debug end here")
  for (num_new_dataset_index in 1:num_new_dataset){
    resampled_dataset_a <- resampled_datasets[[num_new_dataset_index]]
    ret_delta_resampled <- calculate_deltas_gform(resampled_dataset_a, BETA.post, draw_index)
    Variance_matrix_delta_C[num_new_dataset_index,(draw_index-burnin)] <- ret_delta_resampled$Delta_C_list
    Variance_matrix_delta_I[num_new_dataset_index,(draw_index-burnin)] <- ret_delta_resampled$Delta_I_list
  }
  
}
row_means_E_deltaC <- apply(Variance_matrix_delta_C, 1, mean)
row_means_E_deltaI <- apply(Variance_matrix_delta_I, 1, mean)

Var_deltaC <- var(row_means_E_deltaC)+var(Delta_C_list)
Var_deltaI <- var(row_means_E_deltaI)+var(Delta_I_list)


############# Performance Metric ###########################
upper_CI <- mean(Delta_C_list) + 1.96*sqrt(Var_deltaC)
lower_CI <- mean(Delta_C_list) - 1.96*sqrt(Var_deltaC)
print("--------------c-ATE---------------")
cat("Mean c-ATE: ", formatC(mean(Delta_C_list), format = "f", digits = 3), "\n")
cat("The credible interval for c-ATE is (", 
    formatC(lower_CI, format = "f", digits = 3), ", ", 
    formatC(upper_CI, format = "f", digits = 3), ")\n", sep = "")

##############################################################
upper_CI <- mean(Delta_I_list) + 1.96*sqrt(Var_deltaI)
lower_CI <- mean(Delta_I_list) - 1.96*sqrt(Var_deltaI)
print("--------------I-ATE---------------")
cat("Mean i-ATE: ", formatC(mean(Delta_I_list), format = "f", digits = 3), "\n")
cat("The credible interval for i-ATE is (", 
    formatC(lower_CI, format = "f", digits = 3), ", ", 
    formatC(upper_CI, format = "f", digits = 3), ")\n", sep = "")
