###################################################################################################
# Basic setting
# -------------------------------------------------------------------------------------------------

# directory ---------------------------------------------------------------------------------------
setwd("C:\\Users\\yangd\\OneDrive\\0. Research\\2. JQT\\Revision\\JQT\\")
# -------------------------------------------------------------------------------------------------

# library -----------------------------------------------------------------------------------------
library(dplyr)
library(reshape2)
library(tidyr)
library(dqrng)
library(LaplacesDemon)
library(fields)
library(foreach)
library(doParallel)
# -------------------------------------------------------------------------------------------------

# source ------------------------------------------------------------------------------------------
source("Functions.R")
# -------------------------------------------------------------------------------------------------

# load data ---------------------------------------------------------------------------------------
load(file="data\\data_wf939.RData")
# -------------------------------------------------------------------------------------------------

# data --------------------------------------------------------------------------------------------
DATA         <- Data__WF939
dat_dim      <- DATA[[1]] %>% dim
DATA_mat     <- do.call('cbind', sapply(seq(DATA), function(ii){ c(DATA[[ii]]) }, simplify = F))
DATA_mat_org <- DATA_mat 

WF_outise_idx <- which(DATA_mat[,1] == 0)
DATA_mat0     <- DATA_mat[-WF_outise_idx,]
DATA_mat_org0 <- DATA_mat_org[-WF_outise_idx,]
# -------------------------------------------------------------------------------------------------

# grid --------------------------------------------------------------------------------------------
WF_center <- ((dat_dim/2 + 0.5)/dat_dim) %>% rev
X_data <- X_data_polar <- NULL 
for(x in 1:dat_dim[2]){
  for(y in 1:dat_dim[1]){
    
    tmp_x   <- x/dat_dim[2]
    tmp_y   <- (dat_dim[1] - y + 1)/dat_dim[1]
    tmp_vec <- (c(tmp_x, tmp_y) - WF_center)  * 2
    
    tmp_r     <- sum(tmp_vec^2) %>% sqrt
    tmp_theta <- atan2(tmp_vec[2], tmp_vec[1])
    
    X_data_polar <- rbind(X_data_polar, c(tmp_r,tmp_theta))
    X_data       <- rbind(X_data, tmp_vec)
  }
}
X_data       <- X_data[-WF_outise_idx,]
X_data_polar <- X_data_polar[-WF_outise_idx,]

## WF grid
x_grid <- X_data[,1] %>% unique %>% sort
y_grid <- X_data[,2] %>% unique %>% sort
# -------------------------------------------------------------------------------------------------

# settings for Zernike polynomial -----------------------------------------------------------------
## setting for zernike polynomial
nMax <- 4

## basis set
nm_list <- zernike_basis_set(nMax)  

## basis matrix
n_data <- nrow(X_data)
basis_count <- length(nm_list)
Bmat <- matrix(0, nrow=n_data, ncol=basis_count)

for(coli in seq_along(nm_list)){
  n_ <- nm_list[[coli]][1]
  m_ <- nm_list[[coli]][2]
  Bmat[, coli] <- zernike_eval_one(n_, m_, X_data_polar[,1], X_data_polar[,2])
}
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
#
###################################################################################################


###################################################################################################
# load Phase I result
# -------------------------------------------------------------------------------------------------

# Sampling ----------------------------------------------------------------------------------------
tmp_idx       <- DATA_mat0 %>% ncol %>% seq
tmp_DATA_mat0 <- DATA_mat0[,tmp_idx]
# -------------------------------------------------------------------------------------------------

# data transformation to continuous data ----------------------------------------------------------
Ph1_DataN <- 50
Ph1_DATA0 <- Ph1_DATA <- list()
for(i in 1:Ph1_DataN){
  
  set.seed(i)
  tmp_DATA_mat   <- data_transform_cont(tmp_DATA_mat0)
  tmp_step1_coef <- get_Zernike(tmp_DATA_mat, Bmat)
  
  Ph1_DATA0[[i]] <- tmp_DATA_mat
  Ph1_DATA[[i]]  <- tmp_step1_coef
}

## data
total_Y <- NULL
for(i in 1:Ph1_DataN){ total_Y <- rbind(total_Y, Ph1_DATA[[i]]) }

## scaling
Step1_scale_center <- total_Y %>% scale %>% attr("scaled:center")
Step1_scale_scale  <- total_Y %>% scale %>% attr("scaled:scale")
total_Yscale       <- total_Y %>% scale # ((t(total_Y) - Step1_scale_center)/Step1_scale_scale) %>% t

## basic setting
each_n  <- nrow(Ph1_DATA[[1]])
K       <- ncol(Ph1_DATA[[1]])
# -------------------------------------------------------------------------------------------------

load("Ph1_8.RData")
# -------------------------------------------------------------------------------------------------
#
###################################################################################################


###################################################################################################
# Competing model 3
# -------------------------------------------------------------------------------------------------

# Summary for Phase I result ----------------------------------------------------------------------
Ph1_Mean <- total_Yscale %>% apply(2, mean)
Ph1_invS <- total_Yscale %>% cov %>% solve
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
#
###################################################################################################


###################################################################################################
# Calculate L for ARL0
# -------------------------------------------------------------------------------------------------
# basic setting
ARL0           <- 100
ARL0_num       <- 1000
Ph2_n          <- 3000
Phase1_ARL0_SS <- 300

# create the cluster
myCluster <- makeCluster(30)

# register the cluster with the foreach package
registerDoParallel(myCluster)

ARL0_output <- foreach(ww = 1:ARL0_num, .packages = c("dplyr","dqrng"), .combine=rbind) %dopar% {
  
  # calculate L for ARL0
  PS_T2 <- NULL
  for(w in 1:Phase1_ARL0_SS){
    
    # DATA for obtaining ARL0 = 100
    dqset.seed(w + ww*1009)
    tmp_arl0_idx    <- dqsample(1:ncol(DATA_mat_org0), Ph2_n, replace=T) %>% sort
    tmp_arl0_DATA0  <- DATA_mat_org0[,tmp_arl0_idx]
    
    set.seed(w + ww*1309)
    tmp_arl0_DATA    <- data_transform_cont(tmp_arl0_DATA0)
    
    # Zernike polynomial coefficients
    tmp_arl0_coef <- get_Zernike(tmp_arl0_DATA, Bmat)
    
    # scaling
    tmp_arl0_Yscale0 <- (tmp_arl0_coef - matrix(Step1_scale_center, nrow=Ph2_n, ncol=K, byrow=T)) / matrix(Step1_scale_scale,  nrow=Ph2_n, ncol=K, byrow=T)
    
    # plotting statistics
    tmp_mean <- apply(tmp_arl0_Yscale0,2,mean)
    PS_T2    <- c(PS_T2, Ph2_n * c(t(tmp_mean - Ph1_Mean) %*% Ph1_invS %*% (tmp_mean - Ph1_Mean)))
  }
  PS_T2
}

# stop the cluster 
stopCluster(myCluster)

# calculate alpha
alpha_list <- seq(0.001, 0.03, 0.0005)
OOC_result <- NULL
for(j in 1:length(alpha_list)){
  
  Hotelling_alpha <- alpha_list[j]
  
  tmp_UCL <- qchisq(1-Hotelling_alpha, K)
  tmp_LCL <- 0
  
  tmp_OOC <- sapply(1:ARL0_num, function(i){ which((ARL0_output[i,] > tmp_UCL) | (ARL0_output[i,] < tmp_LCL)) %>% min(Phase1_ARL0_SS) })
  OOC_result <- rbind(OOC_result, tmp_OOC)
}
RL_result <- OOC_result %>% apply(1, mean)
RL_idx    <- which(RL_result > ARL0) %>% tail(1)
if(length(RL_idx) == 0){ ALPHA <- NULL; cat("Error !! (Competing model 3) \n") }else{ ALPHA <- alpha_list[RL_idx] }
# -------------------------------------------------------------------------------------------------
#
###################################################################################################
