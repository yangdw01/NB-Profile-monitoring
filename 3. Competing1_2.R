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
# Competing model 1 and 2
# -------------------------------------------------------------------------------------------------

# Summary for Phase I result ----------------------------------------------------------------------

## basic setting
ARL0           <- 100
Ph2_n          <- 3000
Phase1_ARL0_SS <- 200

## Phase I result
Ph1_ordering   <- Ph1_result[[2]]
Ph1_RESULT     <- Ph1_result[[1]]
Ph1___Gamma_i  <- Ph1_result[[3]]
Ph1___Gamma_i2 <- Ph1_result[[4]]

Ph1___phi_t    <- Ph1_RESULT[[Ph1_DataN]][[3]]
Ph1___psi_hj   <- Ph1_RESULT[[Ph1_DataN]][[4]]
Ph1___m_hj     <- Ph1_RESULT[[Ph1_DataN]][[5]]
Ph1___a_hj     <- Ph1_RESULT[[Ph1_DataN]][[6]]
Ph1___b_hj     <- Ph1_RESULT[[Ph1_DataN]][[7]]

Ph1_Gamma_table <- Ph1___Gamma_i %>% na.omit %>% table
Ph1_n           <- length(Ph1___Gamma_i)
Ph1_H           <- length(Ph1_Gamma_table)

## priors
a   <- 1
b   <- 1
m   <- 0
psi <- 1
alpha_stars <- c(0.01, 0.05, c(0.1 + 0.2 * (0:19)))
eta_probs   <- rep(1/length(alpha_stars), length(alpha_stars))
TT          <- length(alpha_stars)

## setting for SUGS algorithm
SUGS_iter_num <- 5
bidirec   <- T
back_uNum <- NULL
Onum      <- 1000
ITERmax   <- 5
RePUpd    <- T

## calculate hat pi for Phase I
M <- 100
Ph1_alpha_hat    <- sum(Ph1___phi_t[2,] * alpha_stars)

Ph1_Gamma_parms2 <- sapply(1:M, function(k){ c(sum(Ph1___Gamma_i2 == k, na.rm=T),sum(Ph1___Gamma_i2 >= k, na.rm=T)) })
Ph1_p_h2         <- (Ph1_Gamma_parms2[1,]+1)/(Ph1_Gamma_parms2[2,]+Ph1_alpha_hat+1)
Ph1_pi_h2        <- sapply(1:M, function(k){ Ph1_p_h2[k] * prod(c(1,1-Ph1_p_h2)[1:k]) }) %>% unlist()
Ph1_hat_pi_h2    <- c(table(Ph1___Gamma_i2)/sum(!is.na(Ph1___Gamma_i2)), 
                      rep(0, M-(Ph1___Gamma_i2 %>% na.omit %>% unique %>% length)))

Ph1_Gamma_parms <- sapply(1:M, function(k){ c(sum(Ph1___Gamma_i == k, na.rm=T),sum(Ph1___Gamma_i >= k, na.rm=T)) })
Ph1_p_h         <- (Ph1_Gamma_parms[1,]+1)/(Ph1_Gamma_parms[2,]+Ph1_alpha_hat+1)
Ph1_pi_h        <- sapply(1:M, function(k){ Ph1_p_h[k] * prod(c(1,1-Ph1_p_h)[1:k]) }) %>% unlist()
Ph1_hat_pi_h    <- c(table(Ph1___Gamma_i)/sum(!is.na(Ph1___Gamma_i)), 
                     rep(0, M-(Ph1___Gamma_i %>% na.omit %>% unique %>% length)))

Ph1_phi_t  <- Ph1___phi_t[2,]
Ph1_pi     <- rbind(matrix(table(Ph1___Gamma_i2), nrow=Ph1_H, ncol=TT, byrow=F), alpha_stars) * 
  matrix(1/(length(Ph1___Gamma_i2)+alpha_stars), nrow=Ph1_H+1, ncol=TT, byrow=T)
Ph1_prob   <- (Ph1_pi * matrix(Ph1_phi_t, nrow=Ph1_H+1, ncol=TT, byrow=T)) %>% apply(1,sum)

## log aPML
Ph1_aPMLs  <- NULL
for(i in 1:Ph1_DataN){ Ph1_aPMLs <- c(Ph1_aPMLs, Ph1_RESULT[[i]][[9]] %>% apply(2,sum) %>% log) }
Ph1_xbar   <- mean(Ph1_aPMLs)
Ph1_s      <- sd(Ph1_aPMLs)

## for DP-based SPC (Marginalized)
zero_psi_hj  <- rep(psi, K)
zero_m_hj    <- rep(m, K)
zero_a_hj    <- rep(a, K)
zero_b_hj    <- rep(b, K)
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
#
###################################################################################################



###################################################################################################
# Calculate L for ARL0
# -------------------------------------------------------------------------------------------------
# basic setting
ARL0           <- 100
ARL0_SimNum    <- 1000
Ph2_n          <- 3000
Phase1_ARL0_SS <- 300

# create the cluster
myCluster <- makeCluster(30)

# register the cluster with the foreach package
registerDoParallel(myCluster)

ARL0_output <- foreach(ww = 1:ARL0_SimNum, .packages = c("dplyr", "dqrng"), .combine=rbind) %dopar% {
  
  # calculate L for ARL0
  PS_xbar <- NULL
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
    
    # likelihood for DP-based SPC 1 and 2
    tmp_Lih <- NULL 
    for(x in 1:Ph1_H){
      
      tmp_Lih <- rbind(tmp_Lih,     sapply(1:Ph2_n, function(i){ Calculate_Noncentral_t(y    =tmp_arl0_Yscale0[i,], 
                                                                                        ai1  =Ph1___a_hj[[x]][2,], 
                                                                                        bi1  =Ph1___b_hj[[x]][2,], 
                                                                                        mi1  =Ph1___m_hj[[x]][2,], 
                                                                                        psii1=Ph1___psi_hj[[x]][2,]) %>% prod }) )
    }
    tmp_Lih <- rbind(tmp_Lih, sapply(1:Ph2_n, function(i){   Calculate_Noncentral_t(y    =tmp_arl0_Yscale0[i,], 
                                                                                    ai1  =zero_a_hj, 
                                                                                    bi1  =zero_b_hj, 
                                                                                    mi1  =zero_m_hj, 
                                                                                    psii1=zero_psi_hj) %>% prod }))
    tmp_prob_table <- tmp_Lih * matrix( Ph1_prob, nrow=Ph1_H+1, ncol=Ph2_n )
    tmp_log_aPML   <- tmp_prob_table %>% apply(2,sum) %>% log %>% mean
    
    PS_xbar <- c(PS_xbar, tmp_log_aPML)
  }
  PS_xbar
}

# stop the cluster 
stopCluster(myCluster)


# calculate h for DP-based SPC 1
## plotting statistics
cusum_k <- 1/2
PS_Cip  <- PS_Cim <- NULL
for(i in 1:ARL0_SimNum){
  
  PS_xbar <- ARL0_output[i,]
  PS_CM2y    <- (PS_xbar - Ph1_xbar)/(Ph1_s/sqrt(Ph2_n))
  PS_CM2_Cip <- PS_CM2_Cim <- NULL
  for(w in 1:Phase1_ARL0_SS){
    
    if(w > 1){
      prev_cip <- PS_CM2_Cip %>% tail(1)
      prev_cim <- PS_CM2_Cim %>% tail(1)
    }else{
      prev_cip <- prev_cim <- 0
    }
    
    PS_CM2_Cip <- c(PS_CM2_Cip, max( PS_CM2y[w] - cusum_k + prev_cip, 0))
    PS_CM2_Cim <- c(PS_CM2_Cim, max(-PS_CM2y[w] - cusum_k + prev_cim, 0))
  }
  
  PS_Cip <- rbind(PS_Cip, PS_CM2_Cip)
  PS_Cim <- rbind(PS_Cim, PS_CM2_Cim)
}

## calculate h
h_list1     <- seq(1, 8, 0.01)
OOC_result1 <- NULL
for(j in 1:length(h_list1)){
  
  tmp_h   <- h_list1[j]
  tmp_UCL <- tmp_h 
  
  tmp_OOC_p <- sapply(1:ARL0_SimNum, function(i){ which(PS_Cip[i,] > tmp_UCL) %>% min(Phase1_ARL0_SS) })
  tmp_OOC_m <- sapply(1:ARL0_SimNum, function(i){ which(PS_Cim[i,] > tmp_UCL) %>% min(Phase1_ARL0_SS) })
  tmp_OOC   <- pmin(tmp_OOC_p, tmp_OOC_m)
  
  OOC_result1 <- rbind(OOC_result1, tmp_OOC)
}

RL_result1 <- OOC_result1 %>% apply(1, mean)
RL_idx1    <- which(RL_result1 > ARL0) %>% head(1)
if(length(RL_idx1) == 0){ CUSUM_h <- NULL; cat("Error !! (Competing model 1) \n") }else{ CUSUM_h <- h_list1[RL_idx1] }


# calculate L for DP-based SPC 2
L_list2     <- seq(1, 3, 0.01)
OOC_result2 <- NULL
for(j in 1:length(L_list2)){
  
  tmp_L   <- L_list2[j]
  
  tmp_UCL <- Ph1_xbar + tmp_L * Ph1_s/sqrt(Ph2_n)
  tmp_LCL <- Ph1_xbar - tmp_L * Ph1_s/sqrt(Ph2_n)
  
  tmp_OOC <- sapply(1:ARL0_SimNum, function(i){ which((ARL0_output[i,] > tmp_UCL) | (ARL0_output[i,] < tmp_LCL)) %>% min(Phase1_ARL0_SS) })
  
  OOC_result2 <- rbind(OOC_result2, tmp_OOC)
}
RL_result2 <- OOC_result2 %>% apply(1, mean)
RL_idx2    <- which(RL_result2 > ARL0) %>% head(1)
if(length(RL_idx2) == 0){ Xbar_L <- NULL; cat("Error !! \n") }else{ Xbar_L <- L_list2[RL_idx2] } 
# -------------------------------------------------------------------------------------------------
#
###################################################################################################