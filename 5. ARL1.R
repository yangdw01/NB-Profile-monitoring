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
library(ggplot2)
library(patchwork)
library(gridExtra)  
library(viridis)    
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

## # of WF
WF_chipN <- prod(dat_dim) - length(WF_outise_idx)
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
Final_Total_Result          <- Ph1_result[[1]]
Final_Total_Result_ordering <- Ph1_result[[2]]
Gamma_i  <- Ph1_result[[3]]
Gamma_i2 <- Ph1_result[[4]]

## check cluster proportion
Gamma_i2_list <- Gamma_i2_list0 <- NULL
for(i in 1:Ph1_DataN){
  
  tmp_i_idx      <- 1:each_n + (i-1) * each_n
  tmp_i_Gamma_i2 <- Gamma_i2[tmp_i_idx]
  Gamma_i2_list  <- rbind(Gamma_i2_list, tmp_i_Gamma_i2)
  
  tmp_i_Gamma_i2 <- tmp_i_Gamma_i2[Final_Total_Result_ordering[[i]] %>% order]
  Gamma_i2_list0 <- rbind(Gamma_i2_list0, tmp_i_Gamma_i2)
}
Ph1_cluster <- Gamma_i2_list0 %>% apply(2, function(x){ x %>% table %>% which.max %>% names %>% as.numeric })
# -------------------------------------------------------------------------------------------------
#
###################################################################################################


###################################################################################################
# ARL1 - setting
# -------------------------------------------------------------------------------------------------
## basic setting
ARL0        <- 100
ARL0_SimNum <- 1000
Ph2_n       <- 3000
Phase2_SS   <- 40

## setting for SUGS algorithm
a   <- 1
b   <- 1
m   <- 0
psi <- 1
alpha_stars <- c(0.01, 0.05, c(0.1 + 0.2 * (0:19)))
eta_probs   <- rep(1/length(alpha_stars), length(alpha_stars))
TT          <- length(alpha_stars)

zero_psi_hj  <- rep(psi, K)
zero_m_hj    <- rep(m, K)
zero_a_hj    <- rep(a, K)
zero_b_hj    <- rep(b, K)

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

Ph1_phi_t  <- Ph1___phi_t[2,]
Ph1_pi     <- rbind(matrix(table(Ph1___Gamma_i2), nrow=Ph1_H, ncol=TT, byrow=F), alpha_stars) * 
  matrix(1/(length(Ph1___Gamma_i2)+alpha_stars), nrow=Ph1_H+1, ncol=TT, byrow=T)
Ph1_prob   <- (Ph1_pi * matrix(Ph1_phi_t, nrow=Ph1_H+1, ncol=TT, byrow=T)) %>% apply(1,sum)


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

## log aPML
Ph1_aPMLs  <- NULL
for(i in 1:Ph1_DataN){ Ph1_aPMLs <- c(Ph1_aPMLs, Ph1_RESULT[[i]][[9]] %>% apply(2,sum) %>% log) }
Ph1_xbar   <- mean(Ph1_aPMLs)
Ph1_s      <- sd(Ph1_aPMLs)

## Summary for Phase I data
Ph1_Mean <- total_Yscale %>% apply(2, mean)
Ph1_invS <- total_Yscale %>% cov %>% solve

## setting for the proposed model
SUGS_iter_num <- 5
bidirec   <- T
back_uNum <- NULL
Onum      <- 1000
ITERmax   <- 5
RePUpd    <- T
# -------------------------------------------------------------------------------------------------
#
###################################################################################################



###################################################################################################
# ARL1 - setting
# -------------------------------------------------------------------------------------------------
# basic setting
coreNUM   <- 30
ARL1_simN <- 250
ARL0      <- 100
Ph2_n     <- 3000
Phase2_SS <- 50

# L, alpha, h
M1_alpha0 <- 0.0006110553
M1_alpha  <- 0.00781407
M2_h      <- 3.61
M3_L      <- 2.61
M4_alpha  <- 0.009

# for proposed model
target_Ph1_pis   <- Ph1_pi_h2
lower_prob_thres <- 0.005
# -------------------------------------------------------------------------------------------------
#
###################################################################################################



###################################################################################################
# ARL1 (Scenario 1) - new cluster!
# -------------------------------------------------------------------------------------------------
SC1_ARL1_result <- list()

delta_list <- c(0.002, 0.004, 0.006, 0.008)
for(delta in delta_list){
  
  delta_idx <- which(delta_list %in% delta)
  
  # parallel computing
  myCluster <- makeCluster(coreNUM)
  registerDoParallel(myCluster)
  output <- foreach(www = 1:ARL1_simN, .packages = c("dplyr", "dqrng")) %dopar% {
    
    ## new cluster
    delta_n        <- round(Ph2_n * delta / (1-delta))
    Ph2_n_new      <- Ph2_n + delta_n
    
    # Phase I parameters
    tmp___Gamma_i  <- Ph1___Gamma_i
    tmp___Gamma_i2 <- Ph1___Gamma_i2
    tmp___phi_t    <- Ph1___phi_t
    tmp___psi_hj   <- Ph1___psi_hj
    tmp___m_hj     <- Ph1___m_hj
    tmp___a_hj     <- Ph1___a_hj
    tmp___b_hj     <- Ph1___b_hj
    
    # calculate L for ARL0
    ARL1_result <- list()
    M2_3_xbar   <- M4_T2 <- NULL
    for(w in 1:Phase2_SS){
      
      ### data generation -------------------------------------------------------------------------------------
      ## raw Data
      dqset.seed(w + www*1301)
      tmp_data_add      <- SC1_RawData(delta_n, X_data_polar, WF_chipN)
      tmp_data_mat_all0 <- cbind(DATA_mat_org0, tmp_data_add)
      tmp_arl1_idx      <- dqsample(1:ncol(tmp_data_mat_all0), Ph2_n_new, replace=T) %>% sort
      tmp_arl1_DATA0    <- tmp_data_mat_all0[,tmp_arl1_idx]
      tmp_new_cl_idx0   <- tmp_arl1_idx[tmp_arl1_idx > Ph2_n]
      
      ## continuous Data
      set.seed(w + www*1301)
      tmp_arl1_DATA     <- SC12_data_transform_cont(tmp_arl1_DATA0)
      
      ## Zernike polynomial coefficients
      tmp_arl1_coef <- get_Zernike(tmp_arl1_DATA, Bmat)
      
      ## scaling
      tmp_arl1_Yscale0 <- (tmp_arl1_coef - matrix(Step1_scale_center, nrow=Ph2_n_new, ncol=K, byrow=T)) / 
        matrix(Step1_scale_scale,  nrow=Ph2_n_new, ncol=K, byrow=T)
      ### ------------------------------------------------------------------------------------------------------
      
      ### proposed model ---------------------------------------------------------------------------------------
      ## BSUGS algorithm
      tmp_RESULT <- tmp_RESULT2 <- list()
      for(ww in 1:SUGS_iter_num){
        
        # ordering
        dqset.seed(ww + w*1000 + www*1301)
        tmp_resample    <- dqsample(1:nrow(tmp_arl1_Yscale0))
        tmp_arl1_Yscale <- tmp_arl1_Yscale0[tmp_resample,]
        tmp_new_cl_idx  <- which(tmp_resample %in% tmp_new_cl_idx0)
        
        # SUGS algorithm
        set.seed(ww + w*1000 + www*1301)
        tmp_result   <- BSUGS(Y=tmp_arl1_Yscale, a=a, b=b, m=m, psi=psi, alpha_stars=alpha_stars, eta_probs=eta_probs,
                              bidirectional = bidirec, back_upperNum = back_uNum,
                              OnceNum       = Onum,    IterMax       = ITERmax,   ReParmUpd   = RePUpd,
                              prior_info=FALSE, prior_Gamma_i=NULL,   Phase1=F,
                              prev_Gamma_i  = tmp___Gamma_i, prev_phi_t=tmp___phi_t, prev_psi_hj=tmp___psi_hj,
                              prev_m_hj     =tmp___m_hj,     prev_a_hj =tmp___a_hj,  prev_b_hj  =tmp___b_hj)
        tmp_RESULT[[ww]]  <- tmp_result
        tmp_RESULT2[[ww]] <- tmp_resample
      }
      tmp_aPMLs    <- tmp_RESULT %>% sapply(function(x){ x[[1]] })
      tmp_max_idx  <- tmp_aPMLs %>% which.max
      tmp_result   <- tmp_RESULT[[tmp_max_idx]]
      tmp_ordering <- tmp_RESULT2[[tmp_max_idx]]
      tmp_aPML     <- tmp_aPMLs %>% max
      
      ## update clustering result
      tmp___Gamma_i      <- c(tmp___Gamma_i,  tmp_result[[2]])
      tmp___Gamma_i2     <- c(tmp___Gamma_i2, tmp_result[[12]])
      
      ## update T2 value
      tmp_nu_star   <- 2*(tmp___a_hj %>% sapply(function(x){ x[2,1] }))
      tmp_t_fun <- function(k){
        
        tmp_nk    <- sum(tmp_result[[2]] == k, na.rm=T)
        
        if(tmp_nk > 0){
          
          tmp_denom <- sqrt(Ph1___b_hj[[k]][2,]/Ph1___a_hj[[k]][2,] * (1/tmp_nk + tmp___psi_hj[[k]][2,]))
          tmp_ymean <- tmp_arl1_Yscale0[tmp_ordering[tmp_result[[2]] == k],,drop=FALSE] %>% apply(2, function(x){ mean(x, na.rm=T) })
          tmp_nom   <- (tmp_ymean - Ph1___m_hj[[k]][2,])
          
          tmp_tvals    <- tmp_nom/tmp_denom
          
        }else{ tmp_tvals <- NA }
        
        return(tmp_tvals)
      }
      tmp_tval <- sapply(1:Ph1_H, function(k){ tmp_t_fun(k) }, simplify=F)
      
      ## update parameters
      tmp___phi_t        <- tmp_result[[3]]
      tmp___psi_hj       <- tmp_result[[4]]
      tmp___m_hj         <- tmp_result[[5]]
      tmp___a_hj         <- tmp_result[[6]]
      tmp___b_hj         <- tmp_result[[7]]
      
      ARL1_result[[w]] <- list(tmp_result[[2]], tmp_result[[12]], 
                               Ph2_n_new - sum(is.na(tmp_result[[2]])),
                               tmp_aPML, tmp___phi_t, tmp___psi_hj, tmp___m_hj, tmp___a_hj, tmp___b_hj,
                               tmp_tval, tmp_nu_star)
      ### ------------------------------------------------------------------------------------------------------
      
      ### Competing model 1 and 2 ------------------------------------------------------------------------------
      tmp_Lih <- NULL 
      for(x in 1:Ph1_H){
        
        tmp_Lih <- rbind(tmp_Lih,     sapply(1:Ph2_n_new, function(i){ Calculate_Noncentral_t(y    =tmp_arl1_Yscale0[i,], 
                                                                                              ai1  =Ph1___a_hj[[x]][2,], 
                                                                                              bi1  =Ph1___b_hj[[x]][2,], 
                                                                                              mi1  =Ph1___m_hj[[x]][2,], 
                                                                                              psii1=Ph1___psi_hj[[x]][2,]) %>% prod }) )
      } 
      tmp_Lih <- rbind(tmp_Lih, sapply(1:Ph2_n_new, function(i){   Calculate_Noncentral_t(y    =tmp_arl1_Yscale0[i,], 
                                                                                          ai1  =zero_a_hj, 
                                                                                          bi1  =zero_b_hj, 
                                                                                          mi1  =zero_m_hj, 
                                                                                          psii1=zero_psi_hj) %>% prod }))
      tmp_prob_table <- tmp_Lih * matrix( Ph1_prob, nrow=Ph1_H+1, ncol=Ph2_n_new )
      tmp_log_aPML   <- tmp_prob_table %>% apply(2,sum) %>% log %>% mean
      M2_3_xbar      <- c(M2_3_xbar, tmp_log_aPML)
      ### ------------------------------------------------------------------------------------------------------
      
      ### Competing 3 ------------------------------------------------------------------------------------------
      tmp_mean <- apply(tmp_arl1_Yscale0,2,mean)
      M4_T2    <- c(M4_T2, Ph2_n_new * c(t(tmp_mean - Ph1_Mean) %*% Ph1_invS %*% (tmp_mean - Ph1_Mean)))
      ### ------------------------------------------------------------------------------------------------------
    }
    
    # proposed model
    ## setting
    new_table <- ARL1_result %>% sapply(function(x){ x[[2]] }) %>% c %>% table
    Int2_idx <- which(target_Ph1_pis >= lower_prob_thres)
    Int3_idx <- which(target_Ph1_pis <  lower_prob_thres) %>% intersect(seq_along(new_table))
    
    ## interval types
    alpha0 <- 0.0006110553
    Interval_type2 <- Interval3(Ph1_pi_h, Ph2_n_new, alpha0)
    UCL_type2      <- Interval_type2[[1]]
    LCL_type2      <- Interval_type2[[2]]
    
    alpha <- 0.00781407
    Interval_type3 <- Interval3(Ph1_pi_h, Ph2_n_new, alpha)
    UCL_type3      <- Interval_type3[[1]]
    LCL_type3      <- Interval_type3[[2]]
    
    ## plot 1
    chk_tf <- NULL
    for(k in seq_along(new_table)){
      
      target_ps <- ARL1_result %>% sapply(function(x){ sum(x[[2]] == k, na.rm=T)/sum(!is.na(x[[2]])) })
      target_CL <- target_Ph1_pis[k] 
      
      int2_yn   <- k %in% Int2_idx
      CtrlL_col <- ifelse(int2_yn, 2, 3)
      
      if(int2_yn){
        
        target_UCL <- UCL_type2[k]
        target_LCL <- LCL_type2[k]
        
        tmp_tf <- (target_ps > target_UCL) | (target_ps < target_LCL)
        
      }else{
        
        target_UCL <- UCL_type3[k]
        target_LCL <- NA
        
        tmp_tf <- (target_ps > target_UCL)
      }
      chk_tf <- rbind(chk_tf, tmp_tf)
    }
    RL <- chk_tf %>% apply(2, function(x){ any(x, na.rm=T) }) %>% which
    RL <- ifelse(length(RL) >= 1, RL[1], Phase2_SS)
    
    # competing model 1 and 2 
    ## competing model 1
    cusum_k    <- 1/2
    PS_CM2y    <- (M2_3_xbar - Ph1_xbar)/(Ph1_s/sqrt(Ph2_n))
    PS_CM2_Cip <- PS_CM2_Cim <- NULL
    for(w in 1:Phase2_SS){
      
      if(w > 1){
        prev_cip <- PS_CM2_Cip %>% tail(1)
        prev_cim <- PS_CM2_Cim %>% tail(1)
      }else{
        prev_cip <- prev_cim <- 0
      }
      
      PS_CM2_Cip <- c(PS_CM2_Cip, max( PS_CM2y[w] - cusum_k + prev_cip, 0))
      PS_CM2_Cim <- c(PS_CM2_Cim, max(-PS_CM2y[w] - cusum_k + prev_cim, 0))
    }
    ARL1_result_cm1 <- rbind(PS_CM2_Cip, PS_CM2_Cim)
    RL_cm1 <- which((PS_CM2_Cip >= M2_h) | (PS_CM2_Cim >= M2_h))
    RL_cm1 <- ifelse(length(RL_cm1) >= 1, RL_cm1[1], Phase2_SS)
    
    ## competing model 2
    ARL1_result_cm2 <- M2_3_xbar
    tmp_ucl <- Ph1_xbar + M3_L * Ph1_s/sqrt(Ph2_n_new)
    tmp_lcl <- Ph1_xbar - M3_L * Ph1_s/sqrt(Ph2_n_new)
    RL_cm2 <- which((ARL1_result_cm2 > tmp_ucl) | (ARL1_result_cm2 < tmp_lcl))
    RL_cm2 <- ifelse(length(RL_cm2) >= 1, RL_cm2[1], Phase2_SS)
    
    ## competing model 3
    ARL1_result_cm3 <- M4_T2
    RL_cm3 <- which(ARL1_result_cm3 > qchisq(1-M4_alpha, K))
    RL_cm3 <- ifelse(length(RL_cm3) >= 1, RL_cm3[1], Phase2_SS)
    
    ## save result 
    Raw_PSs <- list(ARL1_result, ARL1_result_cm1, ARL1_result_cm2, ARL1_result_cm3)
    RLs     <- list(RL, RL_cm1, RL_cm2, RL_cm3)
    
    tmp_RL_result <- list(Raw_PSs, RLs)
    tmp_RL_result
  }
  stopCluster(myCluster)
  
  # save
  SC1_ARL1_result[[delta_idx]] <- output
}
# -------------------------------------------------------------------------------------------------
#
####################################################################################################
timer1  <- proc.time()[3]
Contime <- (timer1 - timer0)/60
cat(Contime, "min! \n")
# save(SC1_ARL1_result, file="SC1_ARL1.RData")


timer2  <- proc.time()[3]
###################################################################################################
# ARL1 (Scenario 2) - two new cluster!
# -------------------------------------------------------------------------------------------------
SC2_ARL1_result <- list()

delta_list <- c(0.0015, 0.003, 0.0045, 0.006)
for(delta in delta_list){
  
  delta_idx <- which(delta_list %in% delta)
  
  # parallel computing
  myCluster <- makeCluster(coreNUM)
  registerDoParallel(myCluster)
  output <- foreach(www = 1:ARL1_simN, .packages = c("dplyr", "dqrng")) %dopar% {
    
    ## new cluster
    delta_n        <- round(Ph2_n * delta / (1-delta))
    Ph2_n_new      <- Ph2_n + 2 * delta_n
    
    # Phase I parameters
    tmp___Gamma_i  <- Ph1___Gamma_i
    tmp___Gamma_i2 <- Ph1___Gamma_i2
    tmp___phi_t    <- Ph1___phi_t
    tmp___psi_hj   <- Ph1___psi_hj
    tmp___m_hj     <- Ph1___m_hj
    tmp___a_hj     <- Ph1___a_hj
    tmp___b_hj     <- Ph1___b_hj
    
    # calculate L for ARL0
    ARL1_result <- list()
    M2_3_xbar   <- M4_T2 <- NULL
    for(w in 1:Phase2_SS){
      
      ### data generation -------------------------------------------------------------------------------------
      ## raw Data
      dqset.seed(w + www*1305)
      tmp_data_add1     <- SC1_RawData(delta_n, X_data_polar, WF_chipN)
      tmp_data_add2     <- SC2_RawData(delta_n, X_data_polar, WF_chipN)
      tmp_data_mat_all0 <- cbind(DATA_mat_org0, tmp_data_add1, tmp_data_add2)
      
      tmp_arl1_idx      <- dqsample(1:ncol(tmp_data_mat_all0), Ph2_n_new, replace=T) %>% sort
      tmp_arl1_DATA0    <- tmp_data_mat_all0[,tmp_arl1_idx]
      tmp_new_cl_idx1   <- tmp_arl1_idx[(tmp_arl1_idx > Ph2_n) & (tmp_arl1_idx <= Ph2_n + delta_n)]
      tmp_new_cl_idx2   <- tmp_arl1_idx[tmp_arl1_idx > Ph2_n + delta_n]
      
      ## continuous Data
      set.seed(w + www*1305)
      tmp_arl1_DATA     <- SC12_data_transform_cont(tmp_arl1_DATA0)
      
      ## Zernike polynomial coefficients
      tmp_arl1_coef <- get_Zernike(tmp_arl1_DATA, Bmat)
      
      ## scaling
      tmp_arl1_Yscale0 <- (tmp_arl1_coef - matrix(Step1_scale_center, nrow=Ph2_n_new, ncol=K, byrow=T)) / 
        matrix(Step1_scale_scale,  nrow=Ph2_n_new, ncol=K, byrow=T)
      ### ------------------------------------------------------------------------------------------------------
      
      ### proposed model ---------------------------------------------------------------------------------------
      ## BSUGS algorithm
      tmp_RESULT <- tmp_RESULT2 <- list()
      for(ww in 1:SUGS_iter_num){
        
        # ordering
        dqset.seed(ww + w*1000 + www*1305)
        tmp_resample    <- dqsample(1:nrow(tmp_arl1_Yscale0))
        tmp_arl1_Yscale <- tmp_arl1_Yscale0[tmp_resample,]
        
        # SUGS algorithm
        set.seed(ww + w*1000 + www*1305)
        tmp_result   <- BSUGS(Y=tmp_arl1_Yscale, a=a, b=b, m=m, psi=psi, alpha_stars=alpha_stars, eta_probs=eta_probs,
                              bidirectional = bidirec, back_upperNum = back_uNum,
                              OnceNum       = Onum,    IterMax       = ITERmax,   ReParmUpd   = RePUpd,
                              prior_info=FALSE, prior_Gamma_i=NULL,   Phase1=F,
                              prev_Gamma_i  = tmp___Gamma_i, prev_phi_t=tmp___phi_t, prev_psi_hj=tmp___psi_hj,
                              prev_m_hj     =tmp___m_hj,     prev_a_hj =tmp___a_hj,  prev_b_hj  =tmp___b_hj)
        tmp_RESULT[[ww]]  <- tmp_result
        tmp_RESULT2[[ww]] <- tmp_resample
      }
      tmp_aPMLs    <- tmp_RESULT %>% sapply(function(x){ x[[1]] })
      tmp_max_idx  <- tmp_aPMLs %>% which.max
      tmp_result   <- tmp_RESULT[[tmp_max_idx]]
      tmp_ordering <- tmp_RESULT2[[tmp_max_idx]]
      tmp_aPML     <- tmp_aPMLs %>% max
      
      ## update clustering result
      tmp___Gamma_i      <- c(tmp___Gamma_i,  tmp_result[[2]])
      tmp___Gamma_i2     <- c(tmp___Gamma_i2, tmp_result[[12]])
      
      ## update T2 value
      tmp_nu_star   <- 2*(tmp___a_hj %>% sapply(function(x){ x[2,1] }))
      tmp_t_fun <- function(k){
        
        tmp_nk    <- sum(tmp_result[[2]] == k, na.rm=T)
        
        if(tmp_nk > 0){
          
          tmp_denom <- sqrt(Ph1___b_hj[[k]][2,]/Ph1___a_hj[[k]][2,] * (1/tmp_nk + tmp___psi_hj[[k]][2,]))
          tmp_ymean <- tmp_arl1_Yscale0[tmp_ordering[tmp_result[[2]] == k],,drop=FALSE] %>% apply(2, function(x){ mean(x, na.rm=T) })
          tmp_nom   <- (tmp_ymean - Ph1___m_hj[[k]][2,])
          
          tmp_tvals    <- tmp_nom/tmp_denom
          
        }else{ tmp_tvals <- NA }
        
        return(tmp_tvals)
      }
      tmp_tval <- sapply(1:Ph1_H, function(k){ tmp_t_fun(k) }, simplify=F)
      
      ## update parameters
      tmp___phi_t        <- tmp_result[[3]]
      tmp___psi_hj       <- tmp_result[[4]]
      tmp___m_hj         <- tmp_result[[5]]
      tmp___a_hj         <- tmp_result[[6]]
      tmp___b_hj         <- tmp_result[[7]]
      
      ARL1_result[[w]] <- list(tmp_result[[2]], tmp_result[[12]], 
                               Ph2_n_new - sum(is.na(tmp_result[[2]])),
                               tmp_aPML, tmp___phi_t, tmp___psi_hj, tmp___m_hj, tmp___a_hj, tmp___b_hj,
                               tmp_tval, tmp_nu_star)
      ### ------------------------------------------------------------------------------------------------------
      
      ### Competing model 1 and 2 ------------------------------------------------------------------------------
      tmp_Lih <- NULL 
      for(x in 1:Ph1_H){
        
        tmp_Lih <- rbind(tmp_Lih,     sapply(1:Ph2_n_new, function(i){ Calculate_Noncentral_t(y    =tmp_arl1_Yscale0[i,], 
                                                                                              ai1  =Ph1___a_hj[[x]][2,], 
                                                                                              bi1  =Ph1___b_hj[[x]][2,], 
                                                                                              mi1  =Ph1___m_hj[[x]][2,], 
                                                                                              psii1=Ph1___psi_hj[[x]][2,]) %>% prod }) )
      } 
      tmp_Lih <- rbind(tmp_Lih, sapply(1:Ph2_n_new, function(i){   Calculate_Noncentral_t(y    =tmp_arl1_Yscale0[i,], 
                                                                                          ai1  =zero_a_hj, 
                                                                                          bi1  =zero_b_hj, 
                                                                                          mi1  =zero_m_hj, 
                                                                                          psii1=zero_psi_hj) %>% prod }))
      tmp_prob_table <- tmp_Lih * matrix( Ph1_prob, nrow=Ph1_H+1, ncol=Ph2_n_new )
      tmp_log_aPML   <- tmp_prob_table %>% apply(2,sum) %>% log %>% mean
      M2_3_xbar      <- c(M2_3_xbar, tmp_log_aPML)
      ### ------------------------------------------------------------------------------------------------------
      
      ### Competing 3 ------------------------------------------------------------------------------------------
      tmp_mean <- apply(tmp_arl1_Yscale0,2,mean)
      M4_T2    <- c(M4_T2, Ph2_n_new * c(t(tmp_mean - Ph1_Mean) %*% Ph1_invS %*% (tmp_mean - Ph1_Mean)))
      ### ------------------------------------------------------------------------------------------------------
    }
    
    # proposed model
    ## setting
    new_table <- ARL1_result %>% sapply(function(x){ x[[2]] }) %>% c %>% table
    Int2_idx <- which(target_Ph1_pis >= lower_prob_thres)
    Int3_idx <- which(target_Ph1_pis <  lower_prob_thres) %>% intersect(seq_along(new_table))
    
    ## interval types
    alpha0 <- 0.0006110553
    Interval_type2 <- Interval3(Ph1_pi_h, Ph2_n_new, alpha0)
    UCL_type2      <- Interval_type2[[1]]
    LCL_type2      <- Interval_type2[[2]]
    
    alpha <- 0.00781407
    Interval_type3 <- Interval3(Ph1_pi_h, Ph2_n_new, alpha)
    UCL_type3      <- Interval_type3[[1]]
    LCL_type3      <- Interval_type3[[2]]
    
    ## plot 1
    chk_tf <- NULL
    for(k in seq_along(new_table)){
      
      target_ps <- ARL1_result %>% sapply(function(x){ sum(x[[2]] == k, na.rm=T)/sum(!is.na(x[[2]])) })
      target_CL <- target_Ph1_pis[k] 
      
      int2_yn   <- k %in% Int2_idx
      CtrlL_col <- ifelse(int2_yn, 2, 3)
      
      if(int2_yn){
        
        target_UCL <- UCL_type2[k]
        target_LCL <- LCL_type2[k]
        
        tmp_tf <- (target_ps > target_UCL) | (target_ps < target_LCL)
        
      }else{
        
        target_UCL <- UCL_type3[k]
        target_LCL <- NA
        
        tmp_tf <- (target_ps > target_UCL)
      }
      chk_tf <- rbind(chk_tf, tmp_tf)
    }
    RL <- chk_tf %>% apply(2, function(x){ any(x, na.rm=T) }) %>% which
    RL <- ifelse(length(RL) >= 1, RL[1], Phase2_SS)
    
    # competing model 1 and 2 
    ## competing model 1
    cusum_k    <- 1/2
    PS_CM2y    <- (M2_3_xbar - Ph1_xbar)/(Ph1_s/sqrt(Ph2_n))
    PS_CM2_Cip <- PS_CM2_Cim <- NULL
    for(w in 1:Phase2_SS){
      
      if(w > 1){
        prev_cip <- PS_CM2_Cip %>% tail(1)
        prev_cim <- PS_CM2_Cim %>% tail(1)
      }else{
        prev_cip <- prev_cim <- 0
      }
      
      PS_CM2_Cip <- c(PS_CM2_Cip, max( PS_CM2y[w] - cusum_k + prev_cip, 0))
      PS_CM2_Cim <- c(PS_CM2_Cim, max(-PS_CM2y[w] - cusum_k + prev_cim, 0))
    }
    ARL1_result_cm1 <- rbind(PS_CM2_Cip, PS_CM2_Cim)
    RL_cm1 <- which((PS_CM2_Cip >= M2_h) | (PS_CM2_Cim >= M2_h))
    RL_cm1 <- ifelse(length(RL_cm1) >= 1, RL_cm1[1], Phase2_SS)
    
    ## competing model 2
    ARL1_result_cm2 <- M2_3_xbar
    tmp_ucl <- Ph1_xbar + M3_L * Ph1_s/sqrt(Ph2_n_new)
    tmp_lcl <- Ph1_xbar - M3_L * Ph1_s/sqrt(Ph2_n_new)
    RL_cm2 <- which((ARL1_result_cm2 > tmp_ucl) | (ARL1_result_cm2 < tmp_lcl))
    RL_cm2 <- ifelse(length(RL_cm2) >= 1, RL_cm2[1], Phase2_SS)
    
    ## competing model 3
    ARL1_result_cm3 <- M4_T2
    RL_cm3 <- which(ARL1_result_cm3 > qchisq(1-M4_alpha, K))
    RL_cm3 <- ifelse(length(RL_cm3) >= 1, RL_cm3[1], Phase2_SS)
    
    ## save result 
    Raw_PSs <- list(ARL1_result, ARL1_result_cm1, ARL1_result_cm2, ARL1_result_cm3)
    RLs     <- list(RL, RL_cm1, RL_cm2, RL_cm3)
    
    tmp_RL_result <- list(Raw_PSs, RLs)
    tmp_RL_result
  }
  stopCluster(myCluster)
  
  # save
  SC2_ARL1_result[[delta_idx]] <- output
}
# -------------------------------------------------------------------------------------------------
#
####################################################################################################
timer3  <- proc.time()[3]
Contime <- (timer3 - timer2)/60
cat(Contime, "min! \n")
# save(SC2_ARL1_result, file="SC2_ARL1.RData")


timer4  <- proc.time()[3]
###################################################################################################
# ARL1 (Scenario 3) - existing cluster!
# -------------------------------------------------------------------------------------------------
SC3_ARL1_result <- list()

delta_list <- c(0.005, 0.01, 0.015, 0.02)
for(delta in delta_list){
  
  delta_idx <- which(delta_list %in% delta)
  
  # parallel computing
  myCluster <- makeCluster(coreNUM)
  registerDoParallel(myCluster)
  output <- foreach(www = 1:ARL1_simN, .packages = c("dplyr", "dqrng")) %dopar% {
    
    ## new cluster
    delta_n        <- round(each_n * delta)
    Ph2_n_new      <- Ph2_n
    target_cl1     <- 2
    target_cl2     <- 20
    
    # Phase I parameters
    tmp___Gamma_i  <- Ph1___Gamma_i
    tmp___Gamma_i2 <- Ph1___Gamma_i2
    tmp___phi_t    <- Ph1___phi_t
    tmp___psi_hj   <- Ph1___psi_hj
    tmp___m_hj     <- Ph1___m_hj
    tmp___a_hj     <- Ph1___a_hj
    tmp___b_hj     <- Ph1___b_hj
    
    # calculate L for ARL0
    ARL1_result <- list()
    M2_3_xbar   <- M4_T2 <- NULL
    for(w in 1:Phase2_SS){
      
      ### data generation -------------------------------------------------------------------------------------
      ## raw Data
      dqset.seed(w + www*1311)
      
      tmp_idx1          <- which(Ph1_cluster != target_cl1)
      tmp_idx2          <- which(Ph1_cluster == target_cl2) %>% dqsample(delta_n) %>% sort 
      tmp_idx3          <- which(Ph1_cluster == target_cl1) %>% dqsample(sum(Ph1_cluster == target_cl1) - delta_n) %>% sort
      tmp_data_mat_all0 <- DATA_mat_org0[,c(tmp_idx1, tmp_idx2, tmp_idx3)]
      
      tmp_arl1_idx      <- dqsample(1:ncol(tmp_data_mat_all0), Ph2_n_new, replace=T) %>% sort
      tmp_arl1_DATA0    <- tmp_data_mat_all0[,tmp_arl1_idx]
      
      ## continuous Data
      set.seed(w + www*1311)
      tmp_arl1_DATA     <- SC12_data_transform_cont(tmp_arl1_DATA0)
      
      ## Zernike polynomial coefficients
      tmp_arl1_coef <- get_Zernike(tmp_arl1_DATA, Bmat)
      
      ## scaling
      tmp_arl1_Yscale0 <- (tmp_arl1_coef - matrix(Step1_scale_center, nrow=Ph2_n_new, ncol=K, byrow=T)) / 
        matrix(Step1_scale_scale,  nrow=Ph2_n_new, ncol=K, byrow=T)
      ### ------------------------------------------------------------------------------------------------------
      
      ### proposed model ---------------------------------------------------------------------------------------
      ## BSUGS algorithm
      tmp_RESULT <- tmp_RESULT2 <- list()
      for(ww in 1:SUGS_iter_num){
        
        # ordering
        dqset.seed(ww + w*1000 + www*1311)
        tmp_resample    <- dqsample(1:nrow(tmp_arl1_Yscale0))
        tmp_arl1_Yscale <- tmp_arl1_Yscale0[tmp_resample,]
        
        # SUGS algorithm
        set.seed(ww + w*1000 + www*1311)
        tmp_result   <- BSUGS(Y=tmp_arl1_Yscale, a=a, b=b, m=m, psi=psi, alpha_stars=alpha_stars, eta_probs=eta_probs,
                              bidirectional = bidirec, back_upperNum = back_uNum,
                              OnceNum       = Onum,    IterMax       = ITERmax,   ReParmUpd   = RePUpd,
                              prior_info=FALSE, prior_Gamma_i=NULL,   Phase1=F,
                              prev_Gamma_i  = tmp___Gamma_i, prev_phi_t=tmp___phi_t, prev_psi_hj=tmp___psi_hj,
                              prev_m_hj     =tmp___m_hj,     prev_a_hj =tmp___a_hj,  prev_b_hj  =tmp___b_hj)
        tmp_RESULT[[ww]]  <- tmp_result
        tmp_RESULT2[[ww]] <- tmp_resample
      }
      tmp_aPMLs    <- tmp_RESULT %>% sapply(function(x){ x[[1]] })
      tmp_max_idx  <- tmp_aPMLs %>% which.max
      tmp_result   <- tmp_RESULT[[tmp_max_idx]]
      tmp_ordering <- tmp_RESULT2[[tmp_max_idx]]
      tmp_aPML     <- tmp_aPMLs %>% max
      
      ## update clustering result
      tmp___Gamma_i      <- c(tmp___Gamma_i,  tmp_result[[2]])
      tmp___Gamma_i2     <- c(tmp___Gamma_i2, tmp_result[[12]])
      
      ## update T2 value
      tmp_nu_star   <- 2*(tmp___a_hj %>% sapply(function(x){ x[2,1] }))
      tmp_t_fun <- function(k){
        
        tmp_nk    <- sum(tmp_result[[2]] == k, na.rm=T)
        
        if(tmp_nk > 0){
          
          tmp_denom <- sqrt(Ph1___b_hj[[k]][2,]/Ph1___a_hj[[k]][2,] * (1/tmp_nk + tmp___psi_hj[[k]][2,]))
          tmp_ymean <- tmp_arl1_Yscale0[tmp_ordering[tmp_result[[2]] == k],,drop=FALSE] %>% apply(2, function(x){ mean(x, na.rm=T) })
          tmp_nom   <- (tmp_ymean - Ph1___m_hj[[k]][2,])
          
          tmp_tvals    <- tmp_nom/tmp_denom
          
        }else{ tmp_tvals <- NA }
        
        return(tmp_tvals)
      }
      tmp_tval <- sapply(1:Ph1_H, function(k){ tmp_t_fun(k) }, simplify=F)
      
      ## update parameters
      tmp___phi_t        <- tmp_result[[3]]
      tmp___psi_hj       <- tmp_result[[4]]
      tmp___m_hj         <- tmp_result[[5]]
      tmp___a_hj         <- tmp_result[[6]]
      tmp___b_hj         <- tmp_result[[7]]
      
      ARL1_result[[w]] <- list(tmp_result[[2]], tmp_result[[12]], 
                               Ph2_n_new - sum(is.na(tmp_result[[2]])),
                               tmp_aPML, tmp___phi_t, tmp___psi_hj, tmp___m_hj, tmp___a_hj, tmp___b_hj,
                               tmp_tval, tmp_nu_star)
      ### ------------------------------------------------------------------------------------------------------
      
      ### Competing model 1 and 2 ------------------------------------------------------------------------------
      tmp_Lih <- NULL 
      for(x in 1:Ph1_H){
        
        tmp_Lih <- rbind(tmp_Lih,     sapply(1:Ph2_n_new, function(i){ Calculate_Noncentral_t(y    =tmp_arl1_Yscale0[i,], 
                                                                                              ai1  =Ph1___a_hj[[x]][2,], 
                                                                                              bi1  =Ph1___b_hj[[x]][2,], 
                                                                                              mi1  =Ph1___m_hj[[x]][2,], 
                                                                                              psii1=Ph1___psi_hj[[x]][2,]) %>% prod }) )
      } 
      tmp_Lih <- rbind(tmp_Lih, sapply(1:Ph2_n_new, function(i){   Calculate_Noncentral_t(y    =tmp_arl1_Yscale0[i,], 
                                                                                          ai1  =zero_a_hj, 
                                                                                          bi1  =zero_b_hj, 
                                                                                          mi1  =zero_m_hj, 
                                                                                          psii1=zero_psi_hj) %>% prod }))
      tmp_prob_table <- tmp_Lih * matrix( Ph1_prob, nrow=Ph1_H+1, ncol=Ph2_n_new )
      tmp_log_aPML   <- tmp_prob_table %>% apply(2,sum) %>% log %>% mean
      M2_3_xbar      <- c(M2_3_xbar, tmp_log_aPML)
      ### ------------------------------------------------------------------------------------------------------
      
      ### Competing 3 ------------------------------------------------------------------------------------------
      tmp_mean <- apply(tmp_arl1_Yscale0,2,mean)
      M4_T2    <- c(M4_T2, Ph2_n_new * c(t(tmp_mean - Ph1_Mean) %*% Ph1_invS %*% (tmp_mean - Ph1_Mean)))
      ### ------------------------------------------------------------------------------------------------------
    }
    
    # proposed model
    ## setting
    new_table <- ARL1_result %>% sapply(function(x){ x[[2]] }) %>% c %>% table
    Int2_idx <- which(target_Ph1_pis >= lower_prob_thres)
    Int3_idx <- which(target_Ph1_pis <  lower_prob_thres) %>% intersect(seq_along(new_table))
    
    ## interval types
    alpha0 <- 0.0006110553
    Interval_type2 <- Interval3(Ph1_pi_h, Ph2_n_new, alpha0)
    UCL_type2      <- Interval_type2[[1]]
    LCL_type2      <- Interval_type2[[2]]
    
    alpha <- 0.00781407
    Interval_type3 <- Interval3(Ph1_pi_h, Ph2_n_new, alpha)
    UCL_type3      <- Interval_type3[[1]]
    LCL_type3      <- Interval_type3[[2]]
    
    ## plot 1
    chk_tf <- NULL
    for(k in seq_along(new_table)){
      
      target_ps <- ARL1_result %>% sapply(function(x){ sum(x[[2]] == k, na.rm=T)/sum(!is.na(x[[2]])) })
      target_CL <- target_Ph1_pis[k] 
      
      int2_yn   <- k %in% Int2_idx
      CtrlL_col <- ifelse(int2_yn, 2, 3)
      
      if(int2_yn){
        
        target_UCL <- UCL_type2[k]
        target_LCL <- LCL_type2[k]
        
        tmp_tf <- (target_ps > target_UCL) | (target_ps < target_LCL)
        
      }else{
        
        target_UCL <- UCL_type3[k]
        target_LCL <- NA
        
        tmp_tf <- (target_ps > target_UCL)
      }
      chk_tf <- rbind(chk_tf, tmp_tf)
    }
    RL <- chk_tf %>% apply(2, function(x){ any(x, na.rm=T) }) %>% which
    RL <- ifelse(length(RL) >= 1, RL[1], Phase2_SS)
    
    # competing model 1 and 2 
    ## competing model 1
    cusum_k    <- 1/2
    PS_CM2y    <- (M2_3_xbar - Ph1_xbar)/(Ph1_s/sqrt(Ph2_n))
    PS_CM2_Cip <- PS_CM2_Cim <- NULL
    for(w in 1:Phase2_SS){
      
      if(w > 1){
        prev_cip <- PS_CM2_Cip %>% tail(1)
        prev_cim <- PS_CM2_Cim %>% tail(1)
      }else{
        prev_cip <- prev_cim <- 0
      }
      
      PS_CM2_Cip <- c(PS_CM2_Cip, max( PS_CM2y[w] - cusum_k + prev_cip, 0))
      PS_CM2_Cim <- c(PS_CM2_Cim, max(-PS_CM2y[w] - cusum_k + prev_cim, 0))
    }
    ARL1_result_cm1 <- rbind(PS_CM2_Cip, PS_CM2_Cim)
    RL_cm1 <- which((PS_CM2_Cip >= M2_h) | (PS_CM2_Cim >= M2_h))
    RL_cm1 <- ifelse(length(RL_cm1) >= 1, RL_cm1[1], Phase2_SS)
    
    ## competing model 2
    ARL1_result_cm2 <- M2_3_xbar
    tmp_ucl <- Ph1_xbar + M3_L * Ph1_s/sqrt(Ph2_n_new)
    tmp_lcl <- Ph1_xbar - M3_L * Ph1_s/sqrt(Ph2_n_new)
    RL_cm2 <- which((ARL1_result_cm2 > tmp_ucl) | (ARL1_result_cm2 < tmp_lcl))
    RL_cm2 <- ifelse(length(RL_cm2) >= 1, RL_cm2[1], Phase2_SS)
    
    ## competing model 3
    ARL1_result_cm3 <- M4_T2
    RL_cm3 <- which(ARL1_result_cm3 > qchisq(1-M4_alpha, K))
    RL_cm3 <- ifelse(length(RL_cm3) >= 1, RL_cm3[1], Phase2_SS)
    
    ## save result 
    Raw_PSs <- list(ARL1_result, ARL1_result_cm1, ARL1_result_cm2, ARL1_result_cm3)
    RLs     <- list(RL, RL_cm1, RL_cm2, RL_cm3)
    
    tmp_RL_result <- list(Raw_PSs, RLs)
    tmp_RL_result
  }
  stopCluster(myCluster)
  
  # save
  SC3_ARL1_result[[delta_idx]] <- output
}
# -------------------------------------------------------------------------------------------------
#
####################################################################################################
timer5  <- proc.time()[3]
Contime <- (timer5 - timer4)/60
cat(Contime, "min! \n")
# save(SC3_ARL1_result, file="SC3_ARL1.RData")



###################################################################################################
# Plot - figure 13, 14, 15, 16, 17
# -------------------------------------------------------------------------------------------------

# figure 13 ---------------------------------------------------------------------------------------
## data
load(file="SC1_ARL1.RData")
ARL1_result <- SC1_ARL1_result[[4]][[2]][[1]][[1]]
Ph2_n_new   <- Ph2_n + round(0.008 * each_n)
tmp_order <- c(1,2,20,6,17,14,24,5,15,18,4,11,3,27,25,8,
               16,7,13,29,23,19,22,12,21,26,9,10,28,30)

## setting
new_table <- ARL1_result %>% sapply(function(x){ x[[2]] }) %>% c %>% table

target_Ph1_pis   <- Ph1_pi_h2
lower_prob_thres <- 0.005
Int2_idx <- which(target_Ph1_pis >= lower_prob_thres)
Int3_idx <- which(target_Ph1_pis <  lower_prob_thres) %>% intersect(seq_along(new_table))

## interval types
alpha1 <- M1_alpha0
Interval_type2 <- Interval3(Ph1_pi_h, Ph2_n, alpha1)
UCL_type2      <- Interval_type2[[1]]
LCL_type2      <- Interval_type2[[2]]

alpha <- M1_alpha
Interval_type3 <- Interval3(Ph1_pi_h, Ph2_n_new, alpha)
UCL_type3      <- Interval_type3[[1]]
LCL_type3      <- Interval_type3[[2]]

# # of subplots
n_plots <- length(tmp_order)

# title for each sbuplot
if(n_plots <= 26){
  plot_labels <- paste0("(", letters[1:n_plots], ")")
} else {
  plot_labels <- c(paste0("(", letters[1:26], ")"),
                   sapply(27:n_plots, function(i) {
                     prefix <- letters[floor((i - 1) / 26)]
                     suffix <- letters[((i - 1) %% 26) + 1]
                     paste0("(", prefix, suffix, ")")
                   }))
}

# plot
plot_list <- list()
for(i in seq_along(tmp_order)) {
  k <- tmp_order[i]
  
  # plotting statistics
  target_ps <- sapply(ARL1_result, function(x) { 
    sum(x[[2]] == k, na.rm = TRUE) / sum(!is.na(x[[2]])) 
  })
  
  # centerline and CLs
  target_CL <- target_Ph1_pis[k]
  int2_yn <- k %in% Int2_idx
  ctrl_col <- if (int2_yn) "#E41A1C" else "#4DAF4A"
  
  if(int2_yn){
    target_UCL <- UCL_type2[k]
    target_LCL <- LCL_type2[k]
  } else {
    target_UCL <- UCL_type3[k]
    target_LCL <- NA
  }
  
  # ylim
  ylim_val <- range(c(target_ps, target_UCL, target_LCL, target_CL), na.rm = TRUE)
  
  # to df
  df <- data.frame(x = seq_along(target_ps), y = target_ps)
  
  # ggplot 
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line(size = 1, color = "#377EB8") +
    geom_hline(yintercept = target_UCL, linetype = "dashed", 
               color = ctrl_col, size = 1) +
    geom_hline(yintercept = target_CL, linetype = "solid", 
               color = "#000000", size = 1) +
    labs(title = plot_labels[i], x = "t", y = "Proportion") +
    scale_y_continuous(limits = ylim_val) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10)
    )
  
  if (!is.na(target_LCL)) {
    p <- p + geom_hline(yintercept = target_LCL, linetype = "dashed", 
                        color = ctrl_col, size = 1)
  }
  
  plot_list[[i]] <- p
}

# 4 x 8 grid
combined_plot <- do.call(grid.arrange, c(plot_list, nrow = 4, ncol = 8))

# # save
# ggsave("figure13.png", combined_plot, width = 20, height = 10, dpi = 1000, bg = "white")
# -------------------------------------------------------------------------------------------------

# figure 14 ---------------------------------------------------------------------------------------
## data
load(file="SC2_ARL1.RData")
ARL1_result <- SC2_ARL1_result[[4]][[1]][[1]][[1]]
Ph2_n_new   <- Ph2_n + 2 * round(0.006 * each_n)
tmp_order <- c(1,2,20,6,17,14,24,5,15,18,4,11,3,27,25,8,
               16,7,13,29,23,19,22,12,21,26,9,10,28,30,31)

## setting
new_table <- ARL1_result %>% sapply(function(x){ x[[2]] }) %>% c %>% table

target_Ph1_pis   <- Ph1_pi_h2
lower_prob_thres <- 0.005
Int2_idx <- which(target_Ph1_pis >= lower_prob_thres)
Int3_idx <- which(target_Ph1_pis <  lower_prob_thres) %>% intersect(seq_along(new_table))

## interval types
alpha1 <- M1_alpha0
Interval_type2 <- Interval3(Ph1_pi_h, Ph2_n, alpha1)
UCL_type2      <- Interval_type2[[1]]
LCL_type2      <- Interval_type2[[2]]

alpha <- M1_alpha
Interval_type3 <- Interval3(Ph1_pi_h, Ph2_n_new, alpha)
UCL_type3      <- Interval_type3[[1]]
LCL_type3      <- Interval_type3[[2]]

# # of subplots
n_plots <- length(tmp_order)

# title for each sbuplot
if(n_plots <= 26){
  plot_labels <- paste0("(", letters[1:n_plots], ")")
} else {
  plot_labels <- c(paste0("(", letters[1:26], ")"),
                   sapply(27:n_plots, function(i) {
                     prefix <- letters[floor((i - 1) / 26)]
                     suffix <- letters[((i - 1) %% 26) + 1]
                     paste0("(", prefix, suffix, ")")
                   }))
}

# plot
plot_list <- list()
for(i in seq_along(tmp_order)) {
  k <- tmp_order[i]
  
  # plotting statistics
  target_ps <- sapply(ARL1_result, function(x) { 
    sum(x[[2]] == k, na.rm = TRUE) / sum(!is.na(x[[2]])) 
  })
  
  # centerline and CLs
  target_CL <- target_Ph1_pis[k]
  int2_yn <- k %in% Int2_idx
  ctrl_col <- if (int2_yn) "#E41A1C" else "#4DAF4A"
  
  if(int2_yn){
    target_UCL <- UCL_type2[k]
    target_LCL <- LCL_type2[k]
  } else {
    target_UCL <- UCL_type3[k]
    target_LCL <- NA
  }
  
  # ylim
  ylim_val <- range(c(target_ps, target_UCL, target_LCL, target_CL), na.rm = TRUE)
  
  # to df
  df <- data.frame(x = seq_along(target_ps), y = target_ps)
  
  # ggplot 
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line(size = 1, color = "#377EB8") +
    geom_hline(yintercept = target_UCL, linetype = "dashed", 
               color = ctrl_col, size = 1) +
    geom_hline(yintercept = target_CL, linetype = "solid", 
               color = "#000000", size = 1) +
    labs(title = plot_labels[i], x = "t", y = "Proportion") +
    scale_y_continuous(limits = ylim_val) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10)
    )
  
  if (!is.na(target_LCL)) {
    p <- p + geom_hline(yintercept = target_LCL, linetype = "dashed", 
                        color = ctrl_col, size = 1)
  }
  
  plot_list[[i]] <- p
}

# 4 x 8 grid
combined_plot <- do.call(grid.arrange, c(plot_list, nrow = 4, ncol = 8))

# # save
# ggsave("figure14.png", combined_plot, width = 20, height = 10, dpi = 1000, bg = "white")
# -------------------------------------------------------------------------------------------------

# figure 15 ---------------------------------------------------------------------------------------
## data
load(file="SC3_ARL1.RData")
ARL1_result <- SC3_ARL1_result[[4]][[1]][[1]][[1]]
Ph2_n_new   <- Ph2_n + round(0.02 * each_n)
tmp_order <- c(1,2,20,6,17,14,24,5,15,18,4,11,3,27,25,8,
               16,7,13,29,23,19,22,12,21,26,9,10,28)

## setting
new_table <- ARL1_result %>% sapply(function(x){ x[[2]] }) %>% c %>% table

target_Ph1_pis   <- Ph1_pi_h2
lower_prob_thres <- 0.005
Int2_idx <- which(target_Ph1_pis >= lower_prob_thres)
Int3_idx <- which(target_Ph1_pis <  lower_prob_thres) %>% intersect(seq_along(new_table))

## interval types
alpha1 <- M1_alpha0
Interval_type2 <- Interval3(Ph1_pi_h, Ph2_n, alpha1)
UCL_type2      <- Interval_type2[[1]]
LCL_type2      <- Interval_type2[[2]]

alpha <- M1_alpha
Interval_type3 <- Interval3(Ph1_pi_h, Ph2_n_new, alpha)
UCL_type3      <- Interval_type3[[1]]
LCL_type3      <- Interval_type3[[2]]

# # of subplots
n_plots <- length(tmp_order)

# title for each sbuplot
if(n_plots <= 26){
  plot_labels <- paste0("(", letters[1:n_plots], ")")
} else {
  plot_labels <- c(paste0("(", letters[1:26], ")"),
                   sapply(27:n_plots, function(i) {
                     prefix <- letters[floor((i - 1) / 26)]
                     suffix <- letters[((i - 1) %% 26) + 1]
                     paste0("(", prefix, suffix, ")")
                   }))
}

# plot
plot_list <- list()
for(i in seq_along(tmp_order)) {
  k <- tmp_order[i]
  
  # plotting statistics
  target_ps <- sapply(ARL1_result, function(x) { 
    sum(x[[2]] == k, na.rm = TRUE) / sum(!is.na(x[[2]])) 
  })
  
  # centerline and CLs
  target_CL <- target_Ph1_pis[k]
  int2_yn <- k %in% Int2_idx
  ctrl_col <- if (int2_yn) "#E41A1C" else "#4DAF4A"
  
  if(int2_yn){
    target_UCL <- UCL_type2[k]
    target_LCL <- LCL_type2[k]
  } else {
    target_UCL <- UCL_type3[k]
    target_LCL <- NA
  }
  
  # ylim
  ylim_val <- range(c(target_ps, target_UCL, target_LCL, target_CL), na.rm = TRUE)
  
  # to df
  df <- data.frame(x = seq_along(target_ps), y = target_ps)
  
  # ggplot 
  p <- ggplot(df, aes(x = x, y = y)) +
    geom_line(size = 1, color = "#377EB8") +
    geom_hline(yintercept = target_UCL, linetype = "dashed", 
               color = ctrl_col, size = 1) +
    geom_hline(yintercept = target_CL, linetype = "solid", 
               color = "#000000", size = 1) +
    labs(title = plot_labels[i], x = "t", y = "Proportion") +
    scale_y_continuous(limits = ylim_val) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text  = element_text(size = 10)
    )
  
  if (!is.na(target_LCL)) {
    p <- p + geom_hline(yintercept = target_LCL, linetype = "dashed", 
                        color = ctrl_col, size = 1)
  }
  
  plot_list[[i]] <- p
}

# 4 x 8 grid
combined_plot <- do.call(grid.arrange, c(plot_list, nrow = 4, ncol = 8))

# # save
# ggsave("figure15.png", combined_plot, width = 20, height = 10, dpi = 1000, bg = "white")
# -------------------------------------------------------------------------------------------------

# figure 16 ---------------------------------------------------------------------------------------
ARL1_result <- SC2_ARL1_result[[4]][[30]][[1]][[1]]

target_cls   <- c(30,31)
target_times <- c(1,50)
tmp_mat      <- NULL

for(ttt in target_times){
  
  tmp_arl1 <- ARL1_result[[ttt]]
  
  # MC - cluster specific mean surface
  I <- 10000
  target_parameters <- list( sapply(tmp_arl1[[7]],   function(x){ x %>% tail(1) %>% c }),
                             sapply(tmp_arl1[[6]],   function(x){ x %>% tail(1) %>% c }),
                             sapply(tmp_arl1[[8]],   function(x){ x %>% tail(1) %>% c }),
                             sapply(tmp_arl1[[9]],   function(x){ x %>% tail(1) %>% c }) )
  tmp_m   <- target_parameters[[1]]
  tmp_psi <- target_parameters[[2]]
  tmp_a   <- target_parameters[[3]]
  tmp_b   <- target_parameters[[4]]
  
  tmp_clusters <- 1:ncol(tmp_m)
  tmp_mean_function <- NULL
  for(cluster in tmp_clusters){
    
    # MC simulation
    tmp_fun12 <- function(i){
      
      tmp_tau <- rgamma(nrow(tmp_a), shape=tmp_a[,cluster], rate=tmp_b[,cluster])
      tmp_mu  <- rnorm(nrow(tmp_m), tmp_m[,cluster], sqrt(tmp_psi[,cluster]/tmp_tau))
      return(tmp_mu)
    }
    tmp_MCsample <- sapply(1:I, function(i){ tmp_fun12(i) }) %>% t
    
    tmp_pararmeter <- tmp_MCsample %>% apply(2,mean)
    tmp_pararmeter <- tmp_pararmeter * Step1_scale_scale + Step1_scale_center
    tmp_mean_vals  <- c(Bmat %*% tmp_pararmeter)
    
    # Remove_zero
    tmp_mean_full                  <- rep(NA, prod(dat_dim))
    tmp_mean_full[WF_outise_idx]   <- NA
    tmp_mean_full[-WF_outise_idx]  <- tmp_mean_vals
    tmp_mean_vals                  <- tmp_mean_full
    
    # save
    tmp_mean_function <- rbind(tmp_mean_function, tmp_mean_vals)
  }
  tmp_mat <- rbind(tmp_mat, tmp_mean_function[target_cls,])
}
tmp_mat <- tmp_mat + 0.5
tmp_mat %>% quantile(c(0,1), na.rm=T)

zlim     <- c(2, 7.1)
tmp_mat2 <- sapply(1:4, function(i){ tmp_mat[i,] %>% matrix(nrow=dat_dim[1], ncol=dat_dim[2]) %>% t }, simplify=F)

# zlim 
zlim <- c(2, 7.1)

# data
tmp_mat2 <- sapply(1:4, function(i){
  tmp_mat[i,] %>% matrix(nrow = dat_dim[1], ncol = dat_dim[2]) %>% t
}, simplify = FALSE)

# titles
plot_titles <- c("(a)", "(b)", "(c)", "(d)")

# plot
plot_list <- list()

for(i in 1:4){
  
  df <- melt(tmp_mat2[[i]] %>% t)
  df$x <- x_grid[df$Var2]
  df$y <- y_grid[df$Var1]
  
  # ggplot2
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(option = "H",
                         limits = zlim, 
                         na.value = "white") +
    labs(title = plot_titles[i], x = "x", y = "y") +
    coord_fixed(expand = FALSE) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
      axis.text  = element_text(size = 10),
      axis.title = element_text(size = 12),
      plot.margin = unit(c(0.5, 0.5, 0.5, 0.5), "cm")
    )
  
  plot_list[[i]] <- p
}

# 1 x 4 grid
combined_plot <- grid.arrange(grobs = plot_list, nrow = 1, ncol = 4)

# # save
# ggsave("figure16.png", combined_plot, width = 16, height = 4, dpi = 1000, bg = "white")
# -------------------------------------------------------------------------------------------------

# figure 17 ---------------------------------------------------------------------------------------
## data for scenario 1
delta_list <- c(0.002, 0.004, 0.006, 0.008)
ARL1_m1 <- c(SC1_ARL1_result[[1]] %>% sapply(function(x){ x[[2]][[1]] }) %>% mean,
             SC1_ARL1_result[[2]] %>% sapply(function(x){ x[[2]][[1]] }) %>% mean,
             SC1_ARL1_result[[3]] %>% sapply(function(x){ x[[2]][[1]] }) %>% mean,
             SC1_ARL1_result[[4]] %>% sapply(function(x){ x[[2]][[1]] }) %>% mean)

ARL1_m2 <- c(SC1_ARL1_result[[1]] %>% sapply(function(x){ x[[2]][[2]] }) %>% mean,
             SC1_ARL1_result[[2]] %>% sapply(function(x){ x[[2]][[2]] }) %>% mean,
             SC1_ARL1_result[[3]] %>% sapply(function(x){ x[[2]][[2]] }) %>% mean,
             SC1_ARL1_result[[4]] %>% sapply(function(x){ x[[2]][[2]] }) %>% mean)

ARL1_m3 <- c(SC1_ARL1_result[[1]] %>% sapply(function(x){ x[[2]][[3]] }) %>% mean,
             SC1_ARL1_result[[2]] %>% sapply(function(x){ x[[2]][[3]] }) %>% mean,
             SC1_ARL1_result[[3]] %>% sapply(function(x){ x[[2]][[3]] }) %>% mean,
             SC1_ARL1_result[[4]] %>% sapply(function(x){ x[[2]][[3]] }) %>% mean)

ARL1_m4 <- c(SC1_ARL1_result[[1]] %>% sapply(function(x){ x[[2]][[4]] }) %>% mean,
             SC1_ARL1_result[[2]] %>% sapply(function(x){ x[[2]][[4]] }) %>% mean,
             SC1_ARL1_result[[3]] %>% sapply(function(x){ x[[2]][[4]] }) %>% mean,
             SC1_ARL1_result[[4]] %>% sapply(function(x){ x[[2]][[4]] }) %>% mean)
ARL1_sc1    <- rbind(ARL1_m1,ARL1_m2,ARL1_m3,ARL1_m4)
delta_list1 <- delta_list

## data for scenario 2
delta_list <- c(0.005, 0.01, 0.015, 0.02)
ARL1_m1 <- c(SC3_ARL1_result[[1]] %>% sapply(function(x){ x[[2]][[1]] }) %>% mean,
             SC3_ARL1_result[[2]] %>% sapply(function(x){ x[[2]][[1]] }) %>% mean,
             SC3_ARL1_result[[3]] %>% sapply(function(x){ x[[2]][[1]] }) %>% mean,
             SC3_ARL1_result[[4]] %>% sapply(function(x){ x[[2]][[1]] }) %>% mean)

ARL1_m2 <- c(SC3_ARL1_result[[1]] %>% sapply(function(x){ x[[2]][[2]] }) %>% mean,
             SC3_ARL1_result[[2]] %>% sapply(function(x){ x[[2]][[2]] }) %>% mean,
             SC3_ARL1_result[[3]] %>% sapply(function(x){ x[[2]][[2]] }) %>% mean,
             SC3_ARL1_result[[4]] %>% sapply(function(x){ x[[2]][[2]] }) %>% mean)

ARL1_m3 <- c(SC3_ARL1_result[[1]] %>% sapply(function(x){ x[[2]][[3]] }) %>% mean,
             SC3_ARL1_result[[2]] %>% sapply(function(x){ x[[2]][[3]] }) %>% mean,
             SC3_ARL1_result[[3]] %>% sapply(function(x){ x[[2]][[3]] }) %>% mean,
             SC3_ARL1_result[[4]] %>% sapply(function(x){ x[[2]][[3]] }) %>% mean)

ARL1_m4 <- c(SC3_ARL1_result[[1]] %>% sapply(function(x){ x[[2]][[4]] }) %>% mean,
             SC3_ARL1_result[[2]] %>% sapply(function(x){ x[[2]][[4]] }) %>% mean,
             SC3_ARL1_result[[3]] %>% sapply(function(x){ x[[2]][[4]] }) %>% mean,
             SC3_ARL1_result[[4]] %>% sapply(function(x){ x[[2]][[4]] }) %>% mean)
ARL1_sc2    <- rbind(ARL1_m1,ARL1_m2,ARL1_m3,ARL1_m4)
delta_list2 <- delta_list

## data for scenario 3
delta_list <- c(0.0015, 0.003, 0.0045, 0.006)
ARL1_m1 <- c(SC2_ARL1_result[[1]] %>% sapply(function(x){ x[[2]][[1]] }) %>% mean,
             SC2_ARL1_result[[2]] %>% sapply(function(x){ x[[2]][[1]] }) %>% mean,
             SC2_ARL1_result[[3]] %>% sapply(function(x){ x[[2]][[1]] }) %>% mean,
             SC2_ARL1_result[[4]] %>% sapply(function(x){ x[[2]][[1]] }) %>% mean)

ARL1_m2 <- c(SC2_ARL1_result[[1]] %>% sapply(function(x){ x[[2]][[2]] }) %>% mean,
             SC2_ARL1_result[[2]] %>% sapply(function(x){ x[[2]][[2]] }) %>% mean,
             SC2_ARL1_result[[3]] %>% sapply(function(x){ x[[2]][[2]] }) %>% mean,
             SC2_ARL1_result[[4]] %>% sapply(function(x){ x[[2]][[2]] }) %>% mean)

ARL1_m3 <- c(SC2_ARL1_result[[1]] %>% sapply(function(x){ x[[2]][[3]] }) %>% mean,
             SC2_ARL1_result[[2]] %>% sapply(function(x){ x[[2]][[3]] }) %>% mean,
             SC2_ARL1_result[[3]] %>% sapply(function(x){ x[[2]][[3]] }) %>% mean,
             SC2_ARL1_result[[4]] %>% sapply(function(x){ x[[2]][[3]] }) %>% mean)

ARL1_m4 <- c(SC2_ARL1_result[[1]] %>% sapply(function(x){ x[[2]][[4]] }) %>% mean,
             SC2_ARL1_result[[2]] %>% sapply(function(x){ x[[2]][[4]] }) %>% mean,
             SC2_ARL1_result[[3]] %>% sapply(function(x){ x[[2]][[4]] }) %>% mean,
             SC2_ARL1_result[[4]] %>% sapply(function(x){ x[[2]][[4]] }) %>% mean)
ARL1_sc3    <- rbind(ARL1_m1,ARL1_m2,ARL1_m3,ARL1_m4)
delta_list3 <- delta_list


# ----- Panel 1 -----
n1 <- length(delta_list1)
data1 <- data.frame(
  delta = rep(delta_list1, times = 4),
  ARL1 = c(ARL1_sc1[1, ], ARL1_sc1[2, ], ARL1_sc1[3, ], ARL1_sc1[4, ]),
  Group = factor(rep(1:4, each = n1), 
                 levels = 1:4,
                 labels = c("Proposed model", "DP-based SPC1", "DP-based SPC2", "T2 chart"))
)
ylim1 <- range(as.vector(ARL1_sc1), na.rm = TRUE)

p1 <- ggplot(data1, aes(x = delta, y = ARL1, color = Group, group = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "(a)", x = "delta", y = "ARL1") +
  scale_y_continuous(limits = ylim1) +
  scale_color_manual(name = "Model", 
                     values = c("black", "#F8766D", "#7CAE00", "#00BFC4")) +
  theme_minimal() +
  theme(plot.title   = element_text(hjust = 0.5, size = 28, face = "bold"),
        axis.text    = element_text(size = 14),
        axis.title   = element_text(size = 18),
        legend.position = "bottom")

# ----- Panel 2 -----
n2 <- length(delta_list2)
data2 <- data.frame(
  delta = rep(delta_list2, times = 4),
  ARL1 = c(ARL1_sc2[1, ], ARL1_sc2[2, ], ARL1_sc2[3, ], ARL1_sc2[4, ]),
  Group = factor(rep(1:4, each = n2), 
                 levels = 1:4,
                 labels = c("Proposed model", "DP-based SPC1", "DP-based SPC2", "T2 chart"))
)
ylim2 <- range(as.vector(ARL1_sc2), na.rm = TRUE)

p2 <- ggplot(data2, aes(x = delta, y = ARL1, color = Group, group = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "(b)", x = "delta", y = "ARL1") +
  scale_y_continuous(limits = ylim2) +
  scale_color_manual(name = "Model", 
                     values = c("black", "#F8766D", "#7CAE00", "#00BFC4")) +
  theme_minimal() +
  theme(plot.title   = element_text(hjust = 0.5, size = 28, face = "bold"),
        axis.text    = element_text(size = 14),
        axis.title   = element_text(size = 18),
        legend.position = "bottom")

# ----- Panel 3 -----
n3 <- length(delta_list3)
data3 <- data.frame(
  delta = rep(delta_list3, times = 4),
  ARL1 = c(ARL1_sc3[1, ], ARL1_sc3[2, ], ARL1_sc3[3, ], ARL1_sc3[4, ]),
  Group = factor(rep(1:4, each = n3), 
                 levels = 1:4,
                 labels = c("Proposed model", "DP-based SPC1", "DP-based SPC2", "T2 chart"))
)
ylim3 <- range(as.vector(ARL1_sc3), na.rm = TRUE)

p3 <- ggplot(data3, aes(x = delta, y = ARL1, color = Group, group = Group)) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  labs(title = "(c)", x = "delta", y = "ARL1") +
  scale_y_continuous(limits = ylim3) +
  scale_color_manual(name = "Model", 
                     values = c("black", "#F8766D", "#7CAE00", "#00BFC4")) +
  theme_minimal() +
  theme(plot.title   = element_text(hjust = 0.5, size = 28, face = "bold"),
        axis.text    = element_text(size = 14),
        axis.title   = element_text(size = 18),
        legend.position = "bottom")

# combine plots
combined_plot <- (p1 + p2 + p3) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# legend
combined_plot <- (p1 + p2 + p3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.title = element_text(size = 24),
        legend.text = element_text(size = 20),
        legend.key.size = unit(1, "cm"))

# # save
# ggsave("figure17.png", combined_plot, width = 20, height = 8, dpi = 1000, bg = "white")
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
#
###################################################################################################