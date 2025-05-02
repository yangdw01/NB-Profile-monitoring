###################################################################################################
# Basic setting
# -------------------------------------------------------------------------------------------------

# directory ---------------------------------------------------------------------------------------
setwd("C:\\Users\\yangd\\OneDrive\\0. Research\\2. JQT\\Revision\\Code\\")
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
library(mvtnorm)
library(mclust)
library(ggplot2)
library(patchwork)
# -------------------------------------------------------------------------------------------------

# source ------------------------------------------------------------------------------------------
source("Functions.R")
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
#
###################################################################################################


###################################################################################################
# Data
# -------------------------------------------------------------------------------------------------
## basic setting 1 for data generation
data_gen <- function(seed){
  
  set.seed(seed)
  
  K <- 2
  cl_num1       <- 2
  each_cl_num1  <- 1000
  each_cl_m1    <- rbind(c(2,-3.5),c(2,2))
  each_cl_sig1  <- c(2,0.5)^2
  
  ## large-size cluster
  True_Cluster <- Raw_Data <- NULL
  for(k in seq(cl_num1)){
    
    tmp_m       <- each_cl_m1[,k]
    tmp_sig     <- each_cl_sig1[k]
    tmp_raw_dat <- rmvnorm(each_cl_num1, mean=tmp_m, sigma=diag(rep(tmp_sig,K)))
    
    Raw_Data     <- rbind(Raw_Data, tmp_raw_dat)
    True_Cluster <- c(True_Cluster, rep(k, each_cl_num1))
  }
  
  ## basic setting 2 for data generation
  cl_num2       <- 1
  each_cl_num2s <- 30
  each_cl_m2    <- rbind(c(-2,-2),c(-1,-3))
  each_cl_sig2  <- c(0.2^2, 0.05^2)
  
  ## small-size cluster
  for(k in seq(cl_num2)){
    
    tmp_n       <- each_cl_num2s[k]
    tmp_m       <- each_cl_m2[,k]
    tmp_sig     <- each_cl_sig2[k]
    tmp_raw_dat <- rmvnorm(tmp_n, mean=tmp_m, sigma=diag(rep(tmp_sig,K)))
    
    Raw_Data     <- rbind(Raw_Data, tmp_raw_dat)
    True_Cluster <- c(True_Cluster, rep(k+cl_num1, tmp_n))
  }
  
  return(list(Raw_Data, True_Cluster))
}



tmp_data <- data_gen(1)
Raw_Data     <- tmp_data[[1]]
True_Cluster <- tmp_data[[2]]

cl_num1 <- 2
cl_num2 <- 1
## plot
xlim <- Raw_Data[,1] %>% quantile(c(0,1))
ylim <- Raw_Data[,2] %>% quantile(c(0,1))
plot(-100,100, xlim=xlim, ylim=ylim)
for(k in seq(cl_num1+cl_num2)){
  
  tmp_pch   <- ifelse(k<=cl_num1, 1, 3)
  tmp_k_idx <- which(True_Cluster == k)
  
  points(Raw_Data[tmp_k_idx,1], Raw_Data[tmp_k_idx,2], col=k, pch=tmp_pch)
}
# -------------------------------------------------------------------------------------------------
#
###################################################################################################



###################################################################################################
# (Step 2-0) Setting for SUGS algorithm
# -------------------------------------------------------------------------------------------------

## scaling
Raw_Data_center <- Raw_Data %>% scale %>% attr("scaled:center")
Raw_Data_scale  <- Raw_Data %>% scale %>% attr("scaled:scale")
Yscale          <- Raw_Data %>% scale # ((t(total_Y) - Step1_scale_center)/Step1_scale_scale) %>% t

## basic setting
n  <- nrow(Yscale)

## priors
a   <- 1
b   <- 1
m   <- 0
psi <- 1
alpha_stars <- c(0.01, 0.05, c(0.1 + 0.2 * (0:19)))
eta_probs   <- rep(1/length(alpha_stars), length(alpha_stars))
TT          <- length(alpha_stars)

## plot (only 2-dim)
xlim <- Yscale[,1] %>% quantile(c(0,1))
ylim <- Yscale[,2] %>% quantile(c(0,1))
plot(-100,100, xlim=xlim, ylim=ylim)
for(k in seq(cl_num1+cl_num2)){
  
  tmp_pch   <- ifelse(k<=cl_num1, 1, 3)
  tmp_k_idx <- which(True_Cluster == k)
  
  points(Yscale[tmp_k_idx,1], Yscale[tmp_k_idx,2], col=k, pch=tmp_pch)
}
# -------------------------------------------------------------------------------------------------
#
###################################################################################################


###################################################################################################
# Simulation
# -------------------------------------------------------------------------------------------------
SimN <- 500

### SUGS algorithm - parallel
myCluster <- makeCluster(30)
registerDoParallel(myCluster)
output <- foreach(w = 1:SimN, .packages = c("dplyr", "dqrng", "mvtnorm", "mclust")) %dopar% {
  
  ## seed
  set.seed(w)
  
  ## data generation
  tmp_data     <- data_gen(w)
  raw_data     <- tmp_data[[1]]
  true_cluster <- tmp_data[[2]]
  
  ## scaling
  Yscale       <- raw_data %>% scale
  
  ## basic setting
  n  <- nrow(raw_data)
  K  <- ncol(raw_data)
  
  ## priors
  a   <- 1
  b   <- 1
  m   <- 0
  psi <- 1
  alpha_stars <- c(0.01, 0.05, c(0.1 + 0.2 * (0:19)))
  eta_probs   <- rep(1/length(alpha_stars), length(alpha_stars))
  TT          <- length(alpha_stars)
  
  ## BSUGS algorithm --------------------------------------------------------------------------------
  ### modeling setting
  bidirec   <- T
  back_uNum <- NULL
  Onum      <- 200
  ITERmax   <- 20
  RePUpd    <- T
  
  ### BSUGS algorithm
  SUGS_ordering_num1 <- 20
  
  tmp_RESULT1 <- tmp_RESULT2 <- list()
  for(ww in 1:SUGS_ordering_num1){
    
    # ordering
    dqset.seed(ww  + w*103)
    tmp_resample <- dqsample(1:n)
    tmp_Y        <- Yscale[tmp_resample,]
    
    # SUGS algorithm
    set.seed(ww + w*103)
    tmp_result   <- BSUGS(Y=tmp_Y, a=a, b=b, m=m, psi=psi, alpha_stars=alpha_stars, eta_probs=eta_probs,
                          bidirectional = bidirec, back_upperNum = back_uNum,
                          OnceNum       = Onum,    IterMax       = ITERmax,   ReParmUpd   = RePUpd)
    tmp_RESULT1[[ww]] <- tmp_result
    tmp_RESULT2[[ww]] <- tmp_resample
  }
  tmp_aPMLs1    <- tmp_RESULT1 %>% sapply(function(x){ x[[1]] })
  tmp_max_idx1  <- tmp_aPMLs1 %>% which.max
  tmp_result1   <- tmp_RESULT1[[tmp_max_idx1]]
  tmp_ordering1 <- tmp_RESULT2[[tmp_max_idx1]]
  tmp_aPML1     <- tmp_aPMLs1 %>% max
  tmp_gamma1    <- tmp_result1[[12]][tmp_ordering1 %>% order]
  tmp_adjRI1    <- adjustedRandIndex(tmp_gamma1, True_Cluster)
  tmp_adjRI1s   <- seq(tmp_RESULT1) %>% sapply(function(x){ tmp_RESULT1[[x]][[12]][tmp_RESULT2[[x]] %>% order] %>% adjustedRandIndex(True_Cluster) })
  
  tmp_BSUGS <- list(tmp_gamma1, tmp_adjRI1, tmp_aPML1, tmp_aPMLs1, tmp_adjRI1s)
  # -------------------------------------------------------------------------------------------------
  
  ## BSUGS algorithm --------------------------------------------------------------------------------
  ### modeling setting
  bidirec   <- T
  back_uNum <- NULL
  Onum      <- 200
  ITERmax   <- 20
  RePUpd    <- T
  
  ### BSUGS algorithm
  SUGS_ordering_num1 <- 20
  
  tmp_RESULT1 <- tmp_RESULT2 <- list()
  for(ww in 1:SUGS_ordering_num1){
    
    # ordering
    dqset.seed(ww  + w*103)
    tmp_resample <- dqsample(1:n)
    tmp_Y        <- Yscale[tmp_resample,]
    
    # SUGS algorithm
    set.seed(ww + w*103)
    tmp_result   <- BSUGS(Y=tmp_Y, a=a, b=1/2, m=m, psi=psi, alpha_stars=alpha_stars, eta_probs=eta_probs,
                          bidirectional = bidirec, back_upperNum = back_uNum,
                          OnceNum       = Onum,    IterMax       = ITERmax,   ReParmUpd   = RePUpd)
    tmp_RESULT1[[ww]] <- tmp_result
    tmp_RESULT2[[ww]] <- tmp_resample
  }
  tmp_aPMLs1    <- tmp_RESULT1 %>% sapply(function(x){ x[[1]] })
  tmp_max_idx1  <- tmp_aPMLs1 %>% which.max
  tmp_result1   <- tmp_RESULT1[[tmp_max_idx1]]
  tmp_ordering1 <- tmp_RESULT2[[tmp_max_idx1]]
  tmp_aPML1     <- tmp_aPMLs1 %>% max
  tmp_gamma1    <- tmp_result1[[12]][tmp_ordering1 %>% order]
  tmp_adjRI1    <- adjustedRandIndex(tmp_gamma1, True_Cluster)
  tmp_adjRI1s   <- seq(tmp_RESULT1) %>% sapply(function(x){ tmp_RESULT1[[x]][[12]][tmp_RESULT2[[x]] %>% order] %>% adjustedRandIndex(True_Cluster) })
  
  tmp_BSUGS2 <- list(tmp_gamma1, tmp_adjRI1, tmp_aPML1, tmp_aPMLs1, tmp_adjRI1s)
  # -------------------------------------------------------------------------------------------------
  
  ## BSUGS algorithm --------------------------------------------------------------------------------
  ### modeling setting
  bidirec   <- T
  back_uNum <- NULL
  Onum      <- 200
  ITERmax   <- 20
  RePUpd    <- T
  
  ### BSUGS algorithm
  SUGS_ordering_num1 <- 20
  
  tmp_RESULT1 <- tmp_RESULT2 <- list()
  for(ww in 1:SUGS_ordering_num1){
    
    # ordering
    dqset.seed(ww  + w*103)
    tmp_resample <- dqsample(1:n)
    tmp_Y        <- Yscale[tmp_resample,]
    
    # SUGS algorithm
    set.seed(ww + w*103)
    tmp_result   <- BSUGS(Y=tmp_Y, a=a, b=1/5, m=m, psi=psi, alpha_stars=alpha_stars, eta_probs=eta_probs,
                          bidirectional = bidirec, back_upperNum = back_uNum,
                          OnceNum       = Onum,    IterMax       = ITERmax,   ReParmUpd   = RePUpd)
    tmp_RESULT1[[ww]] <- tmp_result
    tmp_RESULT2[[ww]] <- tmp_resample
  }
  tmp_aPMLs1    <- tmp_RESULT1 %>% sapply(function(x){ x[[1]] })
  tmp_max_idx1  <- tmp_aPMLs1 %>% which.max
  tmp_result1   <- tmp_RESULT1[[tmp_max_idx1]]
  tmp_ordering1 <- tmp_RESULT2[[tmp_max_idx1]]
  tmp_aPML1     <- tmp_aPMLs1 %>% max
  tmp_gamma1    <- tmp_result1[[12]][tmp_ordering1 %>% order]
  tmp_adjRI1    <- adjustedRandIndex(tmp_gamma1, True_Cluster)
  tmp_adjRI1s   <- seq(tmp_RESULT1) %>% sapply(function(x){ tmp_RESULT1[[x]][[12]][tmp_RESULT2[[x]] %>% order] %>% adjustedRandIndex(True_Cluster) })
  
  tmp_BSUGS3 <- list(tmp_gamma1, tmp_adjRI1, tmp_aPML1, tmp_aPMLs1, tmp_adjRI1s)
  # -------------------------------------------------------------------------------------------------
  
  ## BSUGS algorithm --------------------------------------------------------------------------------
  ### modeling setting
  bidirec   <- T
  back_uNum <- NULL
  Onum      <- 200
  ITERmax   <- 20
  RePUpd    <- T
  
  ### BSUGS algorithm
  SUGS_ordering_num1 <- 20
  
  tmp_RESULT1 <- tmp_RESULT2 <- list()
  for(ww in 1:SUGS_ordering_num1){
    
    # ordering
    dqset.seed(ww  + w*103)
    tmp_resample <- dqsample(1:n)
    tmp_Y        <- Yscale[tmp_resample,]
    
    # SUGS algorithm
    set.seed(ww + w*103)
    tmp_result   <- BSUGS(Y=tmp_Y, a=a, b=1/10, m=m, psi=psi, alpha_stars=alpha_stars, eta_probs=eta_probs,
                          bidirectional = bidirec, back_upperNum = back_uNum,
                          OnceNum       = Onum,    IterMax       = ITERmax,   ReParmUpd   = RePUpd)
    tmp_RESULT1[[ww]] <- tmp_result
    tmp_RESULT2[[ww]] <- tmp_resample
  }
  tmp_aPMLs1    <- tmp_RESULT1 %>% sapply(function(x){ x[[1]] })
  tmp_max_idx1  <- tmp_aPMLs1 %>% which.max
  tmp_result1   <- tmp_RESULT1[[tmp_max_idx1]]
  tmp_ordering1 <- tmp_RESULT2[[tmp_max_idx1]]
  tmp_aPML1     <- tmp_aPMLs1 %>% max
  tmp_gamma1    <- tmp_result1[[12]][tmp_ordering1 %>% order]
  tmp_adjRI1    <- adjustedRandIndex(tmp_gamma1, True_Cluster)
  tmp_adjRI1s   <- seq(tmp_RESULT1) %>% sapply(function(x){ tmp_RESULT1[[x]][[12]][tmp_RESULT2[[x]] %>% order] %>% adjustedRandIndex(True_Cluster) })
  
  tmp_BSUGS4 <- list(tmp_gamma1, tmp_adjRI1, tmp_aPML1, tmp_aPMLs1, tmp_adjRI1s)
  # -------------------------------------------------------------------------------------------------
  
  ## BSUGS algorithm --------------------------------------------------------------------------------
  ### modeling setting
  bidirec   <- T
  back_uNum <- NULL
  Onum      <- 200
  ITERmax   <- 20
  RePUpd    <- F
  
  ### BSUGS algorithm
  SUGS_ordering_num1 <- 20
  
  tmp_RESULT1 <- tmp_RESULT2 <- list()
  for(ww in 1:SUGS_ordering_num1){
    
    # ordering
    dqset.seed(ww  + w*103)
    tmp_resample <- dqsample(1:n)
    tmp_Y        <- Yscale[tmp_resample,]
    
    # SUGS algorithm
    set.seed(ww + w*103)
    tmp_result   <- BSUGS(Y=tmp_Y, a=a, b=1/10, m=m, psi=psi, alpha_stars=alpha_stars, eta_probs=eta_probs,
                          bidirectional = bidirec, back_upperNum = back_uNum,
                          OnceNum       = Onum,    IterMax       = ITERmax,   ReParmUpd   = RePUpd)
    tmp_RESULT1[[ww]] <- tmp_result
    tmp_RESULT2[[ww]] <- tmp_resample
  }
  tmp_aPMLs1    <- tmp_RESULT1 %>% sapply(function(x){ x[[1]] })
  tmp_max_idx1  <- tmp_aPMLs1 %>% which.max
  tmp_result1   <- tmp_RESULT1[[tmp_max_idx1]]
  tmp_ordering1 <- tmp_RESULT2[[tmp_max_idx1]]
  tmp_aPML1     <- tmp_aPMLs1 %>% max
  tmp_gamma1    <- tmp_result1[[12]][tmp_ordering1 %>% order]
  tmp_adjRI1    <- adjustedRandIndex(tmp_gamma1, True_Cluster)
  tmp_adjRI1s   <- seq(tmp_RESULT1) %>% sapply(function(x){ tmp_RESULT1[[x]][[12]][tmp_RESULT2[[x]] %>% order] %>% adjustedRandIndex(True_Cluster) })
  
  tmp_BSUGS0 <- list(tmp_gamma1, tmp_adjRI1, tmp_aPML1, tmp_aPMLs1, tmp_adjRI1s)
  # -------------------------------------------------------------------------------------------------
  
  ## SUGS algorithm ---------------------------------------------------------------------------------
  ### modeling setting
  c <- 1
  d <- 10
  
  ### BSUGS algorithm
  SUGS_ordering_num2 <- 50
  
  tmp_RESULT3 <- tmp_RESULT4 <- list()
  for(ww in 1:SUGS_ordering_num2){
    
    cat(ww, "  ")
    
    # ordering
    dqset.seed(ww  + w*103)
    tmp_resample <- dqsample(1:n)
    tmp_Y        <- Yscale[tmp_resample,]
    
    # SUGS algorithm
    set.seed(ww  + w*103)
    tmp_result <- SUGS(Y=tmp_Y, a=a, b=NULL, m=m, psi=psi, c=c, d=d, alpha_stars=alpha_stars, eta_probs=eta_probs,
                       preliminary=TRUE, prior_info=FALSE, prior_Gamma_i=NULL, first_iter=TRUE, previous_Nt1=0)
    tmp_result <- SUGS(Y=tmp_Y, a=a, b=tmp_result[[2]][2,], m=m, psi=psi, c=c, d=d, alpha_stars=alpha_stars, eta_probs=eta_probs,
                       preliminary=FALSE, prior_info=FALSE, prior_Gamma_i=NULL, first_iter=TRUE, previous_Nt1=0)
    tmp_RESULT3[[ww]] <- tmp_result
    tmp_RESULT4[[ww]] <- tmp_resample
  }
  tmp_aPMLs2    <- tmp_RESULT3 %>% sapply(function(x){ x[[1]] })
  tmp_max_idx2  <- tmp_aPMLs2 %>% which.max
  tmp_result2   <- tmp_RESULT3[[tmp_max_idx2]]
  tmp_ordering2 <- tmp_RESULT4[[tmp_max_idx2]]
  tmp_aPML2     <- tmp_aPMLs2 %>% max
  tmp_gamma3    <- tmp_result2[[3]][tmp_ordering2 %>% order]
  tmp_adjRI3    <- adjustedRandIndex(tmp_gamma3, True_Cluster)
  tmp_adjRI3s   <- seq(tmp_RESULT3) %>% sapply(function(x){ tmp_RESULT3[[x]][[3]][tmp_RESULT4[[x]] %>% order] %>% adjustedRandIndex(True_Cluster) })
  
  tmp_SUGS0 <- list(tmp_gamma3, tmp_adjRI3, tmp_aPML2, tmp_aPMLs2, tmp_adjRI3s)
  # -------------------------------------------------------------------------------------------------
  
  ## SUGS algorithm ---------------------------------------------------------------------------------
  ### modeling setting
  c <- 1
  d <- 10
  
  ### BSUGS algorithm
  SUGS_ordering_num2 <- 50
  
  tmp_RESULT3 <- tmp_RESULT4 <- list()
  for(ww in 1:SUGS_ordering_num2){
    
    # ordering
    dqset.seed(ww  + w*103)
    tmp_resample <- dqsample(1:n)
    tmp_Y        <- Yscale[tmp_resample,]
    
    # SUGS algorithm
    set.seed(ww  + w*103)
    tmp_result <- SUGS(Y=tmp_Y, a=a, b=c(1,1), m=m, psi=psi, c=c, d=d, alpha_stars=alpha_stars, eta_probs=eta_probs,
                       preliminary=FALSE, prior_info=FALSE, prior_Gamma_i=NULL, first_iter=TRUE, previous_Nt1=0)
    tmp_RESULT3[[ww]] <- tmp_result
    tmp_RESULT4[[ww]] <- tmp_resample
  }
  tmp_aPMLs2    <- tmp_RESULT3 %>% sapply(function(x){ x[[1]] })
  tmp_max_idx2  <- tmp_aPMLs2 %>% which.max
  tmp_result2   <- tmp_RESULT3[[tmp_max_idx2]]
  tmp_ordering2 <- tmp_RESULT4[[tmp_max_idx2]]
  tmp_aPML2     <- tmp_aPMLs2 %>% max
  tmp_gamma3    <- tmp_result2[[3]][tmp_ordering2 %>% order]
  tmp_adjRI3    <- adjustedRandIndex(tmp_gamma3, True_Cluster)
  tmp_adjRI3s   <- seq(tmp_RESULT3) %>% sapply(function(x){ tmp_RESULT3[[x]][[3]][tmp_RESULT4[[x]] %>% order] %>% adjustedRandIndex(True_Cluster) })
  
  tmp_SUGS <- list(tmp_gamma3, tmp_adjRI3, tmp_aPML2, tmp_aPMLs2, tmp_adjRI3s)
  # -------------------------------------------------------------------------------------------------
  
  ## SUGS algorithm ---------------------------------------------------------------------------------
  ### modeling setting
  c <- 1
  d <- 10
  
  ### BSUGS algorithm
  SUGS_ordering_num2 <- 50
  
  tmp_RESULT3 <- tmp_RESULT4 <- list()
  for(ww in 1:SUGS_ordering_num2){
    
    # ordering
    dqset.seed(ww  + w*103)
    tmp_resample <- dqsample(1:n)
    tmp_Y        <- Yscale[tmp_resample,]
    
    # SUGS algorithm
    set.seed(ww  + w*103)
    tmp_result <- SUGS(Y=tmp_Y, a=a, b=c(1/2,1/2), m=m, psi=psi, c=c, d=d, alpha_stars=alpha_stars, eta_probs=eta_probs,
                       preliminary=FALSE, prior_info=FALSE, prior_Gamma_i=NULL, first_iter=TRUE, previous_Nt1=0)
    tmp_RESULT3[[ww]] <- tmp_result
    tmp_RESULT4[[ww]] <- tmp_resample
  }
  tmp_aPMLs2    <- tmp_RESULT3 %>% sapply(function(x){ x[[1]] })
  tmp_max_idx2  <- tmp_aPMLs2 %>% which.max
  tmp_result2   <- tmp_RESULT3[[tmp_max_idx2]]
  tmp_ordering2 <- tmp_RESULT4[[tmp_max_idx2]]
  tmp_aPML2     <- tmp_aPMLs2 %>% max
  tmp_gamma3    <- tmp_result2[[3]][tmp_ordering2 %>% order]
  tmp_adjRI3    <- adjustedRandIndex(tmp_gamma3, True_Cluster)
  tmp_adjRI3s   <- seq(tmp_RESULT3) %>% sapply(function(x){ tmp_RESULT3[[x]][[3]][tmp_RESULT4[[x]] %>% order] %>% adjustedRandIndex(True_Cluster) })
  
  tmp_SUGS2 <- list(tmp_gamma3, tmp_adjRI3, tmp_aPML2, tmp_aPMLs2, tmp_adjRI3s)
  # -------------------------------------------------------------------------------------------------
  
  
  ## SUGS algorithm ---------------------------------------------------------------------------------
  ### modeling setting
  c <- 1
  d <- 10
  
  ### BSUGS algorithm
  SUGS_ordering_num2 <- 50
  
  tmp_RESULT3 <- tmp_RESULT4 <- list()
  for(ww in 1:SUGS_ordering_num2){
    
    # ordering
    dqset.seed(ww  + w*103)
    tmp_resample <- dqsample(1:n)
    tmp_Y        <- Yscale[tmp_resample,]
    
    # SUGS algorithm
    set.seed(ww  + w*103)
    tmp_result <- SUGS(Y=tmp_Y, a=a, b=c(1/5,1/5), m=m, psi=psi, c=c, d=d, alpha_stars=alpha_stars, eta_probs=eta_probs,
                       preliminary=FALSE, prior_info=FALSE, prior_Gamma_i=NULL, first_iter=TRUE, previous_Nt1=0)
    tmp_RESULT3[[ww]] <- tmp_result
    tmp_RESULT4[[ww]] <- tmp_resample
  }
  tmp_aPMLs2    <- tmp_RESULT3 %>% sapply(function(x){ x[[1]] })
  tmp_max_idx2  <- tmp_aPMLs2 %>% which.max
  tmp_result2   <- tmp_RESULT3[[tmp_max_idx2]]
  tmp_ordering2 <- tmp_RESULT4[[tmp_max_idx2]]
  tmp_aPML2     <- tmp_aPMLs2 %>% max
  tmp_gamma3    <- tmp_result2[[3]][tmp_ordering2 %>% order]
  tmp_adjRI3    <- adjustedRandIndex(tmp_gamma3, True_Cluster)
  tmp_adjRI3s   <- seq(tmp_RESULT3) %>% sapply(function(x){ tmp_RESULT3[[x]][[3]][tmp_RESULT4[[x]] %>% order] %>% adjustedRandIndex(True_Cluster) })
  
  tmp_SUGS3 <- list(tmp_gamma3, tmp_adjRI3, tmp_aPML2, tmp_aPMLs2, tmp_adjRI3s)
  # -------------------------------------------------------------------------------------------------
  
  
  ## SUGS algorithm ---------------------------------------------------------------------------------
  ### modeling setting
  c <- 1
  d <- 10
  
  ### BSUGS algorithm
  SUGS_ordering_num2 <- 50
  
  tmp_RESULT3 <- tmp_RESULT4 <- list()
  for(ww in 1:SUGS_ordering_num2){
    
    # ordering
    dqset.seed(ww  + w*103)
    tmp_resample <- dqsample(1:n)
    tmp_Y        <- Yscale[tmp_resample,]
    
    # SUGS algorithm
    set.seed(ww  + w*103)
    tmp_result <- SUGS(Y=tmp_Y, a=a, b=c(1/10,1/10), m=m, psi=psi, c=c, d=d, alpha_stars=alpha_stars, eta_probs=eta_probs,
                       preliminary=FALSE, prior_info=FALSE, prior_Gamma_i=NULL, first_iter=TRUE, previous_Nt1=0)
    tmp_RESULT3[[ww]] <- tmp_result
    tmp_RESULT4[[ww]] <- tmp_resample
  }
  tmp_aPMLs2    <- tmp_RESULT3 %>% sapply(function(x){ x[[1]] })
  tmp_max_idx2  <- tmp_aPMLs2 %>% which.max
  tmp_result2   <- tmp_RESULT3[[tmp_max_idx2]]
  tmp_ordering2 <- tmp_RESULT4[[tmp_max_idx2]]
  tmp_aPML2     <- tmp_aPMLs2 %>% max
  tmp_gamma3    <- tmp_result2[[3]][tmp_ordering2 %>% order]
  tmp_adjRI3    <- adjustedRandIndex(tmp_gamma3, True_Cluster)
  tmp_adjRI3s   <- seq(tmp_RESULT3) %>% sapply(function(x){ tmp_RESULT3[[x]][[3]][tmp_RESULT4[[x]] %>% order] %>% adjustedRandIndex(True_Cluster) })
  
  tmp_SUGS4 <- list(tmp_gamma3, tmp_adjRI3, tmp_aPML2, tmp_aPMLs2, tmp_adjRI3s)
  # -------------------------------------------------------------------------------------------------
  
  Result_1iter <- list(tmp_BSUGS, tmp_BSUGS2, tmp_BSUGS3, tmp_BSUGS4, tmp_BSUGS0,
                       tmp_SUGS0, tmp_SUGS,   tmp_SUGS2,  tmp_SUGS3,  tmp_SUGS4)
  Result_1iter
}
stopCluster(myCluster)
# save(output, file="simulation6.RData")
# -------------------------------------------------------------------------------------------------
#
###################################################################################################
# vec1 <- sapply(1:SimN, function(x){ output[[x]][[1]][[2]] })
# vec2 <- sapply(1:SimN, function(x){ output[[x]][[2]][[2]] })
# vec3 <- sapply(1:SimN, function(x){ output[[x]][[3]][[2]] })
# vec4 <- sapply(1:SimN, function(x){ output[[x]][[4]][[2]] })
# vec5 <- sapply(1:SimN, function(x){ output[[x]][[5]][[2]] })
# vec6 <- sapply(1:SimN, function(x){ output[[x]][[6]][[2]] })
# vec7 <- sapply(1:SimN, function(x){ output[[x]][[7]][[2]] })
# vec8 <- sapply(1:SimN, function(x){ output[[x]][[8]][[2]] })
# vec9 <- sapply(1:SimN, function(x){ output[[x]][[9]][[2]] })
# vec10 <- sapply(1:SimN, function(x){ output[[x]][[10]][[2]] })
# boxplot(vec3,vec5,vec8,vec9)
# boxplot(vec1,vec2,vec3,vec4,vec6,vec7,vec8,vec9,vec10, ylim=c(0.85,0.99),
#         names=c("BSUGS\nb=1","BSUGS\nb=1/2","BSUGS\nb=1/5","BSUGS\nb=1/10","Empirical\nSUGS","SUGS\nb=1","SUGS\nb=1/2","SUGS\nb=1/5","SUGS\nb=1/10"))




###################################################################################################
# Figures 5, 6
# -------------------------------------------------------------------------------------------------

# Figure 5 ----------------------------------------------------------------------------------------
# example of clustering structure
idx <- 1
cl0 <- True_Cluster
cl1 <- output[[idx]][[4]][[1]]
cl2 <- output[[idx]][[10]][[1]]

cl00 <- cl0
cl0[cl00 == 1] <- 2
cl0[cl00 == 2] <- 1

# to df 
tmp_data <- data_gen(idx)
Raw_Data <- tmp_data[[1]]

df <- as.data.frame(Raw_Data)
colnames(df)[1:2] <- c("x", "y")
df$cl0 <- cl0
df$cl1 <- cl1  
df$cl2 <- cl2 

# color
colors_general <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#17becf",
                    "#9467bd", "#8c564b", "#e377c2", "#7f7f7f",
                    "#bcbd22", "#d62728")

# raw data
unique_cl0 <- sort(unique(df$cl0))
num_cl0 <- length(unique_cl0)
colors_cl0 <- colors_general[1:num_cl0]
names(colors_cl0) <- unique_cl0  

p0 <- ggplot(df, aes(x = x, y = y, fill = factor(cl0))) +
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.7) +
  scale_fill_manual(values = colors_cl0, name = "Cluster") +
  coord_cartesian(xlim = range(df$x), ylim = range(df$y)) +
  labs(x = "X", y = "Y", title = "(a)") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text  = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)
  ) +
  guides(fill = guide_legend(nrow = 1, 
                             override.aes = list(shape = 21, color = NA, stroke = 0)))

# result 1
unique_cl1 <- sort(unique(df$cl1))
num_cl1 <- length(unique_cl1)
colors_cl1 <- colors_general[1:num_cl1]
names(colors_cl1) <- unique_cl1  

p1 <- ggplot(df, aes(x = x, y = y, fill = factor(cl1))) +
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.7) +
  scale_fill_manual(values = colors_cl1, name = "Cluster") +
  coord_cartesian(xlim = range(df$x), ylim = range(df$y)) +
  labs(x = "X", y = "Y", title = "(b)") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text  = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)
  ) +
  guides(fill = guide_legend(nrow = 1, 
                             override.aes = list(shape = 21, color = NA, stroke = 0)))

# result 2
unique_cl2 <- sort(unique(df$cl2))
num_cl2 <- length(unique_cl2)
colors_cl2 <- colors_general[1:num_cl2]  
names(colors_cl2) <- unique_cl2

p2 <- ggplot(df, aes(x = x, y = y, fill = factor(cl2))) +
  geom_point(size = 4, shape = 21, color = "black", stroke = 0.7) +
  scale_fill_manual(values = colors_cl2, name = "Cluster") +
  coord_cartesian(xlim = range(df$x), ylim = range(df$y)) +
  labs(x = "X", y = "Y", title = "(c)") +
  theme_minimal() +
  theme(
    axis.title = element_text(size = 24, face = "bold"),
    axis.text  = element_text(size = 20),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text  = element_text(size = 18),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.title = element_text(size = 30, face = "bold", hjust = 0.5)
  ) +
  guides(fill = guide_legend(nrow = 1, 
                             override.aes = list(shape = 21, color = NA, stroke = 0)))

# 1 x 2 grid
combined_plot <- p0 + p1 + p2 + plot_layout(ncol = 3)

# print
print(combined_plot)

# # save
# ggsave(filename = "figure5.png", plot = combined_plot, dpi = 1000, width = 20, height = 8, bg = "white")
# # -------------------------------------------------------------------------------------------------



# Figure 6 ----------------------------------------------------------------------------------------
# data
df <- tibble(
  BSUGS_b1        = vec1,
  BSUGS_b05       = vec2,
  BSUGS_b015      = vec3,
  BSUGS_b010      = vec4,
  Empirical_SUGS  = vec6,
  SUGS_b1         = vec7,
  SUGS_b05        = vec8,
  SUGS_b015       = vec9,
  SUGS_b010       = vec10
) %>%
  pivot_longer(
    cols = everything(),
    names_to  = "Method",
    values_to = "Value"
  )

# label
df$Method <- factor(
  df$Method,
  levels = c("BSUGS_b1","BSUGS_b05","BSUGS_b015","BSUGS_b010",
             "Empirical_SUGS","SUGS_b1","SUGS_b05","SUGS_b015","SUGS_b010"),
  labels = c("BSUGS\nb=1","BSUGS\nb=1/2","BSUGS\nb=1/5","BSUGS\nb=1/10",
             "Empirical\nSUGS","SUGS\nb=1","SUGS\nb=1/2","SUGS\nb=1/5","SUGS\nb=1/10")
)

# BSUGS vs SUGS
df$Block <- ifelse(as.integer(df$Method) <= 4, "Type1", "Type2")

# color
cols <- c(Type1 = "#1B9E77", Type2 = "#D95F02")

# plot
p <- ggplot(df, aes(x = Method, y = Value, fill = Block)) +
  geom_boxplot(width = 0.8, outlier.size = 1, fatten = 1) +
  coord_cartesian(ylim = c(0.85, 0.99)) +
  scale_fill_manual(values = cols) +
  scale_x_discrete(expand = expansion(add = 0.5)) +
  labs(x = NULL, y = "Adj. Rand index") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 16, lineheight = .9),
    axis.text.y = element_text(size = 16),
    axis.title.y       = element_text(size = 24)  )

print(p)

# # save
# ggsave("figure6.png", plot = p, width = 10, height = 6, dpi = 1000, bg = "white")
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
#
###################################################################################################
