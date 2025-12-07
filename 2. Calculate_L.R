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
library(gridExtra)
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

## load the Phase I result 
load("Ph1_8.RData")
# -------------------------------------------------------------------------------------------------
#
###################################################################################################



###################################################################################################
# BSUGS algorithm
# -------------------------------------------------------------------------------------------------
# basic setting
ARL0_simN      <- 1000
ARL0           <- 100
Ph2_n          <- 3000
Phase1_ARL0_SS <- 200

# Phase I result
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

# priors
a   <- 1
b   <- 1
m   <- 0
psi <- 1
alpha_stars <- c(0.01, 0.05, c(0.1 + 0.2 * (0:19)))
eta_probs   <- rep(1/length(alpha_stars), length(alpha_stars))
TT          <- length(alpha_stars)

# setting for SUGS algorithm
SUGS_iter_num <- 5
bidirec   <- T
back_uNum <- NULL
Onum      <- 1000
ITERmax   <- 5
RePUpd    <- T

# calculate hat pi for Phase I
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

myCluster <- makeCluster(30)
registerDoParallel(myCluster)
output <- foreach(www = 1:ARL0_simN, .packages = c("dplyr", "dqrng")) %dopar% {
  
  # Phase I parameters
  tmp___Gamma_i  <- Ph1___Gamma_i
  tmp___Gamma_i2 <- Ph1___Gamma_i2
  tmp___phi_t    <- Ph1___phi_t
  tmp___psi_hj   <- Ph1___psi_hj
  tmp___m_hj     <- Ph1___m_hj
  tmp___a_hj     <- Ph1___a_hj
  tmp___b_hj     <- Ph1___b_hj
  
  # calculate L for ARL0
  ARL0_result <- list()
  for(w in 1:Phase1_ARL0_SS){
    
    # DATA for obtaining ARL0 = 100
    dqset.seed(w + (www+(rep_idx-1)*ARL0_simN)*1009)
    tmp_arl0_idx    <- dqsample(1:ncol(DATA_mat_org0), Ph2_n, replace=T) %>% sort
    tmp_arl0_DATA0  <- DATA_mat_org0[,tmp_arl0_idx]
    
    set.seed(w + (www+(rep_idx-1)*ARL0_simN)*1309)
    tmp_arl0_DATA    <- data_transform_cont(tmp_arl0_DATA0)
    
    tmp_org_gamma_mat1 <- tmp_org_gamma_mat2 <- NULL
    for(ix in 1:Ph1_DataN){
      
      tmp_i_idx <- 1:each_n + (ix-1) * each_n
      
      tmp_gam1 <- Ph1___Gamma_i[tmp_i_idx]
      tmp_gam2 <- Ph1___Gamma_i2[tmp_i_idx]
      
      tmp_org1 <- seq_along(tmp_arl0_idx) %>% sapply( function(x){tmp_gam1[which(Ph1_ordering[[ix]] == tmp_arl0_idx[x])]} )
      tmp_org2 <- seq_along(tmp_arl0_idx) %>% sapply( function(x){tmp_gam2[which(Ph1_ordering[[ix]] == tmp_arl0_idx[x])]} )
      
      tmp_org_gamma_mat1 <- rbind(tmp_org_gamma_mat1,tmp_org1)
      tmp_org_gamma_mat2 <- rbind(tmp_org_gamma_mat2,tmp_org2)
    }
    
    tmp_org_gamma_i1  <- tmp_org_gamma_mat1 %>% apply(2, function(x){ x %>% table %>% which.max %>% names %>% as.numeric })
    tmp_org_gamma_i2 <- tmp_org_gamma_mat2 %>% apply(2, function(x){ x %>% table %>% which.max %>% names %>% as.numeric })
    
    # Zernike polynomial coefficients
    tmp_arl0_coef <- get_Zernike(tmp_arl0_DATA, Bmat)
    
    # scaling
    tmp_arl0_Yscale0 <- (tmp_arl0_coef - matrix(Step1_scale_center, nrow=Ph2_n, ncol=K, byrow=T)) / matrix(Step1_scale_scale,  nrow=Ph2_n, ncol=K, byrow=T)
    
    tmp_RESULT <- tmp_RESULT2 <- list()
    for(ww in 1:SUGS_iter_num){
      
      # ordering
      dqset.seed(ww + w*1000 + (www+(rep_idx-1)*ARL0_simN)*2009)
      tmp_resample    <- dqsample(1:nrow(tmp_arl0_Yscale0))
      tmp_arl0_Yscale <- tmp_arl0_Yscale0[tmp_resample,]
      
      # SUGS algorithm
      set.seed(ww + w*1000 + (www+(rep_idx-1)*ARL0_simN)*2901)
      tmp_result   <- BSUGS(Y=tmp_arl0_Yscale, a=a, b=b, m=m, psi=psi, alpha_stars=alpha_stars, eta_probs=eta_probs,
                            bidirectional = bidirec, back_upperNum = back_uNum,
                            OnceNum       = Onum,    IterMax       = ITERmax,   ReParmUpd   = RePUpd,
                            prior_info=FALSE, prior_Gamma_i=NULL,   Phase1=F,
                            prev_Gamma_i = tmp___Gamma_i, prev_phi_t=tmp___phi_t, prev_psi_hj=tmp___psi_hj,
                            prev_m_hj    =tmp___m_hj,     prev_a_hj =tmp___a_hj,  prev_b_hj  =tmp___b_hj)
      tmp_RESULT[[ww]]  <- tmp_result
      tmp_RESULT2[[ww]] <- tmp_resample
    }
    tmp_aPMLs    <- tmp_RESULT %>% sapply(function(x){ x[[1]] })
    tmp_max_idx  <- tmp_aPMLs %>% which.max
    tmp_result   <- tmp_RESULT[[tmp_max_idx]]
    tmp_ordering <- tmp_RESULT2[[tmp_max_idx]]
    tmp_aPML     <- tmp_aPMLs %>% max
    
    # update clustering result
    tmp___Gamma_i      <- c(tmp___Gamma_i,  tmp_result[[2]])
    tmp___Gamma_i2     <- c(tmp___Gamma_i2, tmp_result[[12]])
    tmp___org_Gamma_i  <- tmp_org_gamma_i1[tmp_ordering]
    tmp___org_Gamma_i2 <- tmp_org_gamma_i2[tmp_ordering]
    
    # update T2 value
    tmp_nu_star   <- 2*(tmp___a_hj %>% sapply(function(x){ x[2,1] }))
    tmp_t_fun <- function(k){
      
      tmp_nk    <- sum(tmp_result[[2]] == k, na.rm=T)
      
      if(tmp_nk > 0){
        
        tmp_denom <- sqrt(Ph1___b_hj[[k]][2,]/Ph1___a_hj[[k]][2,] * (1/tmp_nk + tmp___psi_hj[[k]][2,]))
        tmp_ymean <- tmp_arl0_Yscale0[tmp_ordering[tmp_result[[2]] == k],,drop=FALSE] %>% apply(2, function(x){ mean(x, na.rm=T) })
        tmp_nom   <- (tmp_ymean - Ph1___m_hj[[k]][2,])
        
        tmp_tvals    <- tmp_nom/tmp_denom
        
      }else{ tmp_tvals <- NA }
      
      return(tmp_tvals)
    }
    tmp_tval <- sapply(1:Ph1_H, function(k){ tmp_t_fun(k) }, simplify=F)
    
    # update parameters
    tmp___phi_t        <- tmp_result[[3]]
    tmp___psi_hj       <- tmp_result[[4]]
    tmp___m_hj         <- tmp_result[[5]]
    tmp___a_hj         <- tmp_result[[6]]
    tmp___b_hj         <- tmp_result[[7]]
    
    ARL0_result[[w]] <- list(tmp_result[[2]], tmp_result[[12]], 
                             Ph2_n - sum(is.na(tmp_result[[2]])),
                             tmp___org_Gamma_i, tmp___org_Gamma_i2,
                             mean(tmp_result[[2]] == tmp___org_Gamma_i,  na.rm=T),
                             mean(tmp_result[[12]] == tmp___org_Gamma_i2, na.rm=T),
                             tmp_aPML, tmp___phi_t, tmp___psi_hj, tmp___m_hj, tmp___a_hj, tmp___b_hj,
                             tmp_tval, tmp_nu_star)
  }
  
  ARL0_result
}
stopCluster(myCluster)
# save(output, file=paste0("ARL0.RData"))
# -------------------------------------------------------------------------------------------------
#
###################################################################################################


###################################################################################################
# Calculate alpha s.t. ARL0 = 100
# -------------------------------------------------------------------------------------------------
# load("ARL0.RData")
target_Ph1_pis   <- Ph1_pi_h2
lower_prob_thres <- 0.005
Int2_idx <- which(target_Ph1_pis >= lower_prob_thres)
Int3_idx <- which(target_Ph1_pis <  lower_prob_thres) %>% intersect(1:Ph1_H)

Lnum <- 200
alpha_list1 <- seq(0.0001, 0.001, length.out=Lnum)
alpha_list2 <- seq(0.005,  0.01, length.out=Lnum)

# large proportion
ARL0_alpha1 <- NULL
if(length(Int2_idx) > 0){
  
  for(j in 1:Lnum){
    
    cat(j,"  ")
    
    # interval types
    alpha          <- alpha_list1[j]
    Interval_type2 <- Interval3(target_Ph1_pis, Ph2_n, alpha)
    UCL_type2      <- Interval_type2[[1]]
    LCL_type2      <- Interval_type2[[2]]
    
    # calculate RLs
    tmp_funx2 <- function(xx){
      
      tmp_props         <- cl_props[[xx]][,Int2_idx,drop=F]
      tmp_result_mat2_0 <- tmp_props %>%
        apply(1, function(x){ ((x > UCL_type2[Int2_idx]) | (x < LCL_type2[Int2_idx])) %>% any })
      tmp_RL2_0         <- ifelse(sum(tmp_result_mat2_0) == 0, Phase1_ARL0_SS, tmp_result_mat2_0 %>% which %>% .[1])
      
      return(tmp_RL2_0)
    }
    
    tmp_RL_j     <- seq_along(cl_props) %>% sapply(tmp_funx2)
    tmp_val      <- mean(tmp_RL_j)
    ARL0_alpha1  <- c(ARL0_alpha1, tmp_val)
  }
}

# very small proprotion
ARL0_alpha2 <- NULL
if(length(Int3_idx) > 0){
  
  for(j in 1:Lnum){
    
    cat(j,"  ")
    
    alpha          <- alpha_list2[j]
    Interval_type3 <- Interval3(target_Ph1_pis, Ph2_n, alpha)
    UCL_type3      <- Interval_type3[[1]]
    LCL_type3      <- Interval_type3[[2]]
    
    # calculate RLs
    tmp_funx3 <- function(xx){
      
      tmp_props         <- cl_props[[xx]][,Int3_idx,drop=F]
      tmp_result_mat3_0 <- tmp_props %>%
        apply(1, function(x){ (x > UCL_type3[Int3_idx]) %>% any })
      tmp_RL3_0         <- ifelse(sum(tmp_result_mat3_0) == 0, Phase1_ARL0_SS, tmp_result_mat3_0 %>% which %>% .[1])
      
      return(tmp_RL3_0)
    }
    
    tmp_RL_j    <- seq_along(cl_props) %>% sapply(tmp_funx3)
    tmp_val     <- mean(tmp_RL_j)
    ARL0_alpha2 <- c(ARL0_alpha2, tmp_val)
  }
}

alpha_idx1 <- which(ARL0_alpha1 >= ARL0) %>% tail(1)
alpha_idx2 <- which(ARL0_alpha2 >= ARL0) %>% tail(1)

ALPHA1 <- alpha_list1[alpha_idx1] # 0.0006110553
ALPHA2 <- alpha_list2[alpha_idx2] # 0.00781407
# -------------------------------------------------------------------------------------------------
#
####################################################################################################


###################################################################################################
# figure 11
# -------------------------------------------------------------------------------------------------
x           <- 2
ARL0_result <- output[[x]]

# interval types
alpha1 <- ALPHA1
Interval_type2 <- Interval3(Ph1_pi_h, Ph2_n, alpha1)
UCL_type2      <- Interval_type2[[1]]
LCL_type2      <- Interval_type2[[2]]

alpha2 <- ALPHA2
Interval_type3 <- Interval3(Ph1_pi_h, Ph2_n, alpha2)
UCL_type3      <- Interval_type3[[1]]
LCL_type3      <- Interval_type3[[2]]

# cluster reordering
tmp_order <- c(1,2,20,6,17,14,24,5,15,18,4,11,3,27,25,8,
               16,7,13,29,23,19,22,12,21,26,9,10,28)
chk_clusters <- NULL 

# # of plots
n_plots <- length(tmp_order)

# title for each subplots
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

# list for plots for each cluster
plot_list <- list()

for(i in seq_along(tmp_order)) {
  k <- tmp_order[i]
  
  # plotting statistics
  target_ps <- sapply(ARL0_result, function(x) { 
    sum(x[[2]] == k, na.rm = TRUE) / sum(!is.na(x[[2]]))
  })
  
  # centerline
  target_CL <- target_Ph1_pis[k]
  
  # color
  int2_yn <- k %in% Int2_idx
  ctrl_col <- ifelse(int2_yn, "#E41A1C", "#4DAF4A")
  
  if(int2_yn){
    target_UCL <- UCL_type2[k]
    target_LCL <- LCL_type2[k]
  } else {
    target_UCL <- UCL_type3[k]
    target_LCL <- NA
  }
  
  # ylim
  ylim_val <- range(c(target_ps, target_UCL, target_LCL), na.rm = TRUE)
  
  # data frame
  df <- data.frame(Time = seq_along(target_ps), Proportion = target_ps)
  
  # ggplot
  p <- ggplot(df, aes(x = Time, y = Proportion)) +
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
      axis.text = element_text(size = 10)
    )
  
  # LCL
  if(!is.na(target_LCL)) {
    p <- p + geom_hline(yintercept = target_LCL, linetype = "dashed",
                        color = ctrl_col, size = 1)
  }
  
  plot_list[[i]] <- p
}

# 4 x 8 grid
combined_plot <- do.call(grid.arrange, c(plot_list, nrow = 4, ncol = 8))

# # save
# ggsave("figure11.png", combined_plot, width = 20, height = 10, dpi = 1000, bg = "white")
# -----------------------------------------------------------------------------------------------
#
####################################################################################################



###################################################################################################
# Bayesian t-chart & figure 12
# -------------------------------------------------------------------------------------------------
## plot
k <- 2          # which cluster???
i <- 100
ARL1_result <- output[[i]]
alpha_t     <- 0.01

tmp_tvalues <- ARL1_result %>% sapply(function(x){ if(is.na(x[[14]][[k]][1])) rep(NA, K) else x[[14]][[k]] })
tmp_dfs     <- cbind(Ph1___a_hj[[k]][2,],
                     ARL1_result %>% sapply(function(x){ x[[15]][[k]] %>% rep(K) }))

## plot
k <- 2          # which cluster???
i <- 100
ARL1_result <- output[[i]]
alpha_t     <- 0.01

tmp_tvalues <- ARL1_result %>% sapply(function(x){ if(is.na(x[[14]][[k]][1])) rep(NA, K) else x[[14]][[k]] })
tmp_dfs     <- cbind(Ph1___a_hj[[k]][2,],
                     ARL1_result %>% sapply(function(x){ x[[15]][[k]] %>% rep(K) }))
par(mfrow=c(3,5))
for(j in 1:K){
  
  tmp_x   <- tmp_tvalues[j,]
  tmp_ucl <- qt(1-alpha_t/(2*Ph1_H*K), df=tmp_dfs[j,] %>% head(-1))
  tmp_lcl <- qt(alpha_t/(2*Ph1_H*K), df=tmp_dfs[j,] %>% head(-1))
  
  ylim <- c(tmp_x,tmp_ucl,tmp_lcl) %>% quantile(c(0,1), na.rm=T)
  plot(tmp_x, type='b', ylim=ylim, xlab='t', ylab='')
  lines(tmp_ucl, col=2, lwd=2)
  lines(tmp_lcl, col=2, lwd=2)
}


# data
tmp_tvalues <- ARL1_result %>% sapply(function(x) {
  if (is.na(x[[14]][[k]][1])) rep(NA, K) else x[[14]][[k]]
})
tmp_dfs <- cbind(Ph1___a_hj[[k]][2,],
                 ARL1_result %>% sapply(function(x) { rep(x[[15]][[k]], K) }))

# list for each subplot
plot_list <- list()

# title for each subplot
if(K <= 26){
  sub_titles <- paste0("(", letters[1:K], ")")
} else {
  sub_titles <- c(paste0("(", letters[1:26], ")"),
                  sapply(27:K, function(i) {
                    prefix <- letters[floor((i - 1) / 26)]
                    suffix <- letters[((i - 1) %% 26) + 1]
                    paste0("(", prefix, suffix, ")")
                  }))
}

# for each subplot
for(j in 1:K) {
  
  # data
  tmp_x   <- tmp_tvalues[j, ]            
  tmp_ucl <- qt(1 - alpha_t/(2 * Ph1_H * K), df = head(tmp_dfs[j, ], -1))
  tmp_lcl <- qt(alpha_t/(2 * Ph1_H * K), df = head(tmp_dfs[j, ], -1))
  
  # y limit
  ylim_val <- range(c(tmp_x, tmp_ucl, tmp_lcl), na.rm = TRUE)
  
  # x axis
  t_vals <- seq_along(tmp_x)
  
  # df
  df_est   <- data.frame(t = t_vals, value = tmp_x)
  df_ucl   <- data.frame(t = t_vals, value = tmp_ucl)
  df_lcl   <- data.frame(t = t_vals, value = tmp_lcl)
  
  # color
  main_col <- "#1B9E77"
  bound_col <- "#D95F02"
  
  # ggplot2
  p <- ggplot() +
    geom_line(data = df_est, aes(x = t, y = value), color = main_col, size = 1) +
    geom_point(data = df_est, aes(x = t, y = value), color = main_col, size = 2) +
    geom_line(data = df_ucl, aes(x = t, y = value), color = bound_col, linetype = "dashed", size = 1) +
    geom_line(data = df_lcl, aes(x = t, y = value), color = bound_col, linetype = "dashed", size = 1) +
    labs(title = sub_titles[j], x = "t", y = "") +
    scale_x_continuous(breaks = pretty(t_vals)) +
    scale_y_continuous(limits = ylim_val) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  plot_list[[j]] <- p
}

# 3 x 5 grid
combined_plot <- do.call(grid.arrange, c(plot_list, nrow = 3, ncol = 5))

# # save
# ggsave("figure12.png", combined_plot, width = 20, height = 10, dpi = 1000, bg = "white")
# --------------------------------------------------------------------------------------------------
#

####################################################################################################
