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
# (Step 1) Zernike polynomial
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
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
#
###################################################################################################



###################################################################################################
# (Step 2-0) Setting for SUGS algorithm
# -------------------------------------------------------------------------------------------------
## data
total_Y <- NULL
for(i in 1:Ph1_DataN){ total_Y <- rbind(total_Y, Ph1_DATA[[i]]) }

## scaling
Step1_scale_center <- total_Y %>% scale %>% attr("scaled:center")
Step1_scale_scale  <- total_Y %>% scale %>% attr("scaled:scale")
total_Yscale       <- total_Y %>% scale 

## basic setting
each_n  <- nrow(Ph1_DATA[[1]])
K       <- ncol(Ph1_DATA[[1]])
M       <- 100

## priors
a   <- 1
b   <- 1
m   <- 0
psi <- 1
alpha_stars <- c(0.01, 0.05, c(0.1 + 0.2 * (0:19)))
eta_probs   <- rep(1/length(alpha_stars), length(alpha_stars))
TT          <- length(alpha_stars)
# -------------------------------------------------------------------------------------------------
#
###################################################################################################



###################################################################################################
# (Step 2) 2. SUGS algorithm
# -------------------------------------------------------------------------------------------------

## initial value for SUGS algorithm
tmp___Gamma_i  <- NULL
tmp___Gamma_i2 <- NULL
tmp___phi_t    <- NULL
tmp___psi_hj   <- NULL
tmp___m_hj     <- NULL
tmp___a_hj     <- NULL
tmp___b_hj     <- NULL

## BSUGS algorithm
Total_Result <- Total_Result_ordering <- list()
SUGS_ordering_num_list <- c(200, rep(15, Ph1_DataN-1))
for(i in 1:Ph1_DataN){

  cat(i,"th dataset!!!\n")

  ### data at i
  tmp_i_idx <- 1:each_n + (i-1) * each_n
  Yscale    <- total_Yscale[tmp_i_idx,]

  ### setting for BSUGS
  SUGS_ordering_num <- SUGS_ordering_num_list[i]
  num_cl            <- ifelse(i==1, 30, SUGS_ordering_num)

  ### modeling setting
  bidirec   <- T
  back_uNum <- NULL
  Onum      <- 1000
  ITERmax   <- ifelse(i==1, max(round(nrow(Yscale)/Onum) * 2, 40), 20)
  RePUpd    <- T

  ### SUGS algorithm - parallel
  myCluster <- makeCluster(num_cl)
  registerDoParallel(myCluster)
  output <- foreach(ww = 1:SUGS_ordering_num, .packages = c("dplyr", "dqrng")) %dopar% {

    cat(ww,"th iteration !!! \n")
    dqset.seed(ww*23 + i*10019)

    #### ordering
    tmp_resample <- dqsample(1:nrow(Yscale))
    tmp_Y        <- Yscale[tmp_resample,]

    #### SUGS algorithm
    tmp_result   <- BSUGS(Y=tmp_Y, a=a, b=b, m=m, psi=psi, alpha_stars=alpha_stars, eta_probs=eta_probs,
                          bidirectional = bidirec, back_upperNum = back_uNum,
                          OnceNum       = Onum,    IterMax       = ITERmax,   ReParmUpd   = RePUpd,
                          prior_info=FALSE, prior_Gamma_i=NULL, Phase1=ifelse(i==1,T,F),
                          prev_Gamma_i=tmp___Gamma_i, prev_phi_t=tmp___phi_t,
                          prev_psi_hj =tmp___psi_hj,  prev_m_hj =tmp___m_hj,
                          prev_a_hj   =tmp___a_hj,    prev_b_hj =tmp___b_hj)

    Result_1iter <- list(Result = tmp_result, Ordering = tmp_resample)
    Result_1iter
  }
  stopCluster(myCluster)

  ### select
  target_output    <- output
  max_idx          <- sapply(target_output, function(x){ x[[1]][[1]] }) %>% which.max
  RESULT           <- target_output[[max_idx]]
  tmp__ordering    <- RESULT[[2]]
  tmp__result      <- RESULT[[1]]

  ## initial value for SUGS algorithm
  tmp___Gamma_i  <- c(tmp___Gamma_i,  tmp__result[[2]])
  tmp___Gamma_i2 <- c(tmp___Gamma_i2, tmp__result[[12]])
  tmp___phi_t    <- tmp__result[[3]]
  tmp___psi_hj   <- tmp__result[[4]]
  tmp___m_hj     <- tmp__result[[5]]
  tmp___a_hj     <- tmp__result[[6]]
  tmp___b_hj     <- tmp__result[[7]]

  ### save
  Total_Result[[i]]          <- tmp__result
  Total_Result_ordering[[i]] <- tmp__ordering
}

## save the result
Final_result___Gamma_i  <- tmp___Gamma_i
Final_result___Gamma_i2 <- tmp___Gamma_i2
Final_result___phi_t    <- tmp___phi_t
Final_result___psi_hj   <- tmp___psi_hj
Final_result___m_hj     <- tmp___m_hj
Final_result___a_hj     <- tmp___a_hj
Final_result___b_hj     <- tmp___b_hj
Gamma_i                 <- Final_result___Gamma_i
Gamma_i2                <- Final_result___Gamma_i2

## check cluster proportion
Gamma_i2_list <- Gamma_i2_list0 <- NULL
for(i in 1:Ph1_DataN){
  
  tmp_i_idx      <- 1:each_n + (i-1) * each_n
  tmp_i_Gamma_i2 <- Gamma_i2[tmp_i_idx]
  Gamma_i2_list  <- rbind(Gamma_i2_list, tmp_i_Gamma_i2)
  
  tmp_i_Gamma_i2 <- tmp_i_Gamma_i2[Total_Result_ordering[[i]] %>% order]
  Gamma_i2_list0 <- rbind(Gamma_i2_list0, tmp_i_Gamma_i2)
}

total_props      <- table(Gamma_i2)/(each_n*Ph1_DataN)
total_categories <- Gamma_i2 %>% unique %>% sort
each_props       <- sapply(1:Ph1_DataN, function(i){ table(factor(Gamma_i2_list[i,], levels=total_categories))/each_n }) %>% t

## plot - check whether cluster estimation is stable?
par(mfrow=c(5,10))
for(k in total_categories){
  
  total_mean <- total_props[k]
  each_means <-each_props[,k]
  
  plot(each_means, ylim=quantile(c(total_mean, each_means), c(0,1)), ylab="prop", xlab="dataset", main=paste("Cluster",k))
  abline(h=total_mean, lwd=2, col=4)
}
# -------------------------------------------------------------------------------------------------

# Repetition ------------------------------------------------------------------------------------
RepNum <- 2

tmp___Gamma_i              <- Gamma_i
Int1_Total_Result          <- Total_Result
Int1_Total_Result_ordering <- Total_Result_ordering
for(repn in 1:RepNum){
  
  ## initial value for SUGS algorithm
  tmp___Gamma_i2 <- NULL
  tmp___phi_t    <- Int1_Total_Result[[Ph1_DataN]][[3]]
  tmp___psi_hj   <- Int1_Total_Result[[Ph1_DataN]][[4]]
  tmp___m_hj     <- Int1_Total_Result[[Ph1_DataN]][[5]]
  tmp___a_hj     <- Int1_Total_Result[[Ph1_DataN]][[6]]
  tmp___b_hj     <- Int1_Total_Result[[Ph1_DataN]][[7]]
  
  ## BSUGS algorithm
  Int2_Total_Result      <- Int2_Total_Result_ordering <- list()
  SUGS_ordering_num_list <- rep(10, Ph1_DataN)
  for(i in 1:Ph1_DataN){
    
    cat(i,"th dataset!!!\n")
    
    ### data at i
    tmp_i_idx <- 1:each_n + (i-1) * each_n
    Yscale    <- total_Yscale[tmp_i_idx,]
    tmp_Y     <- Yscale[Int1_Total_Result_ordering[[i]],]
    
    ### backward elimination
    tmp_result   <- Repeat_Backward(Y=tmp_Y, a=a, b=b, m=m, psi=psi, alpha_stars=alpha_stars, eta_probs=eta_probs, target_i_idx=tmp_i_idx, 
                                    prev_Gamma_i=tmp___Gamma_i, prev_phi_t=tmp___phi_t, prev_psi_hj=tmp___psi_hj, 
                                    prev_m_hj   =tmp___m_hj,    prev_a_hj =tmp___a_hj,  prev_b_hj  =tmp___b_hj)
    tmp___Gamma_i[tmp_i_idx] <- NA
    tmp___phi_t    <- tmp_result[[1]]
    tmp___psi_hj   <- tmp_result[[2]]
    tmp___m_hj     <- tmp_result[[3]]
    tmp___a_hj     <- tmp_result[[4]]
    tmp___b_hj     <- tmp_result[[5]]
    
    ### setting for BSUGS
    SUGS_ordering_num <- SUGS_ordering_num_list[i]
    num_cl            <- SUGS_ordering_num
    
    ### modeling setting
    bidirec   <- T
    back_uNum <- NULL
    Onum      <- 1000
    ITERmax   <- 20
    RePUpd    <- T
    
    ### SUGS algorithm - parallel
    myCluster <- makeCluster(num_cl)
    registerDoParallel(myCluster)
    output <- foreach(ww = 1:SUGS_ordering_num, .packages = c("dplyr", "dqrng")) %dopar% {
      
      cat(ww,"th iteration !!! \n")
      dqset.seed(ww*17 + i*10000 + repn * 50)
      
      #### ordering
      tmp_resample <- dqsample(1:nrow(Yscale))
      tmp_Y        <- Yscale[tmp_resample,]
      
      #### SUGS algorithm
      tmp_result   <- BSUGS(Y=tmp_Y, a=a, b=b, m=m, psi=psi, alpha_stars=alpha_stars, eta_probs=eta_probs,
                            bidirectional = bidirec, back_upperNum = back_uNum,
                            OnceNum       = Onum,    IterMax       = ITERmax,   ReParmUpd   = RePUpd,
                            prior_info=FALSE, prior_Gamma_i=NULL,   Phase1=F,
                            prev_Gamma_i=tmp___Gamma_i[-tmp_i_idx], prev_phi_t=tmp___phi_t,
                            prev_psi_hj =tmp___psi_hj,              prev_m_hj =tmp___m_hj,
                            prev_a_hj   =tmp___a_hj,                prev_b_hj =tmp___b_hj)
      
      Result_1iter <- list(Result = tmp_result, Ordering = tmp_resample)
      Result_1iter
    }
    stopCluster(myCluster)
    
    ### select
    target_output    <- output
    max_idx          <- sapply(target_output, function(x){ x[[1]][[1]] }) %>% which.max
    RESULT           <- target_output[[max_idx]]
    tmp__ordering    <- RESULT[[2]]
    tmp__result      <- RESULT[[1]]
    
    ## initial value for SUGS algorithm
    tmp___Gamma_i[tmp_i_idx]  <- tmp__result[[2]]
    tmp___Gamma_i2            <- c(tmp___Gamma_i2, tmp__result[[12]])
    tmp___phi_t    <- tmp__result[[3]]
    tmp___psi_hj   <- tmp__result[[4]]
    tmp___m_hj     <- tmp__result[[5]]
    tmp___a_hj     <- tmp__result[[6]]
    tmp___b_hj     <- tmp__result[[7]]
    
    ### save
    Int2_Total_Result[[i]]          <- tmp__result
    Int2_Total_Result_ordering[[i]] <- tmp__ordering
  }
  
  ## save the result
  Int1_Total_Result          <- Int2_Total_Result
  Int1_Total_Result_ordering <- Int2_Total_Result_ordering
}

### Final result
Final_Total_Result          <- Int1_Total_Result
Final_Total_Result_ordering <- Int1_Total_Result_ordering

### For summary
Final_result___Gamma_i  <- tmp___Gamma_i
Final_result___Gamma_i2 <- tmp___Gamma_i2
Final_result___phi_t    <- tmp___phi_t
Final_result___psi_hj   <- tmp___psi_hj
Final_result___m_hj     <- tmp___m_hj
Final_result___a_hj     <- tmp___a_hj
Final_result___b_hj     <- tmp___b_hj
Gamma_i                 <- Final_result___Gamma_i
Gamma_i2                <- Final_result___Gamma_i2

Ph1_result <- list(Final_Total_Result, Final_Total_Result_ordering, Gamma_i, Gamma_i2)
save(Ph1_result, file="Ph1_8.RData")
# -------------------------------------------------------------------------------------------------
#
###################################################################################################


###################################################################################################
# Check Step 2
# -------------------------------------------------------------------------------------------------
# load("Ph1_8.RData")
# Final_Total_Result          <- Ph1_result[[1]]
# Final_Total_Result_ordering <- Ph1_result[[2]]
# Gamma_i  <- Ph1_result[[3]]
# Gamma_i2 <- Ph1_result[[4]]

## check cluster proportion
Gamma_i2_list <- Gamma_i2_list0 <- NULL
for(i in 1:Ph1_DataN){
  
  tmp_i_idx      <- 1:each_n + (i-1) * each_n
  tmp_i_Gamma_i2 <- Gamma_i2[tmp_i_idx]
  Gamma_i2_list  <- rbind(Gamma_i2_list, tmp_i_Gamma_i2)
  
  tmp_i_Gamma_i2 <- tmp_i_Gamma_i2[Final_Total_Result_ordering[[i]] %>% order]
  Gamma_i2_list0 <- rbind(Gamma_i2_list0, tmp_i_Gamma_i2)
}

total_props      <- table(Gamma_i2)/(each_n*Ph1_DataN)
total_categories <- Gamma_i2 %>% unique %>% sort
each_props       <- sapply(1:Ph1_DataN, function(i){ table(factor(Gamma_i2_list[i,], levels=total_categories))/each_n }) %>% t

## plot - check whether cluster estimation is stable?
par(mfrow=c(5,10))
for(k in total_categories){
  
  total_mean <- total_props[k]
  each_means <- each_props[,k]
  
  plot(each_means, ylim=quantile(c(total_mean, each_means), c(0,1)), ylab="prop", xlab="dataset", main=paste("Cluster",k))
  abline(h=total_mean, lwd=2, col=4)
}
# -------------------------------------------------------------------------------------------------
#
###################################################################################################


###################################################################################################
# (Step 2) 3. Posterior inference
# -------------------------------------------------------------------------------------------------
#######################################################
target_i     <- Ph1_DataN                         #####
target_i_idx <- 1:each_n + (target_i-1) * each_n  #####
#######################################################
# MC - cluster specific mean surface
tmp_mean_function_list <- list()
for(iii in 1:Ph1_DataN){
  
  I <- 10000
  target_parameters <- list( sapply(Final_Total_Result[[iii]][[5]],   function(x){ x %>% tail(1) %>% c }),
                             sapply(Final_Total_Result[[iii]][[4]], function(x){ x %>% tail(1) %>% c }),
                             sapply(Final_Total_Result[[iii]][[6]],   function(x){ x %>% tail(1) %>% c }),
                             sapply(Final_Total_Result[[iii]][[7]],   function(x){ x %>% tail(1) %>% c }) )
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
  tmp_mean_function_list[[iii]] <- tmp_mean_function
}

zlim  <- tmp_mean_function_list %>% sapply(function(x){ quantile(x,c(0,1), na.rm=T) }) %>% apply(1, mean) * c(0.8, 1.2)
zlim2 <- tmp_DATA_mat %>% unlist %>% quantile(c(0,1), na.rm=T)
# -------------------------------------------------------------------------------------------------

# cluster specific surface ------------------------------------------------------------------------
gam_unq      <- Gamma_i2_list[target_i,] %>% na.omit %>% unique %>% sort
tmp_mean_mat <- sapply(tmp_clusters, function(i){ tmp_mean_function_list[[target_i]][i,] %>% matrix(nrow=dat_dim[1], ncol=dat_dim[2]) }, simplify=F)

par(mfrow=c(4,8))
for(i in seq(tmp_mean_mat)){
  
  tmp_gam <- gam_unq[i]
  
  txt_title <- paste0("Cl ", i, " (n=", sum(Gamma_i2_list[target_i,] == tmp_gam, na.rm=T), ")")
  image.plot(x_grid, y_grid, t(tmp_mean_mat[[i]]), xlab="x", ylab="y", zlim=zlim)
}
# -------------------------------------------------------------------------------------------------

# plot 5-------------------------------------------------------------------------------------------
target_cluster <- 29
target_idx     <- Final_Total_Result_ordering[[target_i]][Gamma_i2_list[target_i,] == target_cluster]

tmp_seq <- 1:length(target_idx)
if(length(target_idx) > 60){ tmp_seq <- 1:60 }

par(mfrow=c(6,10))
for(i in tmp_seq){ 

  tmp_i                          <- target_idx[i]
  tmp_data_full                  <- rep(NA, prod(dat_dim))
  tmp_data_full[-WF_outise_idx]  <- Ph1_DATA0[[target_i]][,tmp_i]
  tmp_data_mat                   <- matrix(tmp_data_full, nrow=dat_dim[1], ncol=dat_dim[2]) %>% t

  image.plot(x_grid, y_grid, tmp_data_mat, xlab="x", ylab="y", zlim=zlim2, 
             main=paste(tmp_i, "(Cluster",Gamma_i2_list[target_i,][Final_Total_Result_ordering[[target_i]] == tmp_i],")"))
}
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
#
###################################################################################################






###################################################################################################
# plots - figure 1, 8, 9, 10
# -------------------------------------------------------------------------------------------------

# figure 1 ----------------------------------------------------------------------------------------
tmp_idx <- c(2,1,15,8,18,7)

txt_titles <- c("(a)","(b)","(c)","(d)","(e)","(f)")
tmp_mean_mat <- sapply(tmp_clusters, function(i){ tmp_mean_function_list[[target_i]][i,] %>% matrix(nrow=dat_dim[1], ncol=dat_dim[2]) }, simplify=F)

# surface to data frame
matrix_to_df <- function(x_grid, y_grid, mat, panel_label) {
  df <- expand.grid(x = x_grid, y = y_grid)
  df$z <- as.vector(t(mat))
  df$panel <- panel_label
  return(df)
}

# surface to data frame with tmp_idx
combined_df <- do.call(rbind, lapply(seq_along(tmp_idx), function(j) {
  i <- tmp_idx[j]
  matrix_to_df(x_grid, y_grid, tmp_mean_mat[[i]], txt_titles[j])
}))

# ggplot2
p <- ggplot(combined_df, aes(x = x, y = y, fill = z)) +
  geom_tile() +  
  scale_fill_viridis_c(option = "H", limits = c(1.9,6.1), name = "z") +
  facet_wrap(~ panel, nrow = 2) +  
  labs(x = "x", y = "y") +
  theme_minimal() +
  theme(
    strip.text = element_text(color = "black", face = "bold", size = 24),
    panel.grid = element_blank()
  )
print(p)

# # save
# ggsave("figure1.png", p, width = 10, height = 7, dpi = 1000, bg = "white")
# -------------------------------------------------------------------------------------------------

# figure 8 ----------------------------------------------------------------------------------------
# title
titles <- c("(a)", "(b)", "(c)", "(d)", "(e)",
            "(f)", "(g)", "(h)", "(i)", "(j)")
zlim  <- c(-1, 9)
zlim2 <- c(1, 3)

target_i  <- 50
target_is <- c(11857, 7107, 8870, 6966, 9923) # 5202 12165 11857


# list for plot
plot_list <- list()

# first loop
for(idx in seq_along(target_is)) {
  i <- target_is[idx]
  
  # data
  tmp_org <- DATA_mat_org0[, i]
  tmp_data_full0 <- rep(NA, prod(dat_dim))
  tmp_data_full0[-WF_outise_idx] <- tmp_org
  tmp_data_mat0 <- matrix(tmp_data_full0, nrow = dat_dim[1], ncol = dat_dim[2])
  
  # to df
  df0 <- melt(tmp_data_mat0)
  df0$x <- x_grid[df0$Var2]
  df0$y <- y_grid[df0$Var1]
  
  # ggplot
  p <- ggplot(df0, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis(name = "Value", limits = zlim2, option = "H", na.value = "white") +
    labs(title = titles[idx], x = "x", y = "y") +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),  
          axis.title = element_text(size = 18),                              
          axis.text = element_text(size = 12))                              
  
  plot_list[[idx]] <- p
}

# second loop
for(idx in seq_along(target_is)) {
  i <- target_is[idx]
  
  # data
  tmp_con <- Ph1_DATA0[[target_i]][, i]
  tmp_data_full <- rep(NA, prod(dat_dim))
  tmp_data_full[-WF_outise_idx] <- tmp_con
  tmp_data_mat <- matrix(tmp_data_full, nrow = dat_dim[1], ncol = dat_dim[2])
  
  # to df
  df <- melt(tmp_data_mat)
  df$x <- x_grid[df$Var2]
  df$y <- y_grid[df$Var1]
  
  # ggplot
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis(name = "Value", limits = zlim, option = "H", na.value = "white") +
    labs(title = titles[length(target_is) + idx], x = "x", y = "y") +
    theme_minimal() +
    theme(aspect.ratio = 1,
          plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 12))
  
  plot_list[[length(target_is) + idx]] <- p
}

# 2 x 5 grid
combined_plot <- arrangeGrob(grobs = plot_list, nrow = 2, ncol = 5)

# print
grid.arrange(combined_plot)

# # save
# ggsave("figure8.png", combined_plot, width = 15, height = 6, dpi = 1000, bg = "white")
# -------------------------------------------------------------------------------------------------


# figure 9 ----------------------------------------------------------------------------------------
WF_chipN       <- prod(dat_dim) - length(WF_outise_idx)
tmp_data_add   <- SC1_RawData(10, X_data_polar, WF_chipN)
tmp_data_add2  <- SC2_RawData(10, X_data_polar, WF_chipN)
tmp_arl1_DATA  <- SC12_data_transform_cont(tmp_data_add)
tmp_arl1_DATA2 <- SC12_data_transform_cont(tmp_data_add2)

sc1_org <- tmp_data_add[,2]
sc1_con <- tmp_arl1_DATA[,2]

sc2_org <- tmp_data_add2[,2]
sc2_con <- tmp_arl1_DATA2[,2]

# plot setting
titles <- c("(a)", "(b)", "(c)", "(d)")
zlim  <- c(-1.5, 9.5)
zlim2 <- c(1, 3)

# list for ggplot object
plot_list <- list()

## (a)
# data
tmp_data_full0 <- rep(NA, prod(dat_dim))
tmp_data_full0[-WF_outise_idx] <- sc1_org
tmp_data_mat0  <- matrix(tmp_data_full0, nrow = dat_dim[1], ncol = dat_dim[2])

# to df
df1 <- melt(tmp_data_mat0)
df1$x <- x_grid[df1$Var2]
df1$y <- y_grid[df1$Var1]

# ggplot
p1 <- ggplot(df1, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis(name = "Value", limits = zlim2, option = "H", na.value = "white") +
  labs(title = titles[1], x = "x", y = "y") +
  theme_minimal() +
  theme(aspect.ratio = 1,
        plot.title  = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title  = element_text(size = 14),
        axis.text   = element_text(size = 12))

plot_list[[1]] <- p1

## (b) 
tmp_data_full0 <- rep(NA, prod(dat_dim))
tmp_data_full0[-WF_outise_idx] <- sc1_con
tmp_data_mat0  <- matrix(tmp_data_full0, nrow = dat_dim[1], ncol = dat_dim[2])

df2 <- melt(tmp_data_mat0)
df2$x <- x_grid[df2$Var2]
df2$y <- y_grid[df2$Var1]

p2 <- ggplot(df2, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis(name = "Value", limits = zlim, option = "H", na.value = "white") +
  labs(title = titles[2], x = "x", y = "y") +
  theme_minimal() +
  theme(aspect.ratio = 1,
        plot.title  = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title  = element_text(size = 14),
        axis.text   = element_text(size = 12))

plot_list[[2]] <- p2

## (c) 
tmp_data_full0 <- rep(NA, prod(dat_dim))
tmp_data_full0[-WF_outise_idx] <- sc2_org
tmp_data_mat0  <- matrix(tmp_data_full0, nrow = dat_dim[1], ncol = dat_dim[2])

df3 <- melt(tmp_data_mat0)
df3$x <- x_grid[df3$Var2]
df3$y <- y_grid[df3$Var1]

p3 <- ggplot(df3, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis(name = "Value", limits = zlim2, option = "H", na.value = "white") +
  labs(title = titles[3], x = "x", y = "y") +
  theme_minimal() +
  theme(aspect.ratio = 1,
        plot.title  = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title  = element_text(size = 14),
        axis.text   = element_text(size = 12))

plot_list[[3]] <- p3

## (d) 
tmp_data_full0 <- rep(NA, prod(dat_dim))
tmp_data_full0[-WF_outise_idx] <- sc2_con
tmp_data_mat0  <- matrix(tmp_data_full0, nrow = dat_dim[1], ncol = dat_dim[2])

df4 <- melt(tmp_data_mat0)
df4$x <- x_grid[df4$Var2]
df4$y <- y_grid[df4$Var1]

p4 <- ggplot(df4, aes(x = x, y = y, fill = value)) +
  geom_tile() +
  scale_fill_viridis(name = "Value", limits = zlim, option = "H", na.value = "white") +
  labs(title = titles[4], x = "x", y = "y") +
  theme_minimal() +
  theme(aspect.ratio = 1,
        plot.title  = element_text(hjust = 0.5, size = 24, face = "bold"),
        axis.title  = element_text(size = 14),
        axis.text   = element_text(size = 12))

plot_list[[4]] <- p4

# 1 x 4 grid
combined_plot <- arrangeGrob(grobs = plot_list, nrow = 1, ncol = 4)

# print
grid.arrange(combined_plot)

# # save
# ggsave("figure9.png", combined_plot, width = 15, height = 4, dpi = 1000, bg = "white")
# -------------------------------------------------------------------------------------------------

# figure 10 ---------------------------------------------------------------------------------------
gam_unq      <- Gamma_i2_list[target_i,] %>% na.omit %>% unique %>% sort
tmp_mean_mat <- sapply(tmp_clusters, function(i){ tmp_mean_function_list[[target_i]][i,] %>% matrix(nrow=dat_dim[1], ncol=dat_dim[2]) %>% t }, simplify=F)

tmp_order         <- Gamma_i2_list %>% c %>% table %>% prop.table %>% order(decreasing = T)
tmp_order[c(1,2)] <- c(1,2)

tmp_prop <- Gamma_i2_list %>% c %>% table %>% prop.table %>% sort(decreasing = T)
zlim     <- c(1.9, 7.1)

# data reordering
tmp_mean_mat <- tmp_mean_mat[tmp_order]
n_plots <- length(tmp_mean_mat)

# titles for each subplot
if(n_plots <= 26){
  plot_labels <- paste0("(", letters[1:n_plots], ")")
} else {
  plot_labels <- c(paste0("(", letters[1:26], ")"),
                   sapply(27:n_plots, function(i) {
                     prefix <- letters[floor((i-1)/26)]
                     suffix <- letters[((i-1) %% 26) + 1]
                     paste0("(", prefix, suffix, ")")
                   }))
}

# plot
plot_list <- list()
for(i in seq_along(tmp_mean_mat)) {
  
  # to df 
  mat <- tmp_mean_mat[[i]] %>% t
  df <- melt(mat)
  
  # x and y
  df$x <- x_grid[df$Var2]
  df$y <- y_grid[df$Var1]
  
  # ggplot2
  p <- ggplot(df, aes(x = x, y = y, fill = value)) +
    geom_tile() +
    scale_fill_viridis_c(
      option = "H",
      limits = zlim,
      na.value = "white",
      guide = guide_colorbar(
        title = "Value",
        barwidth = 0.5,  
        barheight = 3   
      )
    ) +
    labs(title = plot_labels[i], x = "x", y = "y") +
    theme_minimal() +
    theme(
      aspect.ratio = 1,
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "right",
      legend.key.size = unit(0.6, "cm"),     
      legend.title = element_text(size = 8), 
      legend.text = element_text(size = 6)   
    )
  
  plot_list[[i]] <- p
}

# 4 x 8 grid
combined_plot <- arrangeGrob(grobs = plot_list, nrow = 4, ncol = 8)

# print
grid.arrange(combined_plot)

# # save
# ggsave("figure10.png", combined_plot, width = 20, height = 10, dpi = 1000, bg = "white")
# -------------------------------------------------------------------------------------------------

# -------------------------------------------------------------------------------------------------
#
###################################################################################################



















