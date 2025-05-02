###################################################################################################
# Interval
# -------------------------------------------------------------------------------------------------
Interval1 <- function(hat_pi, Nt, L){
  
  tmp_sd  <- sqrt(hat_pi * (1-hat_pi) / Nt)
  tmp_UCL <- hat_pi + L * tmp_sd
  tmp_LCL <- hat_pi - L * tmp_sd
  tmp_val <- list(tmp_UCL, tmp_LCL)
  
  return(tmp_val)
}


Interval2 <- function(hat_pi, Nt, L){
  
  bar_Nt  <- Nt + L^2
  bar_pi  <- (Nt * hat_pi + L^2/2)/bar_Nt
  
  tmp_sd  <- sqrt(bar_pi * (1-bar_pi) / bar_Nt)
  tmp_UCL <- bar_pi + L * tmp_sd
  tmp_LCL <- bar_pi - L * tmp_sd
  tmp_val <- list(tmp_UCL, tmp_LCL)
  
  return(tmp_val)
}

Interval3 <- function(hat_pi, Nt, alpha){
  
  tmp_parm1 <- Nt * hat_pi + 1/2
  tmp_parm2 <- Nt - Nt * hat_pi + 1/2
  
  tmp_UCL <- qbeta(1-alpha/2, tmp_parm1, tmp_parm2)
  tmp_LCL <- qbeta(alpha/2,   tmp_parm1, tmp_parm2)
  tmp_val <- list(tmp_UCL, tmp_LCL)
  
  return(tmp_val)
}
# -------------------------------------------------------------------------------------------------
#
###################################################################################################


###################################################################################################
# efficient calculation for table
# -------------------------------------------------------------------------------------------------
table_cum <- function(org_table, add_table){
  
  if(length(add_table) > 0){
    
    # all categories for two tables
    all_levels <- union(names(org_table), names(add_table))
    
    # all category name
    org_table_aligned <- setNames(rep(0, length(all_levels)), all_levels)
    add_table_aligned <- setNames(rep(0, length(all_levels)), all_levels)
    
    # fill table
    org_table_aligned[names(org_table)] <- org_table
    add_table_aligned[names(add_table)] <- add_table
    
    # combine two table
    org_table <- org_table_aligned + add_table_aligned
  }
  
  return(org_table)
}
# -------------------------------------------------------------------------------------------------
#
###################################################################################################



###################################################################################################
# Data transformation
# -------------------------------------------------------------------------------------------------

# WBM to WF map pattern ---------------------------------------------------------------------------
data_transform_cont <- function(Data){
  
  n1 <- sum(Data == 1)
  n2 <- sum(Data == 2)
  
  Data[Data == 1] <- rnorm(n1, mean=2, sd=0.8)
  Data[Data == 2] <- rnorm(n2, mean=4, sd=0.8)
  
  return(Data)
}

SC1_RawData <- function(Num, polar_coord, Cnum){
  
  raw_data <- NULL
  for(i in 1:Num){
    
    SC1_F_idx           <- which(polar_coord[,1] > 0.9)
    test_mat            <- rep(1,Cnum)
    test_mat[SC1_F_idx] <- 3
    SC1_R_idx           <- which(test_mat != 3) %>% dqsample(dqsample(10,1))
    test_mat[SC1_R_idx] <- 2
    
    raw_data <- cbind(raw_data,test_mat)
  }
  if(Num == 1){ raw_data <- c(raw_data) }
  
  return(raw_data)
}

SC2_RawData <- function(Num, x_coord, Cnum){
  
  raw_data <- NULL
  for(i in 1:Num){
    
    SC2_F_idx           <- (t(x_coord) - c(0.35,-0.40)) %>% apply(2, function(x){ sum(x^2) < 0.1 }) %>% which
    test_mat            <- rep(1,Cnum)
    test_mat[SC2_F_idx] <- 2
    SC2_R_idx           <- which(test_mat != 2) %>% dqsample(dqsample(10,1))
    test_mat[SC2_R_idx] <- 2
    
    raw_data <- cbind(raw_data,test_mat)
  }
  if(Num == 1){ raw_data <- c(raw_data) }
  
  return(raw_data)
}


SC12_data_transform_cont <- function(Data){
  
  n1 <- sum(Data == 1)
  n2 <- sum(Data == 2)
  n3 <- sum(Data == 3)
  
  Data[Data == 1] <- rnorm(n1, mean=2, sd=0.8)
  Data[Data == 2] <- rnorm(n2, mean=4, sd=0.8)
  Data[Data == 3] <- rnorm(n3, mean=6, sd=0.8)
  
  return(Data)
}
# -------------------------------------------------------------------------------------------------

# zernike polynomial ------------------------------------------------------------------------------
zernike_radial <- function(n, m, r) {
  
  ## (n-m) should be even
  if( (n - m) %% 2 != 0 ) { return(rep(0, length(r))) }
  
  ## radial polynomial
  val <- rep(0, length(r))
  nm2 <- (n - m) / 2
  for(k in 0:floor(nm2)){
    cfac <- ((-1)^k) * factorial(n-k) / ( factorial(k) * factorial((n+m)/2 - k) * factorial((n-m)/2 - k) )
    val  <- val + cfac * r^(n - 2*k)
  }
  return(val)
}

# Zernike (real)
zernike_real <- function(n, m, r, theta) {
  
  mm <- abs(m)
  Rnm <- zernike_radial(n, mm, r)
  
  # normalizing factor
  norm_factor <- 0
  if(m == 0) {
    norm_factor <- sqrt(n+1)
  }else{
    norm_factor <- sqrt(2*(n+1))
  }
  
  # angle part
  if(m == 0) {
    # no angle variation
    return( norm_factor * Rnm )
  } else if(m > 0) {
    # cos(m * theta)
    return( norm_factor * Rnm * cos(m*theta) )
  } else {
    # m<0 => sin(|m| * theta)
    return( norm_factor * Rnm * sin(mm*theta) )
  }
}

zernike_basis_set <- function(nMax=4) {
  combos <- list()
  idx <- 1
  for(n in 0:nMax){
    for(m in (-n):n){
      if( (n - abs(m))%%2 == 0 ) {
        combos[[idx]] <- c(n,m)
        idx<-idx+1
      }
    }
  }
  return(combos)
}

zernike_eval_one <- function(n, m, r, theta){ zernike_real(n,m,r,theta) }

get_Zernike <- function(Data, BMAT){
  
  tmp_n    <- Data %>% ncol
  tmp_coef <- NULL
  for(i in 1:tmp_n){
    
    tmp_data      <- Data[,i]
    alpha_hat     <- solve(t(BMAT)%*%BMAT, t(BMAT)%*%tmp_data)
    hat_parameter <- c(alpha_hat)
    
    tmp_coef      <- rbind(tmp_coef, hat_parameter)
  }
  
  return(tmp_coef)
}
# -------------------------------------------------------------------------------------------------
#
###################################################################################################

###################################################################################################
# Bayesian approach based on SUGS algorithm
# -------------------------------------------------------------------------------------------------
cal_entropy <- function(xvec){
  
  tmp_prob <- xvec/sum(xvec)
  tmp_val  <- -sum(tmp_prob * log(tmp_prob), na.rm=T)
  return(tmp_val)
}

Noncentral_t <- function(y, others){
  
  nu <- others[[1]]
  mu <- others[[2]]
  sig <- others[[3]]
  
  tmp_val1 <- gamma(0.5)/beta(nu/2, 0.5) / (((pi*nu)^(1/2)) * sig)
  tmp_val2 <- (1 + 1/((sig^2)*nu) * ((y-mu)^2))^(-(nu+1)/2)
  result <- tmp_val1 * tmp_val2
  return(result)
}

Calculate_Noncentral_t <- function(y, ai1, bi1, mi1, psii1, w=1){
  
  tmp_nu   <- 2*ai1
  tmp_mu   <- mi1
  tmp_sig  <- sqrt( bi1/ai1 * (w+psii1) )
  tmp_prms <- list(tmp_nu, tmp_mu, tmp_sig)
  
  result   <- Noncentral_t(y=y, others=tmp_prms)
  return(result)
}

BSUGS <- function(Y, a=1, b=1, m=0, psi=1, psi0=100, alpha_stars, eta_probs=NULL,
                  bidirectional = T, back_upperNum = NULL, OnceNum = 1000, IterMax=NULL, ReParmUpd = TRUE, cprob=NULL,
                  prior_info=FALSE, prior_Gamma_i=NULL, Phase1=TRUE,
                  prev_Gamma_i=NULL, prev_phi_t=NULL, prev_psi_hj=NULL, 
                  prev_m_hj=NULL, prev_a_hj=NULL, prev_b_hj=NULL)
{
  # hyperparameters
  if(is.null(eta_probs)){ eta_probs  <- rep(1/length(alpha_stars), length(alpha_stars))}
  if(is.null(IterMax)){   IterMax    <- round(nrow(Y)/OnceNum) * 3 }
  if(is.null(cprob)){       cprob    <- -1 }
  
  # basic setting
  n      <- nrow(Y)
  K      <- ncol(Y)
  TT     <- length(alpha_stars)
  prev_n <- sum(!is.na(prev_Gamma_i))
  
  # 0th value
  zero_psi_hj  <- rep(psi, K)
  zero_psi0_hj <- rep(psi0, K)
  zero_m_hj    <- rep(m, K)
  zero_a_hj    <- rep(a, K)
  zero_b_hj    <- rep(b, K)
  
  # parameters
  if(Phase1 == TRUE){
    
    prev_Gamma_table <- NULL
    Gamma_i          <- rep(NA, n)
    phi_t            <- NULL
    psi_hj           <- m_hj <- a_hj <- b_hj <- list()
    
  }else{
    
    prev_Gamma_table <- prev_Gamma_i %>% na.omit %>% table
    Gamma_i          <- rep(NA, n)
    phi_t            <- prev_phi_t
    psi_hj           <- prev_psi_hj
    m_hj             <- prev_m_hj
    a_hj             <- prev_a_hj
    b_hj             <- prev_b_hj
  }
  
  # algorithm
  for(Iter in 1:IterMax){
    
    last_iter <- Iter == IterMax
    iter_idx  <- Gamma_i %>% is.na %>% which
    if((length(iter_idx)>1) & (Iter>1)){ iter_idx <- dqsample(iter_idx) }
    certain_p <- cprob
    
    ## sequential update
    if(length(iter_idx) > 0){
      
      ### only OnceNum subjects
      tmp_OnceNum <- ifelse(Iter == 1, max(OnceNum, length(prior_Gamma_i)), OnceNum)
      if(!last_iter){ iter_idx  <- iter_idx[1:min(tmp_OnceNum, length(iter_idx))] }
      
      for(i in iter_idx){
        
        # k_i1 <- c(prev_Gamma_i, Gamma_i) %>% na.omit %>% unique %>% length
        unique_prev <- names(prev_Gamma_table)
        unique_pres <- Gamma_i %>% na.omit %>% unique %>% as.character
        k_i1        <- unique_total_count <- length(union(unique_prev, unique_pres))
        
        #### choose gamma
        if((i==1)&Phase1){ #### first iteraton
          
          new_gamma_i <- 1
          Gamma_i[i]  <- new_gamma_i
          
        }else{             #### not first
          
          tmp_Gamma_table  <- prev_Gamma_table
          PRES_Gamma_table <- Gamma_i %>% na.omit %>% table
          tmp_Gamma_table  <- table_cum(tmp_Gamma_table, PRES_Gamma_table)
          
          old_phi_t <- phi_t[2,]
          pi_iht <- rbind(matrix(tmp_Gamma_table, nrow=k_i1, ncol=TT, byrow=F), alpha_stars) *
            matrix(1/(sum(!is.na(Gamma_i)) + prev_n + alpha_stars), nrow=k_i1+1, ncol=TT, byrow=T)
          
          if(prior_info & (i <= length(prior_Gamma_i))){
            
            new_gamma_i <- prior_Gamma_i[i]
            Gamma_i[i]  <- new_gamma_i
            
          }else{
            
            tmp_probs <- pi_iht * matrix(old_phi_t, nrow=k_i1+1, ncol=TT, byrow=T)
            Lih <- sapply(1:k_i1, function(x){  Calculate_Noncentral_t(y    =Y[i,],
                                                                       ai1  =a_hj[[x]][2,],
                                                                       bi1  =b_hj[[x]][2,],
                                                                       mi1  =m_hj[[x]][2,],
                                                                       psii1=psi_hj[[x]][2,]) }) %>% apply(2, prod)
            Lih <- c(Lih, Calculate_Noncentral_t(y    =Y[i,],
                                                 ai1  =zero_a_hj,
                                                 bi1  =zero_b_hj,
                                                 mi1  =zero_m_hj,
                                                 psii1=zero_psi_hj) %>% prod)
            prob_gammai  <- (tmp_probs*matrix(Lih, nrow=k_i1+1, ncol=TT)) %>% apply(1, sum)
            prob_gammai  <- prob_gammai/sum(prob_gammai)
            max_prob     <- prob_gammai %>% max
            max_prob_idx <- prob_gammai %>% which.max
            new_gamma_i  <- max_prob_idx
            Gamma_i[i]   <- new_gamma_i
          }
        }
        
        if(!is.na(new_gamma_i)){
          
          ### update alpha
          if( (Phase1 & (i<=length(prior_Gamma_i)) & prior_info) | (Phase1 & (i==1) & (!prior_info)) ){ # ((i==1)&Phase1) #### first iteraton
            
            old_phi_t <- eta_probs
            new_phi_t <- eta_probs
            phi_t     <- rbind(old_phi_t, new_phi_t)
            
          }else{    #### not first
            
            old_phi_t <- phi_t[2,]
            new_phi_t <- old_phi_t * pi_iht[new_gamma_i,]
            new_phi_t <- new_phi_t/sum(new_phi_t)
            phi_t     <- rbind(old_phi_t, new_phi_t)
          }
          
          ### update theta=
          if(k_i1>0){ #### not first
            
            for(ii in 1:k_i1){
              
              if(ii == new_gamma_i){
                
                old_psi_hj <- psi_hj[[ii]][2,]
                new_psi_hj <- 1/(1/old_psi_hj + 1)
                
                old_m_hj   <- m_hj[[ii]][2,]
                new_m_hj   <- new_psi_hj * ((1/old_psi_hj)*old_m_hj + Y[i,])
                
                old_a_hj   <- a_hj[[ii]][2,]
                new_a_hj   <- old_a_hj + 1/2
                
                old_b_hj   <- b_hj[[ii]][2,]
                new_b_hj   <- old_b_hj + 1/2 * 1/( old_psi_hj + 1 ) * (Y[i,] - old_m_hj)^2
                
                psi_hj[[new_gamma_i]] <- rbind(old_psi_hj, new_psi_hj)
                m_hj[[new_gamma_i]]   <- rbind(old_m_hj,   new_m_hj)
                a_hj[[new_gamma_i]]   <- rbind(old_a_hj,   new_a_hj)
                b_hj[[new_gamma_i]]   <- rbind(old_b_hj,   new_b_hj)
              }
            }
          }
          
          # new component added or first iteration
          count_in_prev <- ifelse(new_gamma_i %in% names(prev_Gamma_table), prev_Gamma_table[[new_gamma_i]], 0L)
          count_in_pres <- sum(Gamma_i == new_gamma_i, na.rm = TRUE)
          new_tf        <- (count_in_prev + count_in_pres) == 1
          # new_tf0       <- sum(c(prev_Gamma_i, Gamma_i) == new_gamma_i, na.rm=T)==1
          
          if(new_tf)
          {
            if(prior_info & (i <= length(prior_Gamma_i))){
              old_psi_hj <- zero_psi0_hj
            }else{
              old_psi_hj <- zero_psi_hj
            }
            # old_psi_hj <- zero_psi_hj
            new_psi_hj <- 1/(1/old_psi_hj + 1)
            
            old_m_hj   <- zero_m_hj
            new_m_hj   <- new_psi_hj * ((1/old_psi_hj)*old_m_hj + Y[i,])
            
            old_a_hj   <- zero_a_hj
            new_a_hj   <- old_a_hj + 1/2
            
            old_b_hj   <- zero_b_hj
            new_b_hj   <- old_b_hj + 1/2 * 1/(old_psi_hj + 1) * (Y[i,] - old_m_hj)^2
            
            psi_hj[[new_gamma_i]] <- rbind(old_psi_hj, new_psi_hj)
            m_hj[[new_gamma_i]]   <- rbind(old_m_hj,   new_m_hj)
            a_hj[[new_gamma_i]]   <- rbind(old_a_hj,   new_a_hj)
            b_hj[[new_gamma_i]]   <- rbind(old_b_hj,   new_b_hj)
          }
        }
      }
    }
    
    ## backward update
    if(prior_info){ fixed_idx <- 1:length(prior_Gamma_i) }else{ fixed_idx <- NULL }
    reconsider_idx <- Gamma_i %>% is.na %>% `!` %>% which
    reconsider_idx <- setdiff(reconsider_idx,fixed_idx)
    if(!bidirectional){ reconsider_idx <- NULL }
    if(length(reconsider_idx)>1){ reconsider_idx <- dqsample(reconsider_idx) }
    if(!is.null(back_upperNum)){
      if(length(reconsider_idx)>1){ 
        reconsider_idx <- dqsample(reconsider_idx, min(length(reconsider_idx),back_upperNum))
        reconsider_idx <- c(reconsider_idx, iter_idx) %>% unique %>% dqsample
      }
    }
    
    if(length(reconsider_idx) > 0){
      
      delete0_Gamma_i <- delete_Gamma_i <- Gamma_i
      delete0_m_hj    <- delete_m_hj    <- m_hj
      delete0_psi_hj  <- delete_psi_hj  <- psi_hj
      delete0_a_hj    <- delete_a_hj    <- a_hj
      delete0_b_hj    <- delete_b_hj    <- b_hj
      delete0_phi_t   <- delete_phi_t   <- phi_t
      delete0_Gamma_i_table <- delete_Gamma_i_table <-
        table_cum(prev_Gamma_table, delete_Gamma_i %>% na.omit %>% table)
      
      for(i in reconsider_idx){
        
        #### delete i
        delete_gam         <- delete_Gamma_i[i]
        delete0_Gamma_i[i] <- NA
        delete0_Gamma_i_table[delete_gam] <- delete0_Gamma_i_table[delete_gam] - 1
        
        ##### if we delete unique element, we need reindexing !
        delete_unq <- delete_Gamma_i_table[delete_gam] == 1
        if(delete_unq){
          
          ### update parameters
          delete0_psi_hj <- delete0_psi_hj[-delete_gam]
          delete0_m_hj   <- delete0_m_hj[-delete_gam]
          delete0_a_hj   <- delete0_a_hj[-delete_gam]
          delete0_b_hj   <- delete0_b_hj[-delete_gam]
          
          ### delete zero
          delete0_Gamma_i_table_nz <- delete0_Gamma_i_table[delete0_Gamma_i_table != 0]
          remaining_orig           <- as.numeric(names(delete0_Gamma_i_table_nz))  # 예: c(1, 2, 3, 5)
          
          ### update table
          delete0_Gamma_i_table_upd        <- delete0_Gamma_i_table_nz
          names(delete0_Gamma_i_table_upd) <- seq_along(delete0_Gamma_i_table_nz)
          
          ### index mapping
          new_mapping <- setNames(seq_along(remaining_orig), remaining_orig)
          
          ### update gamma
          delete0_Gamma_i       <- as.numeric(new_mapping[as.character(delete0_Gamma_i)])
          delete0_Gamma_i_table <- delete0_Gamma_i_table_upd
          
        }else{
          
          ### update parameters
          new_psi_hj <- delete0_psi_hj[[delete_gam]][2,]
          old_psi_hj <- 1/(1/new_psi_hj - 1)
          delete0_psi_hj[[delete_gam]][1,] <- old_psi_hj
          delete0_psi_hj[[delete_gam]][2,] <- old_psi_hj
          
          new_m_hj   <- delete0_m_hj[[delete_gam]][2,]
          old_m_hj   <- old_psi_hj * (new_m_hj/new_psi_hj - Y[i,])
          delete0_m_hj[[delete_gam]][1,]   <- old_m_hj
          delete0_m_hj[[delete_gam]][2,]   <- old_m_hj
          
          new_a_hj   <- delete0_a_hj[[delete_gam]][2,]
          old_a_hj   <- new_a_hj - 1/2
          delete0_a_hj[[delete_gam]][1,]   <- old_a_hj
          delete0_a_hj[[delete_gam]][2,]   <- old_a_hj
          
          new_b_hj   <- delete0_b_hj[[delete_gam]][2,]
          old_b_hj   <- new_b_hj - 1/2 * (Y[i,]-old_m_hj)^2/(old_psi_hj + 1)
          delete0_b_hj[[delete_gam]][1,]   <- old_b_hj
          delete0_b_hj[[delete_gam]][2,]   <- old_b_hj
        }
        
        ### pi_iht
        pi_iht     <- rbind(matrix(delete0_Gamma_i_table, nrow=length(delete0_Gamma_i_table), ncol=TT, byrow=F), alpha_stars) *
          matrix(1/(sum(!is.na(delete0_Gamma_i)) + prev_n + alpha_stars), nrow=length(delete0_Gamma_i_table)+1, ncol=TT, byrow=T)
        
        ### update phi_t
        if(delete_unq){
          tmp_pi <- pi_iht %>% tail(1) %>% c
        }else{
          tmp_pi <- pi_iht[delete_gam,]
        }
        
        new_phi_t  <- delete0_phi_t[2,]
        old_phi_t  <- new_phi_t/tmp_pi
        old_phi_t  <- old_phi_t/sum(old_phi_t)
        delete0_phi_t[1,]  <- old_phi_t
        delete0_phi_t[2,]  <- old_phi_t 
        
        tmp_probs <- pi_iht * matrix(old_phi_t, nrow=length(delete0_Gamma_i_table)+1, ncol=TT, byrow=T)
        Lih <- sapply(1:length(delete0_Gamma_i_table), function(x){  Calculate_Noncentral_t(y    =Y[i,],
                                                                                            ai1  =delete0_a_hj[[x]][2,],
                                                                                            bi1  =delete0_b_hj[[x]][2,],
                                                                                            mi1  =delete0_m_hj[[x]][2,],
                                                                                            psii1=delete0_psi_hj[[x]][2,]) %>% prod })
        Lih <- c(Lih, Calculate_Noncentral_t(y    =Y[i,],
                                             ai1  =zero_a_hj,
                                             bi1  =zero_b_hj,
                                             mi1  =zero_m_hj,
                                             psii1=zero_psi_hj) %>% prod)
        prob_gammai  <- (tmp_probs*matrix(Lih, nrow=length(delete0_Gamma_i_table)+1, ncol=TT)) %>% apply(1, sum)
        prob_gammai  <- prob_gammai/sum(prob_gammai)
        max_prob     <- prob_gammai %>% max
        max_prob_idx <- prob_gammai %>% which.max
        
        remain_tf    <- ifelse(delete_unq, max_prob_idx == length(delete0_Gamma_i_table)+1, max_prob_idx == delete_gam)
        remain_tf    <- remain_tf & (max_prob > certain_p)
        
        if(remain_tf){
          
          delete0_Gamma_i <- delete_Gamma_i
          delete0_m_hj    <- delete_m_hj
          delete0_psi_hj  <- delete_psi_hj
          delete0_a_hj    <- delete_a_hj
          delete0_b_hj    <- delete_b_hj
          delete0_phi_t   <- delete_phi_t
          delete0_Gamma_i_table <- delete_Gamma_i_table
          
        }else{
          
          delete_Gamma_i <- delete0_Gamma_i
          delete_m_hj    <- delete0_m_hj
          delete_psi_hj  <- delete0_psi_hj
          delete_a_hj    <- delete0_a_hj
          delete_b_hj    <- delete0_b_hj
          delete_phi_t   <- delete0_phi_t
          delete_Gamma_i_table <- delete0_Gamma_i_table
        }
      }
      
      ### update after backward process
      Gamma_i <- delete_Gamma_i
      m_hj    <- delete_m_hj
      psi_hj  <- delete_psi_hj
      a_hj    <- delete_a_hj
      b_hj    <- delete_b_hj
      phi_t   <- delete_phi_t
    }
    
    ### update parameters
    if(ReParmUpd){
      
      tmp_fun0 <- function(k){
        
        k_tf <- k %in% names(prev_Gamma_table)
        # k_tf <- sum(prev_Gamma_i == k,na.rm=T) > 0
        if(k_tf){
          
          tmp_psi <- prev_psi_hj[[k]][2,]
          tmp_m   <- prev_m_hj[[k]][2,]
          tmp_a   <- prev_a_hj[[k]][2,]
          tmp_b   <- prev_b_hj[[k]][2,]
          
        }else{
          
          tmp_psi <- rep(psi0,K)
          tmp_m   <- rep(m,K)
          tmp_a   <- rep(a,K)
          tmp_b   <- rep(b,K)
        }
        
        tmp_k_idx <- Gamma_i == k
        tmp_k_idx[is.na(tmp_k_idx)] <- F
        
        if(sum(tmp_k_idx) > 0){
          
          result_psi   <- 1/(1/tmp_psi + sum(tmp_k_idx))
          result_m     <- (tmp_m/tmp_psi + apply(Y[tmp_k_idx,,drop=F],2,sum))/(1/tmp_psi + sum(tmp_k_idx))
          result_a     <- tmp_a + 1/2 * sum(tmp_k_idx)
          result_b     <- tmp_b +
            1/2 * (1/tmp_psi * sum(tmp_k_idx))/(1/tmp_psi + sum(tmp_k_idx)) * (apply(Y[tmp_k_idx,,drop=F],2,mean) - tmp_m)^2
          if(sum(tmp_k_idx) > 1){ result_b <- result_b + 1/2 * (sum(tmp_k_idx)-1) * apply(Y[tmp_k_idx,,drop=F],2,var) }
          
        }else{
          
          result_psi <- tmp_psi
          result_m   <- tmp_m
          result_a   <- tmp_a
          result_b   <- tmp_b
        }
        
        return(rbind(result_psi, result_m, result_a, result_b))
      }
      
      PRES_Gamma_table <- Gamma_i %>% na.omit %>% table
      tmp_Gamma_table  <- table_cum(prev_Gamma_table, PRES_Gamma_table)
      updated_parms    <- sapply(seq_along(tmp_Gamma_table), tmp_fun0, simplify = F)
      
      m_hj <- psi_hj <- a_hj <- b_hj <- list()
      for(k in seq_along(tmp_Gamma_table)){
        
        psi_hj[[k]] <- rbind(updated_parms[[k]][1,], updated_parms[[k]][1,])
        m_hj[[k]]   <- rbind(updated_parms[[k]][2,], updated_parms[[k]][2,])
        a_hj[[k]]   <- rbind(updated_parms[[k]][3,], updated_parms[[k]][3,])
        b_hj[[k]]   <- rbind(updated_parms[[k]][4,], updated_parms[[k]][4,])
      }
    }
    
    cat("Iter :", Iter, "-", Gamma_i %>% is.na %>% sum, "\n")
  }
  
  # approximate PML
  kn <- c(prev_Gamma_i, Gamma_i) %>% na.omit %>% unique %>% length
  
  final_phi_t  <- phi_t[2,]
  final_pi     <- rbind(matrix(c(prev_Gamma_i, Gamma_i) %>% na.omit %>% table, nrow=kn, ncol=TT, byrow=F), alpha_stars) *
    matrix(1/(sum(!is.na(Gamma_i))+prev_n+alpha_stars), nrow=kn+1, ncol=TT, byrow=T)
  final_prob   <- (final_pi * matrix(final_phi_t, nrow=kn+1, ncol=TT, byrow=T)) %>% apply(1,sum)
  final_Lih <- NULL
  for(x in 1:kn){
    
    final_Lih <- rbind(final_Lih, sapply(1:n,
                                         function(ii){   Calculate_Noncentral_t(y    =Y[ii,],
                                                                                ai1  =a_hj[[x]][2,],
                                                                                bi1  =b_hj[[x]][2,],
                                                                                mi1  =m_hj[[x]][2,],
                                                                                psii1=psi_hj[[x]][2,]) }) %>% apply(2, prod) )
  }
  final_Lih <- rbind(final_Lih, sapply(1:n, function(ii){   Calculate_Noncentral_t(y    =Y[ii,],
                                                                                   ai1  =zero_a_hj,
                                                                                   bi1  =zero_b_hj,
                                                                                   mi1  =zero_m_hj,
                                                                                   psii1=zero_psi_hj) }) %>% apply(2, prod) )
  prob_table  <- final_Lih * matrix( final_prob, nrow=kn+1, ncol=n )
  log_aPML    <- prob_table %>% apply(2,sum) %>% log %>% mean
  log_aPML2   <- prob_table %>% apply(2,sum) %>% log %>% .[!is.na(Gamma_i)] %>% mean
  corr_prob   <- mean(apply(prob_table,2,which.max) == Gamma_i)
  EntropyVal  <- apply(prob_table, 2, cal_entropy) %>% mean
  
  Gamma_i2    <- Gamma_i
  tmp_gamma   <- prob_table %>% apply(2, which.max)
  Gamma_i2[is.na(Gamma_i)] <- tmp_gamma[is.na(Gamma_i)]
  
  return(list(log_aPML, Gamma_i, phi_t, psi_hj, m_hj, a_hj, b_hj, corr_prob, prob_table, EntropyVal, log_aPML2, Gamma_i2))
}

Repeat_Backward <- function(Y, a=1, b=1, m=0, psi=1, psi0=100, alpha_stars, eta_probs=NULL, 
                            ReParmUpd=TRUE, target_i_idx=NULL,
                            prev_Gamma_i=NULL, prev_phi_t=NULL, prev_psi_hj=NULL, 
                            prev_m_hj=NULL, prev_a_hj=NULL, prev_b_hj=NULL)
{
  # hyperparameters
  if(is.null(eta_probs)){ eta_probs  <- rep(1/length(alpha_stars), length(alpha_stars))}
  
  # basic setting
  n      <- nrow(Y)
  K      <- ncol(Y)
  TT     <- length(alpha_stars)
  prev_n <- sum(!is.na(prev_Gamma_i)) - sum(!is.na(prev_Gamma_i[target_i_idx]))
  
  # 0th value
  zero_psi_hj  <- rep(psi, K)
  zero_psi0_hj <- rep(psi0, K)
  zero_m_hj    <- rep(m, K)
  zero_a_hj    <- rep(a, K)
  zero_b_hj    <- rep(b, K)
  
  # parameters
  prev_Gamma_table <- prev_Gamma_i[-target_i_idx] %>% na.omit %>% table
  Gamma_i          <- prev_Gamma_i[target_i_idx]
  phi_t            <- prev_phi_t
  psi_hj           <- prev_psi_hj
  m_hj             <- prev_m_hj
  a_hj             <- prev_a_hj
  b_hj             <- prev_b_hj
  
  ## backward update
  delete_idx <- dqsample(which(!is.na(Gamma_i)))
  
  delete0_Gamma_i <- delete_Gamma_i <- Gamma_i
  delete0_m_hj    <- delete_m_hj    <- m_hj
  delete0_psi_hj  <- delete_psi_hj  <- psi_hj
  delete0_a_hj    <- delete_a_hj    <- a_hj
  delete0_b_hj    <- delete_b_hj    <- b_hj
  delete0_phi_t   <- delete_phi_t   <- phi_t
  delete0_Gamma_i_table <- delete_Gamma_i_table <-
    table_cum(prev_Gamma_table, delete_Gamma_i %>% na.omit %>% table)
  
  for(i in delete_idx){
    
    #### delete i
    delete_gam         <- delete_Gamma_i[i]
    delete0_Gamma_i[i] <- NA
    delete0_Gamma_i_table[delete_gam] <- delete0_Gamma_i_table[delete_gam] - 1
    
    ##### if we delete unique element, we need reindexing !
    delete_unq <- delete_Gamma_i_table[delete_gam] == 1
    if(delete_unq){
      
      ### update parameters
      delete0_psi_hj <- delete0_psi_hj[-delete_gam]
      delete0_m_hj   <- delete0_m_hj[-delete_gam]
      delete0_a_hj   <- delete0_a_hj[-delete_gam]
      delete0_b_hj   <- delete0_b_hj[-delete_gam]
      
      ### delete zero
      delete0_Gamma_i_table_nz <- delete0_Gamma_i_table[delete0_Gamma_i_table != 0]
      remaining_orig           <- as.numeric(names(delete0_Gamma_i_table_nz))  # 예: c(1, 2, 3, 5)
      
      ### update table
      delete0_Gamma_i_table_upd        <- delete0_Gamma_i_table_nz
      names(delete0_Gamma_i_table_upd) <- seq_along(delete0_Gamma_i_table_nz)
      
      ### index mapping
      new_mapping <- setNames(seq_along(remaining_orig), remaining_orig)
      
      ### update gamma
      delete0_Gamma_i       <- as.numeric(new_mapping[as.character(delete0_Gamma_i)])
      delete0_Gamma_i_table <- delete0_Gamma_i_table_upd
      
    }else{
      
      ### update parameters
      new_psi_hj <- delete0_psi_hj[[delete_gam]][2,]
      old_psi_hj <- 1/(1/new_psi_hj - 1)
      delete0_psi_hj[[delete_gam]][1,] <- old_psi_hj
      delete0_psi_hj[[delete_gam]][2,] <- old_psi_hj
      
      new_m_hj   <- delete0_m_hj[[delete_gam]][2,]
      old_m_hj   <- old_psi_hj * (new_m_hj/new_psi_hj - Y[i,])
      delete0_m_hj[[delete_gam]][1,]   <- old_m_hj
      delete0_m_hj[[delete_gam]][2,]   <- old_m_hj
      
      new_a_hj   <- delete0_a_hj[[delete_gam]][2,]
      old_a_hj   <- new_a_hj - 1/2
      delete0_a_hj[[delete_gam]][1,]   <- old_a_hj
      delete0_a_hj[[delete_gam]][2,]   <- old_a_hj
      
      new_b_hj   <- delete0_b_hj[[delete_gam]][2,]
      old_b_hj   <- new_b_hj - 1/2 * (Y[i,]-old_m_hj)^2/(old_psi_hj + 1)
      delete0_b_hj[[delete_gam]][1,]   <- old_b_hj
      delete0_b_hj[[delete_gam]][2,]   <- old_b_hj
    }
    
    ### pi_iht
    pi_iht     <- rbind(matrix(delete0_Gamma_i_table, nrow=length(delete0_Gamma_i_table), ncol=TT, byrow=F), alpha_stars) *
      matrix(1/(sum(!is.na(delete0_Gamma_i)) + prev_n + alpha_stars), nrow=length(delete0_Gamma_i_table)+1, ncol=TT, byrow=T)
    
    ### update phi_t
    if(delete_unq){
      tmp_pi <- pi_iht %>% tail(1) %>% c
    }else{
      tmp_pi <- pi_iht[delete_gam,]
    }
    new_phi_t  <- delete0_phi_t[2,]
    old_phi_t  <- new_phi_t/tmp_pi
    old_phi_t  <- old_phi_t/sum(old_phi_t)
    delete0_phi_t[1,]  <- old_phi_t
    delete0_phi_t[2,]  <- old_phi_t
    
    ## update
    delete_Gamma_i <- delete0_Gamma_i
    delete_m_hj    <- delete0_m_hj
    delete_psi_hj  <- delete0_psi_hj
    delete_a_hj    <- delete0_a_hj
    delete_b_hj    <- delete0_b_hj
    delete_phi_t   <- delete0_phi_t
    delete_Gamma_i_table <- delete0_Gamma_i_table
  }
  
  ### update after backward process
  Gamma_i <- delete_Gamma_i
  m_hj    <- delete_m_hj
  psi_hj  <- delete_psi_hj
  a_hj    <- delete_a_hj
  b_hj    <- delete_b_hj
  phi_t   <- delete_phi_t
  
  return(list(phi_t, psi_hj, m_hj, a_hj, b_hj))
}
