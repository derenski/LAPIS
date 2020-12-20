library(glmnet)
library(foreach)
library(doParallel)
library(quadprog)


treat.estimator <- function(Y, L.hat, W){ ## Efficient calculator of treatment estimator
  
  ## Y: Data matrix, L.hat: Estimated untreated matrix, W: Treatment indicator matrix
  
  upper.ones <- upper.tri(matrix(0, nrow=dim(Y)[2], ncol=dim(Y)[2]), diag = T)+0
  
  T.mat <- W %*% upper.ones
  
  treatment.estimator <- rep(NA, max(T.mat))
  
  for (i in 1:length(treatment.estimator)){
    
    treatment.estimator[i] <- mean(Y[T.mat==i]-L.hat[T.mat==i])
    
  } 
  
  return(treatment.estimator)
  
}

make_rho_mat <- function(rho,p){
  
  the_vec <- matrix(NA, nrow = p, ncol = p)
  
  for (i in 1:(p)){
    
    for (j in 1:p){
      
      the_vec[i,j] <- rho^abs(i-j)
      
    }
  }
  
  return(the_vec)
  
}

#### The Effect Size Functions

delta_t_updown <- function(t, arg_max=3, y_max=7, cutoff=3, ...){ ## Max is at exp(mu-1)
  
  implied_mu <- log(arg_max)+1
  
  scaling_factor <- y_max*exp(.5+log(arg_max))
  
  f_t <- scaling_factor*(1/t)*exp(-.5*(log(t)-implied_mu)^2)
  
  f_t[cutoff < cutoff & t > arg_max] <- cutoff
  
  return(f_t)
  
}


delta_t_decay <- function(t, y_max=6, halfway_time=4, cutoff=3, ...){ 
  
  implied_lambda <- -log(2)/(halfway_time-1)
  
  f_t <- y_max*exp(implied_lambda*(t-1))
  
  f_t[f_t < cutoff] <- cutoff
  
  return(f_t)
  
}


delta_t_plateau <- function(t, y_max=5, halfway_time=3, ...){ 
  
  implied_lambda <- -log(2)/halfway_time
  
  f_t <-  y_max*(1-exp(implied_lambda*(t)))
  
  return(f_t)
  
}

delta_t_constant <- function(t, value=10, ...){ 
  
  return(value)
  
}

treated_matrix_creator <- function(x, f_of_t=delta_t_constant,...){
  
  treatment_row <- rep(0, length(x))
  
  treatment_times <- which(x==1)
  
  if (length(treatment_times)!=0){
    
    treatment_row[treatment_times] <- f_of_t(treatment_times-(min(treatment_times)-1),...)
    
  }
  
  return(treatment_row)
  
}


list_of_functions <- list(delta_t_updown, delta_t_decay, delta_t_plateau, delta_t_constant)

names(list_of_functions) <- c("up_down", "decay", "plateau", "constant")

W_maker <- function(N, Time, ones_per_row){
  
  W <- t(sapply(ones_per_row, FUN=function(x) c(rep(0, Time-x),
                                                rep(1, x)) ))
  
  return(W)
  
}

### A Copula Method for Generating 
### Serial correlation in non-gaussian settings

norta <- function(number, corr_mat, desired_mean_matrix, 
                  distribution='gaussian', scalar_sigma=1,
                  df=3){
  
  Zs <- rmvn(n=number, mu=rep(0, dim(corr_mat)[1]), Sigma=corr_mat)
  
  phi_Zs <- pnorm(Zs)
  
  desired_variables <- matrix(NA, nrow=number, ncol=dim(corr_mat)[1])
  
  if (distribution == 'poisson'){
    
    for (row in 1:dim(desired_mean_matrix)[1]){
      
      desired_variables[row, ] <- qpois(phi_Zs[row, ], 
                                  lambda=desired_mean_matrix[row, ])
      
    }
    
    }else if(distribution == 'scaled_gamma'){
      
      for (row in 1:dim(desired_mean_matrix)[1]){
        
        desired_variables[row, ] <- qgamma(phi_Zs[row, ], 
                        shape=desired_mean_matrix[row, ], rate=1)
        
      }
      
    }else if (distribution == 'exponential'){
      
      for (row in 1:dim(desired_mean_matrix)[1]){
        
        desired_variables[row, ] <- qexp(phi_Zs[row, ], 
                      rate=1/desired_mean_matrix[row, ])
        
      }
      
    } else if (distribution=='t'){
      
      if (df < 3){
        
        print("Warning: May Not Have Desired Covariance Structure")
        
      }
      
      for (row in 1:dim(desired_mean_matrix)[1]){
      
      desired_variables[row, ] <- (qt(phi_Zs[row, ], df=df)
                                   +desired_mean_matrix[row, ])
      
      }
      
    }else if (distribution == 'gaussian'){
      
      for (row in 1:dim(desired_mean_matrix)[1]){
        
        desired_variables[row, ] <- qnorm(phi_Zs[row, ], 
        mean=desired_mean_matrix[row, ], sd=scalar_sigma)
        
      }
      
    } 
    
    
    return(desired_variables)
    
  }
  
S_tau <- function(x, tau=0){
  
  return(sign(x)*pmax(abs(x)-tau, 0))
  
}

P_omega <- function(A, W){ ## Helper function for projection
  ## Projects onto {(i,j) : W_{(i,j)} == 0}
  
  A[W!=0] <- 0
  
  return(A)
  
}


D_tau <- function(x, tau=0){
  
  svd_x <- svd(x)
  
  return(svd_x$u %*% S_tau(diag(svd_x$d), tau=tau) %*% t(svd_x$v))
  
}

robust_pca <- function(X, mu=1, lambda=1, delta){
  
  S_k <- matrix(0, nrow=dim(X)[1], ncol=dim(X)[2])
  
  Y_k <- matrix(0, nrow=dim(X)[1], ncol=dim(X)[2])
  
  L_k <- D_tau(x=X, tau=(1/mu))
  
  while (!(norm(X-L_k-S_k, 'F') < delta*norm(X, 'F')) & (norm(L_k, 'F') != 0) ){
    
    L_k <- D_tau(x=X-S_k+(1/mu)*Y_k, tau=(1/mu))
    
    S_k <- S_tau(x=X-L_k+(1/mu)*Y_k, tau=(lambda/mu))
    
    Y_k <- Y_k + mu*(X-L_k-S_k)
    
  }
  
  returning_list <- list(L_k,S_k)
  
  names(returning_list) <- c('L', 'S')
  
  return(returning_list)
  
}


## DID estimator

DID <- function(Y, W){
  
  control_or_treatment_unit <- apply(W, MARGIN=1, FUN = function(x) any(x == 1))
  
  control_units <- which(!control_or_treatment_unit)
  
  treatment_units <- which(control_or_treatment_unit)
  
  control_or_treatment_time <- apply(W, MARGIN=2, FUN = function(x) any(x == 1))
  
  control_times <- which(!control_or_treatment_time)
  
  treatment_times <- which(control_or_treatment_time)

  diff_treatment <- ((mean(Y[treatment_units, treatment_times]) -
                           mean(Y[treatment_units, control_times]))
                     
                     - (mean(Y[control_units, treatment_times]) -
                               mean(Y[control_units, control_times]))
                     
  )
  
  return(diff_treatment)
  
}

## Synthetic Control


synth_cont <- function(Y, W){
  
  control_or_treatment_unit <- apply(W, MARGIN=1, FUN = function(x) any(x == 1))
  
  control_units <- which(!control_or_treatment_unit)
  
  treated_units <- which(control_or_treatment_unit)
  
  control_or_treatment_time <- apply(W, MARGIN=2, FUN = function(x) any(x == 1))
  
  control_times <- which(!control_or_treatment_time)
  
  treatment_times <- which(control_or_treatment_time)
  
  zeta = sd(Y)^2
  
  T_1 <- length(treatment_times)
  
  T_0 <- dim(Y)[2]-T_1
  
  N_0 <- length(control_units)
  
  N_1 <- dim(Y)[1]-N_0
  
  
  ### Estimation of Omegas
  
  train_X_omega <- t(Y[control_units,  control_times])
  
  test_X_omega <- t(Y[control_units,  treatment_times])
  
  if (N_1==1){
    
    train_y_omega <- Y[treated_units,  control_times]
    
  }else{
    
    train_y_omega <- colMeans(Y[treated_units,  control_times])
    
  }
  
  solved_thing_omega <- solve.QP(
    Dmat = ((t(train_X_omega) %*% train_X_omega)/T_0 
            + (zeta/N_1)*diag(rep(1, N_0))),
    dvec=(1/T_0)*(t(train_y_omega) %*% train_X_omega),
    Amat=cbind(rep(1,dim(train_X_omega)[2]), 
               diag(rep(1, dim(train_X_omega)[2]))),
    bvec=c(1, rep(0, dim(train_X_omega)[2])), meq=1)
  
  coeff_vector <- matrix(c(solved_thing_omega$solution))
  
 # coeff_vector[coeff_vector < 1e-10] <- 0
  
  synthetic_counterfact <- t(test_X_omega %*% coeff_vector)
  
  
  if (N_1==1){
    
    synth_cont_est <- Y[treated_units,  treatment_times]-synthetic_counterfact
    
  }else{
    
    Y_treatment_only <- (Y*W)[treated_units, treatment_times]
    
    Y_treatment_only[Y_treatment_only==0] <- NA
    
    average_of_treated <- apply(Y_treatment_only, 2, FUN=mean, na.rm=T)
    
    synth_cont_est <- average_of_treated-synthetic_counterfact
    
  }
  
  return(synth_cont_est[1,])
  
}



# SDID and Alternate Implementation 

SDID_general <- function(Y, W,
                 iterations_for_coord_desc=100){ 
  # Synthetic Differences in Differences Method (Only Unit N Treated at Time T)
  
  control_or_treatment_unit <- apply(W, MARGIN=1, FUN = function(x) any(x == 1))
  
  control_units <- which(!control_or_treatment_unit)
  
  treatment_units <- which(control_or_treatment_unit)
  
  control_or_treatment_time <- apply(W, MARGIN=2, FUN = function(x) any(x == 1))
  
  control_times <- which(!control_or_treatment_time)
  
  treatment_times <- which(control_or_treatment_time)
  
  zeta = sd(Y)^2
  
  T_1 <- length(treatment_times)
  
  T_0 <- dim(Y)[2]-T_1
  
  N_0 <- length(control_units)
  
  N_1 <- dim(Y)[1]-N_0
  
  
  ### Estimation of Omegas
  
  train_X_omega <- t(Y[control_units,  control_times])
  
  if (N_1==1){
    
    train_y_omega <- Y[treatment_units,  control_times]
    
  }else{
    
    train_y_omega <- colMeans(Y[treatment_units,  control_times])
    
  }
  
  solved_thing_omega <- solve.QP(
    Dmat = ((t(train_X_omega) %*% train_X_omega)/T_0 
                           + (zeta/N_1)*diag(rep(1, N_0))),
                           dvec=(1/T_0)*(t(train_y_omega) %*% train_X_omega),
                           Amat=cbind(rep(1,dim(train_X_omega)[2]), 
                                      diag(rep(1, dim(train_X_omega)[2]))),
                           bvec=c(1, rep(0, dim(train_X_omega)[2])), meq=1)
  
  
  
  ### Estimation of Lambdas
  
  train_X_lambda <- Y[control_units,  control_times]
  
  if (T_1==1){
    
    train_y_lambda <- Y[control_units,  treatment_times]
    
  }else{
    
    train_y_lambda <- rowMeans(Y[control_units,  treatment_times])
    
  }
  
  
  

  
  solved_thing_lambda <- solve.QP(
    Dmat = ((t(train_X_lambda) %*% train_X_lambda)/N_0 
                                   + (zeta/T_1)*diag(rep(1, T_0))),
                           dvec=(1/N_0)*(t(train_y_lambda) %*% train_X_lambda),
                           Amat=cbind(rep(1,dim(train_X_lambda)[2]), 
                                      diag(rep(1, dim(train_X_lambda)[2]))),
                           bvec=c(1, rep(0, dim(train_X_lambda)[2])), meq=1)
  
  
  omega_vector <- matrix(c(solved_thing_omega$solution, rep(1/N_1, N_1)),
                         ncol=1)
  
  lambda_vector <- matrix(c(solved_thing_lambda$solution, rep(1/T_1, T_1)),
                          ncol=1)
  
  matrix_of_ones <- matrix(1, nrow=dim(Y)[1], ncol=dim(Y)[2])
  
  A_k <- matrix(0, nrow=dim(Y)[1], ncol=dim(Y)[2])
  
  B_k <- matrix(0, nrow=dim(Y)[1], ncol=dim(Y)[2])
  
  tau_k <- matrix(0, nrow=dim(Y)[2], ncol=1)
  
  mu_k <- 0
  
  for (k in 1:iterations_for_coord_desc){
    
  mu_k <- as.numeric(((t(omega_vector) %*% Y %*% lambda_vector) 
           - (t(omega_vector) %*% (A_k+B_k) %*% lambda_vector)
           -(t(omega_vector) %*% W %*% (lambda_vector*tau_k)))
           /(t(omega_vector) %*% matrix_of_ones %*% lambda_vector))
  
  alpha_k <- (((Y %*% lambda_vector) 
              - (B_k %*% lambda_vector)
              -mu_k*(matrix_of_ones %*% lambda_vector)
              -W %*% (lambda_vector*tau_k))
             /sum(lambda_vector))
  
  
  A_k <- matrix(rep(alpha_k, length.out=prod(dim(Y))), byrow=F,
                nrow=dim(Y)[1], ncol=dim(Y)[2])
  
  
  beta_k <- (((t(omega_vector) %*% Y) 
             - (t(omega_vector) %*% A_k)
             -mu_k*(t(omega_vector) %*% matrix_of_ones)
             -(t(omega_vector)%*% W) * t(tau_k))
            /sum(omega_vector))
  
  
  B_k <- matrix(rep(beta_k, length.out=prod(dim(Y))), byrow=T,
                nrow=dim(Y)[1], ncol=dim(Y)[2])
  
  tau_k <- as.numeric(((t(omega_vector) %*% (Y*W)) 
                       - (t(omega_vector) %*% ((A_k+B_k)*W))
                       -mu_k*(t(omega_vector) %*% W))
                      /(t(omega_vector) %*% W))
  
  tau_k[is.nan(tau_k)] <- 0
    
  }

  return(tau_k[treatment_times])
  
}



























##############################################################################################################
##############################################################################################################

### Helper functions for rank estimation


### Fan 2018 Method


corrected_eigen <- function(z, n, j, eigen_values){
  
  p <- length(eigen_values)
  
  m_nj <- ((p-j)^-1)*(sum((eigen_values[(j+1):p]-z)^-1) 
                           + (.25*(3*eigen_values[j]+eigen_values[j+1])-z)^-1)
  
  rho_j_n_minus_1 <- (p-j)/(n-1)
  
  bar_m_nj <- -(1-rho_j_n_minus_1)*(z^-1)+rho_j_n_minus_1*m_nj 
  
  lambda_C <- -(1/bar_m_nj)
  
  return(lambda_C)
  
}

##### L completion with automatic rank estimation (unknown rank)

completion_with_rank_estimation <- function(Y, W, r_init=40, 
                                            tolerance=1e-04, 
                                            min_iter=10,
                                            max_iter=1000,
                                            mu=1e-10){
  
  shrink_operator <- function(x, mu){ ## A helper function, the shrinkage operator
    
    if (x > 1*mu){
      
      return(x-mu)
      
    }else if (abs(x)< mu){
      
      return(0)
      
    }else{
      
      return(x+mu)
      
    }
    
  }
  
  N <- dim(Y)[1]
  
  Time <- dim(Y)[2]
  
  Us <- rmvn(n=N, mu=rep(0, r_init))
  
  Us <- apply(Us, MARGIN=2, FUN = function(x) x/norm(x, '2'))
  
  W_vec <- as.numeric(rmvn(n=1, mu=rep(0,r_init)))
  
  Vs <- rmvn(n=Time, mu=rep(0,r_init))
  
  Vs<- apply(Vs, MARGIN=2, FUN = function(x) x/norm(x, '2'))
  
  r <- r_init
  
  L_k <- P_omega(A=Y, W=W)
  
  iterating <- TRUE
  
  iter_number <- 1
  
  while (iterating){
    
    L_k_r <- L_k
    
    if(iter_number==1){
      
      parts_to_update <- 1:length(W_vec)
      
    }else{
    
    parts_to_update <- which(W_vec > 0)
    
    }
    
    for (r_number in parts_to_update){ ## Note: weights don't influence update of u,v (because we normalize)
        
        Us[, r_number] <- (L_k_r %*% Vs[, r_number])
        
        Us[, r_number] <- Us[, r_number]/norm(Us[, r_number], '2')
        
        Vs[, r_number] <- (t(L_k_r) %*% Us[, r_number])
        
        Vs[, r_number] <- Vs[, r_number]/norm(Vs[, r_number], '2')
        
        W_vec[r_number] <- max(0, shrink_operator(x=sum((as.matrix(Us[, r_number]) 
                                                         %*%  t(Vs[, r_number])) * L_k_r), mu=mu))
        
        L_k_r <- L_k_r - W_vec[r_number] * (Us[, r_number] %*% t(Vs[, r_number]))
      
    }
    
    L_k_plus_1 <- L_k
    
    Z <- L_k-L_k_r
    
    L_k_plus_1[W!=0] <- Z[W!=0]
    
    first_num <- (norm(P_omega(L_k_plus_1
                               -Z, W=W), 'F')
                  /norm(P_omega(L_k_plus_1, W=W)))
    ## Is rank r approx close to original matrix?
    
    
    second_num <- (norm(L_k_plus_1-L_k, 'F')
                   /norm(L_k_plus_1))
    
    ## Is updating changing much? 
    
    if (((first_num < tolerance) | (second_num  < tolerance) |
         (iter_number >= max_iter)) & !(iter_number < min_iter)){
      
      L_k_final <- Z
      
      iterating <- FALSE
      
      break
      
    }else{
      
      L_k <- L_k_plus_1
      
      iter_number <- iter_number+1
      
      next
      
    }
    
    
  }
  
  threshold_number <- (10^-3)*(sum(W==0)/prod(dim(W)))*sum(W_vec)
  
  R_star <- sum(W_vec > threshold_number)
  
  final_output <- list(Z, R_star)
  
  names(final_output) <- c("L_hat", 'rank_estimate')
  
  return(final_output)
  
}


#### Completion with rank estimation, validating mu

completion_with_rank_estimation_validate_mu <- function(Y, W,  
                                                        initial_rank=20,
                                                        tolerance=1e-03, 
                                                        min_iter=10,
                                                        max_iter=1000,
                                                        mu_grid=10^seq(-2, 2,1),
                                                        K=4){
  if (sum(is.na(Y)) > 0){
    
    Y[is.na(Y)] <- 0
    
  }
  
  
  
  
  doing_validation <- length(mu_grid) > 1
  
  if (doing_validation){
  
    card_O <- sum(W==0)
  
    fold_size <- floor((card_O^2)/prod(dim(Y)))
  
   # mean_squared_errors_mus <- c()
  
    mean_squared_errors_mus <- foreach(potential_mu=mu_grid,
                                       .packages=c('LaplacesDemon'),
                                       .combine = 'c',
                  .export = c('completion_with_rank_estimation',
                              'P_omega')) %dopar% {
    
      mean_squared_errors_this_k <- rep(NA, K)
    
      for (fold in 1:K){ ## Are they doing cv correctly?
      
        untreated_cells <- which(W==0)
      
        train_fold <- sample(untreated_cells, size=fold_size, replace=F)
      
        card_O_k <- length(train_fold)
      
        val_fold <- untreated_cells[!(untreated_cells %in% train_fold)]
      
        W_train <- W 
      
        W_train[val_fold] <- 1 ## Ignores validation cells
      
        completion_this_mu <- completion_with_rank_estimation(Y, W_train, 
                                                    r_init=initial_rank, 
                                                    tolerance=tolerance, 
                                                    min_iter=min_iter,
                                                    max_iter=max_iter,
                                                    mu=potential_mu)
      
        L_k <- completion_this_mu$L_hat
      
        mean_squared_errors_this_k[fold] <- mean((L_k[val_fold]-Y[val_fold])^2, 
                                               na.rm=T)
      
      }
    
      mean(mean_squared_errors_this_k, na.rm=T)
    
    }
  
    chosen_mu <- mu_grid[which.min(mean_squared_errors_mus)]
  
  }else{
    
    chosen_mu <- mu_grid
    
  }
  
  completion_chosen_mu <- completion_with_rank_estimation(Y, W, 
                                                        r_init=initial_rank, 
                                                        tolerance=tolerance, 
                                                        min_iter=min_iter,
                                                        max_iter=max_iter,
                                                        mu=chosen_mu)
  
  final_output <- list(completion_chosen_mu$L_hat, 
                       completion_chosen_mu$rank_estimate,
                       chosen_mu)
  
  names(final_output) <- c("L_hat", 'rank_estimate', 'chosen_mu')
  
  return(final_output)
  
}




##### Athey-based method
### Implementation of Athey et al's matrix recovery method. 

matrix_completion_causal <- function(Y, W, num_iter=1000, K=5, 
                            lambda_grid=c(0, 10^seq(-20, 0, 1)),
                            tol=1e-03){
  
  ### Cross Validation Section 
  
  card_O <- sum(W==0)
  
  fold_size <- (card_O^2)/prod(dim(Y))
  
  mean_squared_errors_lambdas <- foreach(potential_lambda=lambda_grid,
          .packages=c('LaplacesDemon'),
          .combine = 'c',
          .export = c('P_omega', 'D_tau', 'S_tau')) %dopar%{
    
    mean_squared_errors_this_k <- rep(NA, K)
    
    for (fold in 1:K){ ## Are they doing cv correctly?
      
      untreated_cells <- which(W==0)
      
      train_fold <- sample(untreated_cells, size=fold_size, replace=F)
      
      card_O_k <- length(train_fold)
      
      val_fold <- untreated_cells[!(untreated_cells %in% train_fold)]
      
      P_O_Y <- P_omega(Y, W=W)
      
    #  P_O_Y[(W==1)] <- 0
      
      P_O_Y[(untreated_cells %in% val_fold)] <- 0
      
      L_k <- P_O_Y
      
      first_L_k <- L_k
      
      for (iter_num in 1:num_iter){
        
        P_O_perp_L_k <- P_omega(L_k, W=1-W)
        
     #   P_O_perp_L_k[W==0] <- 0
        
        P_O_perp_L_k[(untreated_cells %in% train_fold)] <- 0
        
        L_k_plus_1 <- D_tau(x=P_O_Y+P_O_perp_L_k, tau=.5*potential_lambda*fold_size)
        
        L_difference <- mean((L_k_plus_1[W==1]-L_k[W==1])^2)
        
        if (L_difference > tol){
          
          L_k <- L_k_plus_1
          
          next
          
        }else{
          
          L_k <- L_k_plus_1
          
          break
          
        }
        
        
      }
      
      mean_squared_errors_this_k[fold] <- mean((L_k[val_fold]-Y[val_fold])^2, 
                                               na.rm=T)
      
    }
    
    mean(mean_squared_errors_this_k, na.rm=T)
    
  }
  
  cv_lambda <- lambda_grid[which.min(mean_squared_errors_lambdas)]
  
  card_O <- sum(W==0)
  
  P_O_Y <- P_omega(Y, W=W)
  
 # P_O_Y[W==1] <- 0
  
  L_k <- P_O_Y
  
  for (iter_num in 1:num_iter){
    
    P_O_perp_L_k <- P_omega(L_k, W=1-W)
    
    # P_O_perp_L_k[W==0] <- 0
    
    L_k_plus_1 <- D_tau(x=P_O_Y+P_O_perp_L_k, tau=.5*cv_lambda *card_O)
    
    L_difference <- mean((L_k_plus_1[W==1]-L_k[W==1])^2)
    
    if (L_difference > tol){
      
      L_k <- L_k_plus_1
      
      next
      
    }else{
      
      L_k <- L_k_plus_1
      
      break
      
    }
    
  }
  
  output <- list(L_k, cv_lambda)
  
  names(output) <- c("L_hat", "cv_lambda")
  
  return(output)
  
}


factor_model_rank_estimator <- function(Y){
  
  correlation_Y <- cor(Y)
  
  corr_eigen_values <- eigen(correlation_Y)$values
  
  corrected_eigen_values <- rep(NA, length(corr_eigen_values))
  
  for (k in 1:length(corr_eigen_values)){
    
    corrected_eigen_values[k] <- corrected_eigen(z=corr_eigen_values[k],
                      n=dim(Y)[1], j=k, eigen_values=corr_eigen_values)
    
  }
  
  threshold <- 1 + sqrt(dim(Y)[2]/(dim(Y)[1]-1))
  
  estimated_rank <- max(c(1, which(corrected_eigen_values > threshold)), na.rm=T)
  
  return(estimated_rank)
  
}



#### Rank Estimator Function

rank_estimator <-function(Y, W, num_iter=100, K=5, 
                          lambda_grid=c(0, 10^seq(-20, 3, 2)), 
                          method="threshold"){
  
  if (method == "threshold"){
    
    estimated_rank <- factor_model_rank_estimator(Y)
    
  }else if (method == "completion"){
    
    completion_output <- matrix_completion_causal(Y=Y, W=W, num_iter=num_iter, K=K, 
                                      lambda_grid=lambda_grid)
    
    L_hat <- completion_output$L_hat
    
    estimated_rank <- rankMatrix(L_hat)[1]
    
    return(estimated_rank)
    
  }
  
}


### Our alternative, assuming the rank is known

LAPIS <- function(Y, L_warm_start=NULL, 
                           W, rank_threshold=5, tolerance=1e-03, 
                           min_iter=10,
                           max_iter=1000,
                           method='non_explicit_tau',
                           change_memory_length=30,
                           minimum_change=30){ 
  # An alternative to Synthetic Difference in Differences
  
  treated_units <- as.numeric(which(apply(W, MARGIN=1, FUN = function(x) any(x==1))))
  
  treated_times <- as.numeric(which(apply(W, MARGIN=2, FUN = function(x) any(x==1))))
  
 # if (rank_threshold > rankMatrix(Y)){
    
#    rank_threshold <- rankMatrix(Y)[1]
    
#  }
  
    change_memory <- rep(NA, change_memory_length)
    
    d_k_minus_1 <- matrix(0, nrow=dim(W)[1], ncol=dim(W)[2])
    
    tau_k_minus_1 <- 0
    
    N <- dim(Y)[1]
    
    Time <- dim(Y)[2]
    
    if (is.null(L_warm_start)){
      
      Y_hat_k <- P_omega(A=Y, W=W)
      
    }else{
      
      Y_hat_k <- P_omega(A=L_warm_start, W=W)
      
    }

    iterating <- TRUE
    
    iter_number <- 1
    
    while (iterating){
      
      svd_untreated <- svd(Y_hat_k)
      
      full_diagonal <- svd_untreated$d
      
      truncated_diagonal <- c(full_diagonal[1:rank_threshold], rep(0,
      length(full_diagonal)-rank_threshold))
      
      rank_r_approx <- (svd_untreated$u %*% diag(truncated_diagonal) %*% 
                          t(svd_untreated$v))
      
      Y_hat_k_plus_1 <- P_omega(A=Y, W=W)+P_omega(A=rank_r_approx,W=1-W) 
  
      
      ## Is rank r approx close to original matrix?
      first_num <- (norm(P_omega(Y_hat_k_plus_1
                                 -rank_r_approx, W=W), 'F')
                    /norm(P_omega(Y_hat_k_plus_1, W=W)))
      
      tau_k <- treat.estimator(Y=Y, L.hat=rank_r_approx, W=W)
      
      second_num <- mean(abs(tau_k-tau_k_minus_1))
      
      change_memory[min(which(is.na(change_memory)))] <- second_num
      
      ## Is updating changing much?
      
      if ((((first_num < tolerance) | (second_num  < tolerance)) 
           & (iter_number >= min_iter))  | ((iter_number >= max_iter))
          ){
        
        tau_k_final <- tau_k
        
        iterating <- FALSE
        
        break
        
      }else{
        
        Y_hat_k <- Y_hat_k_plus_1
        
        tau_k_minus_1 <- tau_k
        
        iter_number <- iter_number+1
        
        if (!any(is.na(change_memory))){
          
          if(sum(change_memory) < minimum_change){
            
            if (rank_threshold > 1){
            
              rank_threshold <- rank_threshold-1

            
            }
            
          }
          
          change_memory <- rep(NA, change_memory_length)
        
        }
        
        next
        
      }
    
    }
  
  return(tau_k)
  
  }
  







### Our alternative, with rank estimation

LAPIS_with_rank_estimation <- function(Y, W,
                           initial_rank=20, tolerance=1e-03, 
                           min_iter=10, max_iter=1000, 
                           mu_grid=c(10^seq(-4,2,1),
                           seq(2,5,1)),
                           warm_start=F, 
                           method='non_explicit_tau',change_memory_length=30,
                           minimum_change=30,
                           num_folds=4){ 
  
  ### Could use CV to estimate mu
  
  ## Similar to Athey's Method
  
  if (!(is.null(mu_grid)) | warm_start){
  
  first_step <- completion_with_rank_estimation_validate_mu(Y, W, 
                                                initial_rank=initial_rank, 
                                                tolerance=tolerance, 
                                                min_iter=min_iter,
                                                max_iter=max_iter,
                                                mu_grid=mu_grid,
                                                K=num_folds)
  
  }
  
  if (warm_start){
    
    second_step <- LAPIS(Y=Y, L_warm_start=first_step$L_hat, W=W, 
                                  rank_threshold=first_step$rank_estimate, 
                                  tolerance=tolerance, 
                                  min_iter=min_iter, max_iter=max_iter, 
                                  method=method,
                                  change_memory_length=change_memory_length,
                                  minimum_change = minimum_change)
    
    final_output <- second_step
    
  }else{
    
    if (is.null(mu_grid)){
    
      second_step <- LAPIS(Y=Y, W=W, 
                                    rank_threshold=initial_rank, 
                                    tolerance=tolerance, 
                                    min_iter=min_iter, max_iter=max_iter, 
                                    method=method,change_memory_length=change_memory_length,
                                    minimum_change = minimum_change)
      final_output <- second_step
    
    }else{
      
      second_step <- LAPIS(Y=Y, W=W, 
                                    rank_threshold=first_step$rank_estimate, 
                                    tolerance=tolerance, 
                                    min_iter=min_iter, max_iter=max_iter, 
                                    method=method, change_memory_length=change_memory_length,
                                    minimum_change = minimum_change)
      
      final_output <- second_step
      
    }
    
  }
  
  return(final_output)
  
}


## Improving the rank heuristic

## Integrate Athey's shrinkage method with hard thresholding

# Threshold the singular values, but don't use shrinkage method explicitly

# instead, keep track of shrunken vector, but still use unshrunk singular values

# When shrunk beyond zero, set that singular value to 0

# Can cross-validate for penalty term





