library(ggplot2)
library(LaplacesDemon)
library(glmnet)
library(foreach)
library(doParallel)
library(quadprog)
library(openxlsx)
library(reshape)
library(plyr)
library(dplyr)
library(gsynth)
max_available_clusters <- detectCores()-1
  
desired_clusters <- 20
  
cl <- makeCluster(min(c(max_available_clusters, desired_clusters)))

registerDoParallel(cl)
 
source('causal_inference_methods_code.R')


mse_and_se_of_mse <- function(error_mat){
  
  squared_errors <- error_mat
  
  the_mse <- mean(squared_errors)
  
  se_mse <- sqrt((sd(squared_errors)^2)/prod(dim(squared_errors)))
  
  final_stuff <- c(the_mse, se_mse)
  
  names(final_stuff) <- c("mse", "se_mse")
  
  return(final_stuff)
  
}



params <- read.xlsx("parameters_and_descriptions.xlsx", sheet = "parameter_data", rowNames=T)


quiet <- function(x) { 
  sink(tempfile()) 
  on.exit(sink()) 
  invisible(force(x)) 
} 

### If you want to exclude gsynth from results, as it performs very similarly to your method
EXCLUDE_GSYNTH <- TRUE
# 3729 = GOOD SEED

set.seed(3729)

number_of_L <- as.numeric(params["number_of_Ls", ])

draws_per_L <- as.numeric(params["draws_per_L", ])

N <- as.numeric(params["N", ])

N0 <- N-as.numeric(params["N1", ])

Time <- as.numeric(params["Time", ])

Time0 <- as.numeric(Time)-as.numeric(params["Time1", ])

R <- as.numeric(params["R", ])

rho_parameter <- as.numeric(params['rho_parameter', ])

tau <- as.numeric(params['tau', ])

sigma_squared <- as.numeric(params['sigma_squared', ])

penalized <- as.logical(as.numeric(params['penalized', ]))

exchangable <- as.logical(as.numeric(params['exchangable', ]))

min_iter <- as.numeric(params['min_iter', ])

max_iter <- as.numeric(params['max_iter', ])

tolerance <- as.numeric(params['tolerance', ])

error = params['error', ]

df <- as.numeric(params['df', ])

rank_estimation_method <- params['rank_estimation_method', ]

L_scaling <- as.numeric(params['L_scaling', ])

arg_max <- as.numeric(params['arg_max', ])
  
y_max <- as.numeric(params['y_max', ])
  
halfway_time <- as.numeric(params['halfway_time', ])

cutoff <- as.numeric(params['cutoff', ])

design <- params['design', ]

lag_structure <- params['lag_structure', ]

average_treatment_length <- min(as.numeric(params['average_treatment_length', ]), Time-Time0)

max_lag <- as.numeric(params['max_lag', ])

treatment_function <- list_of_functions[[params['treatment_effect_function', ]]]















if (design=="staggered_adoption"){ ## Come up with a way to vary the lag in the staggered structure
  
  if(lag_structure == "random"){
    
    ones_we_make <- c(rep(0, N0), pmin(rpois(N-N0, 
                                            lambda=average_treatment_length-1)+1, 
                                      min(max_lag*(N-N0), .8*Time)))
    
  }else if (lag_structure=="constant"){ ## Does not control T-T0
    
    ones_we_make <- c(rep(0, N0), pmin(max_lag*seq(1, (N-N0)), floor(.8*Time)))
    
  }

}else if (design=="block_treatment"){
  
  ones_we_make <- c(rep(0, N0), rep(Time-Time0, N-N0))
  
}

W <- W_maker(N=N, Time=Time, ones_per_row = ones_we_make)

tau_matrix <- t(apply(W, MARGIN=1, FUN=treated_matrix_creator, 
                      f_of_t=treatment_function, arg_max=arg_max, 
                      y_max=y_max, halfway_time=halfway_time, cutoff=cutoff))

treated_units <- as.numeric(which(apply(W, MARGIN=1, FUN = function(x) any(x==1))))
  
treatment_times <- as.numeric(which(apply(W, MARGIN=2, FUN = function(x) any(x==1))))

delta_t <- treatment_function(treatment_times-(min(treatment_times)-1),
                              arg_max=arg_max, 
                      y_max=y_max, halfway_time=halfway_time, cutoff=cutoff, value=tau)
  
prediction_error_matrix_did <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

prediction_error_matrix_sc <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

                                          
prediction_error_matrix_gsynth <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)                                      
                                          
prediction_error_matrix_mc_nnm <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

prediction_error_matrix_sdid <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
                                          

prediction_error_matrix_lapis <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

prediction_error_matrix_oracle <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

autocorrelation_matrix <- make_rho_mat(rho=rho_parameter, p=dim(W)[2])

sig_to_noise_ratios <- c()

for (i in 1:number_of_L){
  
 # set.seed(3729)
  
  errors_this_L_did <- rep(NA, draws_per_L)
  
  errors_this_L_sc <- rep(NA, draws_per_L)
  
  errors_this_L_gsynth <- rep(NA, draws_per_L)  
    
  errors_this_L_mc_nnm <- rep(NA, draws_per_L)
  
  errors_this_L_sdid <- rep(NA, draws_per_L)
  
  errors_this_L_lapis <- rep(NA, draws_per_L)
  
  errors_this_L_oracle <- rep(NA, draws_per_L)
  
  if (exchangable){
    
    U_vec <- rexp(n=N*R, rate=1)

    V_vec <- rexp(n=Time*R, rate=1)
  
    U <- matrix(U_vec, nrow=N, ncol=R, byrow=T)

    V <- matrix(V_vec, nrow=Time, ncol=R, byrow=T)
  
  }else{
    
    U <- matrix(NA, nrow=N, ncol=R, byrow=T)

    V <- matrix(NA, nrow=Time, ncol=R, byrow=T)
    
    for (row_unit in 1:N){
      
      U[row_unit,] <- rpois(n=R, lambda=sqrt(row_unit/N))
      
    } 
    
    for (row_time in 1:Time){
      
      V[row_time,] <- rpois(n=R, lambda=sqrt(row_time/Time))
      
    }
    
  }

  L <- L_scaling*(U %*% t(V))
  
  for (j in 1:draws_per_L){
    
    if (error == 'gaussian'){
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W, distribution='gaussian',
                 scalar_sigma=sqrt(sigma_squared))
      
    } else if (error == 't'){
    
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='t', scalar_sigma=sqrt(sigma_squared))
    
    } else if (error == 'poisson'){
      
      if (exchangable == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='poisson', scalar_sigma=sqrt(sigma_squared))
      
    } else if (error == 'scaled_gamma'){
      
      if (exchangable == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='scaled_gamma', scalar_sigma=sqrt(sigma_squared))
      
    }else if (error == 'exponential'){
      
      if (exchangable == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='exponential', scalar_sigma=sqrt(sigma_squared))

    }
    
    #estimated_rank <- rank_estimator(Y, W, num_iter=100, K=5, 
    #                      lambda_grid=c(0, 10^seq(-20, 0, 2)), 
    #                      method="threshold")
    
    #Y_0_LAPIS <- LAPIS(Y, rank_threshold=estimated_rank,
    #                                    min_iter=1, max_iter=max_iter,
    #                                    tolerance=tolerance, W=W)
    
    
    if (N-N0 > 1){
    
    treatment_subjects_averaged <- colMeans(Y[1:dim(Y)[1] > N0,])
    
    W_averaged <- colMeans(W[1:dim(W)[1] > N0,])
    
    new_Y <- rbind(Y[1:dim(Y)[1] <= N0,], treatment_subjects_averaged)
    
    new_W <- rbind(W[1:dim(Y)[1] <= N0,], W_averaged)
    
    } else {
    
      new_Y <- Y
    
      new_W <- W
      
      
    }
      
    meltedD <- melt(W) 

    names(meltedD) <- c('id', 'time', 'treated')

    meltedObservedData <- melt(Y)

    names(meltedObservedData) <- c('id', 'time', 'value')

    joinedDataForGsynth <- quiet(meltedD %>% inner_join(meltedObservedData,
                                                       by=c('id', 'time'))  )
      
      
    
    tau_estimate_did <- DID(Y=Y, W=W)
    
    tau_estimate_sc <- synth_cont(Y=Y, W=W)
    
    tau_estimate_sdid <- SDID_general(Y=Y, W=W,
                 iterations_for_coord_desc=100)
      
      
    meltedD <- melt(W) 

    names(meltedD) <- c('id', 'time', 'treated')

    meltedObservedData <- melt(Y)

    names(meltedObservedData) <- c('id', 'time', 'value')
    
    mc_nnm_info <- matrix_completion_causal(Y=Y, W=W, num_iter=1000, K=5, 
                            lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                            tol=1e-04)
    
    L_mc_nnm <- mc_nnm_info$L_hat
    
    tau_estimate_mc_nnm <- treat.estimator(Y=Y, L.hat=L_mc_nnm, W=W)
    
    estFactors <- rankMatrix(mc_nnm_info$L_hat)[1]  
      
    gsynthInfo <- quiet(gsynth(value~treated, data=joinedDataForGsynth, index=c('id', 'time'), 
                         parallel = TRUE, r=estFactors))

    gsynthContEst <- gsynthInfo$att
    
    if (design=='block_treatment'){
      
    tau_estimate_gsynth <- gsynthContEst[(Time0+1):Time]
        
        }else{
                            
    tau_estimate_gsynth <- gsynthContEst[
        which(names(gsynthContEst)==1):length(gsynthContEst)]
        
    }
    
    tau_estimate_lapis <- LAPIS_with_rank_estimation(Y=Y, 
                           W=W, initial_rank=rankMatrix(mc_nnm_info$L_hat)[1],
                           tolerance=tolerance, 
                           min_iter=min_iter, max_iter=max_iter,   
                           mu_grid=NULL, warm_start=F, method = 'explicit_tau')
    
    tau_estimate_oracle <- treat.estimator(Y=Y, L.hat=L, W=W)
    
    # stopCluster(cl)
    
    # c(10^seq(-4,2,1), seq(2,5,1))


  
    
    
    # tau_estimate_mc_nnm
    
    ## Only oracle in the sense that we know L
    
    error_tau_sc <- mean(abs(tau_estimate_sc-delta_t)^2)
      
    error_tau_gsynth <- mean(abs(tau_estimate_gsynth-delta_t)^2)
    
    error_tau_did <- mean(abs(tau_estimate_did-delta_t)^2)
    
    error_tau_mc_nnm <- mean(abs(tau_estimate_mc_nnm-delta_t)^2)
    
    error_tau_sdid <- mean(abs(tau_estimate_sdid-delta_t)^2)
    
    error_tau_lapis <- mean(abs(tau_estimate_lapis-delta_t)^2)
    
    error_tau_oracle <- mean(abs(tau_estimate_oracle-delta_t)^2)

    errors_this_L_did[j] <- error_tau_did
    
    errors_this_L_sc[j] <- error_tau_sc
      
    errors_this_L_gsynth[j] <- error_tau_gsynth
    
    errors_this_L_mc_nnm[j] <- error_tau_mc_nnm
    
    errors_this_L_sdid[j] <- error_tau_sdid
    
    errors_this_L_lapis[j] <- error_tau_lapis
    
    errors_this_L_oracle[j] <- error_tau_oracle

  }
  
  prediction_error_matrix_did[i,] <- errors_this_L_did
  
  prediction_error_matrix_sc[i,] <- errors_this_L_sc
    
  prediction_error_matrix_gsynth[i,] <- errors_this_L_gsynth
  
  prediction_error_matrix_mc_nnm[i, ] <- errors_this_L_mc_nnm
  
  prediction_error_matrix_sdid[i,] <- errors_this_L_sdid
  
  prediction_error_matrix_lapis[i,] <- errors_this_L_lapis
  
  prediction_error_matrix_oracle[i,] <- errors_this_L_oracle

}


effect_plot <- (ggplot(NULL, aes(x=1:length(delta_t), y=delta_t)) 
                + geom_point() + theme_bw() + xlab("Time") + ylab("Treatment Effect")
                +ggtitle("True Treatment Effect Over Time"))


mse_and_se_of_mse_did <- mse_and_se_of_mse(prediction_error_matrix_did)


mse_and_se_of_mse_sc <- mse_and_se_of_mse(prediction_error_matrix_sc)


mse_and_se_of_mse_gsynth <- mse_and_se_of_mse(prediction_error_matrix_gsynth)

mse_and_se_of_mse_mc_nnm <- mse_and_se_of_mse(prediction_error_matrix_mc_nnm)

mse_and_se_of_mse_sdid <- mse_and_se_of_mse(prediction_error_matrix_sdid)

mse_and_se_of_mse_lapis <- mse_and_se_of_mse(prediction_error_matrix_lapis)

mse_and_se_of_mse_oracle <- mse_and_se_of_mse(prediction_error_matrix_oracle)


fixed_parameter_error_frame <- data.frame(t(cbind.data.frame(mse_and_se_of_mse_did, mse_and_se_of_mse_sc,
                                                mse_and_se_of_mse_gsynth,
                mse_and_se_of_mse_sdid, mse_and_se_of_mse_mc_nnm, mse_and_se_of_mse_lapis,
                mse_and_se_of_mse_oracle)))
                                          
rownames(fixed_parameter_error_frame) <- NULL
                                          
fixed_parameter_error_frame$Method <- c('DID','SC', 'GSYNTH', 'SDID', "MC-NNM", 'LAPIS',
'ORACLE')

colnames(fixed_parameter_error_frame)[1:2] <- c('MSE', 'SE')

tableForPresenting <- (fixed_parameter_error_frame[, c('Method', 'MSE', 'SE')] %>% mutate(MSE = round(MSE, 3),
                                                                                        SE=round(SE, 3),
                                                                                        `MSE (SE)` = paste(MSE
                                                                                                          ,' (', SE, ')',
            sep='')))[, c('Method', 'MSE (SE)')]
                                          
            
if (EXCLUDE_GSYNTH){
    
    tableForPresenting <- tableForPresenting %>% filter(Method != 'GSYNTH')
    
    }


tableForPresenting <- knitr::kable(tableForPresenting, 'latex')


                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
                                          
### WITH N0/N
## Chunk 18

all_N0s <- floor((c(5, seq(10, 90, 10), 98)/100)*N)

mses_N0_did <- c()

se_mses_N0_did <- c()

mses_N0_gsynth <- c()

se_mses_N0_gsynth <- c()

mses_N0_sc<- c()

se_mses_N0_sc <- c()

mses_N0_mc_nnm <- c()

se_mses_N0_mc_nnm <- c()

mses_N0_sdid <- c()

se_mses_N0_sdid <- c()

mses_N0_lapis <- c()

se_mses_N0_lapis <- c()

mses_N0_oracle <- c()

se_mses_N0_oracle <- c()

the_N0_over_Ns <- c()

for (this_N0 in all_N0s){ 
  
  the_N0_over_Ns <- c(the_N0_over_Ns, this_N0/N)
  
  set.seed(3729)
  
  if (design=="staggered_adoption"){ ## Come up with a way to vary the lag in the staggered structure
  
  if(lag_structure == "random"){
    
    ones_we_make <- c(rep(0, this_N0), pmin(rpois(N-this_N0, 
                                            lambda=average_treatment_length-1)+1, 
                                      min(max_lag*(N-this_N0), .8*Time)))
    
  }else if (lag_structure=="constant"){ ## Does not control T-T0
    
    ones_we_make <- c(rep(0, this_N0), pmin(max_lag*seq(1, (N-this_N0)), floor(.8*Time)))
    
  }

}else if (design=="block_treatment"){
  
  ones_we_make <- c(rep(0, this_N0), rep(Time-Time0, N-this_N0))
  
}
  
  W <- W_maker(N=N, Time=Time, ones_per_row = ones_we_make)

  tau_matrix <- t(apply(W, MARGIN=1, FUN=treated_matrix_creator, 
                      f_of_t=treatment_function, arg_max=arg_max, 
                      y_max=y_max, halfway_time=halfway_time, cutoff=cutoff, 
                      value=tau))

  treated_units <- as.numeric(which(apply(W, MARGIN=1, FUN = function(x) any(x==1))))
  
  treatment_times <- as.numeric(which(apply(W, MARGIN=2, FUN = function(x) any(x==1))))

  delta_t <- treatment_function(treatment_times-(min(treatment_times)-1),
                              arg_max=arg_max, 
                      y_max=y_max, halfway_time=halfway_time, cutoff=cutoff, value=tau)
  
  autocorrelation_matrix <- make_rho_mat(rho=rho_parameter, p=dim(W)[2])

  prediction_error_matrix_did <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_sc <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
                                            
  prediction_error_matrix_gsynth <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_mc_nnm <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

  prediction_error_matrix_sdid <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

  prediction_error_matrix_lapis <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_oracle <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
 
  for (i in 1:number_of_L){
  
 # set.seed(3729)
  
  errors_this_L_did <- rep(NA, draws_per_L)
  
  errors_this_L_sc <- rep(NA, draws_per_L)
      
  errors_this_L_gsynth <- rep(NA, draws_per_L)
  
  errors_this_L_mc_nnm <- rep(NA, draws_per_L)
  
  errors_this_L_sdid <- rep(NA, draws_per_L)
  
  errors_this_L_lapis <- rep(NA, draws_per_L)
  
  errors_this_L_oracle <- rep(NA, draws_per_L)
  
  if (exchangable){
    
    U_vec <- rexp(n=N*R, rate=1)

    V_vec <- rexp(n=Time*R, rate=1)
  
    U <- matrix(U_vec, nrow=N, ncol=R, byrow=T)

    V <- matrix(V_vec, nrow=Time, ncol=R, byrow=T)
  
  }else{
    
    U <- matrix(NA, nrow=N, ncol=R, byrow=T)

    V <- matrix(NA, nrow=Time, ncol=R, byrow=T)
    
    for (row_unit in 1:N){
      
      U[row_unit,] <- rpois(n=R, lambda=sqrt(row_unit/N))
      
    } 
    
    for (row_time in 1:Time){
      
      V[row_time,] <- rpois(n=R, lambda=sqrt(row_time/Time))
      
    }
    
  }

  L <- L_scaling*(U %*% t(V))
  
  for (j in 1:draws_per_L){
    
    if (error == 'gaussian'){
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W, distribution='gaussian',
                 scalar_sigma=sqrt(sigma_squared))
      
    } else if (error == 't'){
    
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='t', scalar_sigma=sqrt(sigma_squared))
    
    } else if (error == 'poisson'){
      
      if (exchangable == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='poisson', scalar_sigma=sqrt(sigma_squared))
      
    } else if (error == 'scaled_gamma'){
      
      if (exchangable == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='scaled_gamma', scalar_sigma=sqrt(sigma_squared))
      
    }else if (error == 'exponential'){
      
      if (exchangable == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='exponential', scalar_sigma=sqrt(sigma_squared))

    }
    
    #estimated_rank <- rank_estimator(Y, W, num_iter=100, K=5, 
    #                      lambda_grid=c(0, 10^seq(-20, 0, 2)), 
    #                      method="threshold")
    
    #Y_0_LAPIS <- LAPIS(Y, rank_threshold=estimated_rank,
    #                                    min_iter=1, max_iter=max_iter,
    #                                    tolerance=tolerance, W=W)
    
    
    if (N-this_N0 > 1){
    
    treatment_subjects_averaged <- colMeans(Y[1:dim(Y)[1] > this_N0,])
    
    W_averaged <- colMeans(W[1:dim(W)[1] > this_N0,])
    
    new_Y <- rbind(Y[1:dim(Y)[1] <= this_N0,], treatment_subjects_averaged)
    
    new_W <- rbind(W[1:dim(Y)[1] <= this_N0,], W_averaged)
    
    } else {
    
      new_Y <- Y
    
      new_W <- W
      
      
    }
    
    tau_estimate_did <- DID(Y=Y, W=W)
    
    tau_estimate_sc <- synth_cont(Y=Y, W=W)
      
      
    meltedD <- melt(W) 

    names(meltedD) <- c('id', 'time', 'treated')

    meltedObservedData <- melt(Y)

    names(meltedObservedData) <- c('id', 'time', 'value')

    joinedDataForGsynth <- quiet(meltedD %>% inner_join(meltedObservedData,
                                                       by=c('id', 'time')) )                         

    tau_estimate_sdid <- SDID_general(Y=Y, W=W,
                 iterations_for_coord_desc=100)
    
    mc_nnm_info <- matrix_completion_causal(Y=Y, W=W, num_iter=1000, K=5, 
                            lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                            tol=1e-04)
    
    L_mc_nnm <- mc_nnm_info$L_hat
    
    tau_estimate_mc_nnm <- treat.estimator(Y=Y, L.hat=L_mc_nnm, W=W)
      
      
    estFactors <- rankMatrix(mc_nnm_info$L_hat)[1]  
      
    gsynthInfo <- quiet(gsynth(value~treated, data=joinedDataForGsynth, index=c('id', 'time'), 
                         parallel = TRUE, r=estFactors))

    gsynthContEst <- gsynthInfo$att
    
    if (design=='block_treatment'){
      
    tau_estimate_gsynth <- gsynthContEst[(Time0+1):Time]
        
        }else{
                            
    tau_estimate_gsynth <- gsynthContEst[
        which(names(gsynthContEst)==1):length(gsynthContEst)]
        
    }  
    
    tau_estimate_lapis <- LAPIS_with_rank_estimation(Y=Y, 
                           W=W, initial_rank=rankMatrix(mc_nnm_info$L_hat)[1],
                           tolerance=tolerance, 
                           min_iter=min_iter, max_iter=max_iter,   
                           mu_grid=NULL, warm_start=F, method = 'explicit_tau')
    
    tau_estimate_oracle <- treat.estimator(Y=Y, L.hat=L, W=W)
    

    
    error_tau_sc <- mean(abs(tau_estimate_sc-delta_t)^2)
      
    error_tau_gsynth <- mean(abs(tau_estimate_gsynth-delta_t)^2)
    
    error_tau_did <- mean(abs(tau_estimate_did-delta_t)^2)
    
    error_tau_mc_nnm <- mean(abs(tau_estimate_mc_nnm-delta_t)^2)
    
    error_tau_sdid <- mean(abs(tau_estimate_sdid-delta_t)^2)
    
    error_tau_lapis <- mean(abs(tau_estimate_lapis-delta_t)^2)
    
    error_tau_oracle <- mean(abs(tau_estimate_oracle-delta_t)^2)

    errors_this_L_did[j] <- error_tau_did
    
    errors_this_L_sc[j] <- error_tau_sc
      
    errors_this_L_gsynth[j] <- error_tau_gsynth
    
    errors_this_L_mc_nnm[j] <- error_tau_mc_nnm
    
    errors_this_L_sdid[j] <- error_tau_sdid
    
    errors_this_L_lapis[j] <- error_tau_lapis
    
    errors_this_L_oracle[j] <- error_tau_oracle

  }
  
  prediction_error_matrix_did[i,] <- errors_this_L_did
  
  prediction_error_matrix_sc[i,] <- errors_this_L_sc
      
      
  prediction_error_matrix_gsynth[i,] <- errors_this_L_gsynth
  
  prediction_error_matrix_mc_nnm[i, ] <- errors_this_L_mc_nnm
  
  prediction_error_matrix_sdid[i,] <- errors_this_L_sdid
  
  prediction_error_matrix_lapis[i,] <- errors_this_L_lapis
  
  prediction_error_matrix_oracle[i,] <- errors_this_L_oracle

}
  
  
    mse_and_se_of_mse_did <- mse_and_se_of_mse(prediction_error_matrix_did)
    
    mse_and_se_of_mse_sc <- mse_and_se_of_mse(prediction_error_matrix_sc)
                                            
    mse_and_se_of_mse_gsynth <- mse_and_se_of_mse(prediction_error_matrix_gsynth)
    
    mse_and_se_of_mse_mc_nnm <- mse_and_se_of_mse(prediction_error_matrix_mc_nnm)
    
    mse_and_se_of_mse_sdid <- mse_and_se_of_mse(prediction_error_matrix_sdid)
    
    mse_and_se_of_mse_lapis <- mse_and_se_of_mse(prediction_error_matrix_lapis)
  
    mse_and_se_of_mse_oracle <- mse_and_se_of_mse(prediction_error_matrix_oracle)
    
    mses_N0_did <- c(mses_N0_did, mse_and_se_of_mse_did[1])
    
    se_mses_N0_did <- c(se_mses_N0_did, mse_and_se_of_mse_did[2])
    
    mses_N0_sc <- c(mses_N0_sc, mse_and_se_of_mse_sc[1])
    
    se_mses_N0_sc <- c(se_mses_N0_sc, mse_and_se_of_mse_sc[2])
                                            
                                            
    mses_N0_gsynth <- c(mses_N0_gsynth, mse_and_se_of_mse_gsynth[1])
    
    se_mses_N0_gsynth <- c(se_mses_N0_gsynth, mse_and_se_of_mse_gsynth[2])                                        
    
    mses_N0_mc_nnm <- c(mses_N0_mc_nnm, mse_and_se_of_mse_mc_nnm[1])
    
    se_mses_N0_mc_nnm <- c(se_mses_N0_mc_nnm, mse_and_se_of_mse_mc_nnm[2])
    
    mses_N0_sdid <- c(mses_N0_sdid, mse_and_se_of_mse_sdid[1])
    
    se_mses_N0_sdid <- c(se_mses_N0_sdid, mse_and_se_of_mse_sdid[2])
    
    mses_N0_lapis <- c(mses_N0_lapis, mse_and_se_of_mse_lapis[1])
    
    se_mses_N0_lapis <- c(se_mses_N0_lapis, mse_and_se_of_mse_lapis[2])
    
    mses_N0_oracle <- c(mses_N0_oracle, mse_and_se_of_mse_oracle[1])
    
    se_mses_N0_oracle <- c(se_mses_N0_oracle, mse_and_se_of_mse_oracle[2])
                                            
                                            
    print(paste("Finished N0:", this_N0))

}
  
  # mses_N0_mc_nnm, mses_N0_did, 


N0_data <- cbind(c( mses_N0_did, mses_N0_sc, mses_N0_gsynth, mses_N0_mc_nnm,
              mses_N0_sdid, mses_N0_lapis, mses_N0_oracle),
              
             c(se_mses_N0_did, se_mses_N0_sc, se_mses_N0_gsynth, se_mses_N0_mc_nnm,
              se_mses_N0_sdid, se_mses_N0_lapis, se_mses_N0_oracle)
              
)



## mses_N0_mc_nnm, mses_N0_did,

N0_sensitivity_data <- data.frame(rep(c('DID','SC', 'GSYNTH', "MC-NNM",'SDID', 'LAPIS',
'ORACLE'), 
               each=length(mses_N0_lapis)))

names(N0_sensitivity_data) <- 'Method'

N0_sensitivity_data$N0 <- c(5, seq(10, 90, 10), 98)/100

N0_sensitivity_data$mse <- N0_data[,1]

N0_sensitivity_data$se <- N0_data[,2]
                                            
                                            
                                            
if (EXCLUDE_GSYNTH){                                            
                                            
N0_sensitivity_data <- N0_sensitivity_data %>% filter(Method != 'GSYNTH')   
    
    }
                                            

p_mse_vs_N0 <- (ggplot(N0_sensitivity_data, aes(x=N0, y=sqrt(mse), col=Method)) + geom_line() + 
                   theme_bw()+ ggtitle("rmse as a Function of N0/N") 
                +xlab("N0/N"))

#p_mse_vs_N0 <- (ggplot(N0_sensitivity_data, aes(x=N0, y=mse, col=Method)) + geom_line() + 
#                   geom_ribbon(aes(ymin=mse-1.9*se,
#                   ymax=mse+1.9*se, alpha=.1), fill = "grey70", lty=2) +
#                   theme_bw()+ ggtitle("Mse as a Function of N0"))
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
## Chunk 19

signal_to_noise_ratios <- c()

all_rhos <- seq(0, .95, .05)

mses_rho_did <- c()

se_mses_rho_did <- c()

mses_rho_sc<- c()

se_mses_rho_sc <- c()

mses_rho_gsynth <- c()

se_mses_rho_gsynth <- c()

mses_rho_mc_nnm <- c()

se_mses_rho_mc_nnm <- c()

mses_rho_sdid <- c()

se_mses_rho_sdid <- c()

mses_rho_lapis <- c()

se_mses_rho_lapis <- c()

mses_rho_oracle <- c()

se_mses_rho_oracle <- c()

for (this_rho in all_rhos){ 
  
  set.seed(3729)
  
  this_autocorrelation_matrix <- make_rho_mat(rho=this_rho, p=dim(W)[2])
  
  if (design=="staggered_adoption"){ ## Come up with a way to vary the lag in the staggered structure
  
  if(lag_structure == "random"){
    
    ones_we_make <- c(rep(0, N0), pmin(rpois(N-N0, 
                                            lambda=average_treatment_length-1)+1, 
                                      min(max_lag*(N-N0), .8*Time)))
    
  }else if (lag_structure=="constant"){ ## Does not control T-T0
    
    ones_we_make <- c(rep(0, N0), pmin(max_lag*seq(1, (N-N0)), floor(.8*Time)))
      
      }
    
  

}else if (design=="block_treatment"){
  
  ones_we_make <- c(rep(0, N0), rep(Time-Time0, N-N0))
  
}

  W <- W_maker(N=N, Time=Time, ones_per_row = ones_we_make)

  tau_matrix <- t(apply(W, MARGIN=1, FUN=treated_matrix_creator, 
                      f_of_t=treatment_function, arg_max=arg_max, 
                      y_max=y_max, halfway_time=halfway_time, cutoff=cutoff,
                      value=tau))

  treated_units <- as.numeric(which(apply(W, MARGIN=1, FUN = function(x) any(x==1))))
  
  treatment_times <- as.numeric(which(apply(W, MARGIN=2, FUN = function(x) any(x==1))))

  delta_t <- treatment_function(treatment_times-(min(treatment_times)-1),
                              arg_max=arg_max, 
                      y_max=y_max, halfway_time=halfway_time, cutoff=cutoff, value=tau)
 
  prediction_error_matrix_did <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_sc <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
                                            
  prediction_error_matrix_gsynth <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_mc_nnm <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

  prediction_error_matrix_sdid <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

  prediction_error_matrix_lapis <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_oracle <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  for (i in 1:number_of_L){
  
 # set.seed(3729)
  
  errors_this_L_did <- rep(NA, draws_per_L)
  
  errors_this_L_sc <- rep(NA, draws_per_L)
  
  errors_this_L_mc_nnm <- rep(NA, draws_per_L)
  
  errors_this_L_sdid <- rep(NA, draws_per_L)
  
  errors_this_L_lapis <- rep(NA, draws_per_L)
  
  errors_this_L_oracle <- rep(NA, draws_per_L)
  
  if (exchangable){
    
    U_vec <- rexp(n=N*R, rate=1)

    V_vec <- rexp(n=Time*R, rate=1)
  
    U <- matrix(U_vec, nrow=N, ncol=R, byrow=T)

    V <- matrix(V_vec, nrow=Time, ncol=R, byrow=T)
  
  }else{
    
    U <- matrix(NA, nrow=N, ncol=R, byrow=T)

    V <- matrix(NA, nrow=Time, ncol=R, byrow=T)
    
    for (row_unit in 1:N){
      
      U[row_unit,] <- rpois(n=R, lambda=sqrt(row_unit/N))
      
    } 
    
    for (row_time in 1:Time){
      
      V[row_time,] <- rpois(n=R, lambda=sqrt(row_time/Time))
      
    }
    
  }

  L <- L_scaling*(U %*% t(V))
  
  for (j in 1:draws_per_L){
    
    if (error == 'gaussian'){
      
      Y <- norta(number=N, corr_mat=this_autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W, distribution='gaussian',
                 scalar_sigma=sqrt(sigma_squared))
      
    } else if (error == 't'){
    
      Y <- norta(number=N, corr_mat=this_autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='t', scalar_sigma=sqrt(sigma_squared))
    
    } else if (error == 'poisson'){
      
      if (exchangable == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=this_autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='poisson', scalar_sigma=sqrt(sigma_squared))
      
    } else if (error == 'scaled_gamma'){
      
      if (exchangable == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=this_autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='scaled_gamma', scalar_sigma=sqrt(sigma_squared))
      
    }else if (error == 'exponential'){
      
      if (exchangable == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=this_autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='exponential', scalar_sigma=sqrt(sigma_squared))

    }
    
    #estimated_rank <- rank_estimator(Y, W, num_iter=100, K=5, 
    #                      lambda_grid=c(0, 10^seq(-20, 0, 2)), 
    #                      method="threshold")
    
    #Y_0_LAPIS <- LAPIS(Y, rank_threshold=estimated_rank,
    #                                    min_iter=1, max_iter=max_iter,
    #                                    tolerance=tolerance, W=W)
    
    
    if (N-N0 > 1){
    
    treatment_subjects_averaged <- colMeans(Y[1:dim(Y)[1] > N0,])
    
    W_averaged <- colMeans(W[1:dim(W)[1] > N0,])
    
    new_Y <- rbind(Y[1:dim(Y)[1] <= N0,], treatment_subjects_averaged)
    
    new_W <- rbind(W[1:dim(Y)[1] <= N0,], W_averaged)
    
    } else {
    
      new_Y <- Y
    
      new_W <- W
      
      
    }
    
    tau_estimate_did <- DID(Y=Y, W=W)
    
    tau_estimate_sc <- synth_cont(Y=Y, W=W)
      
      
    meltedD <- melt(W) 

    names(meltedD) <- c('id', 'time', 'treated')

    meltedObservedData <- melt(Y)

    names(meltedObservedData) <- c('id', 'time', 'value')

    joinedDataForGsynth <- quiet(meltedD %>% inner_join(meltedObservedData,
                                                       by=c('id', 'time')))
    
    tau_estimate_sdid <- SDID_general(Y=Y, W=W,
                 iterations_for_coord_desc=100)
    
    mc_nnm_info <- matrix_completion_causal(Y=Y, W=W, num_iter=1000, K=5, 
                            lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                            tol=1e-04)
    
    L_mc_nnm <- mc_nnm_info$L_hat
    
    tau_estimate_mc_nnm <- treat.estimator(Y=Y, L.hat=L_mc_nnm, W=W)
      
      
    estFactors <- rankMatrix(mc_nnm_info$L_hat)[1]  
      
    gsynthInfo <- quiet(gsynth(value~treated, data=joinedDataForGsynth, 
                         index=c('id', 'time'), 
                         parallel = TRUE, r=estFactors))

    gsynthContEst <- gsynthInfo$att
    
    if (design=='block_treatment'){
      
        tau_estimate_gsynth <- gsynthContEst[(Time0+1):Time]
        
        }else{
                            
        tau_estimate_gsynth <- gsynthContEst[
        which(names(gsynthContEst)==1):length(gsynthContEst)]
        
    }
    
    tau_estimate_lapis <- LAPIS_with_rank_estimation(Y=Y, 
                           W=W, initial_rank=rankMatrix(mc_nnm_info$L_hat)[1],
                           tolerance=tolerance, 
                           min_iter=min_iter, max_iter=max_iter,   
                           mu_grid=NULL, warm_start=F, method = 'explicit_tau')
    
    tau_estimate_oracle <- treat.estimator(Y=Y, L.hat=L, W=W)
    
    ## Only oracle in the sense that we know L
    
    error_tau_sc <- mean(abs(tau_estimate_sc-delta_t)^2)
      
    error_tau_gsynth <- mean(abs(tau_estimate_gsynth-delta_t)^2) 
    
    error_tau_did <- mean(abs(tau_estimate_did-delta_t)^2)
    
    error_tau_mc_nnm <- mean(abs(tau_estimate_mc_nnm-delta_t)^2)
    
    error_tau_sdid <- mean(abs(tau_estimate_sdid-delta_t)^2)
    
    error_tau_lapis <- mean(abs(tau_estimate_lapis-delta_t)^2)
    
    error_tau_oracle <- mean(abs(tau_estimate_oracle-delta_t)^2)

    errors_this_L_did[j] <- error_tau_did
    
    errors_this_L_sc[j] <- error_tau_sc
      
    errors_this_L_gsynth[j] <- error_tau_gsynth
      
    errors_this_L_mc_nnm[j] <- error_tau_mc_nnm
    
    errors_this_L_sdid[j] <- error_tau_sdid
    
    errors_this_L_lapis[j] <- error_tau_lapis
    
    errors_this_L_oracle[j] <- error_tau_oracle

  }
  
  prediction_error_matrix_did[i,] <- errors_this_L_did
  
  prediction_error_matrix_sc[i,] <- errors_this_L_sc
      
  prediction_error_matrix_gsynth[i, ] <- errors_this_L_gsynth
  
  prediction_error_matrix_mc_nnm[i, ] <- errors_this_L_mc_nnm
  
  prediction_error_matrix_sdid[i,] <- errors_this_L_sdid
  
  prediction_error_matrix_lapis[i,] <- errors_this_L_lapis
  
  prediction_error_matrix_oracle[i,] <- errors_this_L_oracle
      

}
  
  
    mse_and_se_of_mse_did <- mse_and_se_of_mse(prediction_error_matrix_did)
    
    mse_and_se_of_mse_sc <- mse_and_se_of_mse(prediction_error_matrix_sc)
                                            
    mse_and_se_of_mse_gsynth <- mse_and_se_of_mse(prediction_error_matrix_gsynth)
    
    mse_and_se_of_mse_mc_nnm <- mse_and_se_of_mse(prediction_error_matrix_mc_nnm)
    
    mse_and_se_of_mse_sdid <- mse_and_se_of_mse(prediction_error_matrix_sdid)
    
    mse_and_se_of_mse_lapis <- mse_and_se_of_mse(prediction_error_matrix_lapis)
  
    mse_and_se_of_mse_oracle <- mse_and_se_of_mse(prediction_error_matrix_oracle)
    
    mses_rho_did <- c(mses_rho_did, mse_and_se_of_mse_did[1])
    
    se_mses_rho_did <- c(se_mses_rho_did, mse_and_se_of_mse_did[2])
    
    mses_rho_sc <- c(mses_rho_sc, mse_and_se_of_mse_sc[1])
    
    se_mses_rho_sc <- c(se_mses_rho_sc, mse_and_se_of_mse_sc[2])
                                            
    mses_rho_gsynth <- c(mses_rho_gsynth, mse_and_se_of_mse_gsynth[1])
    
    se_mses_rho_gsynth <- c(se_mses_rho_gsynth, mse_and_se_of_mse_gsynth[2])                                                                          
    mses_rho_mc_nnm <- c(mses_rho_mc_nnm, mse_and_se_of_mse_mc_nnm[1])
    
    se_mses_rho_mc_nnm <- c(se_mses_rho_mc_nnm, mse_and_se_of_mse_mc_nnm[2])
    
    mses_rho_sdid <- c(mses_rho_sdid, mse_and_se_of_mse_sdid[1])
    
    se_mses_rho_sdid <- c(se_mses_rho_sdid, mse_and_se_of_mse_sdid[2])
    
    mses_rho_lapis <- c(mses_rho_lapis, mse_and_se_of_mse_lapis[1])
    
    se_mses_rho_lapis <- c(se_mses_rho_lapis, mse_and_se_of_mse_lapis[2])
    
    mses_rho_oracle <- c(mses_rho_oracle, mse_and_se_of_mse_oracle[1])
    
    se_mses_rho_oracle <- c(se_mses_rho_oracle, mse_and_se_of_mse_oracle[2])
  
  
  
  signal_to_noise_ratios <- c(signal_to_noise_ratios, (svd(L)$d[R]/                             svd(sigma_squared*this_autocorrelation_matrix)$d[1]))
                                            
 print(paste("Finished rho:", this_rho))

  
       
       }
  

# mses_rho_mc_nnm, mses_rho_did,

rho_data <- cbind(c( mses_rho_did, mses_rho_sc, mses_rho_gsynth, mses_rho_sdid, mses_rho_mc_nnm,
               mses_rho_lapis, mses_rho_oracle),
              
             c(se_mses_rho_did, se_mses_rho_sc,  mses_rho_gsynth, se_mses_rho_sdid, se_mses_rho_mc_nnm,
               se_mses_rho_lapis, se_mses_rho_oracle)
              
)


# , 'MC_NNM', 'DID', 
  
rho_sensitivity_data <- data.frame(rep(c('DID','SC', 'GSYNTH', 'SDID', "MC-NNM", 'LAPIS',
'ORACLE'), 
               each=length(mses_rho_lapis)))

names(rho_sensitivity_data) <- 'Method'

rho_sensitivity_data$rho <- all_rhos

rho_sensitivity_data$mse <- rho_data[,1]

rho_sensitivity_data$se <- rho_data[,2]
                                            
                                            
if (EXCLUDE_GSYNTH){
                                            
    
rho_sensitivity_data <- rho_sensitivity_data %>% filter(Method != 'GSYNTH')   
    
    
    }


p_snr_vs_rho <- (ggplot(NULL, aes(x=all_rhos, y=signal_to_noise_ratios)) + geom_line() + theme_bw()+ ggtitle("SNR as a Function of the Correlation Parameter"))

p_mse_vs_rho <- (ggplot(rho_sensitivity_data, aes(x=rho, y=sqrt(mse), col=Method)) + geom_line() + xlab("rho") + theme_bw()+ ggtitle("rmse as a Function of the Correlation Parameter"))

#p_mse_vs_rho <- (ggplot(rho_sensitivity_data, aes(x=rho, y=mse, col=Method)) + geom_line() + 
#                   geom_ribbon(aes(ymin=mse-1.9*se,
#                   ymax=mse+1.9*se, alpha=.1), fill = "grey70", lty=2) +
#                   theme_bw()+ ggtitle("mse as a Function of the Correlation Parameter"))


#p_mse_vs_rho
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
## Chunk 20

mses_tau_did <- c()

se_mses_tau_did <- c()

mses_tau_sc<- c()

se_mses_tau_sc <- c()

mses_tau_gsynth <- c()

se_mses_tau_gsynth <- c()

mses_tau_mc_nnm <- c()

se_mses_tau_mc_nnm <- c()

mses_tau_sdid <- c()

se_mses_tau_sdid <- c()

mses_tau_lapis <- c()

se_mses_tau_lapis <- c()

mses_tau_oracle <- c()

se_mses_tau_oracle <- c()

all_taus <- seq(1, 30, 2)

for (this_tau in all_taus){

  set.seed(3729)
  
  if (design=="staggered_adoption"){ ## Come up with a way to vary the lag in the staggered structure
  
  if(lag_structure == "random"){
    
    ones_we_make <- c(rep(0, N0), pmin(rpois(N-N0, 
                                            lambda=average_treatment_length-1)+1, 
                                      min(max_lag*(N-N0), .8*Time)))
    
  }else if (lag_structure=="constant"){ ## Does not control T-T0
    
    ones_we_make <- c(rep(0, N0), pmin(max_lag*seq(1, (N-N0)), floor(.8*Time)))
    
  }

}else if (design=="block_treatment"){
  
  ones_we_make <- c(rep(0, N0), rep(Time-Time0, N-N0))
  
}

  W <- W_maker(N=N, Time=Time, ones_per_row = ones_we_make)

  this_tau_matrix <- t(apply(W, MARGIN=1, FUN=treated_matrix_creator, f_of_t=delta_t_constant,
                        value=this_tau))

  treated_units <- as.numeric(which(apply(W, MARGIN=1, FUN = function(x) any(x==1))))
  
  treatment_times <- as.numeric(which(apply(W, MARGIN=2, FUN = function(x) any(x==1))))

  delta_t <- delta_t_constant(treatment_times-(min(treatment_times)-1), value=this_tau)
  
  treated_units <- as.numeric(which(apply(W, MARGIN=1, FUN = function(x) any(x==1))))
  
  treatment_times <- as.numeric(which(apply(W, MARGIN=2, FUN = function(x) any(x==1))))
  
  prediction_error_matrix_did <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_sc <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
                                            
  prediction_error_matrix_gsynth <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_mc_nnm <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

  prediction_error_matrix_sdid <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

  prediction_error_matrix_lapis <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_oracle <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  autocorrelation_matrix <- make_rho_mat(rho=rho_parameter, p=dim(W)[2])

  
  for (i in 1:number_of_L){
    
   # set.seed(3729)
    
    errors_this_L_did <- rep(NA, draws_per_L)
    
    errors_this_L_sc <- rep(NA, draws_per_L)
      
    errors_this_L_gsynth <- rep(NA, draws_per_L)
    
    errors_this_L_mc_nnm <- rep(NA, draws_per_L)
    
    errors_this_L_sdid <- rep(NA, draws_per_L)
    
    errors_this_L_lapis <- rep(NA, draws_per_L)
    
    errors_this_L_oracle <- rep(NA, draws_per_L)
    
    if (exchangable){
      
      U_vec <- rexp(n=N*R, rate=1)
  
      V_vec <- rexp(n=Time*R, rate=1)
    
      U <- matrix(U_vec, nrow=N, ncol=R, byrow=T)
  
      V <- matrix(V_vec, nrow=Time, ncol=R, byrow=T)
    
    }else{
      
      U <- matrix(NA, nrow=N, ncol=R, byrow=T)
  
      V <- matrix(NA, nrow=Time, ncol=R, byrow=T)
      
      for (row_unit in 1:N){
        
        U[row_unit,] <- rpois(n=R, lambda=sqrt(row_unit/N))
        
      } 
      
      for (row_time in 1:Time){
        
        V[row_time,] <- rpois(n=R, lambda=sqrt(row_time/Time))
        
      }
      
    }
  
    L <- L_scaling*(U %*% t(V))
    
    for (j in 1:draws_per_L){
      
      if (error == 'gaussian'){
        
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+this_tau_matrix*W, distribution='gaussian',
                   scalar_sigma=sqrt(sigma_squared))
        
      } else if (error == 't'){
      
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+this_tau_matrix*W,
                   distribution='t', scalar_sigma=sqrt(sigma_squared))
      
      } else if (error == 'poisson'){
        
        if (exchangable == F){
        
          L <- abs(L)+1
          
        }
        
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+this_tau_matrix*W,
                   distribution='poisson', scalar_sigma=sqrt(sigma_squared))
        
      } else if (error == 'scaled_gamma'){
        
        if (exchangable == F){
        
          L <- abs(L)+1
          
        }
        
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+this_tau_matrix*W,
                   distribution='scaled_gamma', scalar_sigma=sqrt(sigma_squared))
        
      }else if (error == 'exponential'){
        
        if (exchangable == F){
        
          L <- abs(L)+1
          
        }
        
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+this_tau_matrix*W,
                   distribution='exponential', scalar_sigma=sqrt(sigma_squared))
  
      }
      
      #estimated_rank <- rank_estimator(Y, W, num_iter=100, K=5, 
      #                      lambda_grid=c(0, 10^seq(-20, 0, 2)), 
      #                      method="threshold")
      
      #Y_0_LAPIS <- LAPIS(Y, rank_threshold=estimated_rank,
      #                                    min_iter=1, max_iter=max_iter,
      #                                    tolerance=tolerance, W=W)
      
      
      if (N-N0 > 1){
      
      treatment_subjects_averaged <- colMeans(Y[1:dim(Y)[1] > N0,])
      
      W_averaged <- colMeans(W[1:dim(W)[1] > N0,])
      
      new_Y <- rbind(Y[1:dim(Y)[1] <= N0,], treatment_subjects_averaged)
      
      new_W <- rbind(W[1:dim(Y)[1] <= N0,], W_averaged)
      
      } else {
      
        new_Y <- Y
      
        new_W <- W
        
        
      }
      
      tau_estimate_did <- DID(Y=Y, W=W)
      
      tau_estimate_sc <- synth_cont(Y=Y, W=W)
        
        
     meltedD <- melt(W) 

    names(meltedD) <- c('id', 'time', 'treated')

    meltedObservedData <- melt(Y)

    names(meltedObservedData) <- c('id', 'time', 'value')

    joinedDataForGsynth <- quiet(meltedD %>% inner_join(meltedObservedData,
                                                       by=c('id', 'time')))

      tau_estimate_sdid <- SDID_general(Y=Y, W=W,
                   iterations_for_coord_desc=100)
      
      mc_nnm_info <- matrix_completion_causal(Y=Y, W=W, num_iter=1000, K=5, 
                            lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                            tol=1e-04)
    
    L_mc_nnm <- mc_nnm_info$L_hat
    
    tau_estimate_mc_nnm <- treat.estimator(Y=Y, L.hat=L_mc_nnm, W=W)
        
        
    estFactors <- rankMatrix(mc_nnm_info$L_hat)[1]  
      
    gsynthInfo <- quiet(gsynth(value~treated, data=joinedDataForGsynth, 
                         index=c('id', 'time'), 
                         parallel = TRUE, r=estFactors))

    gsynthContEst <- gsynthInfo$att
    
    if (design=='block_treatment'){
      
    tau_estimate_gsynth <- gsynthContEst[(Time0+1):Time]
        
        }else{
                            
    tau_estimate_gsynth <- gsynthContEst[
        which(names(gsynthContEst)==1):length(gsynthContEst)]
        
    }
    
    tau_estimate_lapis <- LAPIS_with_rank_estimation(Y=Y, 
                           W=W, initial_rank=rankMatrix(mc_nnm_info$L_hat)[1],
                           tolerance=tolerance, 
                           min_iter=min_iter, max_iter=max_iter,   
                           mu_grid=NULL, warm_start=F, method = 'explicit_tau')
    
    
    # tau_estimate_mc_nnm
    
    ## Only oracle in the sense that we know L
    
    tau_estimate_oracle <- treat.estimator(Y=Y, L.hat=L, W=W)
      
      
      
      error_tau_sc <- mean(abs(tau_estimate_sc-this_tau)^2)
        
      error_tau_gsynth <- mean(abs(tau_estimate_gsynth-this_tau)^2)
      
      error_tau_did <- mean(abs(tau_estimate_did-this_tau)^2)
      
      error_tau_mc_nnm <- mean(abs(tau_estimate_mc_nnm-this_tau)^2)
      
      error_tau_sdid <- mean(abs(tau_estimate_sdid-this_tau)^2)
      
      error_tau_lapis <- mean(abs(tau_estimate_lapis-this_tau)^2)
      
      error_tau_oracle <- mean(abs(tau_estimate_oracle-this_tau)^2)
  
      errors_this_L_did[j] <- error_tau_did
      
      errors_this_L_sc[j] <- error_tau_sc
        
      errors_this_L_gsynth[j] <- error_tau_gsynth
      
      errors_this_L_mc_nnm[j] <- error_tau_mc_nnm
      
      errors_this_L_sdid[j] <- error_tau_sdid
      
      errors_this_L_lapis[j] <- error_tau_lapis
      
      errors_this_L_oracle[j] <- error_tau_oracle
  
    }
    
    prediction_error_matrix_did[i,] <- errors_this_L_did
    
    prediction_error_matrix_sc[i,] <- errors_this_L_sc
      
    prediction_error_matrix_gsynth[i,] <- errors_this_L_gsynth
    
    prediction_error_matrix_mc_nnm[i, ] <- errors_this_L_mc_nnm
    
    prediction_error_matrix_sdid[i,] <- errors_this_L_sdid
    
    prediction_error_matrix_lapis[i,] <- errors_this_L_lapis
    
    prediction_error_matrix_oracle[i,] <- errors_this_L_oracle

}

    
    mse_and_se_of_mse_did <- mse_and_se_of_mse(prediction_error_matrix_did)
    
    mse_and_se_of_mse_sc <- mse_and_se_of_mse(prediction_error_matrix_sc)
                                            
    mse_and_se_of_mse_gsynth <- mse_and_se_of_mse(prediction_error_matrix_gsynth)
    
    mse_and_se_of_mse_mc_nnm <- mse_and_se_of_mse(prediction_error_matrix_mc_nnm)
    
    mse_and_se_of_mse_sdid <- mse_and_se_of_mse(prediction_error_matrix_sdid)
    
    mse_and_se_of_mse_lapis <- mse_and_se_of_mse(prediction_error_matrix_lapis)
  
    mse_and_se_of_mse_oracle <- mse_and_se_of_mse(prediction_error_matrix_oracle)
    
    mses_tau_did <- c(mses_tau_did, mse_and_se_of_mse_did[1])
    
    se_mses_tau_did <- c(se_mses_tau_did, mse_and_se_of_mse_did[2])
    
    mses_tau_sc <- c(mses_tau_sc, mse_and_se_of_mse_sc[1])
    
    se_mses_tau_sc <- c(se_mses_tau_sc, mse_and_se_of_mse_sc[2])
                                            
                                            
    mses_tau_gsynth <- c(mses_tau_gsynth, mse_and_se_of_mse_gsynth[1])
    
    se_mses_tau_gsynth <- c(se_mses_tau_gsynth, mse_and_se_of_mse_gsynth[2])
    
    mses_tau_mc_nnm <- c(mses_tau_mc_nnm, mse_and_se_of_mse_mc_nnm[1])
    
    se_mses_tau_mc_nnm <- c(se_mses_tau_mc_nnm, mse_and_se_of_mse_mc_nnm[2])
    
    mses_tau_sdid <- c(mses_tau_sdid, mse_and_se_of_mse_sdid[1])
    
    se_mses_tau_sdid <- c(se_mses_tau_sdid, mse_and_se_of_mse_sdid[2])
    
    mses_tau_lapis <- c(mses_tau_lapis, mse_and_se_of_mse_lapis[1])
    
    se_mses_tau_lapis <- c(se_mses_tau_lapis, mse_and_se_of_mse_lapis[2])
    
    mses_tau_oracle <- c(mses_tau_oracle, mse_and_se_of_mse_oracle[1])
    
    se_mses_tau_oracle <- c(se_mses_tau_oracle, mse_and_se_of_mse_oracle[2])
                                            
                                            
    print(paste("Finished tau:", this_tau))
  
  
  

}



# mses_tau_mc_nnm, mses_tau_did,


tau_data <- cbind(c( mses_tau_did, mses_tau_sc, mses_tau_gsynth, mses_tau_mc_nnm,
              mses_tau_sdid, mses_tau_lapis, mses_tau_oracle),
              
             c(se_mses_tau_did, se_mses_tau_sc, se_mses_tau_gsynth, se_mses_tau_mc_nnm,
              se_mses_tau_sdid, se_mses_tau_lapis, se_mses_tau_oracle)
              
)

tau_sensitivity_data <- data.frame(rep(c('DID','SC', 'GSYNTH', "MC-NNM",'SDID', 'LAPIS',
'ORACLE'), 
               each=length(mses_tau_lapis)))

names(tau_sensitivity_data) <- 'Method'

tau_sensitivity_data$tau <- all_taus

tau_sensitivity_data$mse <- tau_data[,1]

tau_sensitivity_data$se <- tau_data[,2]
                                            
                                            
if (EXCLUDE_GSYNTH){
                                            
tau_sensitivity_data <- tau_sensitivity_data %>% filter(Method != 'GSYNTH')   
    
    }

p_mse_vs_tau <- (ggplot(tau_sensitivity_data, aes(x=tau, y=sqrt(mse), col=Method)) + geom_line() +xlab("tau") + theme_bw()+ ggtitle("rmse as a Function of the Effect Size"))

#p_mse_vs_tau <- (ggplot(tau_sensitivity_data, aes(x=tau, y=mse, col=Method)) + geom_line() + 
#                   geom_ribbon(aes(ymin=mse-1.9*se,
#                   ymax=mse+1.9*se, alpha=.1), fill = "grey70", lty=2) +
#                   theme_bw()+ ggtitle("mse as a Function of the Effect Size"))

                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
## Chunk 21

mses_rank_did <- c()

se_mses_rank_did <- c()

mses_rank_sc<- c()

se_mses_rank_sc <- c()

mses_rank_gsynth <- c()

se_mses_rank_gsynth <- c()

mses_rank_mc_nnm <- c()

se_mses_rank_mc_nnm <- c()

mses_rank_sdid <- c()

se_mses_rank_sdid <- c()

mses_rank_lapis <- c()

se_mses_rank_lapis <- c()

mses_rank_oracle <- c()

se_mses_rank_oracle <- c()

all_ranks <- seq(2, 20, 2)

for (rank in all_ranks){
  
  set.seed(3729)
  
  if (design=="staggered_adoption"){ ## Come up with a way to vary the lag in the staggered structure
  
  if(lag_structure == "random"){
    
    ones_we_make <- c(rep(0, N0), pmin(rpois(N-N0, 
                                            lambda=average_treatment_length-1)+1, 
                                      min(max_lag*(N-N0), .8*Time)))
    
  }else if (lag_structure=="constant"){ ## Does not control T-T0
    
    ones_we_make <- c(rep(0, N0), pmin(max_lag*seq(1, (N-N0)), floor(.8*Time)))
    
  }

}else if (design=="block_treatment"){
  
  ones_we_make <- c(rep(0, N0), rep(Time-Time0, N-N0))
  
}

  W <- W_maker(N=N, Time=Time, ones_per_row = ones_we_make)

  tau_matrix <- t(apply(W, MARGIN=1, FUN=treated_matrix_creator, 
                      f_of_t=treatment_function, arg_max=arg_max, 
                      y_max=y_max, halfway_time=halfway_time, cutoff=cutoff))

  treated_units <- as.numeric(which(apply(W, MARGIN=1, FUN = function(x) any(x==1))))
  
  treatment_times <- as.numeric(which(apply(W, MARGIN=2, FUN = function(x) any(x==1))))

  delta_t <- treatment_function(treatment_times-(min(treatment_times)-1),
                              arg_max=arg_max, 
                      y_max=y_max, halfway_time=halfway_time, cutoff=cutoff, value=tau)
  
  prediction_error_matrix_did <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_sc <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_gsynth <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)                                          
                                            
  prediction_error_matrix_mc_nnm <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

  prediction_error_matrix_sdid <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

  prediction_error_matrix_lapis <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_oracle <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  autocorrelation_matrix <- make_rho_mat(rho=rho_parameter, p=dim(W)[2])
  
for (i in 1:number_of_L){
    
   # set.seed(3729)
    
    errors_this_L_did <- rep(NA, draws_per_L)
    
    errors_this_L_sc <- rep(NA, draws_per_L)
    
    errors_this_L_gsynth <- rep(NA, draws_per_L)
    
    errors_this_L_mc_nnm <- rep(NA, draws_per_L)
    
    errors_this_L_sdid <- rep(NA, draws_per_L)
    
    errors_this_L_lapis <- rep(NA, draws_per_L)
    
    errors_this_L_oracle <- rep(NA, draws_per_L)
    
    if (exchangable){
      
      U_vec <- rexp(n=N*R, rate=1)
  
      V_vec <- rexp(n=Time*R, rate=1)
    
      U <- matrix(U_vec, nrow=N, ncol=rank, byrow=T)
  
      V <- matrix(V_vec, nrow=Time, ncol=rank, byrow=T)
    
    }else{
      
      U <- matrix(NA, nrow=N, ncol=rank, byrow=T)
  
      V <- matrix(NA, nrow=Time, ncol=rank, byrow=T)
      
      for (row_unit in 1:N){
        
        U[row_unit,] <- rpois(n=rank, lambda=sqrt(row_unit/N))
        
      } 
      
      for (row_time in 1:Time){
        
        V[row_time,] <- rpois(n=rank, lambda=sqrt(row_time/Time))
        
      }
      
    }
  
    L <- L_scaling*(U %*% t(V))
    
    for (j in 1:draws_per_L){
      
      if (error == 'gaussian'){
        
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+tau*W, distribution='gaussian',
                   scalar_sigma=sqrt(sigma_squared))
        
      } else if (error == 't'){
      
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+tau*W,
                   distribution='t', scalar_sigma=sqrt(sigma_squared))
      
      } else if (error == 'poisson'){
        
        if (exchangable == F){
        
          L <- abs(L)+1
          
        }
        
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+tau*W,
                   distribution='poisson', scalar_sigma=sqrt(sigma_squared))
        
      } else if (error == 'scaled_gamma'){
        
        if (exchangable == F){
        
          L <- abs(L)+1
          
        }
        
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+tau*W,
                   distribution='scaled_gamma', scalar_sigma=sqrt(sigma_squared))
        
      }else if (error == 'exponential'){
        
        if (exchangable == F){
        
          L <- abs(L)+1
          
        }
        
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+tau*W,
                   distribution='exponential', scalar_sigma=sqrt(sigma_squared))
  
      }
      
      #estimated_rank <- rank_estimator(Y, W, num_iter=100, K=5, 
      #                      lambda_grid=c(0, 10^seq(-20, 0, 2)), 
      #                      method="threshold")
      
      #Y_0_LAPIS <- LAPIS(Y, rank_threshold=estimated_rank,
      #                                    min_iter=1, max_iter=max_iter,
      #                                    tolerance=tolerance, W=W)
      
      
      if (N-N0 > 1){
      
      treatment_subjects_averaged <- colMeans(Y[1:dim(Y)[1] > N0,])
      
      W_averaged <- colMeans(W[1:dim(W)[1] > N0,])
      
      new_Y <- rbind(Y[1:dim(Y)[1] <= N0,], treatment_subjects_averaged)
      
      new_W <- rbind(W[1:dim(Y)[1] <= N0,], W_averaged)
      
      } else {
      
        new_Y <- Y
      
        new_W <- W
        
        
      }
      
      tau_estimate_did <- DID(Y=Y, W=W)
      
      tau_estimate_sc <- synth_cont(Y=Y, W=W)
        
        
        
         meltedD <- melt(W) 

    names(meltedD) <- c('id', 'time', 'treated')

    meltedObservedData <- melt(Y)

    names(meltedObservedData) <- c('id', 'time', 'value')

    joinedDataForGsynth <- quiet(meltedD %>% inner_join(meltedObservedData,
                                                       by=c('id', 'time'))       )

      tau_estimate_sdid <- SDID_general(Y=Y, W=W,
                   iterations_for_coord_desc=100)
      
    mc_nnm_info <- matrix_completion_causal(Y=Y, W=W, num_iter=1000, K=5, 
                            lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                            tol=1e-04)
    
    L_mc_nnm <- mc_nnm_info$L_hat
    
    tau_estimate_mc_nnm <- treat.estimator(Y=Y, L.hat=L_mc_nnm, W=W)
        
    estFactors <- rankMatrix(mc_nnm_info$L_hat)[1]  
      
    gsynthInfo <- quiet(gsynth(value~treated, data=joinedDataForGsynth, index=c('id', 'time'), 
                         parallel = TRUE, r=estFactors))

    gsynthContEst <- gsynthInfo$att
    
    if (design=='block_treatment'){
      
        tau_estimate_gsynth <- gsynthContEst[(Time0+1):Time]
        
        }else{
                            
        tau_estimate_gsynth <- gsynthContEst[
        which(names(gsynthContEst)==1):length(gsynthContEst)]
        
    }
    
    tau_estimate_lapis <- LAPIS_with_rank_estimation(Y=Y, 
                           W=W, initial_rank=rankMatrix(mc_nnm_info$L_hat)[1],
                           tolerance=tolerance, 
                           min_iter=min_iter, max_iter=max_iter,   
                           mu_grid=NULL, warm_start=F, method = 'explicit_tau')
    
    tau_estimate_oracle <- treat.estimator(Y=Y, L.hat=L, W=W)
    
    ## Only oracle in the sense that we know L
      
      
      error_tau_sc <- mean(abs(tau_estimate_sc-delta_t)^2)
        
      error_tau_gsynth <- mean(abs(tau_estimate_gsynth-delta_t)^2)
      
      error_tau_did <- mean(abs(tau_estimate_did-delta_t)^2)
      
      error_tau_mc_nnm <- mean(abs(tau_estimate_mc_nnm-delta_t)^2)
      
      error_tau_sdid <- mean(abs(tau_estimate_sdid-delta_t)^2)
      
      error_tau_lapis <- mean(abs(tau_estimate_lapis-delta_t)^2)
      
      error_tau_oracle <- mean(abs(tau_estimate_oracle-delta_t)^2)
  
      errors_this_L_did[j] <- error_tau_did
      
      errors_this_L_sc[j] <- error_tau_sc
        
      errors_this_L_gsynth[j] <- error_tau_gsynth
      
      errors_this_L_mc_nnm[j] <- error_tau_mc_nnm
      
      errors_this_L_sdid[j] <- error_tau_sdid
      
      errors_this_L_lapis[j] <- error_tau_lapis
      
      errors_this_L_oracle[j] <- error_tau_oracle
  
    }
    
    prediction_error_matrix_did[i,] <- errors_this_L_did
    
    prediction_error_matrix_sc[i,] <- errors_this_L_sc
    
    prediction_error_matrix_gsynth[i,] <- errors_this_L_gsynth
    
    prediction_error_matrix_mc_nnm[i, ] <- errors_this_L_mc_nnm
    
    prediction_error_matrix_sdid[i,] <- errors_this_L_sdid
    
    prediction_error_matrix_lapis[i,] <- errors_this_L_lapis
    
    prediction_error_matrix_oracle[i,] <- errors_this_L_oracle

}

    
    mse_and_se_of_mse_did <- mse_and_se_of_mse(prediction_error_matrix_did)
    
    mse_and_se_of_mse_sc <- mse_and_se_of_mse(prediction_error_matrix_sc)
                                            
    mse_and_se_of_mse_gsynth <- mse_and_se_of_mse(prediction_error_matrix_gsynth)
    
    mse_and_se_of_mse_mc_nnm <- mse_and_se_of_mse(prediction_error_matrix_mc_nnm)
    
    mse_and_se_of_mse_sdid <- mse_and_se_of_mse(prediction_error_matrix_sdid)
    
    mse_and_se_of_mse_lapis <- mse_and_se_of_mse(prediction_error_matrix_lapis)
  
    mse_and_se_of_mse_oracle <- mse_and_se_of_mse(prediction_error_matrix_oracle)
    
    mses_rank_did <- c(mses_rank_did, mse_and_se_of_mse_did[1])
    
    se_mses_rank_did <- c(se_mses_rank_did, mse_and_se_of_mse_did[2])
    
    mses_rank_sc <- c(mses_rank_sc, mse_and_se_of_mse_sc[1])
    
    se_mses_rank_sc <- c(se_mses_rank_sc, mse_and_se_of_mse_sc[2])
                                            
    mses_rank_gsynth <- c(mses_rank_gsynth, mse_and_se_of_mse_gsynth[1])
    
    se_mses_rank_gsynth <- c(se_mses_rank_gsynth, mse_and_se_of_mse_gsynth[2])
    
    mses_rank_mc_nnm <- c(mses_rank_mc_nnm, mse_and_se_of_mse_mc_nnm[1])
    
    se_mses_rank_mc_nnm <- c(se_mses_rank_mc_nnm, mse_and_se_of_mse_mc_nnm[2])
    
    mses_rank_sdid <- c(mses_rank_sdid, mse_and_se_of_mse_sdid[1])
    
    se_mses_rank_sdid <- c(se_mses_rank_sdid, mse_and_se_of_mse_sdid[2])
    
    mses_rank_lapis <- c(mses_rank_lapis, mse_and_se_of_mse_lapis[1])
    
    se_mses_rank_lapis <- c(se_mses_rank_lapis, mse_and_se_of_mse_lapis[2])
    
    mses_rank_oracle <- c(mses_rank_oracle, mse_and_se_of_mse_oracle[1])
    
    se_mses_rank_oracle <- c(se_mses_rank_oracle, mse_and_se_of_mse_oracle[2])
                                            
    print(paste("Finished rank:", rank))
}

rank_data <- cbind(c( mses_rank_did, mses_rank_sc, mses_rank_gsynth, mses_rank_mc_nnm,
              mses_rank_sdid, mses_rank_lapis, mses_rank_oracle),
              
             c(se_mses_rank_did, se_mses_rank_sc, se_mses_rank_gsynth, se_mses_rank_mc_nnm,
              se_mses_rank_sdid, se_mses_rank_lapis, se_mses_rank_oracle)
              
)


rank_sensitivity_data <- data.frame(rep(c('DID','SC', "GSYNTH", "MC-NNM",'SDID', 'LAPIS',
'ORACLE'), 
               each=length(mses_rank_lapis)))

names(rank_sensitivity_data) <- 'Method'

rank_sensitivity_data$rank <- all_ranks

rank_sensitivity_data$mse <- rank_data[,1]

rank_sensitivity_data$se <- rank_data[,2]
  
                                            

if (EXCLUDE_GSYNTH){
                                            
rank_sensitivity_data <- rank_sensitivity_data %>% filter(Method != 'GSYNTH')   
    
    
    }

p_mse_vs_rank <- (ggplot(rank_sensitivity_data, aes(x=rank, y=sqrt(mse), col=Method)) + geom_line() + theme_bw()+ ggtitle("rmse as a Function of the True Rank"))

                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
## Chunk 21

# R+rank_error

mses_rank_error_did <- c()

se_mses_rank_error_did <- c()

mses_rank_error_sc<- c()

se_mses_rank_error_sc <- c()

mses_rank_error_gsynth <- c()

se_mses_rank_error_gsynth <- c()

mses_rank_error_mc_nnm <- c()

se_mses_rank_error_mc_nnm <- c()

mses_rank_error_sdid <- c()

se_mses_rank_error_sdid <- c()

mses_rank_error_lapis <- c()

se_mses_rank_error_lapis <- c()

mses_rank_error_oracle <- c()

se_mses_rank_error_oracle <- c()

all_rank_errors <- seq(max(R-10, R-(R-1)), R+(10), 2)-R

for (rank_error in all_rank_errors){
  
  print(rank_error)
  
  set.seed(3729)
  
  if (design=="staggered_adoption"){ ## Come up with a way to vary the lag in the staggered structure
  
  if(lag_structure == "random"){
    
    ones_we_make <- c(rep(0, N0), pmin(rpois(N-N0, 
                                            lambda=average_treatment_length-1)+1, 
                                      min(max_lag*(N-N0), .8*Time)))
    
  }else if (lag_structure=="constant"){ ## Does not control T-T0
    
    ones_we_make <- c(rep(0, N0), pmin(max_lag*seq(1, (N-N0)), floor(.8*Time)))
    
  }

}else if (design=="block_treatment"){
  
  ones_we_make <- c(rep(0, N0), rep(Time-Time0, N-N0))
  
}

  W <- W_maker(N=N, Time=Time, ones_per_row = ones_we_make)

  tau_matrix <- t(apply(W, MARGIN=1, FUN=treated_matrix_creator, 
                      f_of_t=treatment_function, arg_max=arg_max, 
                      y_max=y_max, halfway_time=halfway_time, cutoff=cutoff))

  treated_units <- as.numeric(which(apply(W, MARGIN=1, FUN = function(x) any(x==1))))
  
  treatment_times <- as.numeric(which(apply(W, MARGIN=2, FUN = function(x) any(x==1))))

  delta_t <- treatment_function(treatment_times-(min(treatment_times)-1),
                              arg_max=arg_max, 
                      y_max=y_max, halfway_time=halfway_time, cutoff=cutoff)
  
  prediction_error_matrix_did <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_sc <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
                                            
  prediction_error_matrix_gsynth <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_mc_nnm <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

  prediction_error_matrix_sdid <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

  prediction_error_matrix_lapis <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_oracle <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  autocorrelation_matrix <- make_rho_mat(rho=rho_parameter, p=dim(W)[2])
  
  for (i in 1:number_of_L){
    
    errors_this_L_did <- rep(NA, draws_per_L)
  
    errors_this_L_sc <- rep(NA, draws_per_L)
      
     errors_this_L_gsynth <- rep(NA, draws_per_L)
  
    errors_this_L_mc_nnm <- rep(NA, draws_per_L)
  
    errors_this_L_sdid <- rep(NA, draws_per_L)
  
    errors_this_L_lapis <- rep(NA, draws_per_L)
    
    errors_this_L_oracle <- rep(NA, draws_per_L)
    
    if (exchangable){
      
      U_vec <- rexp(n=N*R, rate=1)
  
      V_vec <- rexp(n=Time*R, rate=1)
    
      U <- matrix(U_vec, nrow=N, ncol=R, byrow=T)
  
      V <- matrix(V_vec, nrow=Time, ncol=R, byrow=T)
    
    }else{
      
      U <- matrix(NA, nrow=N, ncol=R, byrow=T)
  
      V <- matrix(NA, nrow=Time, ncol=R, byrow=T)
      
      for (row_unit in 1:N){
        
        U[row_unit,] <- rpois(n=R, lambda=sqrt(row_unit/N))
        
      } 
      
      for (row_time in 1:Time){
        
        V[row_time,] <- rpois(n=R, lambda=sqrt(row_time/Time))
        
      }
      
    }
  
    L <- L_scaling*(U %*% t(V))
    
    for (j in 1:draws_per_L){
      
        if (error == 'gaussian'){
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W, distribution='gaussian',
                 scalar_sigma=sqrt(sigma_squared))
      
    } else if (error == 't'){
    
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='t', scalar_sigma=1, df=df)
    
    } else if (error == 'poisson'){
      
      if (exchangable == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='poisson', scalar_sigma=1)
      
    } else if (error == 'scaled_gamma'){
      
      if (exchangable == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='scaled_gamma', scalar_sigma=1)
      
    }else if (error == 'exponential'){
      
      if (exchangable == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='exponential', scalar_sigma=1)

    }
      
   # estimated_rank <- rank_estimator(Y, W, num_iter=100, K=5, 
  #                        lambda_grid=c(0, 10^seq(-20, 0, 2)), 
  #                        method="threshold")
      # estimated_rank <- rank_estimator(Y, W, num_iter=100, K=5, 
    #                      lambda_grid=c(0, 10^seq(-20, 0, 2)), 
    #                      method="threshold")
      
      
    if (N-N0 > 1){
    
      treatment_subjects_averaged <- colMeans(Y[1:dim(Y)[1] > N0,])
    
      W_averaged <- colMeans(W[1:dim(W)[1] > N0,])
    
      new_Y <- rbind(Y[1:dim(Y)[1] <= N0,], treatment_subjects_averaged)
    
      new_W <- rbind(W[1:dim(Y)[1] <= N0,], W_averaged)
    
    } else {
    
      new_Y <- Y
    
      new_W <- W
      
      
    }
    
    tau_estimate_did <- DID(Y=Y, W=W)
    
    tau_estimate_sc <- synth_cont(Y=Y, W=W)
        
     meltedD <- melt(W) 

    names(meltedD) <- c('id', 'time', 'treated')

    meltedObservedData <- melt(Y)

    names(meltedObservedData) <- c('id', 'time', 'value')

    joinedDataForGsynth <- quiet(meltedD %>% inner_join(meltedObservedData,
                                                       by=c('id', 'time')))

    tau_estimate_sdid <- SDID_general(Y=Y, W=W,
                 iterations_for_coord_desc=100)
    
    mc_nnm_info <- matrix_completion_causal(Y=Y, W=W, num_iter=1000, K=5, 
                            lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                            tol=1e-04)
    
    L_mc_nnm <- mc_nnm_info$L_hat
    
    tau_estimate_mc_nnm <- treat.estimator(Y=Y, L.hat=L_mc_nnm, W=W)
    
        
    estFactors <- rankMatrix(mc_nnm_info$L_hat)[1]  
      
    gsynthInfo <- quiet(gsynth(value~treated, data=joinedDataForGsynth, index=c('id', 'time'), 
                         parallel = TRUE, r=estFactors))

    gsynthContEst <- gsynthInfo$att
    
    if (design=='block_treatment'){
      
        tau_estimate_gsynth <- gsynthContEst[(Time0+1):Time]
        
        }else{
                            
        tau_estimate_gsynth <- gsynthContEst[
            which(names(gsynthContEst)==1):length(gsynthContEst)]
        
    }
    
    tau_estimate_lapis <- LAPIS_with_rank_estimation(Y=Y, 
                           W=W, initial_rank=rankMatrix(mc_nnm_info$L_hat)[1],
                           tolerance=tolerance, 
                           min_iter=min_iter, max_iter=max_iter,   
                           mu_grid=NULL, warm_start=F, method = 'explicit_tau')
    
    tau_estimate_oracle <- treat.estimator(Y=Y, L.hat=L, W=W)
    
    
    # tau_estimate_mc_nnm
    
    ## Only oracle in the sense that we know L
    
    
    error_rank_error_sc <- mean(abs(tau_estimate_sc-delta_t)^2)
        
    error_rank_error_gsynth <- mean(abs(tau_estimate_gsynth-delta_t)^2)
    
    error_rank_error_did <- mean(abs(tau_estimate_did-delta_t)^2)
    
    error_rank_error_mc_nnm <- mean(abs(tau_estimate_mc_nnm-delta_t)^2)
    
    error_rank_error_sdid <- mean(abs(tau_estimate_sdid-delta_t)^2)
    
    error_rank_error_lapis <- mean(abs(tau_estimate_lapis-delta_t)^2)
    
    error_rank_error_oracle <- mean(abs(tau_estimate_oracle-delta_t)^2)

    errors_this_L_did[j] <- error_rank_error_did
    
    errors_this_L_sc[j] <- error_rank_error_sc
        
    errors_this_L_gsynth[j] <- error_rank_error_gsynth
    
    errors_this_L_mc_nnm[j] <- error_rank_error_mc_nnm
    
    errors_this_L_sdid[j] <- error_rank_error_sdid
    
    errors_this_L_lapis[j] <- error_rank_error_lapis
    
    errors_this_L_oracle[j] <- error_rank_error_oracle

    }
  
    prediction_error_matrix_did[i,] <- errors_this_L_did
    
    prediction_error_matrix_sc[i,] <- errors_this_L_sc
    
    prediction_error_matrix_gsynth[i,] <- errors_this_L_gsynth
    
    prediction_error_matrix_mc_nnm[i,] <- errors_this_L_mc_nnm
    
    prediction_error_matrix_sdid[i,] <- errors_this_L_sdid
  
    prediction_error_matrix_lapis[i,] <- errors_this_L_lapis
    
    prediction_error_matrix_oracle[i,] <- errors_this_L_oracle

  }

    
    mse_and_se_of_mse_did <- mse_and_se_of_mse(prediction_error_matrix_did)
    
    mse_and_se_of_mse_sc <- mse_and_se_of_mse(prediction_error_matrix_sc)
                                            
    mse_and_se_of_mse_gsynth <- mse_and_se_of_mse(prediction_error_matrix_gsynth)
    
    mse_and_se_of_mse_mc_nnm <- mse_and_se_of_mse(prediction_error_matrix_mc_nnm)
    
    mse_and_se_of_mse_sdid <- mse_and_se_of_mse(prediction_error_matrix_sdid)
    
    mse_and_se_of_mse_lapis <- mse_and_se_of_mse(prediction_error_matrix_lapis)
  
    mse_and_se_of_mse_oracle <- mse_and_se_of_mse(prediction_error_matrix_oracle)
    
    mses_rank_error_did <- c(mses_rank_error_did, mse_and_se_of_mse_did[1])
    
    se_mses_rank_error_did <- c(se_mses_rank_error_did, mse_and_se_of_mse_did[2])
    
    mses_rank_error_sc <- c(mses_rank_error_sc, mse_and_se_of_mse_sc[1])
    
    se_mses_rank_error_sc <- c(se_mses_rank_error_sc, mse_and_se_of_mse_sc[2])
                                            
    mses_rank_error_gsynth <- c(mses_rank_error_gsynth, mse_and_se_of_mse_gsynth[1])
    
    se_mses_rank_error_gsynth <- c(se_mses_rank_error_gsynth, mse_and_se_of_mse_gsynth[2])
    
    mses_rank_error_mc_nnm <- c(mses_rank_error_mc_nnm, mse_and_se_of_mse_mc_nnm[1])
    
    se_mses_rank_error_mc_nnm <- c(se_mses_rank_error_mc_nnm, mse_and_se_of_mse_mc_nnm[2])
    
    mses_rank_error_sdid <- c(mses_rank_error_sdid, mse_and_se_of_mse_sdid[1])
    
    se_mses_rank_error_sdid <- c(se_mses_rank_error_sdid, mse_and_se_of_mse_sdid[2])
    
    mses_rank_error_lapis <- c(mses_rank_error_lapis, mse_and_se_of_mse_lapis[1])
    
    se_mses_rank_error_lapis <- c(se_mses_rank_error_lapis, mse_and_se_of_mse_lapis[2])
    
    mses_rank_error_oracle <- c(mses_rank_error_oracle, mse_and_se_of_mse_oracle[1])
    
    se_mses_rank_error_oracle <- c(se_mses_rank_error_oracle, mse_and_se_of_mse_oracle[2])

    print(paste("Finished rank error:", rank_error))
}


rank_error_data <- cbind(c( mses_rank_error_did, mses_rank_error_sc, mses_rank_error_gsynth,mses_rank_error_mc_nnm,
              mses_rank_error_sdid, mses_rank_error_lapis, mses_rank_error_oracle),
              
             c(se_mses_rank_error_did, se_mses_rank_error_sc, se_mses_rank_error_gsynth, se_mses_rank_error_mc_nnm,
              se_mses_rank_error_sdid, se_mses_rank_error_lapis, se_mses_rank_error_oracle)
              
)



## mses_rank_error_mc_nnm, mses_rank_error_did,

rank_error_sensitivity_data <- data.frame(rep(c('DID','SC', 'GSYNTH', "MC-NNM",'SDID', 'LAPIS',
'ORACLE'), 
               each=length(mses_rank_error_lapis)))

names(rank_error_sensitivity_data) <- 'Method'

rank_error_sensitivity_data$rank_error <- all_rank_errors

rank_error_sensitivity_data$mse <- rank_error_data[,1]

rank_error_sensitivity_data$se <- rank_error_data[,2]
   
                                            
if (EXCLUDE_GSYNTH){
                                            
rank_error_sensitivity_data <- rank_error_sensitivity_data %>% filter(Method != 'GSYNTH')   
                                            
}                                            
                                            
                                            
                                            

p_mse_vs_rank_error <- (ggplot(rank_error_sensitivity_data, aes(x=rank_error, y=sqrt(mse), col=Method)) + geom_line() + 
                   theme_bw()+ ggtitle("rmse as a Function of the Initial Rank Error"))

#p_mse_vs_rank_error <- (ggplot(rank_error_sensitivity_data, aes(x=rank_error, y=mse, #col=Method)) + geom_ribbon(aes(ymin=mse-1.9*se, ymax=mse+1.9*se, alpha=.1), fill = "grey70", #lty=2)+ geom_line()  + theme_bw()+ ggtitle("mse as a Function of the Initial Rank Error"))

                                            
                                            
                                            
server_name <- Sys.info()['nodename']                                        
                                            
                                            
sim_error_date_directory <- paste("../reports/", server_name,'/',
error, "_simulations","/" , Sys.Date(), "_simulations", sep='')

if (!dir.exists(sim_error_date_directory)){
  
  dir.create(path=sim_error_date_directory, recursive = TRUE)
  
}                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
full_output_directory <- sim_error_date_directory

image_directory <- paste(full_output_directory, "/simulation_plots", sep="")

if (!dir.exists(image_directory)){

dir.create(image_directory)
  
}                                            
                                            
                                            
                                            
                                            
if (exists("effect_plot")){

ggsave(paste(image_directory,'/effect_size_plot.pdf' , sep=''), effect_plot,
       width=5, height=3)
}

if (exists("p_mse_vs_tau")){

ggsave(paste(image_directory,'/varied_tau.pdf' , sep=''), p_mse_vs_tau,
       width=5, height=3)
}                                            
                                            
                                            
                                            
                                            
                                            
                                            
 if (exists("p_mse_vs_rho")){
  
ggsave(paste(image_directory,'/varied_rho.pdf' , sep=''), p_mse_vs_rho,
       width=5, height=3)
  


}


if (exists("p_snr_vs_rho")){
  
ggsave(paste(image_directory,'/snr_vs_rho.pdf' , sep=''), p_snr_vs_rho,
       width=5, height=3)
  


}
                                           
                                            
                                            
                                            
                                            
                                            
                                            
                                            
                                            
if (exists("p_mse_vs_block_size")){

ggsave(paste(image_directory,'/varied_block_size.pdf' , sep=''), p_mse_vs_block_size, width=5, height=3)
  
} 

if (exists("p_mse_vs_N0")){

ggsave(paste(image_directory,'/varied_N0.pdf' , sep=''), p_mse_vs_N0, width=5, height=3)
  
}

if (exists("p_mse_vs_rank_error")){

ggsave(paste(image_directory,'/varied_rank_error.pdf' , sep=''), p_mse_vs_rank_error, width=5, height=3)
  
}

if (exists("p_mse_vs_rank")){

ggsave(paste(image_directory,'/varied_rank.pdf' , sep=''), p_mse_vs_rank, width=5, height=3)
  
}                                            
                                            
                                            
makeFilename <- function(directory, fileName){
    
    paste(directory, fileName, sep='/')
}                                            
                                            
fileConn<-file(makeFilename(image_directory, fileName='fixed_parameter_error_table.txt'))
writeLines(tableForPresenting , fileConn)
close(fileConn)                                            
                                            
                                           
                                            