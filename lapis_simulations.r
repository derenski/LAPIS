library(ggplot2)
library(tidyverse)
library(LaplacesDemon)
library(glmnet)
library(foreach)
library(doParallel)
library(quadprog)
library(openxlsx)
library(plyr)
library(dplyr)
library(cubelyr)

max_available_clusters <- detectCores()-1
  
desired_clusters <- 4
  
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

makeFilename <- function(directory, fileName){
    
    paste(directory, fileName, sep='/')
}                            



options(dplyr.summarise.inform = FALSE)

methodNames <- c('DID','SC', 'SDID', "MC-NNM", 'LAPIS', 'ORACLE')

params <- read.xlsx("parameters_and_descriptions.xlsx", sheet = "parameter_data", rowNames=T)

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

balanced <- as.logical(params['balanced', ])

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





msesFixedParameterData <- array(NA, dim=c(length(methodNames), draws_per_L, number_of_L ),
                   dimnames=c(list(Method=methodNames), list(Iteration=1:draws_per_L), list(L_number=1:number_of_L)))

if (design=="staggered_adoption"){ ## Come up with a way to vary the lag in the staggered structure
  
  if(lag_structure == "random"){
    
    ones_we_make <- c(rep(0, N0), pmin(rpois(N-N0, 
                                            lambda=average_treatment_length-1)+1, 
                                      min(max_lag*(N-N0), .8*Time)))
    
  }else if (lag_structure=="constant"){ ## Does not control T-T0
    
    ones_we_make <- c(rep(0, N0), pmin(max_lag*seq(1, (N-N0)), floor(.8*Time)))
    
  }

}else if (design=="simultaneous_adoption"){
  
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

autocorrelation_matrix <- make_rho_mat(rho=rho_parameter, p=dim(W)[2])

sig_to_noise_ratios <- c()

for (i in 1:number_of_L){
  
  if (balanced){
    
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
      
    startTime <- Sys.time()
    
    if (error == 'gaussian'){
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W, distribution='gaussian',
                 scalar_sigma=sqrt(sigma_squared))
      
    } else if (error == 't'){
    
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='t', scalar_sigma=sqrt(sigma_squared))
    
    } else if (error == 'poisson'){
      
      if (balanced == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='poisson', scalar_sigma=sqrt(sigma_squared))
      
    } else if (error == 'scaled_gamma'){
      
      if (balanced == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='scaled_gamma', scalar_sigma=sqrt(sigma_squared))
      
    }else if (error == 'exponential'){
      
      if (balanced == F){
      
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

      
    
    tau_estimate_did <- DID(Y=Y, W=W)
    
    tau_estimate_sc <- synth_cont(Y=Y, W=W)
    
    tau_estimate_sdid <- SDID_general(Y=Y, W=W,
                 iterations_for_coord_desc=100)
    
    mc_nnm_info <- matrix_completion_causal(Y=Y, W=W, num_iter=1000, K=5, 
                            lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                            tol=1e-04)
    
    L_mc_nnm <- mc_nnm_info$L_hat
    
    tau_estimate_mc_nnm <- treat.estimator(Y=Y, L.hat=L_mc_nnm, W=W)
    
    estFactors <- rankMatrix(mc_nnm_info$L_hat)[1]  
    
    lapis_info <- LAPIS_with_rank_estimation(Y=Y, 
                           W=W, initial_rank=rankMatrix(mc_nnm_info$L_hat)[1],
                           tolerance=tolerance, 
                           min_iter=min_iter, max_iter=max_iter,   
                           mu_grid=NULL, warm_start=F, method = 'explicit_tau')
      
      
    tau_estimate_lapis <- treat.estimator(Y, lapis_info, W)
    
    tau_estimate_oracle <- treat.estimator(Y=Y, L.hat=L, W=W)

    msesFixedParameterData['DID', j, i] <- mean(abs(tau_estimate_did-delta_t)^2)
    
    msesFixedParameterData['SC', j, i] <- mean(abs(tau_estimate_sc-delta_t)^2)
    
    msesFixedParameterData['MC-NNM', j, i] <- mean(abs(tau_estimate_mc_nnm-delta_t)^2)
    
    msesFixedParameterData['SDID', j, i] <- mean(abs(tau_estimate_sdid-delta_t)^2)
    
    msesFixedParameterData['LAPIS', j, i] <- mean(abs(tau_estimate_lapis-delta_t)^2)
    
    msesFixedParameterData['ORACLE', j, i] <- mean(abs(tau_estimate_oracle-delta_t)^2)
      
    outMessage <- paste('Iteration: ', j, sep='')
      
    print(outMessage)
      
    endTime <- Sys.time()
      
    print(difftime(endTime, startTime, 'mins'))

  }


}


effect_plot <- (ggplot(NULL, aes(x=1:length(delta_t), y=delta_t)) 
                + geom_point() + theme_bw() + xlab("Time") + ylab("Treatment Effect")
                +ggtitle("True Treatment Effect Over Time"))


meltedFixedParameterData <- as_tibble(as.tbl_cube(msesFixedParameterData, met='MSE'))

                                          
meltedFixedParameterData$Method <- factor(meltedFixedParameterData$Method, levels=methodNames)

fixed_parameter_error_frame <- meltedFixedParameterData %>% filter(Method !='DID')                                          
                                          
tableForPresenting <- (fixed_parameter_error_frame %>% group_by(Method) %>% dplyr::summarize(RMSE = round(mean(sqrt(MSE)), 3),
                        SE=round(sd(sqrt(MSE), na.rm=T)/sqrt(n()), 3),
            `RMSE (SE)` = paste(RMSE
                    ,' (', SE, ')',
            sep='')))[, c('Method', 'RMSE (SE)')]
                                          

tableForPresenting <- knitr::kable(tableForPresenting, 'latex')


### WITH N0/N
## Chunk 18

### floor((c(5, seq(10, 90, 10), 98)/100)*N)

all_N0s <- floor(c(.1, .3, .5, .7, .9)*N)

msesN0Data <- array(NA, dim=c(length(methodNames), draws_per_L, length(all_N0s) ),
                   dimnames=c(list(Method=methodNames), list(Iteration=1:draws_per_L), list(N0=all_N0s )))

for (this_N0 in all_N0s){ 
    
  N0ForIndex <- as.character(this_N0 )
  
  set.seed(3729)
  
  if (design=="staggered_adoption"){ ## Come up with a way to vary the lag in the staggered structure
  
  if(lag_structure == "random"){
    
    ones_we_make <- c(rep(0, this_N0), pmin(rpois(N-this_N0, 
                                            lambda=average_treatment_length-1)+1, 
                                      min(max_lag*(N-this_N0), .8*Time)))
    
  }else if (lag_structure=="constant"){ ## Does not control T-T0
    
    ones_we_make <- c(rep(0, this_N0), pmin(max_lag*seq(1, (N-this_N0)), floor(.8*Time)))
    
  }

}else if (design=="simultaneous_adoption"){
  
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
 
  for (i in 1:number_of_L){
  
 # set.seed(3729)

  if (balanced){
    
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
      
      startTime <- Sys.time()
    
    if (error == 'gaussian'){
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W, distribution='gaussian',
                 scalar_sigma=sqrt(sigma_squared))
      
    } else if (error == 't'){
    
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='t', scalar_sigma=sqrt(sigma_squared))
    
    } else if (error == 'poisson'){
      
      if (balanced == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='poisson', scalar_sigma=sqrt(sigma_squared))
      
    } else if (error == 'scaled_gamma'){
      
      if (balanced == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='scaled_gamma', scalar_sigma=sqrt(sigma_squared))
      
    }else if (error == 'exponential'){
      
      if (balanced == F){
      
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

    tau_estimate_sdid <- SDID_general(Y=Y, W=W,
                 iterations_for_coord_desc=100)
    
    mc_nnm_info <- matrix_completion_causal(Y=Y, W=W, num_iter=1000, K=5, 
                            lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                            tol=1e-04)
    
    L_mc_nnm <- mc_nnm_info$L_hat
    
    tau_estimate_mc_nnm <- treat.estimator(Y=Y, L.hat=L_mc_nnm, W=W)
      
      
    estFactors <- rankMatrix(mc_nnm_info$L_hat)[1]  
      
    lapis_info <- LAPIS_with_rank_estimation(Y=Y, 
                           W=W, initial_rank=rankMatrix(mc_nnm_info$L_hat)[1],
                           tolerance=tolerance, 
                           min_iter=min_iter, max_iter=max_iter,   
                           mu_grid=NULL, warm_start=F, method = 'explicit_tau')
      
      
    tau_estimate_lapis <- treat.estimator(Y, lapis_info, W)
    
    tau_estimate_oracle <- treat.estimator(Y=Y, L.hat=L, W=W)
    

    msesN0Data['DID', j, N0ForIndex] <- mean(abs(tau_estimate_did-delta_t)^2)
    
    msesN0Data['SC', j, N0ForIndex] <- mean(abs(tau_estimate_sc-delta_t)^2)
    
    msesN0Data['MC-NNM', j, N0ForIndex] <- mean(abs(tau_estimate_mc_nnm-delta_t)^2)
    
    msesN0Data['SDID', j, N0ForIndex] <- mean(abs(tau_estimate_sdid-delta_t)^2)
    
    msesN0Data['LAPIS', j, N0ForIndex] <- mean(abs(tau_estimate_lapis-delta_t)^2)
    
    msesN0Data['ORACLE', j, N0ForIndex] <- mean(abs(tau_estimate_oracle-delta_t)^2)
      
    endTime <- Sys.time()
      
    outMessage <- paste('N0/N: ', this_N0/N, '; Iteration: ', j, sep='')
      
    print(outMessage)
      
    print(difftime(endTime, startTime, 'mins'))

  }


}

}
  
meltedN0Data <- as_tibble(as.tbl_cube(msesN0Data, met='MSE'))

meltedN0Data$Method <- factor(meltedN0Data$Method, levels=methodNames)                                            
                                            
meltedN0Data  <- meltedN0Data  %>% filter(Method !='DID') 
 
aggN0Data <- meltedN0Data %>% group_by(Method, N0) %>% dplyr::summarize(RMSE = mean(sqrt(MSE)),
                                                                SE = sd(sqrt(MSE))/sqrt(n()))

p_mse_vs_N0 <- (ggplot(aggN0Data, aes(x=N0/N, y=RMSE, col=Method)) + geom_line(lwd=1.5) + 
                   theme_bw(base_size=20)+ ggtitle("RMSE as a Function of N0/N") 
                +xlab("N0/N") + ylab('RMSE'))
                                            
p_mse_vs_N0_uncertainty <- (ggplot(meltedN0Data, aes(x=N0/N, y=sqrt(MSE), col=Method)) + geom_smooth(lwd=1.5) + 
                   theme_bw(base_size=20)+ ggtitle("RMSE as a Function of N0/N") 
                +xlab("N0/N") + ylab('RMSE'))                                      
                                        

N0RMSEDataForLatex <- aggN0Data %>% mutate(N0_over_N = round(N0/N, 3),
    `RMSE (SE)` = paste(round(RMSE, 3), ' (', round(SE, 3), ')', sep='')) 

N0RMSEDataForLatex <- N0RMSEDataForLatex[c('Method', 'RMSE (SE)', 'N0_over_N')]


N0RMSEDataForLatex <- spread(N0RMSEDataForLatex[, c('Method', "N0_over_N", "RMSE (SE)")], key=N0_over_N, value=`RMSE (SE)`)


N0RMSEDataForLatex<- knitr::kable(N0RMSEDataForLatex, 'latex')

## Chunk 19

signal_to_noise_ratios <- c()

all_rhos <- c(.1, .3, .5, .7, .9)

msesRhoData <- array(NA, dim=c(length(methodNames), draws_per_L, length(all_rhos) ),
                   dimnames=c(list(Method=methodNames), list(Iteration=1:draws_per_L), list(Rho=all_rhos )))




for (this_rho in all_rhos){ 
    
    
  rhoForIndex <- as.character(this_rho)
  
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
    
  

}else if (design=="simultaneous_adoption"){
  
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
  
  for (i in 1:number_of_L){
  
 # set.seed(3729)
  
  errors_this_L_did <- rep(NA, draws_per_L)
  
  errors_this_L_sc <- rep(NA, draws_per_L)
  
  errors_this_L_mc_nnm <- rep(NA, draws_per_L)
  
  errors_this_L_sdid <- rep(NA, draws_per_L)
  
  errors_this_L_lapis <- rep(NA, draws_per_L)
  
  errors_this_L_oracle <- rep(NA, draws_per_L)
  
  if (balanced){
    
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
      
    startTime <- Sys.time()
    
    if (error == 'gaussian'){
      
      Y <- norta(number=N, corr_mat=this_autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W, distribution='gaussian',
                 scalar_sigma=sqrt(sigma_squared))
      
    } else if (error == 't'){
    
      Y <- norta(number=N, corr_mat=this_autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='t', scalar_sigma=sqrt(sigma_squared))
    
    } else if (error == 'poisson'){
      
      if (balanced == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=this_autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='poisson', scalar_sigma=sqrt(sigma_squared))
      
    } else if (error == 'scaled_gamma'){
      
      if (balanced == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=this_autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='scaled_gamma', scalar_sigma=sqrt(sigma_squared))
      
    }else if (error == 'exponential'){
      
      if (balanced == F){
      
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
    
    tau_estimate_sdid <- SDID_general(Y=Y, W=W,
                 iterations_for_coord_desc=100)
    
    mc_nnm_info <- matrix_completion_causal(Y=Y, W=W, num_iter=1000, K=5, 
                            lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                            tol=1e-04)
    
    L_mc_nnm <- mc_nnm_info$L_hat
    
    tau_estimate_mc_nnm <- treat.estimator(Y=Y, L.hat=L_mc_nnm, W=W)
      
      
    estFactors <- rankMatrix(mc_nnm_info$L_hat)[1]  
    
    lapis_info <- LAPIS_with_rank_estimation(Y=Y, 
                           W=W, initial_rank=rankMatrix(mc_nnm_info$L_hat)[1],
                           tolerance=tolerance, 
                           min_iter=min_iter, max_iter=max_iter,   
                           mu_grid=NULL, warm_start=F, method = 'explicit_tau')
      
      
    tau_estimate_lapis <- treat.estimator(Y, lapis_info, W)
    
    tau_estimate_oracle <- treat.estimator(Y=Y, L.hat=L, W=W)
    
    ## Only oracle in the sense that we know L
    
    error_tau_sc <- mean(abs(tau_estimate_sc-delta_t)^2)
    
    error_tau_did <- mean(abs(tau_estimate_did-delta_t)^2)
    
    error_tau_mc_nnm <- mean(abs(tau_estimate_mc_nnm-delta_t)^2)
    
    error_tau_sdid <- mean(abs(tau_estimate_sdid-delta_t)^2)
    
    error_tau_lapis <- mean(abs(tau_estimate_lapis-delta_t)^2)
    
    error_tau_oracle <- mean(abs(tau_estimate_oracle-delta_t)^2)

    errors_this_L_did[j] <- error_tau_did
    
    errors_this_L_sc[j] <- error_tau_sc
      
    errors_this_L_mc_nnm[j] <- error_tau_mc_nnm
    
    errors_this_L_sdid[j] <- error_tau_sdid
    
    errors_this_L_lapis[j] <- error_tau_lapis
    
    errors_this_L_oracle[j] <- error_tau_oracle

    msesRhoData['DID', j, rhoForIndex] <- mean(abs(tau_estimate_did-delta_t)^2)
    
    msesRhoData['SC', j, rhoForIndex] <- mean(abs(tau_estimate_sc-delta_t)^2)
    
    msesRhoData['MC-NNM', j, rhoForIndex] <- mean(abs(tau_estimate_mc_nnm-delta_t)^2)
    
    msesRhoData['SDID', j, rhoForIndex] <- mean(abs(tau_estimate_sdid-delta_t)^2)
    
    msesRhoData['LAPIS', j, rhoForIndex] <- mean(abs(tau_estimate_lapis-delta_t)^2)
    
    msesRhoData['ORACLE', j, rhoForIndex] <- mean(abs(tau_estimate_oracle-delta_t)^2)
      
    endTime <- Sys.time()
      
    outMessage <- paste('rho: ', this_rho, '; Iteration: ', j, sep='')
      
    print(outMessage)
      
    print(difftime(endTime, startTime, 'mins'))

  }   

}
  
  
  signal_to_noise_ratios <- c(signal_to_noise_ratios, (svd(L)$d[R]/svd(sigma_squared*this_autocorrelation_matrix)$d[1]))
       
       }
  

meltedRhoData <- as_tibble(as.tbl_cube(msesRhoData, met='MSE'))

meltedRhoData$Method <- factor(meltedRhoData$Method, levels=methodNames)                                                  
                                            
meltedRhoData  <- meltedRhoData  %>% filter(Method !='DID') 
 
aggRhoData <- meltedRhoData %>% group_by(Method, Rho) %>% dplyr::summarize(RMSE = mean(sqrt(MSE)),
                                                                SE = sd(sqrt(MSE))/sqrt(n()))

p_mse_vs_rho <- (ggplot(aggRhoData, aes(x=Rho, y=RMSE, col=Method)) + geom_line(lwd=1.5) + 
                   theme_bw(base_size=20)+ ggtitle("RMSE as a Function of the Correlation Parameter") 
                +xlab("rho") + ylab('RMSE'))
                                            
p_mse_vs_rho_uncertainty <- (ggplot(meltedRhoData, aes(x=Rho, y=sqrt(MSE), col=Method)) + geom_smooth(lwd=1.5) + 
                   theme_bw(base_size=20)+ ggtitle("RMSE as a Function of Correlation Parameter") 
                +xlab("rho") + ylab('RMSE'))

rhoRMSEDataForLatex <- aggRhoData %>% mutate(
    `RMSE (SE)` = paste(round(RMSE, 3), ' (', round(SE, 3), ')', sep='')) 

rhoRMSEDataForLatex <- rhoRMSEDataForLatex[c('Method', 'RMSE (SE)', 'Rho')]


rhoRMSEDataForLatex <- spread(rhoRMSEDataForLatex, key=Rho, value=`RMSE (SE)`)


rhoRMSEDataForLatex<- knitr::kable(rhoRMSEDataForLatex, 'latex')


all_taus <- c(1, 11, 16, 21, 26)

msesTauData <- array(NA, dim=c(length(methodNames), draws_per_L, length(all_taus) ),
                   dimnames=c(list(Method=methodNames), list(Iteration=1:draws_per_L), list(Tau=all_taus )))

for (this_tau in all_taus){

  set.seed(3729)
    
  tauForIndex <- as.character(this_tau)
  
  if (design=="staggered_adoption"){ ## Come up with a way to vary the lag in the staggered structure
  
  if(lag_structure == "random"){
    
    ones_we_make <- c(rep(0, N0), pmin(rpois(N-N0, 
                                            lambda=average_treatment_length-1)+1, 
                                      min(max_lag*(N-N0), .8*Time)))
    
  }else if (lag_structure=="constant"){ ## Does not control T-T0
    
    ones_we_make <- c(rep(0, N0), pmin(max_lag*seq(1, (N-N0)), floor(.8*Time)))
    
  }

}else if (design=="simultaneous_adoption"){
  
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
  
  autocorrelation_matrix <- make_rho_mat(rho=rho_parameter, p=dim(W)[2])

  
  for (i in 1:number_of_L){
    
   # set.seed(3729)
    
    if (balanced){
      
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
        
      startTime <- Sys.time()
      
      if (error == 'gaussian'){
        
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+this_tau_matrix*W, distribution='gaussian',
                   scalar_sigma=sqrt(sigma_squared))
        
      } else if (error == 't'){
      
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+this_tau_matrix*W,
                   distribution='t', scalar_sigma=sqrt(sigma_squared))
      
      } else if (error == 'poisson'){
        
        if (balanced == F){
        
          L <- abs(L)+1
          
        }
        
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+this_tau_matrix*W,
                   distribution='poisson', scalar_sigma=sqrt(sigma_squared))
        
      } else if (error == 'scaled_gamma'){
        
        if (balanced == F){
        
          L <- abs(L)+1
          
        }
        
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+this_tau_matrix*W,
                   distribution='scaled_gamma', scalar_sigma=sqrt(sigma_squared))
        
      }else if (error == 'exponential'){
        
        if (balanced == F){
        
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

      tau_estimate_sdid <- SDID_general(Y=Y, W=W,
                   iterations_for_coord_desc=100)
      
      mc_nnm_info <- matrix_completion_causal(Y=Y, W=W, num_iter=1000, K=5, 
                            lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                            tol=1e-04)
    
    L_mc_nnm <- mc_nnm_info$L_hat
    
    tau_estimate_mc_nnm <- treat.estimator(Y=Y, L.hat=L_mc_nnm, W=W)
        
        
    estFactors <- rankMatrix(mc_nnm_info$L_hat)[1]  
    
    lapis_info <- LAPIS_with_rank_estimation(Y=Y, 
                           W=W, initial_rank=rankMatrix(mc_nnm_info$L_hat)[1],
                           tolerance=tolerance, 
                           min_iter=min_iter, max_iter=max_iter,   
                           mu_grid=NULL, warm_start=F, method = 'explicit_tau')
      
      
    tau_estimate_lapis <- treat.estimator(Y, lapis_info, W)
    
    
    # tau_estimate_mc_nnm
    
    ## Only oracle in the sense that we know L
    
    tau_estimate_oracle <- treat.estimator(Y=Y, L.hat=L, W=W)
  
    msesTauData['DID', j, tauForIndex] <- mean(abs(tau_estimate_did-this_tau)^2)
    
    msesTauData['SC', j, tauForIndex] <- mean(abs(tau_estimate_sc-this_tau)^2)
    
    msesTauData['MC-NNM', j, tauForIndex] <- mean(abs(tau_estimate_mc_nnm-this_tau)^2)
    
    msesTauData['SDID', j, tauForIndex] <- mean(abs(tau_estimate_sdid-this_tau)^2)
    
    msesTauData['LAPIS', j, tauForIndex] <- mean(abs(tau_estimate_lapis-this_tau)^2)
    
    msesTauData['ORACLE', j, tauForIndex] <- mean(abs(tau_estimate_oracle-this_tau)^2)
        
    endTime <- Sys.time()
        
    outMessage <- paste('tau: ', this_tau, '; Iteration: ', j, sep='')
      
    print(outMessage)
      
    print(difftime(endTime, startTime, 'mins'))
  
    }


}

}
                                            
                                            
 meltedTauData <- as_tibble(as.tbl_cube(msesTauData, met='MSE'))

meltedTauData$Method <- factor(meltedTauData$Method, levels=methodNames)                                                  
                                            
meltedTauData  <- meltedTauData  %>% filter(Method !='DID') 
 
aggTauData <- meltedTauData %>% group_by(Method, Tau) %>% dplyr::summarize(RMSE = mean(sqrt(MSE)),
                                                                SE = sd(sqrt(MSE))/sqrt(n()))

p_mse_vs_tau <- (ggplot(aggTauData, aes(x=Tau, y=RMSE, col=Method)) + geom_line(lwd=1.5) + 
                   theme_bw(base_size=20)+ ggtitle("RMSE as a Function of Tau") )
                                            
p_mse_vs_tau_uncertainty <- (ggplot(meltedTauData, aes(x=Tau, y=sqrt(MSE), col=Method)) + geom_smooth(lwd=1.5) + 
                   theme_bw(base_size=20)+ ggtitle("RMSE as a Function of Tau") )                                           
                                                         
                                            
                                            
                                            


tauRMSEDataForLatex <- aggTauData %>% mutate(
    `RMSE (SE)` = paste(round(RMSE, 3), ' (', round(SE, 3), ')', sep='')) 

tauRMSEDataForLatex <- tauRMSEDataForLatex[c('Method', 'RMSE (SE)', 'Tau')]

tauRMSEDataForLatex <- spread(tauRMSEDataForLatex, key=Tau, value=`RMSE (SE)`) 

tauRMSEDataForLatex<- knitr::kable(tauRMSEDataForLatex, 'latex')

#all_ranks <- seq(2, 20, 2)



#all_ranks <- c(5, 8, 12, 16, 20)

all_ranks <- c(3,5,7,9,11)

msesRankData <- array(NA, dim=c(length(methodNames), draws_per_L, length(all_ranks) ),
                   dimnames=c(list(Method=methodNames), list(Iteration=1:draws_per_L), list(Rank=all_ranks)))

for (rank in all_ranks){
    
  rankForIndex <- as.character(rank)
  
  set.seed(3729)
  
  if (design=="staggered_adoption"){ ## Come up with a way to vary the lag in the staggered structure
  
  if(lag_structure == "random"){
    
    ones_we_make <- c(rep(0, N0), pmin(rpois(N-N0, 
                                            lambda=average_treatment_length-1)+1, 
                                      min(max_lag*(N-N0), .8*Time)))
    
  }else if (lag_structure=="constant"){ ## Does not control T-T0
    
    ones_we_make <- c(rep(0, N0), pmin(max_lag*seq(1, (N-N0)), floor(.8*Time)))
    
  }

}else if (design=="simultaneous_adoption"){
  
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
                                            
  prediction_error_matrix_mc_nnm <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

  prediction_error_matrix_sdid <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)

  prediction_error_matrix_lapis <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  prediction_error_matrix_oracle <- matrix(NA, nrow=number_of_L, ncol=draws_per_L)
  
  autocorrelation_matrix <- make_rho_mat(rho=rho_parameter, p=dim(W)[2])
  
for (i in 1:number_of_L){
    
    if (balanced){
      
      U_vec <- rexp(n=N*rank, rate=1)
  
      V_vec <- rexp(n=Time*rank, rate=1)
    
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
        
      startTime <- Sys.time()
      
      if (error == 'gaussian'){
        
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+tau*W, distribution='gaussian',
                   scalar_sigma=sqrt(sigma_squared))
        
      } else if (error == 't'){
      
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+tau*W,
                   distribution='t', scalar_sigma=sqrt(sigma_squared))
      
      } else if (error == 'poisson'){
        
        if (balanced == F){
        
          L <- abs(L)+1
          
        }
        
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+tau*W,
                   distribution='poisson', scalar_sigma=sqrt(sigma_squared))
        
      } else if (error == 'scaled_gamma'){
        
        if (balanced == F){
        
          L <- abs(L)+1
          
        }
        
        Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                   desired_mean_matrix= L+tau*W,
                   distribution='scaled_gamma', scalar_sigma=sqrt(sigma_squared))
        
      }else if (error == 'exponential'){
        
        if (balanced == F){
        
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


      tau_estimate_sdid <- SDID_general(Y=Y, W=W,
                   iterations_for_coord_desc=100)
      
    mc_nnm_info <- matrix_completion_causal(Y=Y, W=W, num_iter=1000, K=5, 
                            lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                            tol=1e-04)
    
    L_mc_nnm <- mc_nnm_info$L_hat
    
    tau_estimate_mc_nnm <- treat.estimator(Y=Y, L.hat=L_mc_nnm, W=W)
        
    estFactors <- rankMatrix(mc_nnm_info$L_hat)[1]  
    
    lapis_info <- LAPIS_with_rank_estimation(Y=Y, 
                           W=W, initial_rank=rankMatrix(mc_nnm_info$L_hat)[1],
                           tolerance=tolerance, 
                           min_iter=min_iter, max_iter=max_iter,   
                           mu_grid=NULL, warm_start=F, method = 'explicit_tau')
      
    tau_estimate_lapis <- treat.estimator(Y, lapis_info, W)
    
    tau_estimate_oracle <- treat.estimator(Y=Y, L.hat=L, W=W)
    
    ## Only oracle in the sense that we know L

  
      msesRankData['DID', j, rankForIndex] <- mean(abs(tau_estimate_did-delta_t)^2)
    
      msesRankData['SC', j, rankForIndex] <- mean(abs(tau_estimate_sc-delta_t)^2)
    
      msesRankData['MC-NNM', j, rankForIndex] <- mean(abs(tau_estimate_mc_nnm-delta_t)^2)
    
      msesRankData['SDID', j, rankForIndex] <- mean(abs(tau_estimate_sdid-delta_t)^2)
    
      msesRankData['LAPIS', j, rankForIndex] <- mean(abs(tau_estimate_lapis-delta_t)^2)
    
      msesRankData['ORACLE', j, rankForIndex] <- mean(abs(tau_estimate_oracle-delta_t)^2)
    
      endTime <- Sys.time()
        
      outMessage <- paste('rank: ', rank, '; Iteration: ', j, sep='')
      
      print(outMessage)
      
      print(difftime(endTime, startTime, 'mins'))
  
    }

}

}
                                          
meltedRankData <- as_tibble(as.tbl_cube(msesRankData, met='MSE'))
                                            
meltedRankData$Method <- factor(meltedRankData$Method, levels=methodNames)      
                                            
meltedRankData <- meltedRankData  %>% filter(Method !='DID') 
 
aggRankData <- meltedRankData %>% group_by(Method, Rank) %>% dplyr::summarize(RMSE = mean(sqrt(MSE)),
                                                                SE = sd(sqrt(MSE))/sqrt(n()))

p_mse_vs_rank <- (ggplot(aggRankData, aes(x=Rank, y=RMSE, col=Method)) + geom_line(lwd=1.5) + 
                   theme_bw(base_size=20)+ ggtitle("RMSE as a Function of the Rank of L") + xlab("Rank of L") )
                                            
p_mse_vs_rank_uncertainty <- (ggplot(meltedRankData, aes(x=Rank, y=sqrt(MSE), col=Method)) + geom_smooth(lwd=1.5) + 
                   theme_bw(base_size=20)+ ggtitle("RMSE as a Function of the Rank of L")  + xlab("Rank of L") )                                              
                                            
                                            
                                            
                                            


rankRMSEDataForLatex <- aggRankData %>% mutate(
    `RMSE (SE)` = paste(round(RMSE, 3), ' (', round(SE, 3), ')', sep='')) 

rankRMSEDataForLatex <- rankRMSEDataForLatex[c('Method', 'RMSE (SE)', 'Rank')]

rankRMSEDataForLatex <- spread(rankRMSEDataForLatex, key=Rank, value=`RMSE (SE)`)

rankRMSEDataForLatex <- knitr::kable(rankRMSEDataForLatex, 'latex')

all_rank_errors <- c(-9, -5, 1, 3, 7)


msesRankErrorData <- array(NA, dim=c(length(methodNames), draws_per_L, length(all_rank_errors) ),
                   dimnames=c(list(Method=methodNames), list(Iteration=1:draws_per_L), list(Rank_error=all_rank_errors)))


for (rank_error in all_rank_errors){
  
  rankErrorForIndex <- as.character(rank_error)
  
  set.seed(3729)
  
  if (design=="staggered_adoption"){ ## Come up with a way to vary the lag in the staggered structure
  
  if(lag_structure == "random"){
    
    ones_we_make <- c(rep(0, N0), pmin(rpois(N-N0, 
                                            lambda=average_treatment_length-1)+1, 
                                      min(max_lag*(N-N0), .8*Time)))
    
  }else if (lag_structure=="constant"){ ## Does not control T-T0
    
    ones_we_make <- c(rep(0, N0), pmin(max_lag*seq(1, (N-N0)), floor(.8*Time)))
    
  }

}else if (design=="simultaneous_adoption"){
  
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
  
  autocorrelation_matrix <- make_rho_mat(rho=rho_parameter, p=dim(W)[2])
  
  for (i in 1:number_of_L){
    
    if (balanced){
      
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
        
      startTime <- Sys.time()
      
        if (error == 'gaussian'){
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W, distribution='gaussian',
                 scalar_sigma=sqrt(sigma_squared))
      
    } else if (error == 't'){
    
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='t', scalar_sigma=1, df=df)
    
    } else if (error == 'poisson'){
      
      if (balanced == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='poisson', scalar_sigma=1)
      
    } else if (error == 'scaled_gamma'){
      
      if (balanced == F){
      
        L <- abs(L)+1
        
      }
      
      Y <- norta(number=N, corr_mat=autocorrelation_matrix,
                 desired_mean_matrix= L+tau_matrix*W,
                 distribution='scaled_gamma', scalar_sigma=1)
      
    }else if (error == 'exponential'){
      
      if (balanced == F){
      
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

    tau_estimate_sdid <- SDID_general(Y=Y, W=W,
                 iterations_for_coord_desc=100)
    
    mc_nnm_info <- matrix_completion_causal(Y=Y, W=W, num_iter=1000, K=5, 
                            lambda_grid=c(10^seq(-4,2,1), seq(2,5,1)),
                            tol=1e-04)
    
    L_mc_nnm <- mc_nnm_info$L_hat
    
    tau_estimate_mc_nnm <- treat.estimator(Y=Y, L.hat=L_mc_nnm, W=W)
    
        
    estFactors <- rankMatrix(mc_nnm_info$L_hat)[1]  
    
    lapis_info <- LAPIS_with_rank_estimation(Y=Y, 
                           W=W, initial_rank=R+rank_error,
                           tolerance=tolerance, 
                           min_iter=min_iter, max_iter=max_iter,   
                           mu_grid=NULL, warm_start=F, method = 'explicit_tau')
      
      
    tau_estimate_lapis <- treat.estimator(Y, lapis_info, W)
    
    tau_estimate_oracle <- treat.estimator(Y=Y, L.hat=L, W=W)
    
    
    # tau_estimate_mc_nnm
    
    ## Only oracle in the sense that we know L
    
      msesRankErrorData['DID', j, rankErrorForIndex] <- mean(abs(tau_estimate_did-delta_t)^2)
    
      msesRankErrorData['SC', j, rankErrorForIndex] <-  mean(abs(tau_estimate_sc-delta_t)^2)
    
      msesRankErrorData['MC-NNM', j, rankErrorForIndex] <- mean(abs(tau_estimate_mc_nnm-delta_t)^2)
    
      msesRankErrorData['SDID', j, rankErrorForIndex] <- mean(abs(tau_estimate_sdid-delta_t)^2)
    
     msesRankErrorData['LAPIS', j, rankErrorForIndex] <- mean(abs(tau_estimate_lapis-delta_t)^2)
    
      msesRankErrorData['ORACLE', j, rankErrorForIndex] <- mean(abs(tau_estimate_oracle-delta_t)^2)
        
      endTime <- Sys.time()
        
      outMessage <- paste('initial rank error: ', rank_error, '; Iteration: ', j, sep='')
      
      print(outMessage)
      
      print(difftime(endTime, startTime, 'mins'))

    }

  } 

    print(paste("Finished rank error:", rank_error))
}

                                            
meltedRankErrorData <- as_tibble(as.tbl_cube(msesRankErrorData, met="MSE"))
                                            
meltedRankErrorData$Method <- factor(meltedRankErrorData$Method, levels=methodNames)      

                                            
meltedRankErrorData <- meltedRankErrorData  %>% filter(Method !='DID') 
 
aggRankErrorData <- meltedRankErrorData %>% group_by(Method, Rank_error) %>% dplyr::summarize(RMSE = mean(sqrt(MSE)),
                                                                SE = sd(sqrt(MSE))/sqrt(n()))

p_mse_vs_rankError <- (ggplot(aggRankErrorData, aes(x=Rank_error, y=RMSE, col=Method)) + geom_line(lwd=1.5) + 
                   theme_bw(base_size=20)+ ggtitle("RMSE as a Function of the Rank of the Initial Rank Error") 
                +xlab("Initial Rank Error") + ylab('RMSE'))
                                            
p_mse_vs_rankError_uncertainty <- (ggplot(meltedRankErrorData, aes(x=Rank_error, y=sqrt(MSE), col=Method)) + geom_smooth(lwd=1.5) + 
                   theme_bw(base_size=20)+ ggtitle("RMSE as a Function of the Rank of the Initial Rank Error") 
                +xlab("Initial Rank Error") + ylab('RMSE'))   

rankErrorRMSEDataForLatex <- aggRankErrorData %>% mutate(
    `RMSE (SE)` = paste(round(RMSE, 3), ' (', round(SE, 3), ')', sep='')) 

rankErrorRMSEDataForLatex <- rankErrorRMSEDataForLatex[c('Method', 'RMSE (SE)', 'Rank_error')]


rankErrorRMSEDataForLatex <-spread(rankErrorRMSEDataForLatex, key=Rank_error, value=`RMSE (SE)`)


rankErrorRMSEDataForLatex <- knitr::kable(rankErrorRMSEDataForLatex, 'latex')

server_name <- Sys.info()['nodename']                                        
                                            
                                            
sim_error_date_directory <- paste("../reports/", server_name,'/',
design, "_simulations","/" ,c('unbalanced', 'balanced')[balanced+1], '/' , Sys.Date(), "_simulations", sep='')

if (!dir.exists(sim_error_date_directory)){
  
  dir.create(path=sim_error_date_directory, recursive = TRUE)
  
}                                            
                                            
                                            
runsAlreadyDone <- length(list.files(sim_error_date_directory))                                            
                                            
                                            
                        
                                            
full_output_directory <- paste(sim_error_date_directory, paste('run', runsAlreadyDone+1, sep='_'),sep='/')

image_directory <- paste(full_output_directory, "/simulation_plots", sep="")

table_directory <- paste(full_output_directory, "/simulation_tables", sep="")

if (!dir.exists(image_directory)){

    dir.create(image_directory, recursive=TRUE)
  
}          


if (!dir.exists(table_directory)){

    dir.create(table_directory, recursive=TRUE)
  
}       
                                            
          


image_directory

if (exists("effect_plot")){

ggsave(paste(image_directory,'/effect_size_plot.pdf' , sep=''), effect_plot,
       width=11, height=8.5)
}

if (exists("tableForPresenting")){
                                            
    fileConn<-file(makeFilename(table_directory, 
                                fileName='fixed_parameter_error_table.txt'))
    writeLines(tableForPresenting , fileConn)
    close(fileConn)    

    }


 
if (exists("p_mse_vs_N0")){

    ggsave(paste(image_directory,'/varied_N0.pdf' , sep=''), p_mse_vs_N0, width=11, height=8.5)

    ggsave(paste(image_directory,'/varied_N0_with_uncertainty.pdf' , sep=''), p_mse_vs_N0_uncertainty, 
           width=11, height=8.5)    
    
    fileConn<-file(makeFilename(table_directory, 
                                fileName='N0_rmse_table.txt'))
    writeLines(N0RMSEDataForLatex , fileConn)
    close(fileConn) 

}
                                            
                                                                  
if (exists("p_mse_vs_rho")){
  
    ggsave(paste(image_directory,'/varied_rho.pdf' , sep=''), p_mse_vs_rho,
           width=11, height=8.5)

    ggsave(paste(image_directory,'/varied_rho_with_uncertainty.pdf' , sep=''), p_mse_vs_rho_uncertainty,
           width=11, height=8.5)     
    
    
    fileConn<-file(makeFilename(table_directory, 
                                fileName='rho_rmse_table.txt'))
    writeLines(rhoRMSEDataForLatex , fileConn)
    close(fileConn)  
     
     
}


if (exists("p_mse_vs_tau")){

    ggsave(paste(image_directory,'/varied_tau.pdf' , sep=''), p_mse_vs_tau,
           width=11, height=8.5)

    ggsave(paste(image_directory,'/varied_tau_with_uncertainty.pdf' , sep=''), p_mse_vs_tau_uncertainty,
           width=11, height=8.5)
    
    
    fileConn<-file(makeFilename(table_directory, 
                                fileName='tau_rmse_table.txt'))
    writeLines(tauRMSEDataForLatex , fileConn)
    close(fileConn)  
} 

if (exists("p_mse_vs_rank")){

    ggsave(paste(image_directory,'/varied_rank.pdf' , sep=''), p_mse_vs_rank, width=11, height=8.5)

    ggsave(paste(image_directory,'/varied_rank_with_uncertainty.pdf' , sep=''), 
           p_mse_vs_rank_uncertainty, width=11, height=8.5)
    
    fileConn<-file(makeFilename(table_directory, 
                                fileName='rank_rmse_table.txt'))
    writeLines(rankRMSEDataForLatex , fileConn)
    close(fileConn) 
  
}   



if (exists("p_mse_vs_rankError")){

    ggsave(paste(image_directory,'/varied_rank_error.pdf' , sep=''), p_mse_vs_rankError, width=11, height=8.5)

    ggsave(paste(image_directory,'/varied_rank_error_with_uncertainty.pdf' , sep=''), 
           p_mse_vs_rankError_uncertainty, width=11, height=8.5)    
    
    
    fileConn<-file(makeFilename(table_directory, 
                                fileName='rank_error_rmse_table.txt'))
    writeLines(rankErrorRMSEDataForLatex , fileConn)
    close(fileConn) 
    
    
}




                                            
                                           
write.csv(params, file=paste(full_output_directory,
                             "simulation_parameters.csv", sep="/"))                                         

# write.table(L, file="fixed_L.txt", row.names=FALSE, col.names=FALSE)

stopCluster(cl)


