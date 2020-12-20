library(knitr)
library(openxlsx)
library(markdown)
library(rmarkdown)

#install.packages("openxlsx", dependencies = TRUE)

server_name <- Sys.info()['nodename']

parameter_data <- read.xlsx("parameters_and_descriptions.xlsx",
                            sheet="parameter_data", rowNames = T)

number_of_Ls <- as.numeric(parameter_data['number_of_Ls', ])
 
draws_per_L <- as.numeric(parameter_data['draws_per_L', ])

N <- as.numeric(parameter_data['N', ])

N1 <- as.numeric(parameter_data['N1', ])

Time <- as.numeric(parameter_data['Time', ])

Time1 <- as.numeric(parameter_data['Time1', ])

matrix_type <- c("wide", "tall")[(N > Time)+1]

R <- as.numeric(parameter_data['R', ])

rho_parameter <- as.numeric(parameter_data['rho_parameter', ])

tau <- as.numeric(parameter_data['tau', ])

sigma_squared <- as.numeric(parameter_data['sigma_squared', ])

penalized <- as.logical(parameter_data['penalized', ])

exchangable <- as.logical(parameter_data['exchangable', ])

min_iter <- as.numeric(parameter_data['min_iter', ])

max_iter <- as.numeric(parameter_data['max_iter', ])

tolerance <- as.numeric(parameter_data['tolerance', ])

error <- parameter_data['error', ]

df <- as.numeric(parameter_data['df', ])

L_scaling <- as.numeric(parameter_data['L_scaling', ])

arg_max <- as.numeric(parameter_data['arg_max', ])

y_max <- as.numeric(parameter_data['y_max', ])

halfway_time <- as.numeric(parameter_data['halfway_time', ])

cutoff <- as.numeric(parameter_data['cutoff', ])

treatment_effect_function <- parameter_data['treatment_effect_function', ]

design <- parameter_data['design', ]

lag_structure <- parameter_data['lag_structure', ]

max_lag <- parameter_data['max_lag', ]

average_treatment_length <- as.numeric(parameter_data['average_treatment_length', ])

sim_error_date_directory <- paste("../reports/", server_name,'/',
error, "_simulations","/" , Sys.Date(), "_simulations", sep='')

if (!dir.exists(sim_error_date_directory)){
  
  dir.create(path=sim_error_date_directory, recursive = TRUE)
  
}

file_name <- paste("report_", Sys.Date(), sep='')

files_in_corresponding_directory <- list.files(sim_error_date_directory)

sim_number <- length(files_in_corresponding_directory)+1

individual_simulation_directory <- paste(sim_error_date_directory, "/",file_name,
                                         "_run_", matrix_type, "_",sim_number, sep="")

dir.create(path=individual_simulation_directory)

file_name <- paste(file_name, "_run_", sim_number, "_",matrix_type,".pdf", sep="")

all_params <- list(number_of_L=number_of_Ls, draws_per_L=draws_per_L, N=N, 
                   N1=N1, Time=Time, Time1=Time1, R=R, rho_parameter=rho_parameter, 
                   tau=tau, sigma_squared=sigma_squared, 
                   penalized=penalized, exchangable=exchangable,
                   min_iter=min_iter, max_iter=max_iter,
                   tolerance=tolerance, error=error, df=df,
                   L_scaling=L_scaling, 
                   arg_max=arg_max, y_max=y_max, halfway_time=halfway_time,
                   cutoff=cutoff, treatment_effect_function=treatment_effect_function,
                   design=design,average_treatment_length=average_treatment_length,
                   output_directory = individual_simulation_directory)

rmarkdown::render("./synth_did_simulations.Rmd",  
output_file=file_name, output_dir=individual_simulation_directory, 
params=all_params)

