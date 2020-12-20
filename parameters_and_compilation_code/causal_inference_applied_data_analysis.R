library(dplyr)
library(ggplot2)
library(openxlsx)

### For Both These Data sets





#### SUICIDE DATA

suicide_data <- read.csv("./suicide data/suicide_per_100000_people.csv", row.names = 1)

years_for_data <- seq(1950, 2016, 1)


first_recording <- apply(suicide_data, 
                         MARGIN=1, FUN=function(x) 
                           min(which(!is.na(x))))

last_recording <- apply(suicide_data, 
                         MARGIN=1, FUN=function(x) 
                           max(which(!is.na(x))))

start_year <- years_for_data[first_recording]

end_year <- years_for_data[last_recording]


countries_we_consider <- which(start_year <= 1981 & end_year >= 1990)


row.names(suicide_data)[countries_we_consider]

table(start_year)

dim(suicide_data)[1]

## Plan of attack: 

plot(as.numeric(suicide_data["United Kingdom", ]) ~ years_for_data, type='l')

## 1. Imputation of missing values NOT INCLUDING JAPAN

## 2. add Japan back in for counterfactual estimation

## 3. All methods will use the imputed data



###### VOTER TURNOUT DATA

voter_data <- read.csv("turnout.csv")

state_example <- voter_data %>% filter(abb=='ID')




