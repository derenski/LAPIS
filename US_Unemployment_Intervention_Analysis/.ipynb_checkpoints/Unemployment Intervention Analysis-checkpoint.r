library(ggplot2)
library(readxl)
library(stringr)
library(dplyr)
library(reshape)
library(lubridate)
library(zoo)
library(MASS)
library(LaplacesDemon)
library(gsynth)
library(doParallel)
library(foreach)
library(forecast)
library(did)
library(data.table)

source('../causal_inference_methods_code.R')

simParams <- read.csv('./unemployment_simulation_parameters.csv', row.names = 1,
                     stringsAsFactors=FALSE)

numberOfIterations <- as.numeric(simParams['number_of_iterations', 1])

timeConsidered = as.numeric(simParams['time_considered', 1])
maxTreatLength = as.numeric(simParams['max_treat_length',1])
balanced = as.logical(simParams['balanced', 1])

## simultaneous_adoption, staggered_adoption
design <- simParams['design', 1]

max_available_clusters <- detectCores()-1
  
desired_clusters <- 4
  
cl <- makeCluster(min(c(max_available_clusters, desired_clusters)))

registerDoParallel(cl)

variableExpander <- function(keysAndVariable, unitKey, timeKey){ ### Expands into unit/time format
  
  orderedData <- keysAndVariable %>% arrange(get(unitKey), get(timeKey))
  
  outcomeName <- names(keysAndVariable)[!(names(keysAndVariable) %in% c(unitKey, timeKey))]
  
  outcomeMatrixForm <- matrix(NA, nrow=length(unique(keysAndVariable[, unitKey])), 
                              ncol=length(unique(keysAndVariable[, timeKey])), byrow=T)
  
  rownames(outcomeMatrixForm) <- unique(keysAndVariable[, unitKey])
  
  colnames(outcomeMatrixForm) <- str_replace_all(unique(keysAndVariable[, timeKey])[order(unique(keysAndVariable[, timeKey]))], pattern="-", replace="")
  
  for (index in 1:length(keysAndVariable[, unitKey])){
    
    outcomeMatrixForm[keysAndVariable[, unitKey][index], str_replace_all(keysAndVariable[, timeKey][index], pattern="-", replace="")] <- keysAndVariable[, outcomeName][index]
    
  }
  
  return(outcomeMatrixForm)
  
}

seriesWithUnemploymentDecay <- function(aSeries, finalProp=.96, decayRate=.1, treatmentIndexer){ ### Creates treatment effect for implementation of drunk driving law
    
    if (length(aSeries) != length(treatmentIndexer)){
        
        stop("Series and Treatment Indexer Must be the Same Length.")
    }
    
    if (all(treatmentIndexer==0)){
        
        return(aSeries)
        
    }else{
    
        dataWhereTransformNeeded <- aSeries[treatmentIndexer==1]
    
        transformedPart <- finalProp*dataWhereTransformNeeded+(1-finalProp)*dataWhereTransformNeeded*exp(-1*decayRate*(1:length(dataWhereTransformNeeded)))
        
        aSeries[treatmentIndexer==1] <- transformedPart
        
        return(aSeries)
        
    }
    
}

probFun <- function(r, xBeta){
      
      exp(xBeta+r)/(1+(exp(xBeta+r)))
      
      
    }



meanProbFun <- function(r, xBeta, maxTreatLength, aveTreatDuration){
    
    initialVals <- exp(xBeta+r)/(1+(exp(xBeta+r)))
    
    return(maxTreatLength*mean(initialVals, na.rm=T)-aveTreatDuration)
    
    
}

rFunBin <- function(r, xBeta, C, lengthCutoff){
      
      vals <- (1-probFun(r=r, xBeta=xBeta))^lengthCutoff
      
      return(mean(vals, na.rm=T)-(1-C))
      
    }

simParams <- read.csv('./unemployment_simulation_parameters.csv', row.names = 1,
                     stringsAsFactors=FALSE)

numberOfIterations <- as.numeric(simParams['number_of_iterations', 1])

timeConsidered = as.numeric(simParams['time_considered', 1])
maxTreatLength = as.numeric(simParams['max_treat_length',1])
balanced = as.logical(simParams['balanced', 1])

## simultaneous_adoption, staggered_adoption
design <- simParams['design', 1]

firstSet <- read_excel('GeoFRED_Unemployment_Rate_by_County_Percent.xls', sheet=1)
colnames(firstSet) <- tolower(str_replace_all(colnames(firstSet), pattern=' ', replace='_'))


secondSet <- read_excel('GeoFRED_Unemployment_Rate_by_County_Percent.xls', sheet=2)
colnames(secondSet) <- tolower(str_replace_all(colnames(secondSet), pattern=' ', replace='_'))

fullSet <- (firstSet %>% inner_join(secondSet))

fullSet <- fullSet[!apply(fullSet, MARGIN=1, FUN=function(x) any(is.na(x))), 1:(min(which(str_detect(colnames(fullSet), '2020')))-1)]

justTimeSeriesData <- as.numeric(t(as.matrix(fullSet[4:dim(fullSet)[2]])))

sum(apply(fullSet, MARGIN=1, FUN=function(x) any(is.na(x))))

allIds <- fullSet[, 1:3] %>% slice(rep(1:n(), each = length(4:dim(fullSet)[2])))

meltedUnemploymentData <- cbind.data.frame(allIds, justTimeSeriesData)
names(meltedUnemploymentData)[4] <- 'unemployment_rate'
meltedUnemploymentData <- meltedUnemploymentData %>% 
mutate(region_designation=str_extract(region_name, pattern='^.*(?=\\,)'), 
       state=str_extract(region_name, pattern='[A-Z]{2}$'))

timePoints <- seq(1990, 2019+(11/12), 1/12)

numericSeries <- as.numeric(firstSet[1,4:256])

bestLambda <- BoxCox.lambda(numericSeries)


testSeries <- seasadj(stl(ts(as.numeric(firstSet[1,4:256]), frequency=12), s.window=21))



seasAdjustGeneral <- function(aSeries, frequency=12){
    
    aSeries <- ts(aSeries, frequency=frequency)
    
    bestLambda <- BoxCox.lambda(aSeries)
    
    transformedSeries <- BoxCox(aSeries, lambda=bestLambda)
    
    stlDecompTransformedSeries <- stl(transformedSeries, s.window=7)
    
    seasAdjustedSeries <- seasadj(stlDecompTransformedSeries)
    
    finalSeries <- InvBoxCox(seasAdjustedSeries, bestLambda)
    
    return(as.numeric(finalSeries))
    
    
}

meltedUnemploymentData <- meltedUnemploymentData %>% group_by(region_code) %>% 
mutate(time=timePoints, 
      seasonal_adjust_unemp = seasAdjustGeneral(unemployment_rate, frequency=12)) %>% ungroup()
meltedUnemploymentDataForModeling <- meltedUnemploymentData %>% filter(!(state %in% c('AK', 'HI')))

countiesToPlot <- 'Ventura|San Bernardino|Riverside|Los Angeles|Imperial|San Diego|Santa Barbara|Orange|San Luis Obispo|Kern'

soCalData <- meltedUnemploymentData %>% filter(str_detect(region_designation, countiesToPlot) &
                                              state=='CA')


alaskaData <- meltedUnemploymentData %>% filter(
                                              state=='AK')

names(soCalData)

socalUnempPlot <- ggplot(soCalData, aes(x=time, y=seasonal_adjust_unemp, col=region_designation)) + geom_line(lwd=1.5) + 
theme_bw(base_size=20) + ylab('Unemployment Rate (Seasonally Adjusted)') + xlab('Time') +
guides(col=guide_legend(title="County")) + ggtitle("Unemployment in Southern California")


alaskaUnempPlot <- ggplot(alaskaData, aes(x=time, y=seasonal_adjust_unemp, col=region_designation)) + geom_line(lwd=1.5) + 
theme_bw(base_size=20) + ylab('Unemployment Rate (Seasonally Adjusted)') + xlab('Time') +
guides(col=guide_legend(title="County")) + ggtitle("Unemployment in Alaska") 

validSmoothedDataOnly <- meltedUnemploymentDataForModeling  %>% filter(!is.na(seasonal_adjust_unemp))

aCountyExample = validSmoothedDataOnly %>% filter(region_name=='Ventura County, CA')
plot(aCountyExample$unemployment_rate, type='l')
lines(aCountyExample$seasonal_adjust_unemp, col='red', lwd=1.5)
#plot(aCountyExample$seasonal_adjust_unemp, type='l')

fullOutcomeData <- variableExpander(data.frame(validSmoothedDataOnly[c('region_code', 'time', 
                    'seasonal_adjust_unemp')]), unitKey='region_code', timeKey='time')

sum(is.na(fullOutcomeData))

unique(table(meltedUnemploymentData$region_code)) ### Verifies all series observed at all time points

aCounty <- validSmoothedDataOnly %>% filter(region_designation=='Ventura County')

recessionDates <- rbind(c("1990-07-01", "1991-03-01"), 
                        c("2001-03-01", "2001-11-01"), 
                        c("2007-12-01", "2009-06-01"))



recessionDates

histAves <- rowMeans(fullOutcomeData)

valueDataFrames <- list()

compTimeDataFrames <- list()

for(propTreat in seq(.1, .9, .1)){

    startTime <- Sys.time()
    
    valueMatrix <- array(NA, dim=c(4, numberOfIterations),
                         dimnames=c(list(c('MC-NNM', 'Synthetic Control', 'Synthetic Difference in Differences',
                                     'LAPIS')), list(paste('Iteration', 1:numberOfIterations, sep='_'))))
    
    compTimeMatrix <- array(NA, dim=c(4, numberOfIterations),
                         dimnames=c(list(c('MC-NNM', 'Synthetic Control', 'Synthetic Difference in Differences',
                                     'LAPIS')), list(paste('Iteration', 1:numberOfIterations, sep='_'))))

    for (simNumber in 1:numberOfIterations){ 
        
        ### Key insight: control group doesn't change, but treated units have duration distribution
        
        unitsToSample <- sample(1:dim(fullOutcomeData)[1], size=600, replace=FALSE)

        fullOutcomeDataSample <- fullOutcomeData[unitsToSample, ]


        numberTreat <- floor(dim(fullOutcomeDataSample)[1]*propTreat)

        numberUntreat <- dim(fullOutcomeDataSample)[1]-numberTreat

        ### For a Balanced Treatment Design

        if (balanced){

         if (design=='staggered_adoption'){

             durations <- c(rep(0, numberUntreat) , rbinom(numberTreat, maxTreatLength, .3))

            }else{

            durations <- 30*rbern(dim(fullOutcomeDataSample)[1], p=propTreat) 

         }

         durations <- sample(durations, size=length(durations), replace=FALSE)

            }else{
            #.3*rbeta(1, 3, 1.5)
            betaScalar <- .7

            xBetas <- betaScalar*histAves[unitsToSample]

            xBetas[xBetas <= quantile(xBetas, 1-propTreat)] <- -Inf

                   # r <- uniroot(rFunBin, xBeta= xBetas[xBetas != -Inf], C=propTreat,
                   #          lengthCutoff=maxTreatLength,
                   #          lower=-500, upper=500)$root


            r <- uniroot(meanProbFun, xBeta= xBetas[!is.infinite(xBetas)], maxTreatLength=100, aveTreatDuration=30,
                             lower=-500, upper=500)$root

            pis <- probFun(r=r, xBetas)

            if (design=='staggered_adoption'){

            durations <- rbinom(length(pis), maxTreatLength, pis)

                }else{

               # durations <- 30*rbern(lengths(pis), p=pis)
                
                durations <- 30*rep(1, length(pis))
                
                durations[is.infinite(xBetas)] <- 0

            }

        }

        D <- array(0, dim=dim(fullOutcomeDataSample))

        dimnames(D) <- dimnames(fullOutcomeDataSample)

        for (i in 1:length(durations)){

            allTimes <- 1:dim(D)[2]

            D[i, dim(D)[2]-durations[i]+1 <= allTimes] <- 1
        }

        timeConsidered <- min(c(timeConsidered, max(rowSums(D))))

        observedData<- t(sapply(1:dim(fullOutcomeDataSample)[1], FUN=function(x) 
            seriesWithUnemploymentDecay(aSeries=fullOutcomeDataSample[x,], finalProp=.6, decayRate=.05, 
                                        treatmentIndexer = D[x,])))

        dimnames(observedData) <- dimnames(fullOutcomeDataSample)  
                                
                                
        startTimeMC_NNM <- Sys.time()                        

        mcnnm_info <- matrix_completion_causal(Y=observedData, W=D)     

        mcnnm_est <- treat.estimator(observedData, mcnnm_info$L_hat, D)
                                
        endTimeMC_NNM <- Sys.time() 
                                
        compTimeMC_NNM <- as.numeric(difftime(endTimeMC_NNM, startTimeMC_NNM))
                                
                                
        startTimeLAPIS <- Sys.time()

        withMu <- LAPIS(Y=observedData, W=D)
                                
        endTimeLAPIS <- Sys.time()                       
                                
        compTimeLAPIS <- as.numeric(difftime(endTimeLAPIS, startTimeLAPIS))
                                     
                                
        startTimeSynthCont <- Sys.time()                           

        synthContEst <- synth_cont(Y=observedData, W=D)
                                
        endTimeSynthCont <- Sys.time()
                                
        compTimeSynthCont <- as.numeric(difftime(endTimeSynthCont, startTimeSynthCont))
                                
                                
        startTimeSDID <- Sys.time()  

        SDIDEst <- SDID_general(Y=observedData, W=D,
                         iterations_for_coord_desc=100)
                                
                                
        endTimeSDID <- Sys.time()  
                                
        compTimeSDID <- as.numeric(difftime(endTimeSDID, startTimeSDID))

        trueEffect <- treat.estimator(observedData, fullOutcomeDataSample, D)


        errorsThisIteration <- c(mean((mcnnm_est-trueEffect)[1:timeConsidered]^2), 
                                 mean((synthContEst-trueEffect)[1:timeConsidered]^2),
                                 mean((SDIDEst-trueEffect)[1:timeConsidered]^2),
                                 mean((withMu-trueEffect)[1:timeConsidered]^2))
                                
        compTimesThisIteration <- c(compTimeMC_NNM, compTimeSynthCont, compTimeSDID,
                                  compTimeLAPIS)
                                
        compTimeMatrix[, simNumber] <- compTimesThisIteration


        valueMatrix[ , simNumber] <- errorsThisIteration

        endTime <- Sys.time()

        print(simNumber)

        print(difftime(endTime, startTime, 'mins'))

        }

    meltedValueMatrix <- reshape::melt(valueMatrix)

    names(meltedValueMatrix) <- c('Method', 'iteration', 'MSE')

    meltedValueMatrix$prop_treat <- propTreat
                                
    valueDataFrames <- c(valueDataFrames, list(meltedValueMatrix))
                                
    meltedCompTimeMatrix <- reshape::melt(compTimeMatrix)

    names(meltedCompTimeMatrix) <- c('Method', 'iteration', 'computation_time_secs')

    meltedCompTimeMatrix$prop_treat <- propTreat
                                
    compTimeDataFrames <- c(compTimeDataFrames, list(meltedCompTimeMatrix))                                    
                                                             
    print(propTreat)
                                
    print(difftime(endTime, startTime, 'mins'))
                            
    }
                                    

errorTable <- rbindlist(valueDataFrames)
aggErr <- errorTable %>% group_by(Method, prop_treat) %>% summarise(mean_RMSE = mean(sqrt(MSE)))
errorFunctionOfPropUntreat <- ggplot(aggErr, aes(x=1-prop_treat, y=mean_RMSE, col=Method)) + geom_line(lwd=1.5) + 
theme_bw(base_size=20) + 
xlab('N0/N') + ylab('RMSE')


errorFunctionOfPropUntreatWithUncertainty <- ggplot(errorTable , 
                aes(x=1-prop_treat, y=sqrt(MSE), col=Method)) + geom_smooth(lwd=1.5) + theme_bw(base_size=20) + 
xlab('N0/N') + ylab('RMSE')






compTimeTable <- rbindlist(compTimeDataFrames)
aggCompTime <- compTimeTable %>% group_by(Method, prop_treat) %>% summarise(mean_compTime = 
                                                                            mean(computation_time_secs))
compTimeFunctionOfPropUntreat <- ggplot(aggCompTime, aes(x=1-prop_treat, y=mean_compTime, col=Method)) + geom_line(lwd=1.5) + 
theme_bw(base_size=20) + 
xlab('N0/N') + ylab('Mean Computation Time (Seconds)')


compTimeFunctionOfPropUntreatWithUncertainty <- ggplot(compTimeTable , 
                aes(x=1-prop_treat, y=computation_time_secs, col=Method)) + geom_smooth(lwd=1.5) + theme_bw(base_size=20) + 
xlab('N0/N') + ylab('Mean Computation Time (Seconds)')

if (balanced){

    oneYearOfAdoption <- rownames(D)[rowSums(D) >=  min(c(max(rowSums(D)), 12))   ]
    CACountiesThatAdopted <- validSmoothedDataOnly %>% filter(region_code %in% oneYearOfAdoption & state=='CA')

    countyOfInterest <- sample(CACountiesThatAdopted$region_code, 1)




    adoptionTime <- as.numeric(colnames(D)[max(which((D[countyOfInterest, ]==0)))])


    unempVals <- as.numeric(c(observedData[countyOfInterest ,], fullOutcomeData[countyOfInterest, ]))
    times <- as.numeric(rep(colnames(D), 2))

    treatVals <- rep(c("Observed", "Untreated Potential Outcome"), each=dim(D)[2])


    adoptionDataCountyExample <- cbind.data.frame(times, treatVals, unempVals)

    names(adoptionDataCountyExample) <- c('Time', 'Outcome', 'Unemployment')


    adoptionExamplePlot <- ggplot(adoptionDataCountyExample %>% filter(Time >= adoptionTime-36/12), 
                                  aes(x=Time, y=Unemployment, col=Outcome, lty=Outcome)) + 
    geom_line(lwd=2) + theme_bw(base_size=20) + ylab('Unemployment Rate (Seasonally Adjusted)') + xlab('Time') +
     ggtitle("A Successful Retraining Program?") + geom_vline(xintercept = adoptionTime, lwd=1.25) 

}

if (!balanced){
    
    pretreatmentPeriod <- max(which(colSums(D)==0))

    pretreatmentData <- observedData[, 1:pretreatmentPeriod]

    untreatedData <- colMeans(pretreatmentData[which(rowSums(D)==0),])
    treatedData <- colMeans(pretreatmentData[which(rowSums(D)> 0),])

    groupData <- cbind.data.frame(as.numeric(names(untreatedData), 2), c(rep("Treated", length(treatedData)), 
                       rep('Untreated', length(untreatedData))) , c(treatedData, untreatedData))

    names(groupData) <- c('Time', 'Group', 'Unemployment Rate (Seasonally Adjusted)')

    unemploymentGroupPlot <- ggplot(groupData, aes(x=Time, y=`Unemployment Rate (Seasonally Adjusted)`, col=Group)) + 
    geom_line(lwd=1.5) + theme_bw(base_size=20) + ggtitle('Unemployment by Treatment Group')
    
    }

meltedErrors <- reshape::melt(valueMatrix)

names(meltedErrors) <- c('Method', 'Itertation', 'MSE')

meltedErrors <- meltedErrors %>% mutate(RMSE = sqrt(MSE))

ggplot(meltedErrors, aes(x=Method, y=RMSE, fill=Method)) + geom_boxplot() + theme_bw(base_size=12)

aggregatedErrors <- errorTable %>% group_by(Method, prop_treat) %>% dplyr::summarize(RMSE = round(mean(sqrt(MSE), na.rm=T),
                                                                                                  3),
                                               SE = round(sd(sqrt(MSE), na.rm=T)/n(), 3))


aggregatedErrors$`RMSE (SE)` <- paste(aggregatedErrors$RMSE, ' (', aggregatedErrors$SE, ')', sep='')

aggregatedErrors$Method <- factor(aggregatedErrors$Method, levels=c("Synthetic Control", 'Generalized Synthetic Control',
                                                            'Synthetic Difference in Differences',
                                                            'MC-NNM', 'LAPIS'))


aggregatedErrors <- aggregatedErrors[order(aggregatedErrors$Method), ]

tableForLatexRMSE <- aggregatedErrors[, c('Method', 'prop_treat','RMSE (SE)')]

tableForLatexRMSE <- reshape2::dcast(tableForLatexRMSE , Method ~ (1-prop_treat))

latexRMSETable <- knitr::kable(tableForLatexRMSE, format='latex')

aggregatedCompTime <- compTimeTable %>% group_by(Method, prop_treat) %>% dplyr::summarize(mean_comp_time = round(mean(computation_time_secs, na.rm=T),
                                                                                                  3),
                                               SE = round(sd(computation_time_secs, na.rm=T)/n(), 3))


aggregatedCompTime$`Mean Computation Time (SE), Seconds` <- paste(aggregatedCompTime$mean_comp_time, ' (', aggregatedErrors$SE, ')', sep='')

aggregatedCompTime$Method <- factor(aggregatedCompTime$Method, levels=c("Synthetic Control", 'Generalized Synthetic Control',
                                                            'Synthetic Difference in Differences',
                                                            'MC-NNM', 'LAPIS'))


aggregatedCompTime <- aggregatedCompTime[order(aggregatedCompTime$Method), ]

tableForLatexCompTime <- aggregatedCompTime[, c('Method', 'prop_treat','Mean Computation Time (SE), Seconds')]

tableForLatexCompTime <- reshape2::dcast(tableForLatexCompTime, Method ~ (1-prop_treat))

latexCompTimeTable <- knitr::kable(tableForLatexCompTime, format='latex')

estData <- c(synthContEst, SDIDEst, mcnnm_est, withMu, trueEffect)

methodValues <- rep(c("Synthetic Control",
                                                            'Synthetic Difference in Differences',
                                                            'MC-NNM', 'LAPIS', 'Truth'), each=length(trueEffect))

relativeTime <- rep(1:length(trueEffect), length.out=length(methodValues))

plotExample <- cbind.data.frame(methodValues, relativeTime ,estData)

names(plotExample) <- c('Method', 'Time', 'Change in Seasonally Adjusted Unemployment Rate')

plotExample <- plotExample %>% filter(Time <= timeConsidered)

effectPlot <- ggplot(plotExample, aes(x=Time, y=`Change in Seasonally Adjusted Unemployment Rate`, col=Method)) + 
geom_line(lwd=2) + theme_bw(base_size=20)

adoptionDurations <- rowSums(D)

adoptionDistribution <- table(adoptionDurations)

adoptDistData <- cbind.data.frame(as.numeric(names(adoptionDistribution)), as.numeric(adoptionDistribution))

names(adoptDistData) <- c('Duration', 'Count')

adoptDistExamp <- ggplot(adoptDistData, aes(x=Duration, y=log(Count))) + geom_bar(stat="identity",width=.1) + 
geom_point() + theme_bw(base_size=20) + theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank()) + xlab("Adoption Duration") + ggtitle("Adoption of a Job Retraining Program")

adoptDistExamp 

plotsAndTablesDir <- paste(paste('./plots_and_tables', c('unbalanced', 'balanced')[balanced+1], sep='_'), design, sep='/')

if (!exists(plotsAndTablesDir)){
    
    
    dir.create(plotsAndTablesDir, recursive = TRUE)
}

ggsave(paste(plotsAndTablesDir,'unemployment_effect.pdf', sep='/'), effectPlot, width=11, height=8.5 )
ggsave(paste(plotsAndTablesDir, 'unemployment_change.pdf', sep='/'), socalUnempPlot, width=11, height=8.5)

if (balanced){
    ggsave(paste(plotsAndTablesDir, 'program_adoption_example.pdf', sep='/'), adoptionExamplePlot, width=11, height=8.5)
    
}

if (!balanced){
    
  ggsave(paste(plotsAndTablesDir, 'unemployment_group_comparison.pdf', sep='/'), 
         unemploymentGroupPlot, width=11, height=8.5)  
    
    
}

ggsave(paste(plotsAndTablesDir,'adoption_curve_example.pdf', sep='/'), adoptDistExamp, width=11, height=8.5)

ggsave(paste(plotsAndTablesDir,'error_vs_prop_treat.pdf', sep='/'), errorFunctionOfPropUntreat, width=11, height=8.5)

ggsave(paste(plotsAndTablesDir,'error_vs_prop_treat_with_se.pdf', sep='/'),
       errorFunctionOfPropUntreatWithUncertainty , width=11, height=8.5)

ggsave(paste(plotsAndTablesDir,'comp_time_vs_prop_treat.pdf', sep='/'), compTimeFunctionOfPropUntreat, width=11, height=8.5)

ggsave(paste(plotsAndTablesDir,'comp_time_vs_prop_treat_with_se.pdf', sep='/'),
       compTimeFunctionOfPropUntreatWithUncertainty, width=11, height=8.5)



fileConn<-file(paste(plotsAndTablesDir, 'unemployment_error_table.txt', sep='/'))
writeLines(latexRMSETable, fileConn)
close(fileConn)

fileConn<-file(paste(plotsAndTablesDir, 'unemployment_comp_time_table.txt', sep='/'))
writeLines(latexCompTimeTable, fileConn)
close(fileConn)

write.table(simParams, paste(plotsAndTablesDir, 'parameter_table.csv', sep='/'))


#ggsave(paste(plotsAndTablesDir, 'unemployment_change_alaska.pdf', sep='/'), alaskaUnempPlot, width=11, height=8.5)


