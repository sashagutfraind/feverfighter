#
# sensitivity analysis for epidemic outcomes using LHS
# 
#

library(truncnorm)

### configure the user directory - this should be the box directory
if(Sys.info()['effective_user'] == 'agutf') {
  if(Sys.info()['sysname'] == 'Windows') {
    setwd("C:\\Users\\agutf\\Sync\\school-flu\\code")
  } else {
    setwd("/home/sasha/projects/schools-flu/")
  }
} else if (Sys.info()['effective_user'] == 'adam') {
  setwd("/home/adam/projects/schools-flu/code")
} else if (Sys.info()['effective_user'] == 'Adam') {
  setwd("C:\\Users\\Adam\\Documents\\schools-flu\\code")
}
if(!file.exists("output")) {
  dir.create("output")
}

source('schools_solver.R')
source('lhslib.R')



compute_shape <-function(x) {
  a <- ((4 + x^2 - x) / (1 - x))
  b <- ((4 + x^2 - x) / (x))
  return(list(shape1=a, shape2=b))
}

param_description_ranges_penn <- list(
  #Known values, look them up from our paper or online
  vac_rate          = 0, #this is in the paper
  vac_efficacy      = 0,
  #Unknown values, use our calibrated data to define the mean, and let the sd = mean/2 (depending on func type)
  #transmissibility_weekend_ratio  = list(func=qgamma, args=list(shape=1.0, rate=1.0)), 
  #transmissibility_closure_ratio  = list(func=qgamma, args=list(shape=1.0, rate=1.0)),
  seasonality       = list(func=qbeta, args=compute_shape(0.8830765)),   #the maximal seasonal effect on transmissibility
  transmissibility  = list(func=qtruncnorm, args=list(a=0, mean=0.003401314, sd=0.003401314/2)), 
  transmissibility_weekend_ratio = list(func=qbeta, args=compute_shape(0.7000992)),  #reduction during weekends
  transmissibility_closure_ratio = list(func=qbeta, args=compute_shape(0.1774011)),
  cross_grade_contact = list(func=qbeta, args=compute_shape(0.3552442)),
  start_date_shift  = list(func=qtruncnorm, args=list(a=-7, mean=0, sd=3.5, b=7)),  #epidemic start day varies by about 30 days each year
  name = "penn_pandemic"
) 


########################################
##SEASONAL FLU STARTS HERE
########################################
param_description_ranges_thames <- list(
  #Known values, look them up from our paper or online
  vac_rate          =list(func=qtruncnorm, args=list(a=0, b=1, mean=0.6, sd=0.1)), #this is in the paper
  vac_efficacy      =list(func=qtruncnorm, args=list(a=0.65, b=0.95, mean=0.67, sd=0.075)),
  #Unknown values, use our calibrated data to define the mean, and let the sd = mean/2 (depending on func type)
  seasonality       = list(func=qbeta, args=compute_shape(0.007458082)),   #the maximal seasonal effect on transmissibility
  transmissibility  = list(func=qtruncnorm, args=list(a=0, mean=0.002163292 , sd=0.002163292 /2)), 
  transmissibility_weekend_ratio = list(func=qbeta, args=compute_shape(0.6250557)),  #reduction during weekends.  
  transmissibility_closure_ratio = list(func=qbeta, args=compute_shape(0.1131701)),
  cross_grade_contact = list(func=qbeta, args=compute_shape(0.2449116)),
  start_date_shift  = list(func=qtruncnorm, args=list(a=-7, mean=0, sd=3.5, b=7)),  #epidemic start day varies by about 30 days each year
  name = "thames_seasonal"
)

param_description_ranges_thames_scenario <- list(
  vac_rate          = list(func=qtruncnorm, args=list(a=0, mean=0.80, sd=0.10, b=1)),  #fraction of children vaccinated
  vac_efficacy      = list(func=qtruncnorm, args=list(a=0, mean=0.50, sd=0.10, b=1)),  #efficiency of the vaccine in preventing illness        
  symptom_attention   = list(func=qtruncnorm, args=list(a=0, mean=0.6806607, sd=0.6806607/2, b=1)),
  compliance        = list(func=qtruncnorm, args=list(a=0, mean=0.7975382, sd=0.7975382/2, b=1)),
  seasonality       = list(func=qtruncnorm, args=list(a=0, mean=0.1174269, sd=0.1174269/2, b=1)),   #the maximal seasonal effect on transmissibility
  transmissibility = list(func=qtruncnorm, args=list(a=0, mean=0.002029424, sd=0.002029424/2, b=1)),  #the transmissibility for the outbreak
  transmissibility_weekend_ratio = list(func=qtruncnorm, args=list(a=0, mean=0.7026213, sd=0.7026213/2, b=1)),  #reduction during weekends.  
  transmissibility_closure_ratio = list(func=qtruncnorm, args=list(a=0, mean=0.1775939, sd=0.1775939/2, b=1)),
  symptom_propensity   = 0.88,
  start_date_shift=list(func=qtruncnorm, args=list(a=-30, mean=0, sd=5, b=30)),  
  cross_grade_contact = list(func=qtruncnorm, args=list(a=0, mean=0.1409936, sd=0.1409936/2, b=1)),
  name = "scenario_seasonal"
)
########################################
##PANDEMIC FLU STARTS HERE
########################################

##DO NOT CHANGE THE FOLLOWING PARAMETERS##
param_description_ranges_flu_scenario <- list(
  #VALUES BELOW
  vac_rate          = 0.0,  #fraction of children vaccinated
  vac_efficacy      = 0.0,  #efficiency of the vaccine in preventing illness        
  symptom_propensity=list(func=qtruncnorm, args=list(a=0.81, mean=0.84, sd=0.015, b=0.87)),
  transmissibility=list(func=qtruncnorm, args=list(a=0.000402, mean=0.000402, sd=0.0000402, b=0.000402)),
  transmissibility_weekend_ratio=list(func=qtruncnorm, args=list(a=0.280141, mean=0.280141, sd=0.008754, b=0.298801)),
  transmissibility_closure_ratio=list(func=qtruncnorm, args=list(a=0.141479, mean=0.141479, sd=0.001822, b=0.161479)),
  symptom_attention=list(func=qtruncnorm, args=list(a=0.698836, mean=0.698836, sd=0.001567, b=0.736156)),
  compliance=list(func=qtruncnorm, args=list(a=0.102405, mean=0.158385, sd=0.199545, b=0.996142)),
  seasonality=list(func=qtruncnorm, args=list(a=0.458733, mean=0.478733, sd=0.001187, b=0.478733)),
  cross_grade_contact=list(func=qtruncnorm, args=list(a=0.589565, mean=0.605565, sd=0.000672, b=0.605565)),
  start_date_shift=list(func=qtruncnorm, args=list(a=-30, mean=0, sd=5, b=30)),  
  name = "flu_baseline"
)
##DO NOT CHANGE THE ABOVE PARAMETERS - computed by calibration##

param_description_ranges_corona_scenario <- list(
  vac_rate          = 0.0,  #fraction of children vaccinated
  vac_efficacy      = 0.0,  #efficiency of the vaccine in preventing illness        
  symptom_propensity=list(func=qtruncnorm, args=list(a=0.81, mean=0.84, sd=0.015, b=0.87)),
  transmissibility=list(func=qtruncnorm, args=list(a=0.000402, mean=0.000402, sd=0, b=0.000402)),
  transmissibility_weekend_ratio=list(func=qtruncnorm, args=list(a=0.280141, mean=0.280141, sd=0.008754, b=0.298801)),
  transmissibility_closure_ratio=list(func=qtruncnorm, args=list(a=0.141479, mean=0.141479, sd=0.001822, b=0.161479)),
  symptom_attention=list(func=qtruncnorm, args=list(a=0.698836, mean=0.698836, sd=0.001567, b=0.736156)),
  compliance=list(func=qtruncnorm, args=list(a=0.102405, mean=0.158385, sd=0.199545, b=0.996142)),
  seasonality=list(func=qtruncnorm, args=list(a=0.458733, mean=0.478733, sd=0.001187, b=0.478733)),
  cross_grade_contact=list(func=qtruncnorm, args=list(a=0.589565, mean=0.605565, sd=0.000672, b=0.605565)),
  start_date_shift=list(func=qtruncnorm, args=list(a=-30, mean=0, sd=5, b=30)),
  name = "scenario_corona"
)


########################################
##COVID STARTS HERE
########################################

##DO NOT CHANGE THE FOLLOWING PARAMETERS##
param_description_ranges_corona <- list(
  vac_rate          = 0.0,  #fraction of children vaccinated
  vac_efficacy      = 0.0,  #efficiency of the vaccine in preventing illness        
  symptom_propensity=list(func=qtruncnorm, args=list(a=0.139, mean=0.181, sd=0.024, b=0.229)),
  transmissibility=list(func=qtruncnorm, args=list(a=0.000511, mean=0.000511, sd=0.0000511, b=0.000511)),
  transmissibility_weekend_ratio=list(func=qtruncnorm, args=list(a=0.280141, mean=0.280141, sd=0.008754, b=0.298801)),
  transmissibility_closure_ratio=list(func=qtruncnorm, args=list(a=0.141479, mean=0.141479, sd=0.001822, b=0.161479)),
  symptom_attention=list(func=qtruncnorm, args=list(a=0.698836, mean=0.698836, sd=0.001567, b=0.736156)),
  compliance=list(func=qtruncnorm, args=list(a=0.102405, mean=0.158385, sd=0.199545, b=0.996142)),
  seasonality=list(func=qtruncnorm, args=list(a=0.458733, mean=0.478733, sd=0.001187, b=0.478733)),
  cross_grade_contact=list(func=qtruncnorm, args=list(a=0.589565, mean=0.605565, sd=0.000672, b=0.605565)),
  start_date_shift=list(func=qtruncnorm, args=list(a=-30, mean=0, sd=5, b=30)),
  name = "static_corona"
)
##DO NOT CHANGE THE ABOVE PARAMETERS##

param_description_ranges_corona_scenario <- list(
  vac_rate          = 0.0,  #fraction of children vaccinated
  vac_efficacy      = 0.0,  #efficiency of the vaccine in preventing illness        
  symptom_propensity=list(func=qtruncnorm, args=list(a=0.139, mean=0.181, sd=0.024, b=0.229)),
  transmissibility=list(func=qtruncnorm, args=list(a=0.000511, mean=0.000511, sd=0.0000511, b=0.000511)),
  transmissibility_weekend_ratio=list(func=qtruncnorm, args=list(a=0.280141, mean=0.280141, sd=0.008754, b=0.298801)),
  transmissibility_closure_ratio=list(func=qtruncnorm, args=list(a=0.141479, mean=0.141479, sd=0.001822, b=0.161479)),
  symptom_attention=list(func=qtruncnorm, args=list(a=0.698836, mean=0.698836, sd=0.001567, b=0.736156)),
  compliance=list(func=qtruncnorm, args=list(a=0.102405, mean=0.158385, sd=0.199545, b=0.996142)),
  seasonality=list(func=qtruncnorm, args=list(a=0.458733, mean=0.478733, sd=0.001187, b=0.478733)),
  cross_grade_contact=list(func=qtruncnorm, args=list(a=0.589565, mean=0.605565, sd=0.000672, b=0.605565)),
  start_date_shift=list(func=qtruncnorm, args=list(a=-30, mean=0, sd=5, b=30)),
  name = "scenario_corona"
)

#test that the shift to higher value is due to truncation of the distribution
#Function below based on https://stat.ethz.ch/R-manual/R-devel/library/stats/html/Beta.html
#Variance is presumed to be (mean/2)^2, and x is the mean of the distribution


compute_epidemic_metrics <- function(params) {
  outbreak_result <- compute_school_outbreak(params = params)
  final_metrics   <- summarize_epidemic_metrics(outbreak_result)
  ret = list(summary_metrics=final_metrics, infected_series=outbreak_result$Infected)
  return(ret)
}

model_func_wrapper <- function(wrapper_params,
                                model_func=compute_epidemic_metrics) {
  outcomeDT <- model_func(wrapper_params$sample)

  return(outcomeDT)
}


#'
#' A function that summarizes the results of lhs runs
#' @param paramset   the parameter set, if NULL loads from fname 
#' @param funcparams parameters directly passed to the lower level function
lhs_summarize_results <- function(run_func, funcparams=NULL, fname=NULL, saveFileID='', paramset=NULL) {
  x.name="summary_metrics"
  y.name="infected_series"
  if(! is.null(paramset)) {
    results_runs <- lhs_run_samples(samples=paramset, 
                                    name="validation_with_data", 
                                    func=run_func, 
                                    funcparams=list(doPlot=FALSE))  
  } else {
    load(file=fname, newenv<-new.env())
    results_pack<-get("results", newenv)
    results_runs <- results_pack[["fullres"]]
  }
  all_runs_data  <- NULL
  all_runs_xdata  <- NULL
  for(sample_num in seq(length(results_runs))) {
    run_data <- results_runs[[sample_num]]
    if(sample_num==1) {
      all_runs_xdata  <- as.data.table(run_data[[x.name]])
      all_runs_data   <- as.data.table(run_data[[y.name]])
    } else {
      all_runs_xdata  <- rbind(all_runs_xdata, as.data.table(run_data[[x.name]]))
      all_runs_data   <- cbind(all_runs_data,  as.data.table(run_data[[y.name]]))
    }
  }
  #helps calculate stats - peak date is calculate as the day of the year relative to the midyear  
  all_runs_xdata$peak_date_New_ILI = as.numeric(difftime(as.POSIXct(all_runs_xdata$peak_date_New_ILI, tz="UTC"), as.POSIXct(all_runs_xdata$start_date, tz="UTC"), units="days"))
  all_runs_xdata$peak_date_New_absent = as.numeric(difftime(as.POSIXct(all_runs_xdata$peak_date_New_absent, tz="UTC"), as.POSIXct(all_runs_xdata$start_date, tz="UTC"), units="days"))
  all_runs_xdata$peak_last_day = as.numeric(difftime(as.POSIXct(all_runs_xdata$peak_last_day, tz="UTC"), as.POSIXct(all_runs_xdata$start_date, tz="UTC"), units="days"))
  all_runs_xdata$peak_New_ILI_last_day = as.numeric(difftime(as.POSIXct(all_runs_xdata$peak_New_ILI_last_day, tz="UTC"), as.POSIXct(all_runs_xdata$start_date, tz="UTC"), units="days"))
  all_runs_xdata$peak_New_absent_last_day = as.numeric(difftime(as.POSIXct(all_runs_xdata$peak_New_absent_last_day, tz="UTC"), as.POSIXct(all_runs_xdata$start_date, tz="UTC"), units="days"))
  all_runs_xdata$peak_date_infected = as.numeric(difftime(as.POSIXct(all_runs_xdata$peak_date_infected, tz="UTC"), as.POSIXct(all_runs_xdata$start_date, tz="UTC"), units="days"))
  all_runs_xdata$start_date = yday(all_runs_xdata$start_date)
  num_failed_cases = sum(!complete.cases(all_runs_data))
  if(num_failed_cases > 0) {
    warning(paste(num_failed_cases, "failed out of ", nrow(all_runs_data)))    
    all_runs_data <- all_runs_data[complete.cases(all_runs_data),]
  }
  if((num_failed_cases/nrow(all_runs_data)) > 0.1) {
    stop("more than 10% of runs failed")
  }
  #drop row of peak dates and make it a difference from start date
  ret<-list()
  ret$summary_metrics <- data.frame(metric=1:length(run_data[[x.name]]))
  ret$summary_metrics[,"metric"] = colnames(all_runs_xdata)
  ret$summary_metrics[,"mean"]   = t(all_runs_xdata[,lapply(.SD, mean)])
  ret$summary_metrics[,"median"]   = t(all_runs_xdata[,lapply(.SD, median)])
  ret$summary_metrics[,"25pc"]   = t(all_runs_xdata[,lapply(.SD, quantile, probs=0.25)])
  ret$summary_metrics[,"75pc"]   = t(all_runs_xdata[,lapply(.SD, quantile, probs=0.75)])
  #add the quantiles here (25, 75) as below
  ret$summary_metrics[,"sd"]     = t(all_runs_xdata[,lapply(.SD, sd)])
  ret$summary_metrics[,"min"]    = t(all_runs_xdata[,lapply(.SD, min)])
  ret$summary_metrics[,"max"]    = t(all_runs_xdata[,lapply(.SD, max)])
  ret$summary_metrics[,"N"]      = rep(dim(all_runs_xdata)[1], dim(all_runs_xdata)[2])
  
  ret$y_series <- data.frame(day=1:length(run_data[[y.name]]))
  ret$y_series[,"mean"] <- apply(all_runs_data, 1, mean)
  ret$y_series[,"sd"]   <- apply(all_runs_data, 1, sd)
  ret$y_series[,"min"]  <- apply(all_runs_data, 1, min)
  ret$y_series[,"max"]  <- apply(all_runs_data, 1, max)
  ret$y_series[,"5pc"]  <- apply(all_runs_data, 1, quantile, 0.05)
  ret$y_series[,"10pc"] <- apply(all_runs_data, 1, quantile, 0.10)
  ret$y_series[,"25pc"] <- apply(all_runs_data, 1, quantile, 0.25)
  ret$y_series[,"75pc"] <- apply(all_runs_data, 1, quantile, 0.75)
  ret$y_series[,"90pc"] <- apply(all_runs_data, 1, quantile, 0.90)
  ret$y_series[,"95pc"] <- apply(all_runs_data, 1, quantile, 0.95)

  saveFname <- time_stamp(paste0('output/summary_metrics', saveFileID), '.csv')
  fwrite(ret$summary_metrics, file=saveFname)
  print(saveFname)
  
  return(ret)
}


z_compute_lhs <- function(num_samples, param_descriptions, baseline_sample=baseline_sample, scenario_name='') {
  print(scenario_name)
  paramset_1  <- lhs_make_samples(num_samples=num_samples, 
                                       param_descriptions=param_descriptions, 
                                       baseline_sample=baseline_sample)
  for(sample_num in seq(length(paramset_1))) {
    paramset_1[[sample_num]]$exclusion_days <- 0
    paramset_1[[sample_num]]$day_start <- paramset_1[[sample_num]]$day_start - round(paramset_1[[sample_num]]$start_date_shift)
    #paramset_1[[sample_num]]$day_end <- paramset_1[[sample_num]]$day_end - round(paramset_1[[sample_num]]$start_date_shift)
    paramset_1[[sample_num]]$day_end <- paramset_1[[sample_num]]$day_start +  + (baseline_sample$day_end - baseline_sample$day_start)
    paramset_1[[sample_num]]$contact_rates <- update_off_diagonal_values(paramset_1[[sample_num]]$contact_rates,
                                                                         paramset_1[[sample_num]]$cross_grade_contact)
  }
  outcomes_1_ranges <- lhs_summarize_results(fname=fname, 
                                              paramset=paramset_1, 
                                              saveFileID=scenario_name, 
                                              run_func=model_func_wrapper)
  fwrite(outcomes_1_ranges$summary_metrics, 
            file=time_stamp(paste0("output/outcomes1_metrics_", scenario_name), ".csv"), 
            row.names = F)
  fwrite(outcomes_1_ranges$y_series, 
            file=time_stamp(paste0("output/outcomes1_infected_", scenario_name), ".csv"), 
            row.names = F)
  return(outcomes_1_ranges)
}

z_compute_lhs_policy <- function(num_samples, param_descriptions, baseline_sample=baseline_sample, scenario_name='', default_policy=0, default_workday=5) {
  print(scenario_name)
  paramset_1  <- lhs_make_samples(num_samples=num_samples, 
                                  param_descriptions=param_descriptions, 
                                  baseline_sample=baseline_sample)
  for(sample_num in seq(length(paramset_1))) {
    paramset_1[[sample_num]]$exclusion_days <- default_policy
    paramset_1[[sample_num]]$workweek_days <- default_workday
    paramset_1[[sample_num]]$day_start <- paramset_1[[sample_num]]$day_start - round(paramset_1[[sample_num]]$start_date_shift)
    #paramset_1[[sample_num]]$day_end <- paramset_1[[sample_num]]$day_end - round(paramset_1[[sample_num]]$start_date_shift)
    paramset_1[[sample_num]]$day_end <- paramset_1[[sample_num]]$day_start + (baseline_sample$day_end - baseline_sample$day_start)
    paramset_1[[sample_num]]$contact_rates <- update_off_diagonal_values(paramset_1[[sample_num]]$contact_rates,
                                                                         paramset_1[[sample_num]]$cross_grade_contact)
  }
  outcomes_1_ranges <- lhs_summarize_results(fname=fname, 
                                             paramset=paramset_1, 
                                             saveFileID=scenario_name, 
                                             run_func=model_func_wrapper)

  return(outcomes_1_ranges)
}


compute_policy_row <- function(data, row) {
  med    <- trunc(data$summary_metrics[row, 3]*1000)/1000 #median
  q1   <- trunc(data$summary_metrics[row, 4]*1000)/1000 #25th percentile
  q3   <- trunc(data$summary_metrics[row, 5]*1000)/1000 #75th percentile
  return(paste(med, ' (', q1, '-', q3, ')', sep=""))
}

compute_policy_table <- function(num_samples, param_descriptions, baseline_sample=baseline_sample, scenario_name='', max_policy) {
  #create empty data table here and rbind below
  policy_table <- data.table(
    policy_day = integer(),
    total_students = character(),
    infected = character(),
    vaccinated = character(),
    persondays_isolated = character(),
    persondays_unisolated = character(),
    peak_New_ILI_cases = character(),
    peak_New_absent = character(),
    peak_infected = character(),
    attack_rate = character(),
    start_date = character(),
    peak_date_New_ILI = character(),
    peak_date_New_absent = character(),
    peak_date_infected = character(),
    peak_duration = character(),
    peak_last_day = character(),
    peak_New_ILI_duration = character(),
    peak_New_ILI_last_day = character(),
    peak_New_absent_duration = character(),
    peak_New_absent_last_day = character()
  )
  for (i in 0:max_policy) {
    print(paste("Currently performing operations for exclusion days = ", i))
    temp_lhs_compute <- z_compute_lhs_policy(num_samples, 
                                      param_descriptions,
                                      baseline_sample,
                                      scenario_name, 
                                      default_policy=i)
    #do for all summary metrics, just just a few
    temp_row <- list()
    temp_row["policy_day"] <- i
    for (row in seq(nrow(temp_lhs_compute$summary_metrics))) {
      temp_row[temp_lhs_compute$summary_metrics[row, 1]] <- compute_policy_row(temp_lhs_compute, row)
    }
    policy_table <- rbind(policy_table, temp_row)
    #Initialise an empty data table and perhaps use rbind to append to it
    
  }

  fname <- paste(time_stamp("output/lhs_output", ".csv"))
  fwrite(policy_table, file = fname)
  print(sprintf("Saved lhs output: %s", fname))
  
  param_descriptions$num_samples = num_samples
  cat(paste(param_descriptions), file=paste0(fname, ".param_descriptions.txt"))
  cat(paste(baseline_sample), file=paste0(fname, ".baseline_sample.txt"))
  
  return(policy_table)
}

compute_workday_table <- function(num_samples, param_descriptions, baseline_sample=baseline_sample, scenario_name='', max_workdays) {
  #create empty data table here and rbind below
  workday_table <- data.table(
    workweek_days = integer(),
    total_students = character(),
    infected = character(),
    vaccinated = character(),
    persondays_isolated = character(),
    persondays_unisolated = character(),
    peak_New_ILI_cases = character(),
    peak_New_absent = character(),
    peak_infected = character(),
    attack_rate = character(),
    start_date = character(),
    peak_date_New_ILI = character(),
    peak_date_New_absent = character(),
    peak_date_infected = character(),
    peak_duration = character(),
    peak_last_day = character(),
    peak_New_ILI_duration = character(),
    peak_New_ILI_last_day = character(),
    peak_New_absent_duration = character(),
    peak_New_absent_last_day = character()
  )
  for (i in 0:max_workdays) {
    print(paste("Currently performing operations for workday = ", i))
    temp_lhs_compute <- z_compute_lhs_policy(num_samples, 
                                             param_descriptions,
                                             baseline_sample,
                                             scenario_name, 
                                             default_workday=i)
    #do for all summary metrics, just just a few
    temp_row <- list()
    temp_row["workweek_days"] <- i
    for (row in seq(nrow(temp_lhs_compute$summary_metrics))) {
      temp_row[temp_lhs_compute$summary_metrics[row, 1]] <- compute_policy_row(temp_lhs_compute, row)
    }
    workday_table <- rbind(workday_table, temp_row)
    #Initialise an empty data table and perhaps use rbind to append to it
    
  }
  
  fname <- paste(time_stamp("output/lhs_output", ".csv"))
  fwrite(workday_table, file = fname)
  print(sprintf("Saved lhs output: %s", fname))
  
  param_descriptions$num_samples = num_samples
  cat(paste(param_descriptions), file=paste0(fname, ".param_descriptions.txt"))
  cat(paste(baseline_sample), file=paste0(fname, ".baseline_sample.txt"))
  
  return(workday_table)
}

#' @desc plot the infected series for multiple levels of workday variable
#' 
z_workday_fig <- function(baseline, workdays_set, savePlots) {
  results = data.table()
  params = baseline
  
  for(d in workdays_set) {
    params$workweek_days = d
    run = compute_school_outbreak(params)
    results <- cbind(results, run$Infected)
  }
  names(results) = as.character(workdays_set)
  results$dates = as.Date(run$dates)
  
  df_long <- reshape2::melt(results, id="dates") # convert to long format
  pl<-ggplot(data=df_long, aes(x=dates, y=value, linetype=variable, color=variable), geom=c("point", "path")) +
    geom_line() + #geom_point() + 
    theme(axis.title.y = element_text(size = 20)) + 
    xlab("Individuals") +
    xlab("") +
    ylab("Currently infected") + 
    scale_x_date() +
    guides(linetype=guide_legend(title = "Weekly\ndays open"),
           color=guide_legend(title =    "Weekly\ndays open"))+
    NULL
  print(pl)
  if(savePlots) {
    basestr = time_stamp("output/infected_with_workdays_", baseline$name)
    ggsave(paste0(basestr, ".png"), width=6, height=4.5, units="in", dpi=300)  
    ggsave(paste0(basestr, ".tiff"), width=6, height=4.5, units="in", dpi=300)  
    ggsave(paste0(basestr, ".pdf"), width=6, height=4.5, units="in")  
  }
  
  return(results)
}


#' @desc plot the infected series for multiple levels of exclusion
#' 
z_policy_combined_fig <- function(baseline, exclusion_set, savePlots) {
  results = data.table()
  params = baseline
  
  for(d in exclusion_set) {
    params$exclusion_days = d
    run = compute_school_outbreak(params)
    results <- cbind(results, run$Infected)
  }
  names(results) = as.character(exclusion_set)
  results$dates = as.Date(run$dates)
  
  df_long <- reshape2::melt(results, id="dates") # convert to long format
  pl<-ggplot(data=df_long, aes(x=dates, y=value, linetype=variable, color=variable), geom=c("point", "path")) +
    geom_line() + #geom_point() + 
    theme(axis.title.y = element_text(size = 20)) + 
    xlab("Individuals") +
    xlab("") +
    ylab("Currently infected") + 
    scale_x_date() +
    guides(linetype=guide_legend(title = "Required\nisolation\ndays"),
           color=guide_legend(title =    "Required\nisolation\ndays"))+
    NULL
  print(pl)
  if(savePlots) {
    basestr = time_stamp("output/infected_with_policy_", baseline$name)
    ggsave(paste0(basestr, ".png"), width=6, height=4.5, units="in", dpi=300)  
    ggsave(paste0(basestr, ".tiff"), width=6, height=4.5, units="in", dpi=300)  
    ggsave(paste0(basestr, ".pdf"), width=6, height=4.5, units="in")  
  }
  
  return(results)
}


## 
## MAIN
##



params_def_penn_partial <- list(
  day_start       = ymd("2009-04-29"),
  day_end         = ymd("2009-06-02"),
  wb_start   = ymd("2008-12-24"), 
  wb_end     = ymd("2009-01-04"),
  sb_start   = ymd("2009-05-14"),
  sb_end     = ymd("2009-05-20"),
  sv_start   = ymd("2009-06-11"),
  sv_end     = ymd("2009-09-02"),
  exclusion_days    = 7,
  workweek_days     = 5,
  symptom_attention   = 0.9899626,
  compliance        = 0.9145737,
  vac_rate          = 0.0,  #fraction of children vaccinated
  vac_efficacy      = 0.0,  #efficiency of the vaccine in preventing illness
  seasonality       = 0.8830765,   #the maximal seasonal effect on transmissibility
  transmissibility = 0.003401314,  #the transmissibility for the outbreak
  transmissibility_weekend_ratio = 0.7000992,  #reduction during weekends.  
  transmissibility_closure_ratio = 0.1774011,
  symptom_propensity   = 0.88,
  start_date_shift = 0,
  focal_symptom        = 'fever',
  cross_grade_contact = 0.3552442
)

params_def_penn <- params_update_school(params = params_def_penn_partial,
                                        students_per_cohort=c(80,95,90,89,100,0,0),
                                        initially_infected_home=c(0,0,0,0,2,0,0),
                                        initially_infected_school=c(0,0,0,0,0,0,0),  #comes second, as in the app
                                        contact_data=default_elementary_contact_data,
                                        symptoms_data=influenza_symptoms)

params_def_thames_partial <- list(
  day_start       = ymd("2012-11-26"),
  day_end         = ymd("2012-12-21"),
  wb_start   = ymd("2008-12-24"), 
  wb_end     = ymd("2009-01-04"),
  sb_start   = ymd("2009-05-14"),
  sb_end     = ymd("2009-05-20"),
  sv_start   = ymd("2009-06-11"),
  sv_end     = ymd("2009-09-02"),
  exclusion_days    = 0,
  workweek_days     = 5,
  symptom_attention   = 0.6597556,
  compliance        = 0.1332141,
  vac_rate          = 0.80,  #fraction of children vaccinated
  vac_efficacy      = 0.50,  #efficiency of the vaccine in preventing illness
  seasonality       = 0.2043411,   #the maximal seasonal effect on transmissibility
  transmissibility = 0.002029424,  #the transmissibility for the outbreak
  transmissibility_weekend_ratio = 0.9038479,  #reduction during weekends.  
  transmissibility_closure_ratio = 0.1233679,
  symptom_propensity   = 0.84,
  focal_symptom        = 'fever',
  start_date_shift = 5,
  cross_grade_contact = 0.05683333
)

params_def_thames <- params_update_school(params = params_def_thames_partial,
                                        students_per_cohort=c(42,60,60,49,45,45,45,45),
                                        initially_infected_home=c(3,13,0,0,0,1,1,2),
                                        initially_infected_school=c(0,0,0,0,0,0,0,0),  #comes second, as in the app
                                        contact_data=default_elementary_contact_data,
                                        symptoms_data=influenza_symptoms)

params_def_boarding_partial <- list(
  day_start       = ymd("2009-05-09"),
  day_end         = ymd("2009-06-03"),
  wb_start   = ymd("2008-12-24"), 
  wb_end     = ymd("2009-01-04"),
  sb_start   = ymd("2009-05-27"),
  sb_end     = ymd("2009-06-07"),
  sv_start   = ymd("2009-06-11"),
  sv_end     = ymd("2009-09-02"),
  exclusion_days    = 0,
  workweek_days     = 5,
  vac_rate          = 0.0,  #fraction of children vaccinated
  vac_efficacy      = 0.0,  #efficiency of the vaccine in preventing illness
  symptom_propensity=0.84,
  transmissibility=0.000875389,
  transmissibility_weekend_ratio=0.5889637,
  transmissibility_closure_ratio=0.14081,
  symptom_attention=0.8863097,
  compliance=0.3663083,
  seasonality=0.4664589,
  cross_grade_contact=0.1614829,
  start_date_shift=-0.6232494,
  focal_symptom = 'fever'
)

params_def_boarding <- params_update_school(params = params_def_boarding_partial,
                                        students_per_cohort=c(258,260,262,268,259,0,0,0),
                                        initially_infected_home=c(0,0,0,0,0,0,0,0),
                                        initially_infected_school=c(1,0,0,0,0,0,0,0),  #comes second, as in the app
                                        contact_data=default_elementary_contact_data,
                                        symptoms_data=influenza_symptoms)


params_baseline_flu_partial <- list(  #based on the boarding school, with adjustments (dates, transmission etc)
  name = 'flu',
  day_start       = ymd("2020-01-01"),
  day_end         = ymd("2020-09-01"),
  wb_start   = ymd("2019-12-23"), 
  wb_end     = ymd("2020-01-05"),
  sb_start   = ymd("2020-04-06"),
  sb_end     = ymd("2020-04-10"),
  sv_start   = ymd("2020-06-17"),
  sv_end     = ymd("2020-09-08"), #SOURCE: https://chicago.suntimes.com/education/2020/2/24/21150863/cps-calendar-2020-2021-school-year-public
  exclusion_days    = 0,
  workweek_days     = 5,
  vac_rate          = 0.0,  #fraction of children vaccinated
  vac_efficacy      = 0.0,  #efficiency of the vaccine in preventing illness
  symptom_propensity=0.84,
  transmissibility=0.000402,
  transmissibility_weekend_ratio=0.280141,
  transmissibility_closure_ratio=0.141479,
  symptom_attention=0.698836,
  compliance=0.158385,
  seasonality=0.478733,
  cross_grade_contact=0.605565,
  start_date_shift=0,
  focal_symptom = 'fever'
)

#params_def_boarding_scenario <- params_update_school(params = params_def_boarding_scenario_partial,
#                                            students_per_cohort=c(70,70,70,70,70,70),
#                                            initially_infected_home=c(0,0,0,0,0,0),
#                                            initially_infected_school=c(1,0,0,0,0,0),  #comes second, as in the app
#                                            contact_data=default_elementary_contact_data,
#                                            symptoms_data=influenza_symptoms)

params_def_thames_scenario <- params_update_school(params = params_def_thames_partial,
                                            students_per_cohort=c(70,70,70,70,70,70,0,0),
                                            initially_infected_home=c(0,0,0,0,0,0,0,0),
                                            initially_infected_school=c(1,0,0,0,0,0,0,0),  #comes second, as in the app
                                            contact_data=default_elementary_contact_data,
                                            symptoms_data=influenza_symptoms)

params_def_corona_partial <- list(
  name = 'corona',
  day_start       = ymd("2020-01-01"), 
  day_end         = ymd("2020-09-01"),
  wb_start   = ymd("2019-12-23"), #SOURCE: https://www.tafths.org/ourpages/auto/2018/6/6/68132177/CPS%20Calendar%202019-20.pdf
  wb_end     = ymd("2020-01-05"),
  sb_start   = ymd("2020-04-06"),
  sb_end     = ymd("2020-04-10"),
  sv_start   = ymd("2020-06-17"),
  sv_end     = ymd("2020-09-08"), #SOURCE: https://chicago.suntimes.com/education/2020/2/24/21150863/cps-calendar-2020-2021-school-year-public
  exclusion_days    = 0,
  workweek_days     = 5, #Israel
  vac_rate          = 0.0,  #fraction of children vaccinated
  vac_efficacy      = 0.0,  #efficiency of the vaccine in preventing illness
  symptom_propensity=0.181,
  transmissibility=0.000511,
  transmissibility_weekend_ratio=0.280141,
  transmissibility_closure_ratio=0.141479,
  symptom_attention=0.698836,
  compliance=0.158385,
  seasonality=0.478733,
  cross_grade_contact=0.605565,
  start_date_shift=0,
  focal_symptom = 'fever' 
)

params_def_flu_scenario <- params_update_school(params = params_baseline_flu_partial,
                                                     students_per_cohort=c(70,70,70,70,70,70),
                                                     initially_infected_home=c(0,0,0,0,0,0),
                                                     initially_infected_school=c(1,0,0,0,0,0),  #comes second, as in the app
                                                     contact_data=default_elementary_contact_data,
                                                      symptoms_data=influenza_symptoms)

params_def_corona_scenario <- params_update_school(params = params_def_corona_partial,
                                                     students_per_cohort=c(70,70,70,70,70,70),
                                                     initially_infected_home=c(0,0,0,0,0,0),
                                                     initially_infected_school=c(1,0,0,0,0,0),  #comes second, as in the app
                                                     contact_data=default_elementary_contact_data,
                                                     symptoms_data=coronavirus_symptoms)

#here - further adjustments available from default
params_flu_seasonal_baseline <- params_def_thames_scenario
params_flu_seasonal_baseline$transmissibility <- 0.001535

#params_flu_seasonal_ranges <- param_description_ranges_thames_scenario
params_flu_seasonal_ranges <- param_description_ranges_flu_scenario
params_flu_seasonal_ranges$compliance         <- list(func=qtruncnorm, args=list(a=params_flu_seasonal_baseline$compliance*0.9, mean=params_flu_seasonal_baseline$compliance, sd=params_flu_seasonal_baseline$compliance/10, b=params_flu_seasonal_baseline$compliance*1.1))
params_flu_seasonal_ranges$symptom_attention  <- list(func=qtruncnorm, args=list(a=params_flu_seasonal_baseline$symptom_attention*0.9, mean=params_flu_seasonal_baseline$symptom_attention, sd=params_flu_seasonal_baseline$symptom_attention/10, b=params_flu_seasonal_baseline$symptom_attention*1.1))
params_flu_seasonal_ranges$transmissibility   <- list(func=qtruncnorm, args=list(a=params_flu_seasonal_baseline$transmissibility*0.9, mean=params_flu_seasonal_baseline$transmissibility, sd=params_flu_seasonal_baseline$transmissibility/10, b=params_flu_seasonal_baseline$transmissibility*1.1))
params_flu_seasonal_ranges$transmissibility_weekend_ratio <- list(func=qtruncnorm, args=list(a=params_flu_seasonal_baseline$transmissibility_weekend_ratio*0.9, mean=params_flu_seasonal_baseline$transmissibility_weekend_ratio, sd=params_flu_seasonal_baseline$transmissibility_weekend_ratio/10, b=params_flu_seasonal_baseline$transmissibility_weekend_ratio*1.1))
params_flu_seasonal_ranges$transmissibility_closure_ratio <- list(func=qtruncnorm, args=list(a=params_flu_seasonal_baseline$transmissibility_closure_ratio*0.9, mean=params_flu_seasonal_baseline$transmissibility_closure_ratio, sd=params_flu_seasonal_baseline$transmissibility_closure_ratio/10, b=params_flu_seasonal_baseline$transmissibility_closure_ratio*1.1))
params_flu_seasonal_ranges$seasonality <- list(func=qtruncnorm, args=list(a=params_flu_seasonal_baseline$seasonality*0.9, mean=params_flu_seasonal_baseline$seasonality, sd=params_flu_seasonal_baseline$seasonality/10, b=params_flu_seasonal_baseline$seasonality*1.1))
params_flu_seasonal_ranges$cross_grade_contact <- list(func=qtruncnorm, args=list(a=params_flu_seasonal_baseline$cross_grade_contact*0.9, mean=params_flu_seasonal_baseline$cross_grade_contact, sd=params_flu_seasonal_baseline$cross_grade_contact/10, b=params_flu_seasonal_baseline$cross_grade_contact*1.1))
params_flu_seasonal_ranges$symptom_propensity <- list(func=qtruncnorm, args=list(a=0.6, mean=0.84, sd=0.1, b=0.95))

#NOTE: we are assuming no vaccine, which doesn't matter, as we know

lhs_output <- compute_policy_table(num_samples=499,
                                     param_descriptions=params_flu_seasonal_ranges,
                                     baseline_sample=params_flu_seasonal_baseline,
                                     scenario_name='flu_seasonal_FI',
                                     max_policy=9
)

lhs_output <- compute_workday_table(num_samples=499,
                                   param_descriptions=params_flu_seasonal_ranges,
                                   baseline_sample=params_flu_seasonal_baseline,
                                   scenario_name='flu_seasonal_SW',
                                   max_workdays=6
)

#Figure 1 of the paper
#z_policy_combined_fig(baseline=params_def_flu_scenario,
#                      exclusion_set=c(0,1,2,3,4,5,6),
#                      savePlots = TRUE)

#z_policy_combined_fig(baseline=params_def_corona_scenario,
#                      exclusion_set=c(0,1,2,4,6,8,10),
#                      savePlots = TRUE)

#z_workday_fig(baseline=params_def_flu_scenario,
#              workdays_set=c(0,1,2,3,4,5),
#              savePlots = TRUE)

#z_workday_fig(baseline=params_def_corona_scenario,
#              workdays_set=c(0,1,2,3,4,5),
#              savePlots = TRUE)

#z_policy_combined_fig(baseline=params_def_scenario,
#                      exclusion_set=c(0,1,2,4,6),
#                      savePlots = TRUE)


