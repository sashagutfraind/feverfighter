#'
#' School Outbreaks Predictor ("FeverFighter")
#' (c) 2016-2017 By Authors.
#' You are free to use this software under CC-BY 4.0 https://creativecommons.org/licenses/by/4.0/
#'

#'
#' Optimization of a school flu outbreak and policies data
#' 
library(lubridate)
library(genalg)

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
  setwd("C:\\Users\\Adam\\Projects\\schools-flu\\code")
}
if(!file.exists("output")) {
  dir.create("output")
}

print('Path configured..')


print('sourcing the model code...')
source('schools_solver.R', local = T);
print('source completed successfully!')

gene_list <- data.frame('run_num' = numeric(),
                    'transmissibility' = numeric(),
                    'transmissibility_weekend_ratio' = numeric(),
                    'transmissibility_closure_ratio' = numeric(),
                    'fever_attention' = numeric(),
                    'compliance' = numeric(),
                    'seasonality' = numeric(),
                    'cross-grade contacts' = numeric(),
                    'start date shift' = numeric(),
                    'retVal' = numeric())

evaluate_penn <- function(soln=c(), return_params=F) {
  returnVal = NA;
  if(length(soln) != 8) {
    stop("Expecting a chromosome of length 8!");
  }
  cal_students_per_grade=c(80,95,90,89,100,0,0)
  cal_initial_school    =c(0,0,0,0,2,0,0)
  cal_initial_home      =c(0,0,0,0,0,0,0)  #comes second, as in the app
  cal_params_partial <- list(
    day_start       = ymd("2009-05-01") + round(soln[8]),
    day_end         = ymd("2009-06-02"),
    wb_start   = ymd("2008-12-24"), 
    wb_end     = ymd("2009-01-04"),
    sb_start   = ymd("2009-05-14"),
    sb_end     = ymd("2009-05-20"),
    sv_start   = ymd("2009-06-11"),
    sv_end     = ymd("2009-09-02"),
    exclusion_days    = 7,
    workweek_days     = 5,
    #fever_attention   = 0.3,
    #compliance        = 0.5,
    vac_rate          = 0.00,  #fraction of children vaccinated
    vac_efficacy      = 0.00,  #efficiency of the vaccine in preventing illness
    #seasonality       = 0.3,   #the maximal seasonal effect on transmissibility
    transmissibility               = soln[1],  #the transmissibility for the outbreak
    transmissibility_weekend_ratio = soln[2],  #reduction during weekends
    transmissibility_closure_ratio = soln[3],
    symptom_attention                = soln[4],
    compliance                     = soln[5],
    seasonality                    = soln[6],
    symptom_propensity   = 0.84
    )
  new_contact_data <- update_off_diagonal_values(default_elementary_contact_data, 
                                     soln[7])
  cal_params <- params_update_school(params = cal_params_partial, 
                       students_per_cohort=cal_students_per_grade, 
                       initially_infected_home=cal_initial_home, 
                       initially_infected_school=cal_initial_school, 
                       contact_data=new_contact_data,
                       symptoms_data=influenza_symptoms)
  
  if(return_params) {
    return(cal_params)
  }
  returnVal = calibration_difference(new_params = cal_params, 
                                     params = params_def_elementary, 
                                     calibration_outbreak = calibration_outbreak_results_pennsylvania, 
                                     combine_weekend_to_monday = FALSE)

  cat(sprintf("R2=%.*f \t:", 2, returnVal), soln, "\n")
  return(returnVal)
}


evaluate_thames <- function(soln=c(), return_params=F) {
  returnVal = NA;
  if(length(soln) != 8) {
    stop("Expecting a chromosome of length 8!");
  }
  #NOTE: List size for school cohort is 8, may need to change some code in schools_solver.R
  cal_students_per_grade=c(42,60,60,49,45,45,45,45)
  #cal_initial_school    =c(3,13,0,0,0,1,1,2)
  cal_initial_school    =c(0,1,0,0,0,0,0,0)
  cal_initial_home      =c(0,0,0,0,0,0,0,0)  #comes second, as in the app
  cal_params_partial <- list(
    day_start       = ymd("2012-11-26") + round(soln[8]),
    day_end         = ymd("2012-12-21"),
    wb_start   = ymd("2008-12-24"), 
    wb_end     = ymd("2009-01-04"),
    sb_start   = ymd("2009-05-27"),
    sb_end     = ymd("2009-06-07"),
    sv_start   = ymd("2009-06-11"),
    sv_end     = ymd("2009-09-02"),
    exclusion_days    = 0,
    workweek_days     = 5,
    vac_rate          = 0.00,  #fraction of children vaccinated
    vac_efficacy      = 0.00,  #efficiency of the vaccine in preventing illness
    transmissibility               = soln[1],  #the transmissibility for the outbreak
    transmissibility_weekend_ratio = soln[2],  #reduction during weekends
    transmissibility_closure_ratio = soln[3],
    symptom_attention                = soln[4],
    compliance                     = soln[5],
    seasonality                    = soln[6],
    focal_symptom                  = 'fever',
    symptom_propensity   = 0.84 #based on Cowling et al.
  )                  
  new_contact_data <- update_off_diagonal_values(default_elementary_contact_data, 
                                                 soln[7])
  cal_params <- params_update_school(params = cal_params_partial, 
                                     students_per_cohort=cal_students_per_grade, 
                                     initially_infected_home=cal_initial_home, 
                                     initially_infected_school=cal_initial_school, 
                                     contact_data=new_contact_data,
                                     symptoms_data=influenza_symptoms)
  
  if(return_params) {
    return(cal_params)
  }
  returnVal = calibration_difference(new_params = cal_params, 
                                     params = params_def_elementary, 
                                     calibration_outbreak = calibration_outbreak_results_thames, 
                                     combine_weekend_to_monday = TRUE)
  
  cat(sprintf("R2=%.*f \t:", 2, returnVal), soln, "\n")
  return(returnVal)
}

evaluate_boarding <- function(soln=c(), return_params=F) {
  returnVal = NA;
  if(length(soln) != 8) {
    stop("Expecting a chromosome of length 8!");
  }
  cal_students_per_grade=c(258,260,262,268,259,0,0)
  cal_initial_school    =c(1,0,0,0,0,0,0)
  cal_initial_home      =c(0,0,0,0,0,0,0)  #comes second, as in the app
  cal_params_partial <- list(
    day_start       = ymd("2009-05-09") + round(soln[8]),
    day_end         = ymd("2009-06-15"),
    wb_start   = ymd("2008-12-24"), 
    wb_end     = ymd("2009-01-04"),
    sb_start   = ymd("2009-05-27"),
    sb_end     = ymd("2009-06-07"),
    sv_start   = ymd("2009-06-11"),
    sv_end     = ymd("2009-09-02"),
    exclusion_days    = 0,
    workweek_days     = 5,
    vac_rate          = 0.00,  #fraction of children vaccinated
    vac_efficacy      = 0.00,  #efficiency of the vaccine in preventing illness
    transmissibility               = soln[1],  #the transmissibility for the outbreak
    transmissibility_weekend_ratio = soln[2],  #reduction during weekends
    transmissibility_closure_ratio = soln[3],
    symptom_attention                = soln[4],
    compliance                     = soln[5],
    seasonality                    = soln[6],
    focal_symptom                  = 'fever',
    symptom_propensity   = 0.84 
  )
  new_contact_data <- update_off_diagonal_values(default_elementary_contact_data, 
                                                 soln[7])
  cal_params <- params_update_school(params = cal_params_partial, 
                                     students_per_cohort=cal_students_per_grade, 
                                     initially_infected_home=cal_initial_home, 
                                     initially_infected_school=cal_initial_school, 
                                     contact_data=new_contact_data,
                                     symptoms_data=influenza_symptoms)
  
  if(return_params) {
    return(cal_params)
  }
  returnVal = calibration_difference(new_params = cal_params, 
                                     params = params_def_elementary, 
                                     calibration_outbreak = calibration_outbreak_results_boarding, 
                                     combine_weekend_to_monday = FALSE)
  gene_list[nrow(gene_list)+1,] <<- append(append((nrow(gene_list)+1), soln), returnVal)
  cat(sprintf("R2=%.*f \t:", 2, returnVal), soln, "\n")
  return(returnVal)
}


##To change calibration data set, swap below function. Search for "To change calibration data set" for others.
##This will produce multiple starts to the genetic algorithm
for (i in seq(10)) {
  rbga.results = rbga(c(0,     0.1, 0.0, 0.1, 0.1, 0.0, 0.0, -14), 
                      c(0.005, 1.0, 0.5, 1.0, 1.0, 1.0, 0.8,  0),  #start day cannot be later than earliest detection
                      monitorFunc = NULL, 
                      evalFunc = evaluate_boarding,
                      verbose = T, 
                      popSize = 80, 
                      iters   = 50
                      )
}

#print(rbga.results)

optimal_solution        <- gene_list[which.min(gene_list$retVal),]
#rbga.results$population[which.min(rbga.results$evaluations), ]
##To change calibration data set, swap below function. Search for "To change calibration data set" for others.
#optimal_solution_params <- evaluate_boarding(optimal_solution, return_params = T)
gene_range <- gene_list[which(gene_list$retVal <= (optimal_solution$retVal * 1.10)),]
##Calculate the 1st and 3rd quartile range
gene_25pct <- apply(gene_range, 2, quantile, 0.25)
gene_75pct <- apply(gene_range, 2, quantile, 0.75)
gene_median <- apply(gene_range, 2, median)
gene_sd <- apply(gene_range, 2, sd, FALSE)
gene_min <- apply(gene_range, 2, min)
gene_max <- apply(gene_range, 2, max)

calibration = data.table(rbind(optimal_solution, gene_min, gene_max, gene_sd, gene_median, gene_25pct, gene_75pct))
row.names(calibration) <- c("opt", "min", "max", "sd", "med", "pct25", "pct75")
fname <- paste(time_stamp("output/calibration", ".csv"))
fwrite(calibration, file = fname, row.names = T)
print(sprintf("Saved calibration output: %s", fname))


print("Optimal solution:")
print(optimal_solution)

print(calibration)

#evaluate_boarding(as.double(unname(optimal_solution[2:9])))
