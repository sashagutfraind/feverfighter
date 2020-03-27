#'
#' Outbreaks Predictor and Policy Analyzer ("FeverFighter")
#' (c) 2016-2020 By Authors.
#' You are free to use this software under CC-BY 4.0 https://creativecommons.org/licenses/by/4.0/
#'
#' Model of a respiratory outbreak and policies, motivated by outbreaks in schools
#' 
#' 

#options(verbose=T)
#you will need:
#install.packages(c('truncnorm', 'data.table', 'lubridate', 'chron', 'ggplot2', 'shiny', 'lhs', 'markdown', 'genalg', 'testthat', 'lhs'))

library(ggplot2)
library(scales)
library(reshape2)
library(chron) #allows sub-day resolution
library(lubridate)
library(data.table)
library(tibble)

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

####################################### Parameters ################################
# TODO: to document; note differences with the paper

if(! ("influenza_symptoms" %in% ls())) {
  influenza_symptoms <- read.csv("influenza_symptoms.csv")
}
if(! ("coronavirus_symptoms" %in% ls())) {
  coronavirus_symptoms <- read.csv("coronavirus_symptoms.csv")
}
if(! ("default_elementary_contact_data" %in% ls())) {
  default_elementary_contact_data <- read.csv('contact_data_default_K-G6.csv', header=T, row.names=1)
}

if(! ("pennsylvania_time_series" %in% ls())) {
  pennsylvania_time_series <- fread('pennsylvania_time_series.csv')
}
if(! ("thames_time_series" %in% ls())) {
  thames_time_series <- fread('thames_time_series.csv')
}
if(! ("thames_time_series_weekend" %in% ls())) {
  thames_time_series_weekend <- fread('thames_time_series_weekend.csv')
}
if(! ("boarding_time_series" %in% ls())) {
  boarding_time_series <- fread('boarding_time_series.csv')
}

get_calibration_dataset <- function(dataset_name) {
  if(dataset_name == "Smith pH1H1 2009") {
    return(boarding_time_series) 
  } else if (dataset_name == "Marchbanks pH1N1 2009") {
    return(pennsylvania_time_series) 
  } else if (dataset_name == "McCann Influenza B 2012") {
    return(thames_time_series) 
  } else {
    return(NULL)
  } 
}

# default params.  
params_outbreak_partial <- list(
  day_start       = ymd("2020-01-01"),# chron() format="Y-m-d"),  #happens to be Monday
  day_end         = ymd("2020-12-31"),# chron() format="Y-m-d"),
  wb_start   = ymd("2020-12-24"), 
  wb_end     = ymd("2021-01-04"),
  sb_start   = ymd("2020-04-15"),
  sb_end     = ymd("2020-04-19"),
  sv_start   = ymd("2020-06-19"),
  sv_end     = ymd("2020-09-02"),
  exclusion_days    = 1,
  symptom_attention = 0.3,   
  compliance        = 0.5,
  vac_rate          = 0.4,  #fraction of children vaccinated
  vac_efficacy      = 0.2,  #efficiency of the vaccine in preventing illness
  seasonality       = 0.4,   #the maximal seasonal effect on transmissibility
  transmissibility = 0.002,  #the transmissibility for the outbreak
  transmissibility_weekend_ratio = 0.1,  #reduction during weekends.  TODO: check literature
  transmissibility_closure_ratio = 0.1,
  focal_symptom = 'fever',
  symptom_propensity   = 0.88 
)



#' @desc update the more complex parameters about the school and outbreak
#' @param params             information about the outbreak
#' @param students_per_cohort vector of students per grade  
#' @param initially_infected_school vector of students that are infected and at school on day 1
#' @param initially_infected_home   vector of students that are infected and at home on day 1
#' @param contact_data       the structure of the school. a square matrix with rows and columns named by grade.  the diagonal is 1.0.
#' @param symptoms_data      a list of vectors, SystemicSymptoms, ViralShedding, RespiratorySymptoms, New_ILI_rate e.g.: c(0.1, 0.4, 0.2) etc.  systemic and respiratory need to be in [0..1] on all days 
params_update_school <- function(params, 
                                 students_per_cohort, 
                                 initially_infected_home, 
                                 initially_infected_school, 
                                 contact_data,
                                 symptoms_data) {
  contact_data           = contact_data[1:length(students_per_cohort), 1:length(students_per_cohort)]
  params$num_stages      = dim(symptoms_data)[1]
  if(params$focal_symptom == 'fever') {
    params$symptom_rates = symptoms_data$SystemicSymptoms
  } else if(params$focal_symptom == 'cough') {
    params$symptom_rates = symptoms_data$RespiratorySymptoms
  } else if(params$focal_symptom == 'sneezing') {
    params$symptom_rates = symptoms_data$NasalSymptoms
  } else { 
    stop("invalid focal symptom")
  }
  stopifnot(all(params$symptom_rates >= 0) & all(params$symptom_rates <= 1))
  
  params$New_ILI_Rates   = symptoms_data$New_ILI_Rates
  params$shedding_rates   = symptoms_data$ViralShedding
  params$SystemicSymptoms = symptoms_data$SystemicSymptoms  #should only be used for calculating daily ILI rates
  #generally, linked to log titers.  multiplies with the transmission
  
  params$cohort_names    = names(contact_data)
  params$num_cohorts     = length(params$num_cohorts)
  stopifnot(row.names(contact_data) == params$cohort_names)
  params$contact_rates   = contact_data
  
  if (!is.null(params$cross_grade_contact)) {
    params$contact_rates <- update_off_diagonal_values(contact_data, params$cross_grade_contact)
    #wishlist - make sure the warning does not interfere with sensitivity analysis
    #warning('Conflicting parameters - contact rates overriden with params$cross_grade_contact')
  }
  
  params$cohort_sizes            = structure(students_per_cohort,      names=params$cohort_names)
  params$initial_infected_home   = structure(initially_infected_home,   names=params$cohort_names) 
  params$initial_infected_school = structure(initially_infected_school, names=params$cohort_names)

  return(params)
}

params_def_elementary <- params_update_school(params = params_outbreak_partial,
                                                    students_per_cohort=c(90, 90, 90, 90, 90, 90, 90, 90),
                                                    initially_infected_home=c(0, 0, 0, 0, 0, 0, 0, 0),
                                                    initially_infected_school=c(1, 0, 0, 0, 0, 0, 0, 0),
                                                    contact_data=default_elementary_contact_data,
                                                    symptoms_data=influenza_symptoms)



######################### ######################### ######################### ###########################
#' @desc MAIN FUNCTION Basic model of school disease outbreak
#' @param params = parameters of the outbreak and of the school
#' @return data.frame containing the day-by-day outbreak sizes
#' 
compute_school_outbreak <- function(params) {
  current_day <- params$day_start

  stages        <- seq(params$num_stages)
  cohort_names  <- params$cohort_names
  cohort_cross_stage <- expand.grid(stages, cohort_names)
  cohort_cross_stage <- sprintf("%s.%s", as.character(cohort_cross_stage$Var2), as.character(cohort_cross_stage$Var1))
  varm <- list()  #for each class, we have a matrix of stages cross levels.  names use the convention "A_l.s" (b/c the paper's "Al,s" creates invalid name)
  varm[["S"]] <- structure(sprintf("S_%s", cohort_names), names=cohort_names)
  # Suseptible
  varm[["VS"]] <- structure(sprintf("VS_%s", cohort_names), names=cohort_names)
  # Vaccinated but ineffectually (remaining susceptible)
  varm[["VR"]] <- structure(sprintf("VR_%s", cohort_names), names=cohort_names)
  # Vaccinated but truly resilient
  varm[["H"]] <- structure(sprintf("H_%s", cohort_cross_stage), names=cohort_cross_stage)
  # Infected && Home
  varm[["I"]] <- structure(sprintf("I_%s", cohort_cross_stage), names=cohort_cross_stage)
  # Infected && In School
  varm[["R"]] <- structure(sprintf("R_%s", cohort_names), names=cohort_names)
  # Recovered

  varnames <- as.character(c(varm[["S"]], varm[["VS"]], varm[["VR"]], varm[["H"]], varm[["I"]], varm[["R"]]))

  p <- function(..., sep='') {paste(..., sep=sep, collapse=sep)}
  
  init_state <- rep(0, length(varnames))
  names(init_state) <- varnames
  
  
  
  for(l in params$cohort_names) {
    init_state[varm[["H"]][[p(l,".",1)]]] <- params$initial_infected_home[[l]]
    init_state[varm[["I"]][[p(l,".",1)]]] <- params$initial_infected_school[[l]]
    uninfected_cohort_size <- max(0, params$cohort_sizes[[l]] - params$initial_infected_home[[l]] - params$initial_infected_school[[l]])
    init_state[varm[["S"]][[l]]]  <- uninfected_cohort_size * (1-params$vac_rate)                          #susceptible, including vaccinated ineffectually
    init_state[varm[["VS"]][[l]]] <- uninfected_cohort_size * (params$vac_rate)*(1-params$vac_efficacy)  #vaccinated but remain susceptible
    init_state[varm[["VR"]][[l]]] <- uninfected_cohort_size * (params$vac_rate)*(  params$vac_efficacy)  #vaccinated and fully removed
  }

  sim_dates<-seq.dates(as.numeric(params$day_start), as.numeric(params$day_end))
  
  newparameters <- params
  newparameters$return_rates <- effective_return_rates(params)

  #compute the constant portion of the rates of infection from grade to grade  
  infection_stencils <- matrix(0, nrow=length(params$cohort_names), ncol=length(varnames))
  dimnames(infection_stencils)<-list(params$cohort_names, varnames)
  for(l1 in params$cohort_names) {
    for(l2 in params$cohort_names) {
      for(s in stages) {
        #only "I" states transmit - no transmission from other states
        #shedding_rates is defined as the log titers
        infection_stencils[l1,varm[["I"]][[p(l2,".",s)]]] <- params$shedding_rates[s]*params$contact_rates[l1,l2]
      }
    }
  }
  

  #run the calculation for days 2, 3, ..., sim_dates
  outdf <- init_state
  state <- init_state
  for(current_day in sim_dates[2:length(sim_dates)]) {
    new_state <- compute_school_helper(current_day, state, newparameters, infection_stencils, varm)
    outdf <- rbind(outdf, new_state)
    state <- new_state
  }
  outdf <- data.frame(outdf, row.names=sim_dates)
  rownames(outdf) <- sim_dates
  outdf$dates <- sim_dates
  
  outdf$InfectedIsolated   <- 0  
  outdf$InfectedUnisolated <- 0 
  for(l in params$cohort_names) {
    outdf[,p("InfectedIsolated", l)]   <- rowSums(outdf[,varm[["H"]][sprintf("%s.%s", l, stages)]])
    outdf[,p("InfectedUnisolated", l)] <- rowSums(outdf[,varm[["I"]][sprintf("%s.%s", l, stages)]])
    outdf$InfectedIsolated   <- outdf$InfectedIsolated   + outdf[,p("InfectedIsolated", l)]
    outdf$InfectedUnisolated <- outdf$InfectedUnisolated + outdf[,p("InfectedUnisolated", l)]
    outdf[,p("Infected", l)] <- outdf[,p("InfectedIsolated", l)] + outdf[,p("InfectedUnisolated", l)]
  }

  newly_absent_rate <- cumprod(shift(newparameters$return_rates, fill=1))*(1-newparameters$return_rates)
  #for those infected on day s, how many will be absent on day s for the first time
  #test: 
  #   never absent == cumprod(newparameters$return_rates)[9] == 1-sum(newly_absent_rate)
  #   sum(outdf$New_absent) == approximately sum(outdf$Infected.1)
  
  outdf$ILI <- 0  #total with ILI.  calculated by summing over stages, weighted by the symptom rate of each stage
  outdf$New_absent <- 0  #newly absent on a given day.  calculated by summing over stages, weighted by the symptom rate of each stage
  outdf$New_ILI <- 0
  for(s in stages) {
    daily_total_state_s       <- rowSums(outdf[,varm[["H"]][sprintf("%s.%s", params$cohort_names, s)]]) + 
                                 rowSums(outdf[,varm[["I"]][sprintf("%s.%s", params$cohort_names, s)]])
    outdf[,p("Infected.", s)] <- daily_total_state_s 
    outdf$ILI                 <- outdf$ILI +  daily_total_state_s*params$SystemicSymptoms[s]*params$symptom_propensity
    outdf$New_absent          <- outdf$New_absent +  newly_absent_rate[s]*outdf[,p("Infected.", s)]
  }

  for(s in stages) {
    outdf$New_ILI <- outdf$New_ILI + params$New_ILI_Rates[s]*params$symptom_propensity*outdf[,p("Infected.", s)]
  }  

  outdf$S     <- rowSums(outdf[,varm[["S"]][sprintf("%s", params$cohort_names)]])
  outdf$R     <- rowSums(outdf[,varm[["R"]][sprintf("%s", params$cohort_names)]])
  outdf$VS     <- rowSums(outdf[,varm[["VS"]][sprintf("%s", params$cohort_names)]])
  outdf$VR     <- rowSums(outdf[,varm[["VR"]][sprintf("%s", params$cohort_names)]])
  outdf$Infected <- outdf$InfectedIsolated + outdf$InfectedUnisolated
  
  return(outdf)
}

#' @desc Runs a single day for the school
#' 
compute_school_helper <- function(current_day, state, parameters, infection_stencils, varm) {
  stopifnot(state >= -1E-5)
  p <- function(..., sep='') {paste(..., sep=sep, collapse=sep)}
  
  #wishlist: precompute this
  #represents heighted contact rates in the winter - peakes on Jan 1 and declines on either side of Jan 1
  seasonality_multiplier <- 1 + parameters$seasonality*cos(2*pi*yday(chron(current_day))/365)  
  base_transmissibility <- parameters$transmissibility*seasonality_multiplier
  
  transmissibility <-       ifelse(is.holiday(current_day,
                                              c(seq.dates(as.numeric(parameters$wb_start), as.numeric(parameters$wb_end)), 
                                                seq.dates(as.numeric(parameters$sb_start), as.numeric(parameters$sb_end)),
                                                seq.dates(as.numeric(parameters$sv_start), as.numeric(parameters$sv_end)))),
                                   parameters$transmissibility_closure_ratio*base_transmissibility, 
                                   ifelse(is.weekend(current_day), 
                                          parameters$transmissibility_weekend_ratio*base_transmissibility,
                                          base_transmissibility)
                                   )
                            #'extended ifelse function: if weekend, calculate transmissibility (as normal)
                            #'                            if not, check if holiday
                            #'                              if holiday, calculate transmissibility, otherwise return base transmissibility
  
  
  new_state <- state
  for (l in parameters$cohort_names) {
    SUSC <- varm[["S"]][[l]]
    VS <- varm[["VS"]][[l]]  #vaccinated who are still susceptible
    new_infections_S  <- transmissibility*state[SUSC]*sum(infection_stencils[l,]*state)
    new_infections_VS <- transmissibility*state[VS]*sum(infection_stencils[l,]*state)
    #sometimes we see more predicted daily infected than remaining susceptibles
    if(new_infections_S > new_state[SUSC]) {
      new_infections_S = new_state[SUSC]
    }
    if(new_infections_VS > new_state[VS]) {
      new_infections_VS = new_state[VS]
    }
    new_state[SUSC] <- new_state[SUSC] - new_infections_S
    new_state[VS]   <- new_state[VS]   - new_infections_VS
    
    new_state[varm[["I"]][[p(l,".","1")]]] <- new_infections_S + new_infections_VS #go to school b/c no symptoms yet
    new_state[varm[["H"]][[p(l,".","1")]]] <- 0

    for(s in 2:parameters$num_stages) {
      infected <- state[varm[["H"]][[p(l,".",s-1)]]] + state[varm[["I"]][[p(l,".",s-1)]]]
      new_state[varm[["H"]][[p(l,".",s)]]] <- infected*(1-parameters$return_rates[[s]])
      new_state[varm[["I"]][[p(l,".",s)]]] <- infected*parameters$return_rates[[s]]
    }
    new_state[varm[["R"]][[l]]] <- state[varm[["R"]][[l]]] + 
                                    state[varm[["H"]][[p(l,".",parameters$num_stages)]]]+ 
                                    state[varm[["I"]][[p(l,".",parameters$num_stages)]]]
  }
  return(new_state)
}

#' @desc thin wrapper for running the school model. allows overwriting of select parameters
compute_school_outbreak_alt <- function(params        = params_def_elementary, 
                                        compute_func  = compute_school_outbreak,
                                        alt_params    = list(), 
                                        doplot=FALSE) {
  for (p in names(alt_params)) {
    stopifnot(p %in% names(params_def_elementary))
    params[p] <- alt_params[p]
  }
  outdf <- compute_func(params = params)
  #outdf<-data.frame(out)
  #outdf$dates <- as.Date(outdf$dates)
  #outdf$Infected <- outdf$S[1] - outdf$S + outdf$I[1]
  #rownames(outdf)<-outdf[,"time"]
  if(doplot) {
    plot_basic(outdf)
  }
  
  return(outdf)
}

#'
#' @desc Compare the differences of outcomes of parameter changes
#' @param baseline output from a simulation run
#' @param alt output from an alternative simulation run
#' 
diff_scenarios <- function(baseline, alt) {
  ret <- list()
  
  base_S_final <-tail(baseline$S,1)
  alt_S_final <-tail(alt$S,1)
  
  base_R_final <-tail(baseline$R,1)
  alt_R_final <-tail(alt$R,1)
  
  base_AR <- 100*base_R_final/head(baseline$S,1)
  alt_AR  <- 100*alt_R_final/head(alt$S,1)
  
  ret$Sbase <- sprintf("%1.0f", base_S_final)
  ret$Salt  <- sprintf("%1.0f", alt_S_final)
  ret$Spcchange <- sprintf("%1.0f", 100*(alt_S_final-base_S_final)/base_S_final)
  ret$Sstring   <- sprintf("%s(%s)", ret$Salt, ret$Spcchange)
  
  ret$Rbase <- sprintf("%1.0f", base_R_final)
  ret$Ralt  <- sprintf("%1.0f", alt_R_final)
  ret$Rpcchange <- sprintf("%1.0f", 100*(alt_R_final-base_R_final)/base_R_final)
  ret$Rstring   <- sprintf("%s(%s)", ret$Ralt, ret$Rpcchange)

  ret$ARbase <- sprintf("%1.0f", base_AR)
  ret$ARalt  <- sprintf("%1.0f", alt_AR)
  ret$ARpcchange <- sprintf("%1.0f", 100*(alt_AR-base_AR)/base_AR)
  ret$ARstring   <- sprintf("%s(%s)", ret$ARalt, ret$ARpcchange)
  
  return(ret)
}

#' @desc Calculate the fraction of children that return based on policies   
#' @return probability, indexed by day of infection (index=2 is 24-48 hrs from exposure)
#' @params params$symptom_rates the underlying rate of the focal symptom
#' @params params$compliance     fraction of kids not allowed to return
#' @params params$exclusion_days     number of days blocked
#' @params params$symptom_propensity fraction of cases with apparent symptoms
#' 
#' 
effective_return_rates <- function(params) {
  if(params$exclusion_days==0) {
    #parents look at previous day's symptoms
    return_rates = 1 - (params$symptom_attention*params$symptom_propensity*shift(params$symptom_rates, n=1, fill=0))
    return(return_rates)
  }
  #the compliance parameter increases the attention from its normal value to a maximum of 1  
  #the result is modified by symptom_propensity - when it's low, the policy fails
  adjusted_attention <- params$symptom_propensity*
                          (params$symptom_attention + (1-params$symptom_attention)*params$compliance)
  
  if(params$exclusion_days>length(params$symptom_rates)) {
    stop('Invalid number of excluded days')
  }
  temp_max <- rep(0, length(params$symptom_rates))
  for (day in 1:params$exclusion_days) {
    temp_max <- pmax(shift(params$symptom_rates, n=day, fill=0), temp_max)
  }
  return_rates <- 1 - (adjusted_attention*temp_max)
  return(return_rates)
}

#' @desc Calculate the effective return rates for different policies
effective_return_rates_table <- function(symptom_rates, symptom_propensity=0.88, symptom_attention=0.832331) {
  
  
  #' @params params$symptom_rates the underlying rate of the focal symptom
  #' @params params$compliance     fraction of kids not allowed to return
  #' @params params$exclusion_days     number of days blocked
  #' @params params$symptom_propensity fraction of cases with apparent symptoms
  
  #Return rates: A = no policy rate and
  p0 = list(symptom_rates=symptom_rates,
            symptom_propensity=symptom_propensity,
            symptom_attention=symptom_attention,
            compliance=0,
            exclusion_days=0
            )
  
  header <- c('')
  rows <- c("Policy", seq(length(symptom_rates)))
  row_open <- c("Open policy", sprintf("%1.3f", effective_return_rates(p0)))
  
  
  #B = 1 day policy rate with 0.5 compliance
  p0 = list(symptom_rates=symptom_rates,
            symptom_propensity=symptom_propensity,
            symptom_attention=symptom_attention,
            compliance=0.5,
            exclusion_days=1
  )
  
  header <- c('')
  rows <- c("Policy", seq(length(symptom_rates)))
  row_B <- c("1 day 0d5 compliance", sprintf("%1.3f", effective_return_rates(p0)))
  
  #C = baseline scenario parameters with 1 day policy rate and 100% compliance.
  #B = 1 day policy rate with 0.5 compliance
  p0 = list(symptom_rates=symptom_rates,
            symptom_propensity=symptom_propensity,
            symptom_attention=symptom_attention,
            compliance=1.0,
            exclusion_days=1
  )
  
  header <- c('')
  rows <- c("Policy", seq(length(symptom_rates)))
  row_C <- c("1 day 1d0 compliance", sprintf("%1.3f", effective_return_rates(p0)))
  
  lxmatrix <- t(rbind(row_open, row_B, row_C))
    
  fname <- paste(time_stamp("output/return_rates", ".csv"))
  fwrite(lxmatrix, file=fname)

  #lxmatrix <- matrix(rows, dim(rows)[1], dim(rows)[2])
  #write.table(lxmatrix, file=fname, sep=",")
  print(sprintf("Saved matrix: %s", fname))
  
  return(lxmatrix)
}



#' @desc plot the outcome of an outbreak
#' @param outdf a data.frame from an outbreak 
plot_base <- function(outdf) {
  df<-data.frame(dates=as.Date(outdf$dates), 
                 #InfectedK=outdf$InfectedK,
                 #InfectedG1=outdf$InfectedG1,
                 #InfectedG2=outdf$InfectedG2,
                 #InfectedG3=outdf$InfectedG3,
                 #InfectedG4=outdf$InfectedG4,
                 #InfectedG5=outdf$InfectedG5,
                 #InfectedG6=outdf$InfectedG6,
                 Susceptible=outdf$S,
                 Infected=outdf$Infected,
                 ILI=outdf$ILI,
                 New_ILI=outdf$New_ILI,
                 Isolated=outdf$InfectedIsolated,
                 Unisolated=outdf$InfectedUnisolated,
                 Recovered=outdf$R
  )
  
  df_long <- reshape2::melt(df, id="dates") # convert to long format
  
  #pl<-ggplot(data=df_long, aes(x=dates, y=value, size=variable, color=variable), geom=c("point", "path")) +
  pl<-ggplot(data=df_long, aes(x=dates))
  pl <- pl +
    geom_line(data=df, aes(y = Isolated,   color="Isolated"), size=1) +
    geom_line(data=df, aes(y = Unisolated,   color="Unisolated"), size=1.5) +
    geom_line(data=df, aes(y = Infected,   color="Infected")) +
    geom_line(data=df, aes(y = ILI,      color="ILI")) +
    geom_line(data=df, aes(y = New_ILI,  color="New_ILI")) +
    geom_line(data=df, aes(y = Susceptible,   color="Susceptible")) + 
    geom_line(data=df, aes(y = Recovered,   color="Recovered"), size=2) +
    #linetype=variable, 
    ylab("") +
    xlab("") + 
    guides(color=guide_legend(title=""))+
    #guides(variable=guide_legend(title=""))+
    NULL
  print(pl)
  ggsave(time_stamp("output/cases_vs_time", ".pdf"), width=6, height=4.5, units="in")  
  ggsave(time_stamp("output/cases_vs_time", ".png"), width=6, height=4.5, units="in")  
}




#' @desc Basic plot of cases
plot_basic <- function(outdf, showNow=F) {
  df<-data.frame(dates=as.Date(outdf$dates), 
                 Susceptible=outdf$S,
                 Infected=outdf$Infected,
                 ILI=outdf$ILI,
                 New_ILI=outdf$New_ILI,
                 Isolated=outdf$InfectedIsolated,
                 Unisolated=outdf$InfectedUnisolated,
                 Recovered=outdf$R
  )
  df_long <- reshape2::melt(df, id="dates") # convert to long format
  pl<-ggplot(data=df_long, aes(x=dates, y=value, linetype=variable, color=variable), geom=c("point", "path")) +
    geom_line() + #geom_point() + 
    ylab("Individuals") +
    xlab("") + 
    guides(variable=guide_legend(title=""))+
    NULL
  if(showNow) {
    print(pl)
    ggsave(time_stamp("output/cases_vs_time", ".pdf"), width=6, height=4.5, units="in")  
  }
  return(pl)
}

savePlots <- F #used for the paper.  generally, False 

#' @desc plot of cases for the app
plot_curves_app <- function(outdf) {
  df<-data.frame(dates=as.Date(outdf$dates), 
                 #Susceptible=outdf$S,
                 Infected=outdf$Infected,
                 ILI=outdf$ILI,
                 New_ILI=outdf$New_ILI,
                 Isolated=outdf$InfectedIsolated,
                 Unisolated=outdf$InfectedUnisolated,
                 Recovered=outdf$R
  )
  df_long <- reshape2::melt(df, id="dates") # convert to long format
  pl<-ggplot(data=df_long, aes(x=dates, y=value, linetype=variable, color=variable), geom=c("point", "path")) +
    geom_line() + #geom_point() + 
    #geom_point(aes(shape=variable, size=2)) +
    ylab("Individuals") +
    xlab("") + 
    guides(variable=guide_legend(title=""))+
    NULL
  if(savePlots) {
    ggsave(time_stamp("output/cases_vs_time", ".eps"), width=6, height=4.5, units="in")  
    ggsave(time_stamp("output/cases_vs_time", ".png"), width=6, height=4.5, units="in", dpi=300)  
    ggsave(time_stamp("output/cases_vs_time", ".tiff"), width=6, height=4.5, units="in", dpi=300)  
  }
  return(pl)
}

savePlots = F
#' @desc Basic plot of calibration date
plot_calibration <- function(outdf, calibdf=boarding_time_series) {
  if (is.null(calibdf)) {
    return(NULL)
  }
  if ('New_ILI' %in% names(calibdf)){
    df <- data.table(date=as.Date(outdf$dates), 
                     model=outdf$New_ILI)
    calib_df       <- data.table(date=as.Date(calibdf$date),
                                 new_ili=calibdf$New_ILI)
                                 
    ylab <- 'New daily ILI cases'
  } else if ('New_absent' %in% names(calibdf)) {
    df <- data.table(date=as.Date(outdf$date), 
                     model=outdf$New_absent)
    calib_df       <- data.table(date=as.Date(calibdf$date),
                                 new_absent=calibdf$New_absent)
    ylab <- 'New daily absent cases'
  } else {
    df <- data.table(date=as.Date(outdf$dates), 
                     model=outdf$New_ILI)
    calib_df       <- data.table(date=as.Date(outdf$dates),
                                 d=0)
    ylab <- 'error: calib data should give New_ILI or New_absent'
  }
  
  df_long <- data.table::melt(df, id="date") # convert to long format
  d_long <- data.table::melt(calib_df, id="date")
  pl<-ggplot() +
    geom_line(data=d_long, aes(x=date, y=value, linetype=variable, color=variable), geom=c("point", "path")) +
    geom_line(data=df_long, aes(x=date, y=value, linetype=variable, color=variable), geom=c("point", "path")) +
    xlab("Individuals") +
    xlab("") +
    ylab(ylab) + 
    guides(variable=guide_legend(title=""))+
    NULL
  if(savePlots) {
    #print(pl)
    ggsave(time_stamp("output/calibration_cases_vs_time", ".png"), width=6, height=4.5, units="in", dpi=300)  
    ggsave(time_stamp("output/calibration_cases_vs_time", ".tiff"), width=6, height=4.5, units="in", dpi=300)  
    ggsave(time_stamp("output/calibration_cases_vs_time", ".pdf"), width=6, height=4.5, units="in")  
  }
  return(pl)
}
  

#' @desc plot the outcome of an outbreak in detail
#' @param outdf a data.frame from an outbreak
plot_detailed <- function(outdf) {
  df<-data.frame(dates=as.Date(outdf$dates), 
                 InfectedK=outdf$InfectedK,
                 InfectedG1=outdf$InfectedG1,
                 InfectedG2=outdf$InfectedG2,
                 InfectedG3=outdf$InfectedG3,
                 InfectedG4=outdf$InfectedG4,
                 InfectedG5=outdf$InfectedG5,
                 InfectedG6=outdf$InfectedG6
                 #Susceptible=outdf$S,
                 #Infected=outdf$Infected,
                 #Isolated=outdf$InfectedIsolated,
                 #Unisolated=outdf$InfectedUnisolated,
                 #Recovered=outdf$R
  )
  
  df_long <- reshape2::melt(df, id="dates") # convert to long format
  
  pl<-ggplot(data=df_long, aes(x=dates, y=value, linetype=variable, color=variable), geom=c("point", "path")) +
    geom_line() + #geom_point() + 
    ylab("") +
    xlab("") + 
    guides(variable=guide_legend(title=""))+
    NULL
  ggsave(time_stamp("output/cases_vs_time_detailed", ".pdf"), width=6, height=4.5, units="in")  
  ggsave(time_stamp("output/cases_vs_time_detailed", ".png"), width=6, height=4.5, units="in")  
  #print(pl)
}

summarize_epidemic_metrics <- function(prediction) {
#report basic measures about an outbreak
  final_metrics <- list()
  final_metrics$total_students <- as.integer(head(prediction$S,1) + head(prediction$VS,1) + head(prediction$VR,1) + head(prediction$Infected,1) + head(prediction$R,1))

  #basic dynamic constant: rowSums(data.table(prediction)[,.(S,VS,VR,Infected,R)])
  
  final_metrics$infected    <- as.integer(tail(prediction$R,1)  + tail(prediction$Infected,1)) 
  final_metrics$vaccinated  <- as.integer(head(prediction$VS,1) + head(prediction$VR,1))
  final_metrics$persondays_isolated   <- as.integer(sum(prediction$InfectedIsolated))
  final_metrics$persondays_unisolated <- as.integer(sum(prediction$InfectedUnisolated))
  final_metrics$peak_new_ILI_cases <- as.integer(max(prediction$New_ILI))
  final_metrics$peak_new_absent <- as.integer(max(prediction$New_absent))
  final_metrics$peak_infected <- as.integer(max(prediction$Infected))
  
  final_metrics$attack_rate <- (tail(prediction$R,1) + tail(prediction$Infected,1))/final_metrics$total_students
  final_metrics$peak_date_new_ILI    <- strftime(prediction[which.max(prediction$New_ILI),"dates"], tz='GMT')
  final_metrics$peak_date_new_absent <- strftime(prediction[which.max(prediction$New_absent),"dates"], tz='GMT')
  final_metrics$peak_date_infected   <- strftime(prediction[which.max(prediction$Infected),"dates"], tz='GMT')
  
  
  #peak_duration: time from first to last day with >=2 Infected cases
  if(sum(prediction$Infected >= 2) > 0) {
    peak_outbreak_times <- which(prediction$Infected >= 2)
    final_metrics$peak_duration <- as.integer(max(peak_outbreak_times) - min(peak_outbreak_times) + 1)
  } else {
    final_metrics$peak_duration <- 0
  }
  
  return(final_metrics)
}

#based on McCann et al. for School A;   
# parameters: vaccination_rate=0
# there are 7 grades of sizes: 42, 60, 60, 49, 45, 45, 45, 45
calibration_outbreak_results_thames <- list(
  infected    = 264,
  attack_rate = 0.675,
  peak_new_absent = 38,
  peak_date_absent   = as.Date('2012-12-06'),  #absenteeism data after weekend depile
  #peak_date_absent   = as.Date('2012-12-10'),  #absenteeism data
  peak_duration   = 30 #Fig. 1: 11/26/2012 - 3 days (probably started on Friday) through 12/21/2012 + 3 days (Fig 1 shows first absence).
)

#wishlist: gently review this - is the peak corresponding to ILI or actual infections?
calibration_outbreak_results_pennsylvania <- list(
  infected           = 109.4,
  attack_rate        = 109.4/456,
  peak_new_ILI_cases = 14.12,
  #peak_date_new_ILI          = as.Date('2009-05-11'), #this is the real peak date
  peak_date_new_ILI          = as.Date('2009-05-10'),
  peak_duration      = 12 #Fig. 2: 5/5/2009 through 5/16/2009.  this corresponds to peak ILI, not necessarily peak infected
)

calibration_outbreak_results_boarding <- list(
  infected           = 101,
  attack_rate        = 101/1307,
  peak_new_ILI_cases = 10,
  peak_date_new_ILI          = as.Date('2009-05-25'),
  peak_duration      = 19#Fig. 1: 5/9/2009 through 5/31/2009.  this corresponds to peak ILI, not necessarily peak infected
)

get_calibration_dataset_results <- function(dataset_name) {
  if(dataset_name == "Smith pH1H1 2009") {
    return(calibration_outbreak_results_boarding) 
  } else if (dataset_name == "Marchbanks pH1N1 2009") {
    return(calibration_outbreak_results_pennsylvania) 
  } else if (dataset_name == "McCann Influenza B 2012") {
    return(calibration_outbreak_results_thames) 
  } else {
    return(NULL)
  } 
}

calibration_difference <- function(new_params    = list(), 
                                   params        = params_def_elementary, 
                                   calibration_outbreak,
                                   combine_weekend_to_monday) {
  #@desc: runs an outbreak and reports the difference to a calibration outbreak 
  #       used in calibration
  #@returns: a float
  #@param new_params: new candidate parameters constructed during optimization
  #@param calibration_outbreak: metrics of the outbreak which we use as the calibration
  for (p in names(new_params)) {
    stopifnot(p %in% names(params_def_elementary))
    params[p] <- new_params[p]
  }
  
  outbreak_result   <- compute_school_outbreak(params = params)
  if(combine_weekend_to_monday) {
    outbreak_result$dates[is.weekend(outbreak_result$dates)] <- (outbreak_result$dates[is.weekend(outbreak_result$dates)] 
                                                                 + (wday(outbreak_result$dates[is.weekend(outbreak_result$dates)])) %% 5)
    outbreak_result <- as.data.table(outbreak_result)[, lapply(.SD, sum), by = .(dates = dates)]
    outbreak_result <- data.frame(outbreak_result, row.names = outbreak_result$dates)
    
  }
  ##Strip weekends here, pass a boolean to determine if a strip occurs.
  final_metrics     <- summarize_epidemic_metrics(outbreak_result)
  result_difference <- calibration_result_difference(calibration_outbreak, final_metrics)
  
  final_list = c(final_metrics$infected,
                    final_metrics$peak_new_ILI_cases,
                    #final_metrics$peak_new_absent,
                    final_metrics$attack_rate,
                    as.Date(final_metrics$peak_date_new_ILI),
                    #as.Date(final_metrics$peak_date_absent),
                    final_metrics$peak_duration
                    )
  calib_list = c(calibration_outbreak$infected,
                    calibration_outbreak$peak_new_ILI_cases,
                    #calibration_outbreak$peak_new_absent,
                    calibration_outbreak$attack_rate,
                    as.Date(calibration_outbreak$peak_date_new_ILI),
                    #as.Date(calibration_outbreak$peak_date_absent),
                    calibration_outbreak$peak_duration
                    )
  

  
  #print('Final')
  #print(final_list)
  #print('Calib')
  #print(calib_list)
  return(result_difference)  
}

calibration_result_difference <- function(calibration_outbreak, final_metrics) {
  #wishlist: allow weight parameters instead of hardcoding
  ret <- sqrt(
         0.5*(final_metrics$infected                      - calibration_outbreak$infected)**2 
    +    3.0*(final_metrics$peak_new_ILI_cases            - calibration_outbreak$peak_new_ILI_cases)**2 
    #+    3.0*(final_metrics$peak_new_absent            - calibration_outbreak$peak_new_absent)**2 
    +    5.0*(as.integer(as.Date(final_metrics$peak_date_new_ILI) - as.Date(calibration_outbreak$peak_date_new_ILI)))**2 
    #   +    5.0*(as.integer(as.Date(final_metrics$peak_date_absent) - as.Date(calibration_outbreak$peak_date_absent)))**2 
    +    0.5*(final_metrics$peak_duration                 - calibration_outbreak$peak_duration)**2)
  return(ret)
}


#' @desc  generates a unique time stamp for every data generated
time_stamp <- function(prefix="", suffix="", outmsg=TRUE) {
  t <- format(Sys.time(), "%Y-%m-%d__%H-%M-%S");
  s <- as.integer(runif(1, max=1000))
  filename <- paste(prefix, t, s, suffix, sep="")
  if (outmsg) {
    print(filename)
  }
  return(filename)
}

update_off_diagonal_values <- function(contact_data, off_diagonal_vals) {
  contact_data[contact_data!=1] <- off_diagonal_vals
  return(contact_data)
}


# ' @desc 
# is_school_break <- function(current_date)

### RUNNING THE CODE ####
### this is normally commented-out when running the app
# Basic usage
# compute_school_outbreak(params = params_def_elementary)

#V1
# new_params = params_def_elementary
# print(new_params$num_stages)
# compute_school_outbreak(params = new_params)
# outbreak_result<-compute_school_outbreak(params = new_params)
#initial in each grade = susceptible + infected + vaccinated; the eqn is true for all times during the outbreak

#Table T1
#covid_rates_table=effective_return_rates_table(coronavirus_symptoms$SystemicSymptoms)
#flu_rates_table=effective_return_rates_table(influenza_symptoms$SystemicSymptoms)

# figure F1
# outdf <- compute_school_outbreak_alt(params = params_def_elementary, doplot=FALSE)
# plot_base(outdf)
# plot_detailed(outdf)

#Table T2
# report_sensitivity_grant(def_params = params_def_elementary)
