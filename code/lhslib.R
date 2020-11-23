#'
#' Library for latin hypercube sampling
#'
#' (c) 2016-2020 By Authors.
#' You are free to use this software under CC-BY 4.0 https://creativecommons.org/licenses/by/4.0/
#'
#' 

#' Usage:
#' 1. lhs_make_samples()
#' 2. lhs_run_samples()
#' 3. User specified lhs_summarize() 
#'

library(lhs)

lhs_make_samples <- function(num_samples=1000, params_override=list(), param_descriptions, baseline_sample) {
  #generates samples. access with ret[[idx]], ret[[idx]]$data, ret[[idx]]$tp, ret[[idx]]$data
  #use "info" for sampling data
  lhs_samples <- list()
  lhs_samples["num_samples"] <- num_samples
  #https://stat.ethz.ch/pipermail/r-help/2011-June/279931.html
  #method: 1. create uniform distributions, 2.  draw from quantiles
  var_descriptors <- param_descriptions
  num_vars <- length(var_descriptors)
  var_names <- names(var_descriptors)
  lhs_samples$num_vars <- num_vars
  lhs_samples$var_names <- var_names
  
  uniform_samples <- maximinLHS(n=num_samples,k=num_vars)
  #reduce variance to just the median
  for(i in seq(num_vars)) {
    descr <- var_descriptors[[i]]
    unif_quantiles <- uniform_samples[,i]
    var_name <- var_names[i]
    #either a constant, a constant with data, or a function (distribution)
    if((! ("func" %in% names(descr))) ) {
      varq <- list(name=var_names[i], tp="const", data=descr)
    } else if(is.character(descr$func) && ( descr$func== "const")) {
      varq <- list(name=var_names[i], tp="const", data=descr$args)     
    } else {
      descr$args$p <- unif_quantiles
      varq <- list(name=var_names[i], tp="varying", data=do.call(descr$func, descr$args))
      stopifnot(length(varq$data) == length(descr$args$p))
    }
    lhs_samples[[var_name]] = varq
  }
  #lhs_samples are arrays, one per variable.  these are now converted to transverse input sets
  clean_samples <- lhs_make_samples_distinct(lhs_samples, params_override, baseline_sample=baseline_sample)
  return(clean_samples)  
}

lhs_make_samples_distinct <- function(samples, params_override, baseline_sample=NULL) {
  #generates actual samples, copying any constants
  num_vars <- samples[["num_vars"]]
  var_names <- samples[["var_names"]]
  #creates a list of finished samples
  clean_samples <- list()
  for(sample_num in seq(samples[["num_samples"]])) {
    if(! is.null(baseline_sample)) {
      sample_params = baseline_sample
    } else {
      sample_params <- list()
    }
    #sample_params <- default_vals
    for(var_name in var_names) {
      var_data <- samples[[var_name]]
      if(is.null(var_data)) {
        next
      }
      if(var_data$tp == "const") {
        var_val <- var_data$data
      } else {
        #var_val <- var_data$data[sample_num]
        var_val <- var_data$data[[sample_num]]
      }
      sample_params[[var_name]] = var_val
    }
    #print(sprintf("Sample %d:", sample_num))
    for (p in names(params_override)) {
      sample_params[p] <- params_override[p]
    }
    clean_samples[[sample_num]] <- sample_params
  }  
  if(! is.null(baseline_sample)) {
    #the baseline sample is inserted as the first sample
    clean_samples <- append(list(NULL), clean_samples)
    clean_samples[[1]] <- baseline_sample
  }
  return(clean_samples)
}

lhs_run_samples <- function(samples, name="set1", func=z_transmission_for_titers, rsavefile=NULL, funcparams=NULL) {
  #run each sample from the samples
  if(is.null(rsavefile)) {
    res1       <- c()
    fullres   <- list()
    #fullres[[name1]] = list()
    sample_num = 1
    print(name)
    for(sample in samples) {
      cat(sample_num)
      args <- copy(funcparams)
      args[["sample"]] <- sample
      #res_sample  <- do.call(func, args=args)
      res_sample <- func(args)
      #fullres[[name1]][[sample_num]] <- res_sample
      fullres[[sample_num]] <- res_sample
      sample_num = sample_num+1
    }
    results=list(samples=samples, name=name, fullres=as.list(fullres))
    save(results, file=time_stamp(paste("output/lhs_runs_", name, "_num_samples=", length(samples), "__", sep=""), ".RSave"))
  } else {
    results <- load(rsavefile)
    fullres <- results$fullres
  }
  
  return(fullres)
}

#' demo function for summarizing results
#' WRITE YOUR OWN
#' assumes func returns a DT with x.name and y.name
#' @param func or paramset need to be specified
#' 
lhs_summarize_demo <- function(func, fname=NULL, paramset=NULL, x.name="logTiters", y.name="p") {
  #(optionally runs and LHS) and return a summary of the outcomes from an LHS
  if(! is.null(paramset)) {
    results_runs <- lhs_run_samples(samples=paramset, name="validation_with_data", func=func)
  } else {
    load(file=fname, newenv<-new.env())
    results_pack<-get("results", newenv)
    results_runs <- results_pack[["fullres"]]
  }
  all_runs_data  <- NULL
  for(sample_num in seq(length(results_runs))) {
    run_data <- results_runs[[sample_num]]
    all_runs_data <- rbind(all_runs_data, run_data[[y.name]])
  }
  ret <- data.frame(num=1:length(run_data[[x.name]]))
  ret[,x.name] <- run_data[[x.name]]
  ret$y_mean  <- apply(all_runs_data, 2, mean)
  #other analysis ...
  
  return(ret)
}


#other helpers
repl <- function(obj, times){lapply(seq(times), FUN=function(i){obj})}

time_stamp <- function(prefix="", suffix="", outmsg=TRUE) {
  #generates a unique time stamp for every data generated
  t <- format(Sys.time(), "%Y-%m-%d__%H-%M-%S");
  s <- as.integer(runif(1, max=1000))
  filename <- paste(prefix, t, s, suffix, sep="")
  if (outmsg) {
    print(filename)
  }
  return(filename)
}
