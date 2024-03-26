
#' Run history matching and emulation
#' 
#' Original code taken from SIR example in hmer documentation.
#' Depends on which wave.
#' 
#' @importFrom lhs maximinLHS
#' 
run_wave <- function(x,
                     n_sim = 100,
                     n_validation = 10) {
  
  wave_name <- paste0("wave", x$wave_no)
  x[[wave_name]] <- list()
  
  # initial wave
  if (x$wave == 0) {
    # latin hypercube sampling
    initial_LHS_training <- maximinLHS(90, 9)
    initial_LHS_validation <- maximinLHS(90, 9)
    initial_LHS <- rbind(initial_LHS_training, initial_LHS_validation)
    
    x$wave0$inputs <-
      setNames(data.frame(t(apply(initial_LHS, 1, 
                                  function(x) x*unlist(lapply(x$ranges, function(y) y[2] - y[1])) + 
                                    unlist(lapply(x$ranges, function(y) y[1]))))), names(x$ranges))
    
    # run model to obtain output sample
    ##TODO: hard coded input and output parameters
    x$wave0$results <-
      setNames(data.frame(t(apply(x$wave0$inputs, 1, x$model, 
                                  c(25, 40, 100, 200, 300, 350),
                                  c('I', 'R')))),
               names(x$targets))
    
    x$wave0$data <- cbind(x$wave0$inputs, x$wave0$results)
    
    # split data set
    x$wave0$training <- x$wave0$data[1:x$n_sim, ]
    x$wave0$validation <- x$wave0$data[(x$n_sim+1):nrow(x$wave0$data), ]
    
    return(x) 
  } else if (x$wave == 1){
    
    # fit emulator model
    x$wave1$ems <-
      emulator_from_data(x$wave0$training,
                         names(x$targets),
                         x$ranges, 
                         specified_priors = list(hyper_p = rep(0.55, length(x$targets))))
    return(x)
  } else {
    pre_wave <- paste0("wave", x$wave_no - 1)
    
    # sample new set of input parameter values
    x[[pre_wave]]$inputs <-
      generate_new_design(x[[wave_name]]$ems, n_points = 180, x$targets, verbose=TRUE)
    
    # run model
    x[[wave_name]]$results <-
      setNames(data.frame(t(apply(x[[pre_wave]]$inputs, 1, 
                                  x$model,
                                  c(25, 40, 100, 200, 300, 350), 
                                  c('I', 'R')))),
               names(x$targets))
    
    x[[wave_name]]$data <- cbind(x[[pre_wave]]$inputs,
                                 x[[wave_name]]$results)
    
    # split data set
    x[[wave_name]]$training <- x[[wave_name]]$data[1:x$n_sim, ]
    x[[wave_name]]$validation <- x[[wave_name]]$data[(x$n_sim+1):nrow(x[[wave_name]]$data), ]
    
    x[[wave_name]]$ems <-
      emulator_from_data(
        input_data = x[[wave_name]]$training,
        output_names = names(x$targets),
        range = x$ranges,
        emulator_type = "deterministic",
        order = 2)
        
    return(x)
  }
}
