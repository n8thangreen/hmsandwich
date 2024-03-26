
#' Run history matching and emulation
#' 
#' Original code taken from SIR example in hmer documentation.
#' Depends on which wave.
#' 
#' @importFrom lhs maximinLHS
#' @import hmer
#' 
run_wave <- function(x) {
  
  wave_name <- paste0("wave", x$wave_no)
  x[[wave_name]] <- list()
  
  # initial wave
  if (x$wave_no == 0) {
    # latin hypercube sampling
    initial_LHS_training <- lhs::maximinLHS(90, 9)
    initial_LHS_validation <- lhs::maximinLHS(90, 9)
    initial_LHS <- rbind(initial_LHS_training, initial_LHS_validation)
    
    x$wave0$inputs <-
      setNames(data.frame(t(apply(initial_LHS, 1, 
                                  function(y) y*unlist(lapply(x$ranges, function(y) y[2] - y[1])) + 
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
    x$wave0$training <- x$wave0$data[1:x$n_train, ]
    x$wave0$validation <- x$wave0$data[(x$n_train + 1):nrow(x$wave0$data), ]
    
    x$wave_no <- 1
    
    return(x) 
  } else if (x$wave_no == 1){
    
    # fit emulator model
    x$wave1$ems <-
      hmer::emulator_from_data(x$wave0$training,
                               names(x$targets),
                               x$ranges, 
                               specified_priors = list(hyper_p = rep(0.55, length(x$targets))))
    x$wave_no <- 2
    
    return(x)
  } else {
    pre_wave <- paste0("wave", x$wave_no - 1)
    
    # sample new set of input parameter values
    x[[pre_wave]]$inputs <-
      hmer::generate_new_design(x[[pre_wave]]$ems, n_points = 180, x$targets, verbose=TRUE)
    
    # run model
    x[[pre_wave]]$results <-
      setNames(data.frame(t(apply(x[[pre_wave]]$inputs, 1, 
                                  x$model,
                                  c(25, 40, 100, 200, 300, 350), 
                                  c('I', 'R')))),
               names(x$targets))
    
    x[[pre_wave]]$data <- cbind(x[[pre_wave]]$inputs,
                                 x[[pre_wave]]$results)
    
    # split data set
    x[[pre_wave]]$training <- x[[pre_wave]]$data[1:x$n_train, ]
    x[[pre_wave]]$validation <- x[[pre_wave]]$data[(x$n_train + 1):nrow(x[[pre_wave]]$data), ]
    
    x[[wave_name]]$ems <-
      hmer::emulator_from_data(
        input_data = x[[pre_wave]]$training,
        output_names = names(x$targets),
        range = x$ranges,
        emulator_type = "deterministic",
        order = 2)
    
    x$wave_no <- x$wave_no + 1
        
    return(x)
  }
}
