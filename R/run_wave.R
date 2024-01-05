
#
run_wave <- function(x,
                     sim = 100,
                     n_validation = 10) {
  
  wave_name <- paste0("wave", x$wave_no)
  x[[wave_name]] <- list()
  
  # initial wave
  if (x$wave == 0) {
    lhs_points <- lhs::maximinLHS(x$n_sim, x$n_grps_in)
    lhs_points_validation <- lhs::maximinLHS(x$n_validation, x$n_grps_in)
    
    # all LHS inputs
    init_points <-
      rbind(lhs_points, lhs_points_validation) |> 
      `colnames<-`(x$groups_in)
    
    x$wave0$inputs <- rescale(x$ranges_in, init_points)
    
    # run model
    x$wave0$results <- t(apply(x$wave0$inputs, 1,
                               x$model,
                               indx_in = x$indx_in,
                               indx_out = x$indx_out))
    
    x$wave0$data <- inp_out_df(x$wave0$inputs, x$wave0$results)
    
    # split data set
    x$wave0$training <- x$wave0$data[1:x$n_sim, ]
    x$wave0$validation <- x$wave0$data[(x$n_sim+1):nrow(x$wave0$data), ]
    
    return(x) 
  } else if (x$wave == 1){
    
    x$wave1$ems <-
      emulator_from_data(
        input_data = x$wave0$training,     # named inputs and outputs from full model
        output_names = names(x$targets),
        range = x$ranges_in,               # min, max inputs
        emulator_type = "deterministic",
        order = 2)
    
    return(x)
  } else {
    pre_wave <- paste0("wave", x$wave_no - 1)
    
    x[[pre_wave]]$inputs <-
      generate_new_design(x[[wave_name]]$ems,
                          n_points = 60,
                          x$targets,
                          verbose = TRUE)
    
    # run model
    x[[wave_name]]$results <- t(apply(x[[pre_wave]]$inputs, 1,
                                      x$model,
                                      indx_in = x$indx_in,
                                      indx_out = x$indx_out))
        
    x[[wave_name]]$data <- inp_out_df(x[[pre_wave]]$inputs,
                                      x[[wave_name]]$results)
    
    # split data set
    x[[wave_name]]$training <- x[[wave_name]]$data[1:x$n_sim, ]
    x[[wave_name]]$validation <- x[[wave_name]]$data[(x$n_sim+1):nrow(x[[wave_name]]$data), ]
    
    x[[wave_name]]$ems <-
      emulator_from_data(
        input_data = xx[[wave_name]]$training,
        output_names = names(x$targets),
        range = x$ranges_in,
        emulator_type = "deterministic",
        order = 2)
        
    return(x)
  }
}
