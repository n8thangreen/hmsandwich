# SIR model history matching
# with emulation workflow
#
# see tidy models for inspiration
# https://www.tmwr.org/workflows.html


##TODO: read from file?

# input data
ranges <- list(
  b = c(1e-5, 1e-4),       # birth rate
  mu = c(1e-5, 1e-4),      # rate of death from other causes
  beta1 = c(0.2, 0.3),     # infection rate at time t=0
  beta2 = c(0.1, 0.2),     # infection rates at time t=100
  beta3 = c(0.3, 0.5),     # infection rates at time t=270
  epsilon = c(0.07, 0.21), # rate of becoming infectious after infection
  alpha = c(0.01, 0.025),  # rate of death from the disease
  gamma = c(0.05, 0.08),   # recovery rate
  omega = c(0.002, 0.004)  # rate at which immunity is lost following recovery
)

targets <- list(
  I25 = list(val = 115.88, sigma = 5.79),
  I40 = list(val = 137.84, sigma = 6.89),
  I100 = list(val = 26.34, sigma = 1.317),
  I200 = list(val = 0.68, sigma = 0.034),
  I300 = list(val = 29.55, sigma = 1.48),
  I350 = list(val = 68.89, sigma = 3.44),
  R25 = list(val = 125.12, sigma = 6.26),
  R40 = list(val = 256.80, sigma = 12.84),
  R100 = list(val = 538.99, sigma = 26.95),
  R200 = list(val = 444.23, sigma = 22.21),
  R300 = list(val = 371.08, sigma = 15.85),
  R350 = list(val = 549.42, sigma = 27.47)
)

# SIR model functions

# provides the solution of the differential equations for a given 
# set of parameters. assumes an initial population of 
# 900 susceptible individuals, 100 exposed individuals, and no infectious 
# or recovered individuals.
ode_results <- function(parms, end_time = 365*2) {
  forcer = matrix(c(0,   parms['beta1'],
                    100, parms['beta2'],
                    180, parms['beta2'],
                    270, parms['beta3']),
                  ncol = 2,
                  byrow = TRUE)
  force_func <- approxfun(
    x = forcer[, 1],
    y = forcer[, 2],
    method = "linear",
    rule = 2)
  
  des <- function(time, state, parms) {
    with(as.list(c(state, parms)), {
      # N <- S + E + I + R
      dS <- b*(S+E+I+R) - force_func(time)*I*S / (S+E+I+R) + omega*R - mu*S
      dE <- force_func(time)*I*S / (S+E+I+R) - epsilon*E - mu*E
      dI <- epsilon*E - alpha*I - gamma*I - mu*I
      dR <- gamma*I - omega*R - mu*R
      return(list(c(dS, dE, dI, dR)))
    })
  }
  yini <- c(S = 900, E = 100, I = 0, R = 0)
  times <- seq(0, end_time, by = 1)
  
  deSolve::ode(yini, times, des, parms)
}

#' Wrapper for `ode_results` to subset which outputs and times should be returned
#'
#' For example, to obtain the number of infected and susceptible individuals 
#' at t=25 and t=50
#' times = c(25,50)
#' outputs = c('I','S')
#' @returns
#'
sir_model <- function(params, times, outputs) {
  t_max <- max(times)
  all_res <- ode_results(params, t_max)
  actual_res <- all_res[all_res[,'time'] %in% times, c('time', outputs)]
  shaped <- reshape2::melt(actual_res[, outputs])
  
  names_shaped <- paste0(shaped$Var2, actual_res[,'time'], sep = "")

  return(setNames(shaped$value, names_shaped))
}

###########
# analysis

# create analysis object
sir_hmer <- 
  hmer_analysis() |> 
  add_input_ranges(ranges) |> 
  add_targets(targets) |> 
  add_model(sir_model)

# run initial wave
sir_hmer <- run_wave(sir_hmer)
# sir_hmer@run_wave


# run first wave
sir_hmer <- run_wave(sir_hmer)
# sir_hmer@run_wave



