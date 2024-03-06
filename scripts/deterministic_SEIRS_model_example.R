# deterministic SEIRS model example
# https://danny-sc.github.io/determ_workshop/index.html

example_params <- c(
  b = 1/(60*365),
  mu = 1/(76*365),
  beta1 = 0.2, beta2 = 0.1, beta3 = 0.3,
  epsilon = 0.13,
  alpha = 0.01,
  gamma = 0.08,
  omega = 0.003
)
solution <- ode_results(example_params)
par(mar = c(2, 2, 2, 2))
plot(solution)

ranges = list(
  b = c(1e-5, 1e-4), # birth rate
  mu = c(1e-5, 1e-4), # rate of death from other causes
  beta1 = c(0.2, 0.3), # infection rate at time t=0
  beta2 = c(0.1, 0.2), # infection rates at time t=100
  beta3 = c(0.3, 0.5), # infection rates at time t=270
  epsilon = c(0.07, 0.21), # rate of becoming infectious after infection
  alpha = c(0.01, 0.025), # rate of death from the disease
  gamma = c(0.05, 0.08), # recovery rate
  omega = c(0.002, 0.004) # rate at which immunity is lost following recovery
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

initial_LHS_training <- maximinLHS(90, 9)
initial_LHS_validation <- maximinLHS(90, 9)
initial_LHS <- rbind(initial_LHS_training, initial_LHS_validation)

initial_points <- setNames(data.frame(t(apply(initial_LHS, 1, 
                  function(x) x*unlist(lapply(ranges, function(x) x[2]-x[1])) + 
                  unlist(lapply(ranges, function(x) x[1]))))), names(ranges))

initial_results <- setNames(data.frame(t(apply(initial_points, 1, get_results, 
                   c(25, 40, 100, 200, 300, 350), c('I', 'R')))), names(targets))

wave0 <- cbind(initial_points, initial_results)

training <- wave0[1:90,]
validation <- wave0[91:180,]

ems_wave1 <- emulator_from_data(training, names(targets), ranges, 
                                specified_priors = list(hyper_p = rep(0.55, length(targets))))

emulator_plot(ems_wave1$R200, params = c('beta1', 'gamma'))
plot_actives(ems_wave1)
emulator_plot(ems_wave1$R200, plot_type = 'var', params = c('beta1', 'gamma'))
summary(ems_wave1$R200$model)$adj.r.squared
emulator_plot(ems_wave1$R200, plot_type = 'imp', 
              targets = targets, params = c('beta1', 'gamma'), cb=TRUE)
emulator_plot(ems_wave1, plot_type = 'imp', 
              targets = targets, params = c('beta1', 'gamma'), cb=TRUE)

vd <- validation_diagnostics(ems_wave1$R200, validation = validation, targets = targets, plt=TRUE)
sigmadoubled_emulator <- ems_wave1$R200$mult_sigma(2)
vd <- validation_diagnostics(sigmadoubled_emulator, 
                             validation = validation, targets = targets, plt=TRUE)
for (j in 1:length(ems)) {
      misclass <- nrow(classification_diag(ems[[j]], targets, validation, plt = FALSE))
      while(misclass > 0) {
        ems[[j]] <- ems[[j]]$mult_sigma(1.1)
        misclass <- nrow(classification_diag(ems[[j]], targets, validation, plt = FALSE))
      }
}
bad.ems <- c()
for (j in 1:length(ems)) {
          bad.model <- nrow(comparison_diag(ems[[j]], targets, validation, plt = FALSE))
          if (bad.model > floor(nrow(validation)/10)) {
            bad.ems <- c(bad.ems, j)
   }
}
ems <- ems[!seq_along(ems) %in% bad.ems]

# proposing new points
restricted_ems <- ems_wave1[c(1,2,3,4,7,8,9,10)]
new_points_restricted <- generate_new_design(restricted_ems, 180, targets, verbose=TRUE)
plot_wrap(new_points_restricted, ranges)
space_removed(ems_wave1, targets, ppd=3) + geom_vline(xintercept = 3, lty = 2) + 
  geom_text(aes(x=3, label="x = 3",y=0.33), colour="black", 
  angle=90, vjust = 1.2, text=element_text(size=11))

# second wave
R_squared_new <- list()
for (i in 1:length(ems_wave2)) {
  R_squared_new[[i]] <- summary(ems_wave2[[i]]$model)$adj.r.squared
}
names(R_squared_new) <- names(ems_wave2)
unlist(R_squared_new)

 ems_wave1_linear <- emulator_from_data(training, names(targets), 
                                        ranges, quadratic=FALSE,
                                        specified_priors = list(hyper_p = rep(0.55, length(targets))))
 R_squared_linear <- list()
 for (i in 1:length(ems_wave1_linear)) {
   R_squared_linear[[i]] <- summary(ems_wave1_linear[[i]]$model)$adj.r.squared
 }
 names(R_squared_linear) <- names(ems_wave1_linear)
 unlist(R_squared_linear)

 emulator_plot(ems_wave1_linear$I200, plot_type = 'var', 
               params = c('beta1', 'gamma'))

emulator_plot(ems_wave1_linear$I200, plot_type = 'imp', targets = targets, 
              params = c('beta1', 'gamma'), cb=TRUE)

ems_wave1_linear$I200 <- ems_wave1_linear$I20$set_hyperparams(
              list(theta=ems_wave1_linear$I200$corr$hyper_p$theta *3))
emulator_plot(ems_wave1_linear$I200, plot_type = 'var', 
              params = c('beta1', 'gamma'))

emulator_plot(ems_wave1_linear$I200, plot_type = 'imp', targets = targets, 
              params = c('beta1', 'gamma'), cb=TRUE)
wave_points(list(initial_points, new_points, new_new_points), input_names = names(ranges), p_size=1)

new_new_initial_results <- setNames(data.frame(t(apply(new_new_points, 1, 
                                    get_results, c(25, 40, 100, 200, 300, 350), 
                                    c('I', 'R')))), names(targets))
wave2 <- cbind(new_new_points, new_new_initial_results)

all_points <- list(wave0, wave1, wave2)
simulator_plot(all_points, targets)
wave_values(all_points, targets, l_wid=1, p_size=1)

