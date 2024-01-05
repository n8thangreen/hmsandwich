
#
add_targets <- function(x, file = "Inputs/Calibration targets.txt") {
  
  # incidence group by year
  targets_dat <-
    read.delim(file, sep = "\t", header = FALSE) |> 
    select_if(~ !any(is.na(.)))
  
  # convert to vector
  # mean and standard deviation
  target_val <- data.frame(
    val = inc_mat_to_vector(targets_dat)) |> 
    mutate(sigma = val/(10*3.92))
  
  # convert to named list
  all_targets <-
    split(target_val, 1:nrow(target_val)) |> 
    setNames(x$full_groups_out) |> 
    map(as.list)
  
  # subset parameters and data for testing
  c(x, targets = all_targets[x$indx_out])
}
