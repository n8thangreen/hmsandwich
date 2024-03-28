
#' Add input ranges
#' @importFrom tibble lst
#' @export
#' 
add_input_ranges <- function(x, ranges = NULL) {
  # ranges <- 
  #   rep(list(c(0, 0.05)), x$n_grps_in) |> 
  #   setNames(x$groups_in)
  
  c(x, ranges = list(ranges))
}
