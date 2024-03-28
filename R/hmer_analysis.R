
#' create an hmer object
#' 
#' @importFrom tibble lst
#' @export
#'
hmer_analysis <- function(n_train = NULL,
                          n_grps_in = NULL,
                          groups_in = NULL,
                          indx_out = NULL,
                          full_groups_out = NULL) {
  lst(
    wave_no = 0,
    n_train = 100,
    n_grps_in,
    groups_in,
    indx_out,
    full_groups_out
  )
}
