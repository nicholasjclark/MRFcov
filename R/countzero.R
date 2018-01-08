#'Count proportions of non-zero coefficients
#'
#'This function counts the proportion of bootstrap estimates
#'in which each tested model coefficient is nonzero. Note that
#'this function is not called directly, but is used in the function
#'\code{\link{bootstrap_MRF}}
#'
#'@importFrom magrittr %>%
#'
#'@param data Matrix of coefficients returned from a fitted \code{MRFcov} model
#'@param x Integer for indexing rows in \code{data}
#'@param y Integer for indexing columns in \code{data}
#'@return Matrix containing proportions of replicates in which each coefficient was retained
#'
#'@author Nicholas Clark: \url{nicholas.j.clark1214@gmail.com}
#'@seealso \code{\link{bootstrap_MRF}}
#'
#'@export
#'
countzero <- function(data, x, y){
  bs.unlist <- data %>% purrr::map('direct_coefs')
  estimatesinxy <- unlist(lapply(bs.unlist, '[', x, y))
  zeron <- length(which(estimatesinxy == 0))

  #correct for finite sampling
  ((.0001 * length(data)) + zeron) / ((.0001 * length(data)) + length(data))
}
