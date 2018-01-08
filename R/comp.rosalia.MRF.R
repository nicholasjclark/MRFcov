#'Compare MRFcov and rosalia outputs
#'
#'This function compares coefficients estimated from binary
#'\code{\link[rosalia]{rosalia}} MRF models (fit by Maximum Likelihood) and
#'\code{\link{MRFcov}} models (fit by separate penalized regressions). Here,
#'\code{\link[rosalia]{rosalia}} coefficients are considered the \emph{true} values,
#' and \code{\link{MRFcov}} coefficients are then compared to assess positive and negative
#' predictive capacity of the separate penalized regression approach.
#'
#'@param MRF A fitted \code{\link{MRFcov}} model
#'@param rosalia A fitted \code{\link[rosalia]{rosalia}} model
#'@return A \code{list} of two objects:
#'\itemize{
#'    \item \code{Positive_predictive_value}: The positive predictive value
#'    \item \code{Negative_predictive_value}: The negative predictive value
#'    }
#'
#'@seealso \code{\link{MRFcov}}, \code{\link[rosalia]{rosalia}}
#'@export
comp.rosalia.MRF = function(MRF,rosalia){

true.pos <- missed.pos <- true.neg <- missed.neg <- matrix(NA, ncol = ncol(MRF),
                                               nrow = nrow(MRF))
for(i in 1:nrow(true.pos)){
  for(j in 1:ncol(true.pos)){
    true.pos[i, j] <- isTRUE(all.equal(sign(MRF[i, j]),
                                      sign(rosalia[i, j]))) &
      isTRUE(sign(rosalia[i,j]) == 1)

    missed.pos[i, j] <- !isTRUE(all.equal(sign(MRF[i, j]),
                                       sign(rosalia[i, j]))) &
      isTRUE(sign(rosalia[i, j]) == 1)

    true.neg[i, j] <- isTRUE(all.equal(sign(MRF[i, j]),
                                      sign(rosalia[i, j]))) &
      isTRUE(sign(rosalia[i, j])==-1)

    missed.neg[i, j] <- !isTRUE(all.equal(sign(MRF[i, j]),
                                       sign(rosalia[i, j]))) &
      isTRUE(sign(rosalia[i, j]) == -1)
  }
}

pos.pred <- sum(true.pos, na.rm = TRUE)/
                (sum(true.pos, na.rm = TRUE) + sum(missed.pos, na.rm = TRUE))

neg.pred <- sum(true.neg, na.rm = TRUE)/
                (sum(true.neg, na.rm = TRUE) + sum(missed.neg, na.rm = TRUE))
return(list(Positive_predictive_value = pos.pred,
            Negative_predictive_value = neg.pred))
}
