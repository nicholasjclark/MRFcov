#' Blood parasite occurrences in New Caledonian birds.
#'
#' A dataset containing binary occurrences of four blood parasite species
#' in New Caledonian birds. The first four variables represent the parasite occurrences
#' and the last variable is a scaled continuous covariate representing host relative
#' abundance.
#'
#' @format A data frame with 449 rows and 5 variables:
#' \describe{
#'   \item{Hzosteropis}{binary occurrence of \emph{Haemoproteus zosteropis}}
#'   \item{Hkillangoi}{binary occurrence of \emph{Haemoproteus killangoi}}
#'   \item{Plas}{binary occurrence of \emph{Plasmdodium} species}
#'   \item{Microfilaria}{binary occurrence of Microfilaria species}
#'   \item{scale.prop.zos}{scaled numeric variable representing relative abundance
#'   of \emph{Zosterops} host species}
#' }
#' @references Clark, N.J., Wells, K., Dimitrov, D. & Clegg, S.M. (2016)
#' Co-infections and environmental conditions drive the distributions of blood
#' parasites in wild birds. Journal of Animal Ecology, 85, 1461-1470.
#' @source \url{http://dx.doi.org/10.5061/dryad.pp6k4}
"Bird.parasites"
