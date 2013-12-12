\name{plfm-package}
\alias{plfm-package}
\docType{package}
\title{Probabilistic latent feature analysis for two-way two-mode frequencies}
\description{
The package can be used for disjunctive and conjunctive probabilistic latent feature analysis on two-way two-mode frequencies
}
\details{
\tabular{ll}{
Package: \tab plfm\cr
Type: \tab Package\cr
Version: \tab 1.1.1\cr
Date: \tab 2013-12-12\cr
License: \tab GPL(>=2)\cr
LazyLoad: \tab yes\cr
}

The \code{\link{plfm}} function can be used to locate the posterior mode(s) of probabilistic feature models, 
 and to compute information criteria for model selection, and measures of statistical and descriptive model fit.
 The \code{\link{stepplfm}} function can be used to fit a series of disjunctive and/or conjunctive probabilistic latent feature models with different number of latent features. 
 In addition, the \code{\link{bayesplfm}} function can be used to compute a sample of the posterior distribution 
 in the neigbourhood of a specific posterior mode.

To see the preferable citation of the package, type citation("plfm").

}
\author{
Michel Meulders

Maintainer: <michel.meulders@kuleuven.be>
}
\references{


Candel, M. J. J. M., and Maris, E. (1997). Perceptual analysis of two-way two-mode
frequency data: probability matrix decomposition and two alternatives.
\emph{International Journal of Research in Marketing, 14}, 321-339.

Gelman, A., Van Mechelen, I., Verbeke, G., Heitjan, D. F., and Meulders, M. (2005). 
Multiple imputation for model checking: Completed-data plots with missing and latent data. 
\emph{Biometrics, 61}, 74-85. 

Maris, E., De Boeck, P., and Van Mechelen, I. (1996). Probability matrix decomposition models. \emph{Psychometrika, 61}, 7-29.

Meulders, M. (2013). An R Package for Probabilistic Latent Feature Analysis of Two-Way Two-Mode Frequencies. \emph{Journal of Statistical Software, 54(14)}, 1-29. 
   URL http://www.jstatsoft.org/v54/i14/.

Meulders, M., De Boeck, P., Kuppens, P., and Van Mechelen, I. (2002). Constrained latent
class analysis of three-way three-mode data. \emph{Journal of Classification, 19}, 277-302.

Meulders, M., De Boeck, P., and Van Mechelen, I. (2001). Probability matrix decomposition
models and main-effects generalized linear models for the analysis of replicated
binary associations. \emph{Computational Statistics and Data Analysis, 38}, 217-233.

Meulders, M., De Boeck, P., and Van Mechelen, I. (2003). A taxonomy of latent structure
assumptions for probability matrix decomposition models. \emph{Psychometrika, 68}, 61-77.

Meulders, M., De Boeck, P., Van Mechelen, I., and Gelman, A. (2005). Probabilistic feature analysis of facial perception of emotions. 
\emph{Applied Statistics, 54}, 781-793.

Meulders, M., De Boeck, P., Van Mechelen, I., Gelman, A., and Maris, E. (2001). Bayesian inference with probability matrix decomposition models. 
\emph{Journal of Educational and Behavioral Statistics, 26}, 153-179.

Meulders, M., Gelman, A., Van Mechelen, I., and De Boeck P. (1998). Generalizing the probability matrix decomposition model: 
An example of Bayesian model checking and model expansion. 
In J. Hox, and E. De Leeuw (Eds.), \emph{Assumptions, robustness, and estimation methods in multivariate modeling} (pp. 1-19). 
TT Publicaties: Amsterdam.


}
\keyword{ package}
