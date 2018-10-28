\name{plfm-package}
\alias{plfm-package}
\docType{package}
\title{Probabilistic Latent Feature Analysis}
\description{
Functions for estimating disjunctive, conjunctive or additive probabilistic latent feature models on (aggregated) binary three-way data
}
\details{
\tabular{ll}{
Package: \tab plfm\cr
Type: \tab Package\cr
Version: \tab 2.2.2\cr
Date: \tab 2018-10-27\cr
License: \tab GPL(>=2)\cr
LazyLoad: \tab yes\cr
}

Probabilistic latent feature models can be used to model three-way three-mode binary observations (e.g. persons who indicate for each of a number of products and for each of 
a set of attributes whether a product has a certain attribute). A basic probabilistic feature model (referred to as \code{\link{plfm}}) uses aggregated three-way three-mode binary 
data as input, namely the two-way two-mode frequency table that is obtained by summing the binary three-way three-mode data across persons. The basic probabilistic feature model 
(Maris, De Boeck and Van Mechelen, 1996) is based
on the assumption that observations are statistically independent and that model parameters are homogeneous across persons. The \code{\link{plfm}} function can be used to locate 
the posterior mode(s) of basic probabilistic feature models, and to compute information criteria for model selection, and measures of statistical and descriptive model fit. 
The \code{\link{stepplfm}} function can be used to fit a series of disjunctive, conjunctive or additive basic probabilistic feature models with different number of latent features. 
In addition, the \code{\link{bayesplfm}} function can be used to compute a sample of the posterior distribution of the basic probabilistic feature model 
in the neigbourhood of a specific posterior mode.
 
Latent class extensions of the probabilistic feature model (referred to as \code{\link{LCplfm}}) take binary three-way three-mode observations as input. In contrast to the basic probabilistic 
feature model, latent class probabilistic feature models allow to model dependencies between (subsets of) observations (Meulders, De Boeck and Van Mechelen, 2003) 
and/or to account for heterogeneity in model parameters across persons (Meulders, Tuerlinckx, and Vanpaemel, 2013).
The \code{\link{LCplfm}} function can be used to compute posterior mode estimates (of different types of) latent class probabilistic 
feature models as well as to compute information criteria for model selection, and measures of descriptive model fit. The \code{\link{stepLCplfm}} function can be used to compute a series of
latent class probabilistic feature models with different numbers of latent features and latent classes.

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

Meulders, M. and De Bruecker, P. (2018). Latent class probabilistic latent feature analysis of three-way three-mode binary data. 
\emph{Journal of Statistical Software, 87(1)}, 1-45.

Meulders, M., Gelman, A., Van Mechelen, I., and De Boeck P. (1998). Generalizing the probability matrix decomposition model: 
An example of Bayesian model checking and model expansion. 
In J. Hox, and E. De Leeuw (Eds.), \emph{Assumptions, robustness, and estimation methods in multivariate modeling} (pp. 1-19). 
TT Publicaties: Amsterdam.

Meulders, M., Tuerlinckx, F., and Vanpaemel, W. (2013). Constrained multilevel latent class models for the analysis of three-way three-mode binary data. 
\emph{Journal of Classification, 30 (3)}, 306-337.

}
\keyword{ package}
