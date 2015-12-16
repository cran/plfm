\name{summary.stepLCplfm}
\alias{summary.stepLCplfm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summary method for \code{\link{stepLCplfm}} objects}
\description{
  The function \code{summary.stepLCplfm} summarizes the fit measures of a series of \code{\link{LCplfm}} objects. 
}
\usage{\method{summary}{stepLCplfm}(object, \dots)}
\arguments{
  \item{object}{Latent class probabilistic latent feature analysis object returned by \code{\link{stepLCplfm}}}
  \item{\dots}{Further arguments are ignored}
}
\details{
  The summary of the \code{\link{stepLCplfm}} generates a table of fit measures (information criteria and descriptive fit measures) 
  for models with different numbers of latent class/features that are involved in the analysis.
}
\author{Michel Meulders}
\seealso{
\code{\link{stepLCplfm}}, \code{\link{print.stepLCplfm}}
}


